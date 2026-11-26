#!/usr/bin/env python3
"""Filter an all-sites VCF (GATK GenotypeGVCFs -all-sites) into callable sites.

Reads VCF text on stdin and writes VCF text on stdout: only records that pass are written,
with FILTER set to PASS. Two kinds of record are kept, everything else is dropped:

  invariant sites (ALT ".", REF one base)
      genotype: missing unless homozygous reference, RGQ (GQ if RGQ is absent) >= --rgq-min
                and min_dp <= DP <= max_dp of that sample
  biallelic SNPs (REF and ALT one base each, ALT not "*")
      site:     FILTER must be PASS or "." (e.g. SnpGap from bcftools filter -g fails);
                QUAL >= --qual-min, QD >= --qd-min, FS <= --fs-max, SOR <= --sor-max,
                MQ >= --mq-min, MQRankSum >= --mqranksum-min,
                ReadPosRankSum >= --readposranksum-min, ExcessHet <= --excesshet-max
                (a missing annotation passes, as in the GATK hard-filter recommendations)
      genotype: missing unless GQ >= --gq-min and min_dp <= DP <= max_dp; a heterozygous
                call also needs ALT reads >= --het-alt-min and
                --het-ab-min <= ALT / (REF + ALT) <= --het-ab-max
      site, after masking: heterozygous fraction of called genotypes <= --het-fraction-max

Both kinds are then dropped when the fraction of missing genotypes exceeds --f-missing-max.
The same depth window and quality threshold apply to reference calls and to variant calls,
so the callable denominator and the variant numerator are filtered alike.

INFO/MLEAC and INFO/MLEAF are removed; AC, AN, ... are not updated here: pipe the output
through bcftools +fill-tags.
Counts of dropped sites and masked genotypes go to --stats (TSV).

Usage:
  bcftools view in.vcf.gz | filter_allsites.py --depth-thresholds depth.tsv [options] \
      --stats stats.tsv | bgzip > out.vcf.gz
depth.tsv: TSV with header and columns sample, min_dp, max_dp (other columns ignored).
"""
import argparse
import sys
from collections import Counter

BASES = {"A", "C", "G", "T"}
INFO_TESTS = (  # (INFO key, option, direction): the site fails when value < min / > max
    ("QD", "qd_min", "min"),
    ("FS", "fs_max", "max"),
    ("SOR", "sor_max", "max"),
    ("MQ", "mq_min", "min"),
    ("MQRankSum", "mqranksum_min", "min"),
    ("ReadPosRankSum", "readposranksum_min", "min"),
    ("ExcessHet", "excesshet_max", "max"),
)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--depth-thresholds", required=True)
    ap.add_argument("--stats", required=True)
    ap.add_argument("--qual-min", type=float, default=30)
    ap.add_argument("--qd-min", type=float, default=2.0)
    ap.add_argument("--fs-max", type=float, default=60.0)
    ap.add_argument("--sor-max", type=float, default=3.0)
    ap.add_argument("--mq-min", type=float, default=50.0)
    ap.add_argument("--mqranksum-min", type=float, default=-12.5)
    ap.add_argument("--readposranksum-min", type=float, default=-8.0)
    ap.add_argument("--excesshet-max", type=float, default=54.69)
    ap.add_argument("--het-fraction-max", type=float, default=0.5)
    ap.add_argument("--gq-min", type=float, default=20)
    ap.add_argument("--rgq-min", type=float, default=20)
    ap.add_argument("--het-alt-min", type=int, default=3)
    ap.add_argument("--het-ab-min", type=float, default=0.2)
    ap.add_argument("--het-ab-max", type=float, default=0.8)
    ap.add_argument("--f-missing-max", type=float, default=0.5)
    return ap.parse_args(argv)


def read_depth_thresholds(path):
    thresholds = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for c in ("sample", "min_dp", "max_dp"):
            if c not in col:
                sys.exit(f"ERROR: {path} lacks column {c}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f == [""]:
                continue
            thresholds[f[col["sample"]]] = (int(f[col["min_dp"]]), int(f[col["max_dp"]]))
    return thresholds


def to_int(x):
    try:
        return int(x)
    except ValueError:
        return None


def to_float(x):
    try:
        return float(x)
    except ValueError:
        return None


def parse_info(info, keys):
    out = {}
    for item in info.split(";"):
        k, _, v = item.partition("=")
        if k in keys:
            out[k] = v
    return out


def missing_gt(gt):
    """Missing genotype with the same ploidy and separator as gt."""
    sep = "|" if "|" in gt else "/"
    n = gt.count("/") + gt.count("|") + 1
    return sep.join(["."] * n)


class SiteFilter:
    def __init__(self, args, samples, thresholds):
        self.a = args
        self.samples = samples
        missing = [s for s in samples if s not in thresholds]
        if missing:
            sys.exit(f"ERROR: no depth thresholds for sample(s): {', '.join(missing)}")
        self.min_dp = [thresholds[s][0] for s in samples]
        self.max_dp = [thresholds[s][1] for s in samples]
        self.n = len(samples)
        self.sites = Counter()      # (kind, outcome) -> count
        self.gt_mask = Counter()    # (sample, kind, reason) -> count
        self.info_keys = {k for k, _, _ in INFO_TESTS}

    def site_fails_info(self, info):
        vals = parse_info(info, self.info_keys)
        for key, opt, direction in INFO_TESTS:
            v = to_float(vals.get(key, "."))
            if v is None:
                continue
            limit = getattr(self.a, opt)
            if (direction == "min" and v < limit) or (direction == "max" and v > limit):
                return key
        return None

    def mask_genotypes(self, f, kind):
        """Mask failing genotypes in place; return (n_called, n_het)."""
        fmt = f[8].split(":")
        idx = {k: i for i, k in enumerate(fmt)}
        i_dp = idx.get("DP")
        i_ad = idx.get("AD")
        i_q = idx.get("GQ") if kind == "snp" else idx.get("RGQ", idx.get("GQ"))
        q_min = self.a.gq_min if kind == "snp" else self.a.rgq_min
        called = het = 0
        for s in range(self.n):
            parts = f[9 + s].split(":")
            gt = parts[0]
            if "." in gt:
                continue
            alleles = gt.replace("|", "/").split("/")
            reason = None
            if kind == "inv" and any(a != "0" for a in alleles):
                reason = "not_ref"
            if reason is None:
                q = to_float(parts[i_q]) if i_q is not None and i_q < len(parts) else None
                if q is None or q < q_min:
                    reason = "low_gq"
            if reason is None:
                dp = to_int(parts[i_dp]) if i_dp is not None and i_dp < len(parts) else None
                if dp is None or dp < self.min_dp[s]:
                    reason = "low_dp"
                elif dp > self.max_dp[s]:
                    reason = "high_dp"
            is_het = len(set(alleles)) > 1
            if reason is None and is_het:
                ad = parts[i_ad].split(",") if i_ad is not None and i_ad < len(parts) else []
                ref_ad = to_int(ad[0]) if len(ad) > 0 else None
                alt_ad = to_int(ad[1]) if len(ad) > 1 else None
                if ref_ad is None or alt_ad is None:
                    reason = "no_ad"
                elif alt_ad < self.a.het_alt_min:
                    reason = "het_low_alt"
                else:
                    ab = alt_ad / (ref_ad + alt_ad) if ref_ad + alt_ad > 0 else 0.0
                    if ab < self.a.het_ab_min or ab > self.a.het_ab_max:
                        reason = "het_ab"
            if reason is not None:
                parts[0] = missing_gt(gt)
                f[9 + s] = ":".join(parts)
                self.gt_mask[(self.samples[s], kind, reason)] += 1
                continue
            called += 1
            het += is_het
        return called, het

    def process(self, line):
        """Return the output line or None."""
        f = line.rstrip("\n").split("\t")
        ref, alt, filt = f[3], f[4], f[6]
        if alt == "." and len(ref) == 1:
            kind = "inv"
        elif len(ref) == 1 and alt in BASES:
            kind = "snp"
        else:
            if "," in alt:
                self.sites[("other", "drop_multiallelic")] += 1
            elif alt == "*" or len(ref) != len(alt):
                self.sites[("other", "drop_indel")] += 1
            else:
                self.sites[("other", "drop_other")] += 1
            return None

        if kind == "snp":
            if filt not in ("PASS", "."):
                self.sites[(kind, f"drop_filter_{filt}")] += 1
                return None
            qual = to_float(f[5])
            if qual is None or qual < self.a.qual_min:
                self.sites[(kind, "drop_QUAL")] += 1
                return None
            failed = self.site_fails_info(f[7])
            if failed:
                self.sites[(kind, f"drop_{failed}")] += 1
                return None

        called, het = self.mask_genotypes(f, kind)
        if called == 0:
            self.sites[(kind, "drop_all_missing")] += 1
            return None
        if 1 - called / self.n > self.a.f_missing_max:
            self.sites[(kind, "drop_F_MISSING")] += 1
            return None
        if kind == "snp" and het / called > self.a.het_fraction_max:
            self.sites[(kind, "drop_het_fraction")] += 1
            return None
        f[6] = "PASS"
        f[7] = strip_info(f[7])
        self.sites[(kind, "kept")] += 1
        return "\t".join(f) + "\n"

    def write_stats(self, path):
        with open(path, "w") as out:
            out.write("level\tsample\tkind\toutcome\tcount\n")
            for (kind, outcome), c in sorted(self.sites.items()):
                out.write(f"site\t.\t{kind}\t{outcome}\t{c}\n")
            for (sample, kind, reason), c in sorted(self.gt_mask.items()):
                out.write(f"genotype\t{sample}\t{kind}\tmask_{reason}\t{c}\n")


## Header lines that name the interval being processed (command lines with -L / file paths).
## They are dropped so that all intervals share one header and can be joined with
## bcftools concat --naive.
## MLEAC/MLEAF (GATK) are stale once genotypes are masked; they are removed from the header
## and from INFO of the written records.
DROP_HEADER_PREFIXES = ("##GATKCommandLine", "##bcftools_",
                        "##INFO=<ID=MLEAC,", "##INFO=<ID=MLEAF,")
DROP_INFO_KEYS = ("MLEAC", "MLEAF")


def strip_info(info):
    kept = [x for x in info.split(";") if x.partition("=")[0] not in DROP_INFO_KEYS]
    return ";".join(kept) if kept else "."


def threshold_header(args):
    keys = [k for k in vars(args) if k not in ("depth_thresholds", "stats")]
    return "##filter_allsites=" + ",".join(f"{k}={getattr(args, k)}" for k in keys) + \
        ",genotype_depth=per-sample min_dp/max_dp\n"


def main(argv=None, stdin=sys.stdin, stdout=sys.stdout):
    args = parse_args(argv)
    thresholds = read_depth_thresholds(args.depth_thresholds)
    sf = None
    for line in stdin:
        if line.startswith("##"):
            if not line.startswith(DROP_HEADER_PREFIXES):
                stdout.write(line)
            continue
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            sf = SiteFilter(args, samples, thresholds)
            stdout.write(threshold_header(args))
            stdout.write(line)
            continue
        if sf is None:
            sys.exit("ERROR: no #CHROM header line before the records")
        out = sf.process(line)
        if out is not None:
            stdout.write(out)
    if sf is None:
        sys.exit("ERROR: empty VCF")
    sf.write_stats(args.stats)


if __name__ == "__main__":
    main()
