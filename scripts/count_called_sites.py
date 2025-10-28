#!/usr/bin/env python3
"""
Count per-sample *called* sites (non-missing genotypes) from an all-sites VCF/VCF.GZ.
No external deps (no pysam). Streams line-by-line; works with bgzip/gzip/plain VCF.

Usage:
  python count_called_sites.py input.vcf.gz > called_counts.tsv

Options:
  --pass-only      Only count records with FILTER == PASS (default: count all).
"""
import sys
import argparse
import gzip

def open_any(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def nth_field(s: str, idx: int, sep: str = ":") -> str:
    """Return field #idx from s split by sep without allocating all fields."""
    if idx == 0:
        i = s.find(sep)
        return s if i == -1 else s[:i]
    k = 0
    start = 0
    for i, ch in enumerate(s):
        if ch == sep:
            if k == idx:
                return s[start:i]
            k += 1
            start = i + 1
    # last or only field
    return s[start:] if k == idx else ""

def main():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("vcf")
    ap.add_argument("--pass-only", action="store_true")
    args = ap.parse_args()

    with open_any(args.vcf) as fh:
        sample_names = None
        gt_idx = None
        n_samples = 0
        called = None
        total_sites = 0  # number of VCF records considered

        for line in fh:
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                sample_names = cols[9:]
                n_samples = len(sample_names)
                called = [0] * n_samples
                continue

            if sample_names is None:
                # malformed VCF (no header)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue  # not a genotype line

            # Optional: only PASS
            if args.pass_only:
                flt = fields[6]
                if flt != "PASS":
                    continue

            fmt = fields[8]

            # Find GT index once; if absent, skip this record
            if gt_idx is None:
                fmt_keys = fmt.split(":")
                gt_idx = fmt_keys.index("GT") if "GT" in fmt_keys else -1
                if gt_idx == -1:
                    # No genotype field; nothing to count in this VCF
                    print("ERROR: FORMAT has no GT field; cannot count calls.", file=sys.stderr)
                    sys.exit(1)

            # Count this record toward totals for every sample
            total_sites += 1

            # Iterate samples and count non-missing GTs
            # Called if GT contains no '.' (e.g., '0/0', '0|1', '1', but NOT './.', '.|.', '0/.', etc.)
            # Faster path when GT is first field
            if gt_idx == 0:
                for i, s in enumerate(fields[9:]):
                    # grab up to first ':' (or whole string if no ':')
                    j = s.find(":")
                    gt = s if j == -1 else s[:j]
                    if "." not in gt and gt != "":
                        called[i] += 1
            else:
                for i, s in enumerate(fields[9:]):
                    gt = nth_field(s, gt_idx, ":")
                    if "." not in gt and gt != "":
                        called[i] += 1

        if sample_names is None:
            print("ERROR: No #CHROM header found; is this a valid VCF?", file=sys.stderr)
            sys.exit(1)

        # Output TSV: sample, called_sites, missing_sites, total_sites, frac_called
        out = sys.stdout
        print("\t".join(["sample", "called_sites", "missing_sites", "total_sites", "frac_called"]), file=out)
        for sname, c in zip(sample_names, called):
            miss = total_sites - c
            frac = (c / total_sites) if total_sites > 0 else 0.0
            print(f"{sname}\t{c}\t{miss}\t{total_sites}\t{frac:.6f}", file=out)

if __name__ == "__main__":
    main()

