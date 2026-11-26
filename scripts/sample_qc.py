#!/usr/bin/env python3
"""Per-sample QC table for one species.

Joins the depth thresholds, the called/heterozygous site counts (genome-wide and, if given,
per region label) and the allele balance (AB = ALT / (REF + ALT)) of heterozygous calls in the
HQSNPs VCF. Writes the table to stdout and an AB histogram (0.05 bins) to --ab-out.

A sample made of two individuals (mixed DNA) shows extra heterozygous calls whose AB is far
from 0.5; frac_het_ab_outer = fraction of its het calls with AB < 0.3 or AB > 0.7.

Usage: sample_qc.py --depth depth.tsv --counts counts.tsv [--counts-by-region by_region.tsv] \
           --hqsnps sp.HQSNPs.vcf.gz --ab-out ab_hist.tsv > sample_qc.tsv
"""
import argparse
import statistics
import subprocess
import sys
from collections import defaultdict

N_BINS = 20


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def ab_by_sample(vcf):
    """AB values of heterozygous calls per sample."""
    cmd = ["bcftools", "query", "-i", 'GT="het"', "-f", "[%SAMPLE\t%GT\t%AD\n]", vcf]
    ab = defaultdict(list)
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as p:
        for line in p.stdout:
            sample, gt, ad = line.rstrip("\n").split("\t")
            if "." in gt or len(set(gt.replace("|", "/").split("/"))) < 2:
                continue
            ad = ad.split(",")
            try:
                ref, alt = int(ad[0]), int(ad[1])
            except (ValueError, IndexError):
                continue
            if ref + alt > 0:
                ab[sample].append(alt / (ref + alt))
    if p.returncode:
        sys.exit(f"ERROR: bcftools query failed on {vcf}")
    return ab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--counts-by-region")
    ap.add_argument("--hqsnps", required=True)
    ap.add_argument("--ab-out", required=True)
    args = ap.parse_args()

    depth = {r["sample"]: r for r in read_tsv(args.depth)}
    counts = {r["sample"]: r for r in read_tsv(args.counts)}
    by_region = defaultdict(dict)
    regions = []
    if args.counts_by_region:
        for r in read_tsv(args.counts_by_region):
            by_region[r["sample"]][r["region"]] = r
            if r["region"] not in regions:
                regions.append(r["region"])
    ab = ab_by_sample(args.hqsnps)

    cols = ["sample", "median_depth", "min_dp", "max_dp", "flag",
            "called_sites", "het_sites", "frac_called", "het_rate",
            "n_het_ab", "median_het_ab", "frac_het_ab_outer"]
    for reg in regions:
        cols += [f"called_sites[{reg}]", f"het_rate[{reg}]"]
    print("\t".join(cols))
    for sample, d in depth.items():
        c = counts.get(sample, {})
        called = int(c.get("called_sites", 0) or 0)
        het = int(c.get("het_sites", 0) or 0)
        vals = ab.get(sample, [])
        outer = sum(1 for v in vals if v < 0.3 or v > 0.7)
        row = [sample, d["median_depth"], d["min_dp"], d["max_dp"], d["flag"],
               str(called), str(het), c.get("frac_called", "NA"),
               f"{het / called:.6g}" if called else "NA",
               str(len(vals)),
               f"{statistics.median(vals):.3f}" if vals else "NA",
               f"{outer / len(vals):.4f}" if vals else "NA"]
        for reg in regions:
            r = by_region[sample].get(reg, {})
            rc = int(r.get("called_sites", 0) or 0)
            rh = int(r.get("het_sites", 0) or 0)
            row += [str(rc), f"{rh / rc:.6g}" if rc else "NA"]
        print("\t".join(row))

    with open(args.ab_out, "w") as out:
        out.write("sample\tab_bin_start\tab_bin_end\tn_het\n")
        for sample in depth:
            hist = [0] * N_BINS
            for v in ab.get(sample, []):
                hist[min(int(v * N_BINS), N_BINS - 1)] += 1
            for i, n in enumerate(hist):
                out.write(f"{sample}\t{i / N_BINS:.2f}\t{(i + 1) / N_BINS:.2f}\t{n}\n")


if __name__ == "__main__":
    main()
