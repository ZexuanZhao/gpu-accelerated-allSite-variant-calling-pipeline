#!/usr/bin/env python3
"""
Fast per-sample call counts (and heterozygous counts) from an indexed VCF/BCF using bcftools.
Counts sites where GT is present and not missing (no '.' in GT). Heterozygous = alleles in GT differ.

Requires:
  - VCF.gz + .tbi (or BCF + .csi)
  - bcftools in PATH

Usage:
  python count_called_sites_bcftools.py file.vcf.gz > called_counts.tsv
Options:
  --pass-only       only count FILTER=PASS sites
  --threads N       (optional) threads for `bcftools view` decompression
"""
import argparse, subprocess, sys

def run_stdout(cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout

def main():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("vcf")
    ap.add_argument("--pass-only", action="store_true")
    ap.add_argument("--threads", type=int, default=0)
    args = ap.parse_args()

    # Sample list
    samples = run_stdout(["bcftools", "query", "-l", args.vcf]).strip().splitlines()
    if not samples:
        print("ERROR: no samples found", file=sys.stderr); sys.exit(1)
    n = len(samples)

    # Build `bcftools view` (filter + optional threads) → uncompressed BCF on stdout
    view_cmd = ["bcftools", "view", "-Ou", args.vcf]
    if args.pass_only:
        view_cmd[2:2] = ["-i", 'FILTER="PASS"']
    if args.threads and args.threads > 0:
        view_cmd[2:2] = ["--threads", str(args.threads)]

    # Pipe into `bcftools query` to stream only GTs
    query_cmd = ["bcftools", "query", "-f", "[%GT\t]\n", "-"]

    p_view  = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)
    p_query = subprocess.Popen(query_cmd, stdin=p_view.stdout,
                               stdout=subprocess.PIPE, text=True, bufsize=1)

    called = [0]*n
    het    = [0]*n
    total  = 0

    for line in p_query.stdout:
        row = line.rstrip("\n").split("\t")
        for i, gt in enumerate(row):
            if not gt or "." in gt:
                continue  # missing genotype like ./., .|1, etc.
            called[i] += 1
            # Heterozygous if any two alleles differ (handles phased/unphased, any ploidy)
            alleles = gt.replace("|", "/").split("/")
            if len(set(alleles)) > 1:
                het[i] += 1
        total += 1

    p_query.stdout.close()
    p_query.wait()
    p_view.wait()

    print("sample\tcalled_sites\thet_sites\tmissing_sites\ttotal_sites\tfrac_called")
    for s, c, h in zip(samples, called, het):
        miss = total - c
        frac = (c/total) if total else 0.0
        print(f"{s}\t{c}\t{h}\t{miss}\t{total}\t{frac:.6f}")

if __name__ == "__main__":
    main()

