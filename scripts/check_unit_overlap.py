#!/usr/bin/env python3
"""Check that the read pairs (units) of one sample are different reads.

Reads the first --n-reads read IDs of each R1 file and reports, for every pair of files, the
fraction of IDs they share. Exits with status 1 when any pair shares more than --max-shared
(the same sequencing run delivered twice, e.g. demultiplexed with two index settings).
Files are compared from their start, which works for FASTQ in instrument (tile) order.

Usage: check_unit_overlap.py --out report.tsv [--n-reads N] [--max-shared F] R1 [R1 ...]
"""
import argparse
import gzip
import itertools
import sys


def read_ids(path, n):
    ids = set()
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                name = line[1:].split()[0]
                if name.endswith("/1") or name.endswith("/2"):
                    name = name[:-2]
                ids.add(name)
                if len(ids) >= n:
                    break
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("r1", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-reads", type=int, default=200000)
    ap.add_argument("--max-shared", type=float, default=0.01)
    args = ap.parse_args()

    ids = {f: read_ids(f, args.n_reads) for f in args.r1}
    failed = []
    with open(args.out, "w") as out:
        out.write("file_a\tfile_b\tn_a\tn_b\tshared\tfrac_shared\tstatus\n")
        if len(args.r1) == 1:
            f = args.r1[0]
            out.write(f"{f}\t.\t{len(ids[f])}\t0\t0\t0\tsingle_unit\n")
        for a, b in itertools.combinations(args.r1, 2):
            shared = len(ids[a] & ids[b])
            frac = shared / max(1, min(len(ids[a]), len(ids[b])))
            status = "FAIL" if frac > args.max_shared else "ok"
            if status == "FAIL":
                failed.append((a, b, frac))
            out.write(f"{a}\t{b}\t{len(ids[a])}\t{len(ids[b])}\t{shared}\t{frac:.4f}\t{status}\n")
    if failed:
        for a, b, frac in failed:
            print(f"ERROR: {a} and {b} share {frac:.1%} of read IDs: the same reads are listed "
                  f"twice. Keep one of them in the sample sheet.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
