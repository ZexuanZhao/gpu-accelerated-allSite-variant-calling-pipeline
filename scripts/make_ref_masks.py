#!/usr/bin/env python3
"""Write BED files of low-complexity regions (from a dustmasker acclist) and N gaps.

dustmasker reports 0-based intervals with an inclusive end (a run at 0-based 100..159 is
printed as "100 159"), so the BED end is end + 1.

Usage: make_ref_masks.py --fasta ref.fa --acclist ref.lcm.acclist --lcm-out lcm.bed --gaps-out N.bed
"""
import argparse


def acclist_to_bed(acclist, out_path):
    with open(acclist) as fh, open(out_path, "w") as out:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            name, start, end = line.split("\t")
            name = name.lstrip(">").split()[0]
            out.write(f"{name}\t{int(start)}\t{int(end) + 1}\n")


def gaps_to_bed(fasta, out_path):
    with open(fasta) as fh, open(out_path, "w") as out:
        name, pos, gap_start = None, 0, None

        def close(gap_start, pos):
            if gap_start is not None:
                out.write(f"{name}\t{gap_start}\t{pos}\n")

        for line in fh:
            if line.startswith(">"):
                close(gap_start, pos)
                name, pos, gap_start = line[1:].split()[0], 0, None
                continue
            for ch in line.rstrip("\n\r"):
                if ch in "Nn":
                    if gap_start is None:
                        gap_start = pos
                elif gap_start is not None:
                    out.write(f"{name}\t{gap_start}\t{pos}\n")
                    gap_start = None
                pos += 1
        close(gap_start, pos)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--acclist", required=True)
    ap.add_argument("--lcm-out", required=True)
    ap.add_argument("--gaps-out", required=True)
    args = ap.parse_args()
    acclist_to_bed(args.acclist, args.lcm_out)
    gaps_to_bed(args.fasta, args.gaps_out)


if __name__ == "__main__":
    main()
