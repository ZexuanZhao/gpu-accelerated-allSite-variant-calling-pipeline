#!/usr/bin/env python3
"""Per-sample genotype depth window from mosdepth region distributions.

m = median depth of the sample over the regions given to mosdepth (the largest depth c with
at least half of the bases at depth >= c, from the "total" rows of *.mosdepth.region.dist.txt).
  min_dp = max(--min-abs, ceil(--min-frac * m))
  max_dp = floor(--max-frac * m)
A value given in the sample sheet ("auto" = none) overrides the computed one.
Samples with m < --low-median-flag get flag "low_median" (they are kept).

Usage: depth_thresholds.py [--min-abs 8 --min-frac 0.5 --max-frac 2 --low-median-flag 10] \
           --sample NAME DIST MIN|auto MAX|auto [--sample ...] > depth_thresholds.tsv
"""
import argparse
import math
import sys


def median_depth(dist_path):
    """Median depth from a mosdepth *.region.dist.txt (cumulative proportions)."""
    best = None
    with open(dist_path) as fh:
        for line in fh:
            label, depth, prop = line.rstrip("\n").split("\t")
            if label != "total":
                continue
            if float(prop) >= 0.5:
                d = int(depth)
                best = d if best is None else max(best, d)
    if best is None:
        sys.exit(f"ERROR: no 'total' rows in {dist_path}")
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-abs", type=int, default=8)
    ap.add_argument("--min-frac", type=float, default=0.5)
    ap.add_argument("--max-frac", type=float, default=2.0)
    ap.add_argument("--low-median-flag", type=float, default=10)
    ap.add_argument("--sample", nargs=4, action="append", required=True,
                    metavar=("NAME", "DIST", "MIN", "MAX"))
    args = ap.parse_args()

    print("sample\tmedian_depth\tmin_dp\tmax_dp\tmin_dp_source\tmax_dp_source\tflag")
    for name, dist, min_override, max_override in args.sample:
        m = median_depth(dist)
        min_dp = max(args.min_abs, math.ceil(args.min_frac * m))
        max_dp = math.floor(args.max_frac * m)
        min_src = max_src = "auto"
        if min_override != "auto":
            min_dp, min_src = int(min_override), "sample_sheet"
        if max_override != "auto":
            max_dp, max_src = int(max_override), "sample_sheet"
        flags = []
        if m < args.low_median_flag:
            flags.append("low_median")
        if max_dp < min_dp:
            flags.append("empty_window")
            print(f"WARNING: {name}: max_dp {max_dp} < min_dp {min_dp}; all its genotypes will "
                  f"be masked", file=sys.stderr)
        print(f"{name}\t{m}\t{min_dp}\t{max_dp}\t{min_src}\t{max_src}\t{','.join(flags) or '.'}")


if __name__ == "__main__":
    main()
