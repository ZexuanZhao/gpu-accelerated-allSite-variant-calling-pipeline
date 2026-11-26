#!/usr/bin/env python3
"""Cross-species QC summary for MultiQC (custom content) and sample renaming.

Per species it reads: <sp>.sample_qc.tsv, <sp>.filter_stats.tsv, plink2 KING table (.kin0),
PCA (.eigenvec, .eigenval), the pixy pi file of the largest window, and samtools stats files.
Writes to --outdir:
  sample_flags_mqc.tsv     one row per sample: alignment, depth, heterozygosity, AB, kinship,
                           PCA, and status PASS/WARN with the flags that fired
  kinship_mqc.tsv          sample pairs with KING kinship >= --kinship-report-min
  pca_mqc.tsv              PC1-PC3 per sample, variance explained, outlier flag
  pca_<species>_mqc.json   PC1 vs PC2 scatter per species
  filter_stats_mqc.tsv     per species: sites kept and dropped per filter
  pi_mqc.tsv               per species: genome-wide pi (sum of differences / sum of comparisons)
  multiqc_sample_names.yaml  MultiQC config renaming <sample>* to <species>_<sample>*

Flags (each relative to the species median unless stated):
  low_mapping       mapped / total reads < --mapping-rate-min (absolute)
  high_mismatch     samtools error rate > --error-rate-rel-max x median
  low_depth         flag from depth_thresholds (median depth below the configured floor)
  het_low           het rate < --het-rate-rel-low x median (haploid male? inbred?)
  het_high          het rate > --het-rate-rel-high x median (mixed DNA?)
  ab_skew           fraction of het calls with AB < 0.3 or > 0.7 > median + --ab-outer-excess
  possible_duplicate  KING kinship >= --kinship-duplicate with another sample (full sisters in
                    haplodiploids are ~0.375, duplicates ~0.5)
  pca_outlier       |PC1| or |PC2| more than --pca-outlier-sd robust SDs (1.4826 x MAD) from
                    the species median

Usage: qc_summary.py --outdir DIR [thresholds] \
    --species SP SAMPLE_QC FILTER_STATS KIN0 EIGENVEC EIGENVAL PI SAMTOOLS_DIR [--species ...]
"""
import argparse
import json
import os
import re
import statistics
from collections import defaultdict


def read_tsv(path, comment_header=False):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if comment_header:
            header[0] = header[0].lstrip("#")
        return [dict(zip(header, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def samtools_sn(path):
    sn = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("SN\t"):
                f = line.rstrip("\n").split("\t")
                sn[f[1].rstrip(":")] = fnum(f[2])
    return sn


def median(values):
    v = [x for x in values if x is not None]
    return statistics.median(v) if v else None


def robust_z(values):
    """{key: robust z} using median and 1.4826 x MAD; empty when MAD is 0."""
    v = [x for x in values.values() if x is not None]
    if len(v) < 4:
        return {}
    med = statistics.median(v)
    mad = statistics.median([abs(x - med) for x in v]) * 1.4826
    if mad == 0:
        return {}
    return {k: (x - med) / mad for k, x in values.items() if x is not None}


def fmt(x, nd=4):
    if x is None:
        return "NA"
    if isinstance(x, float):
        return f"{x:.{nd}g}"
    return str(x)


def mqc_header(id_, name, desc, plot="table"):
    return (f"# id: '{id_}'\n# section_name: '{name}'\n# description: '{desc}'\n"
            f"# plot_type: '{plot}'\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mapping-rate-min", type=float, default=0.90)
    ap.add_argument("--error-rate-rel-max", type=float, default=2.0)
    ap.add_argument("--het-rate-rel-low", type=float, default=0.2)
    ap.add_argument("--het-rate-rel-high", type=float, default=2.0)
    ap.add_argument("--ab-outer-excess", type=float, default=0.15)
    ap.add_argument("--kinship-duplicate", type=float, default=0.45)
    ap.add_argument("--kinship-report-min", type=float, default=0.0884)
    ap.add_argument("--pca-outlier-sd", type=float, default=3.0)
    ap.add_argument("--species", nargs=8, action="append", required=True,
                    metavar=("SP", "SAMPLE_QC", "FILTER_STATS", "KIN0", "EIGENVEC",
                             "EIGENVAL", "PI", "SAMTOOLS_DIR"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    rows, kin_rows, pca_rows, filt_rows, pi_rows = [], [], [], [], []
    rename = {}
    for sp, sqc_path, filt_path, kin_path, vec_path, val_path, pi_path, st_dir in a.species:
        sqc = read_tsv(sqc_path)
        samples = [r["sample"] for r in sqc]
        name = {s: f"{sp}_{s}" for s in samples}
        for s in samples:
            rename[s] = name[s]

        # alignment
        aln = {}
        for s in samples:
            sn = samtools_sn(os.path.join(st_dir, f"{s}.stats"))
            total = sn.get("raw total sequences") or 0
            aln[s] = {
                "reads_M": total / 1e6 if total else None,
                "mapping_rate": (sn.get("reads mapped") or 0) / total if total else None,
                "error_rate": sn.get("error rate"),
                "dup_rate": (sn.get("reads duplicated") or 0) / total if total else None,
                "insert_size": sn.get("insert size average"),
            }

        # kinship
        best = {s: (None, None) for s in samples}
        for r in read_tsv(kin_path, comment_header=True):
            i1, i2 = r.get("IID1"), r.get("IID2")
            k = fnum(r.get("KINSHIP"))
            if k is None:
                continue
            for x, y in ((i1, i2), (i2, i1)):
                if x in best and (best[x][0] is None or k > best[x][0]):
                    best[x] = (k, y)
            if k >= a.kinship_report_min:
                kin_rows.append((f"{sp}: {name.get(i1, i1)} | {name.get(i2, i2)}", sp, k,
                                 r.get("NSNP", "NA"), r.get("IBS0", "NA")))

        # PCA
        vec = read_tsv(vec_path, comment_header=True)
        vals = [fnum(l.strip()) for l in open(val_path) if l.strip()]
        tot = sum(v for v in vals if v)
        var_pct = [100 * v / tot if tot else None for v in vals]
        pcs = {r["IID"]: [fnum(r.get(f"PC{i}")) for i in (1, 2, 3)] for r in vec}
        z1 = robust_z({s: pcs.get(s, [None])[0] for s in samples})
        z2 = robust_z({s: (pcs.get(s, [None, None])[1] if len(pcs.get(s, [])) > 1 else None)
                       for s in samples})
        pca_out = {s: max(abs(z1.get(s, 0)), abs(z2.get(s, 0))) > a.pca_outlier_sd for s in samples}
        scatter = {}
        for s in samples:
            p = pcs.get(s, [None, None, None])
            pca_rows.append((name[s], sp, p[0], p[1], p[2],
                             var_pct[0] if var_pct else None,
                             var_pct[1] if len(var_pct) > 1 else None,
                             "yes" if pca_out[s] else "no"))
            if p[0] is not None and p[1] is not None:
                scatter[name[s]] = {"x": p[0], "y": p[1]}
        safe = re.sub(r"[^A-Za-z0-9._-]", "_", sp)
        with open(os.path.join(a.outdir, f"pca_{safe}_mqc.json"), "w") as out:
            json.dump({"id": f"pca_{safe}", "section_name": f"PCA {sp}",
                       "description": "PC1 vs PC2 of LD-pruned HQSNPs (plink2).",
                       "plot_type": "scatter",
                       "pconfig": {"id": f"pca_{safe}_plot", "title": f"PCA {sp}",
                                   "xlab": f"PC1 ({fmt(var_pct[0] if var_pct else None, 3)}%)",
                                   "ylab": f"PC2 ({fmt(var_pct[1] if len(var_pct) > 1 else None, 3)}%)"},
                       "data": scatter}, out, indent=1)

        # species medians and flags
        het = {r["sample"]: fnum(r["het_rate"]) for r in sqc}
        outer = {r["sample"]: fnum(r["frac_het_ab_outer"]) for r in sqc}
        m_err = median([aln[s]["error_rate"] for s in samples])
        m_het = median(list(het.values()))
        m_outer = median(list(outer.values()))
        for r in sqc:
            s = r["sample"]
            flags = []
            mr = aln[s]["mapping_rate"]
            if mr is not None and mr < a.mapping_rate_min:
                flags.append("low_mapping")
            er = aln[s]["error_rate"]
            if er is not None and m_err and er > a.error_rate_rel_max * m_err:
                flags.append("high_mismatch")
            if "low_median" in r.get("flag", "") or "empty_window" in r.get("flag", ""):
                flags.append("low_depth")
            h = het[s]
            if h is not None and m_het:
                if h < a.het_rate_rel_low * m_het:
                    flags.append("het_low")
                elif h > a.het_rate_rel_high * m_het:
                    flags.append("het_high")
            o = outer[s]
            if o is not None and m_outer is not None and o > m_outer + a.ab_outer_excess:
                flags.append("ab_skew")
            k, partner = best[s]
            if k is not None and k >= a.kinship_duplicate:
                flags.append(f"possible_duplicate:{partner}")
            if pca_out[s]:
                flags.append("pca_outlier")
            p = pcs.get(s, [None, None, None])
            rows.append((name[s], sp, "WARN" if flags else "PASS", ",".join(flags) or ".",
                         aln[s]["reads_M"], mr, er, aln[s]["dup_rate"], aln[s]["insert_size"],
                         fnum(r["median_depth"]), r["min_dp"], r["max_dp"], fnum(r["frac_called"]),
                         h, fnum(r["median_het_ab"]), o, k,
                         name.get(partner, partner) if partner else "NA", p[0], p[1]))

        # filter stats (site level)
        counts = defaultdict(int)
        for r in read_tsv(filt_path):
            if r["level"] == "site":
                counts[f"{r['kind']}:{r['outcome']}"] += int(r["count"])
        filt_rows.append((sp, counts))

        # pi
        diffs = comps = sites = 0.0
        win_pi = []
        for r in read_tsv(pi_path):
            d, c, n = fnum(r.get("count_diffs")), fnum(r.get("count_comparisons")), fnum(r.get("no_sites"))
            if d is not None and c:
                diffs += d
                comps += c
                win_pi.append(d / c)
            if n:
                sites += n
        pi_rows.append((sp, len(samples), diffs / comps if comps else None,
                        median(win_pi), int(sites)))

    cols = ["species", "status", "flags", "reads_M", "mapping_rate", "error_rate", "dup_rate",
            "insert_size", "median_depth", "min_dp", "max_dp", "frac_called", "het_rate",
            "median_het_ab", "frac_het_ab_outer", "max_kinship", "kinship_partner", "PC1", "PC2"]
    with open(os.path.join(a.outdir, "sample_flags_mqc.tsv"), "w") as out:
        out.write(mqc_header("sample_flags", "Sample QC flags",
                             "One row per sample; status WARN lists the flags that fired "
                             "(see scripts/qc_summary.py)."))
        out.write("Sample\t" + "\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(fmt(x) for x in r) + "\n")

    with open(os.path.join(a.outdir, "kinship_mqc.tsv"), "w") as out:
        out.write(mqc_header("kinship", "KING kinship",
                             f"Sample pairs with KING kinship >= {a.kinship_report_min}. In "
                             "haplodiploids full sisters are ~0.375 and duplicates ~0.5; few samples "
                             "and population structure make the values noisy."))
        out.write("Pair\tspecies\tkinship\tn_snp\tibs0\n")
        for r in sorted(kin_rows, key=lambda x: -x[2]):
            out.write("\t".join(fmt(x) for x in r) + "\n")

    with open(os.path.join(a.outdir, "pca_mqc.tsv"), "w") as out:
        out.write(mqc_header("pca_table", "PCA per species",
                             f"PC1-PC3 of LD-pruned HQSNPs; outlier = more than {a.pca_outlier_sd} "
                             "robust SDs from the species median on PC1 or PC2."))
        out.write("Sample\tspecies\tPC1\tPC2\tPC3\tPC1_var_pct\tPC2_var_pct\toutlier\n")
        for r in pca_rows:
            out.write("\t".join(fmt(x) for x in r) + "\n")

    keys = sorted({k for _, c in filt_rows for k in c})
    with open(os.path.join(a.outdir, "filter_stats_mqc.tsv"), "w") as out:
        out.write(mqc_header("filter_stats", "Site filters",
                             "Sites kept and dropped per filter (kind:outcome), summed over intervals."))
        out.write("Species\t" + "\t".join(keys) + "\n")
        for sp, c in filt_rows:
            out.write(sp + "\t" + "\t".join(str(c.get(k, 0)) for k in keys) + "\n")

    with open(os.path.join(a.outdir, "pi_mqc.tsv"), "w") as out:
        out.write(mqc_header("pi_summary", "Nucleotide diversity",
                             "Genome-wide pi = sum of differences / sum of comparisons over pixy "
                             "windows (largest window size); median of window pi; sites with >= 1 "
                             "genotype."))
        out.write("Species\tn_samples\tpi_genome\tpi_median_window\tn_sites\n")
        for r in pi_rows:
            out.write("\t".join(fmt(x, 5) for x in r) + "\n")

    with open(os.path.join(a.outdir, "multiqc_sample_names.yaml"), "w") as out:
        out.write("# Written by qc_summary.py: <sample> and <sample>.<suffix> -> <species>_<sample>...\n")
        out.write("sample_names_replace_regex: true\nsample_names_replace:\n")
        for s, new in sorted(rename.items(), key=lambda x: -len(x[0])):
            pat = "^" + re.escape(s) + r"(?=$|[._])"
            out.write(f"  {json.dumps(pat)}: {json.dumps(new)}\n")


if __name__ == "__main__":
    main()
