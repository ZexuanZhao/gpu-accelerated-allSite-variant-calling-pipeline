import csv
import os
import re

## ---------------------------------------------------------------------
## Configuration, species table and sample sheet
## ---------------------------------------------------------------------
OUTDIR = config["outdir"]
SPECIES_CFG = config["species"]
ENABLED_SPECIES = [s for s, v in SPECIES_CFG.items() if v.get("enabled", True)]
if not ENABLED_SPECIES:
    raise ValueError("No species is enabled in the config.")


def load_sample_sheet(path):
    """Read the sample sheet (one row per read pair) and validate it.

    Columns: species_id, sample, ploidy, r1, r2 and optionally min_dp, max_dp, include
    (yes/no, default yes; all rows of a sample must agree) and note (free text, not used).
    Samples with include = no are left out of the run.
    Returns {species: {sample: {"ploidy": int, "min_dp": str, "max_dp": str,
                                "units": [(r1, r2), ...]}}} for enabled species only.
    """
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if any((v or "").strip() for v in r.values())]
    required = ["species_id", "sample", "ploidy", "r1", "r2"]
    if not rows:
        raise ValueError(f"Sample sheet {path} has no rows.")
    missing = [c for c in required if c not in rows[0]]
    if missing:
        raise ValueError(f"Sample sheet {path} lacks column(s): {', '.join(missing)}")
    sheet = {}
    seen_fastq = set()
    include_flag = {}
    for n, row in enumerate(rows, start=2):
        row = {k: (v or "").strip() for k, v in row.items() if k is not None}
        sp, sm = row["species_id"], row["sample"]
        if sp not in SPECIES_CFG:
            raise ValueError(f"Sample sheet line {n}: species '{sp}' is not in the config species table.")
        inc = (row.get("include", "") or "yes").lower()
        if inc not in ("yes", "no"):
            raise ValueError(f"Sample sheet line {n}: include must be yes, no or empty.")
        if include_flag.setdefault((sp, sm), inc) != inc:
            raise ValueError(f"Sample sheet line {n}: include differs between rows of {sp}/{sm}.")
        if inc == "no":
            continue
        for f in (row["r1"], row["r2"]):
            if f in seen_fastq:
                raise ValueError(f"Sample sheet line {n}: FASTQ listed twice: {f}")
            seen_fastq.add(f)
        if not re.fullmatch(r"[1-9][0-9]*", row["ploidy"]):
            raise ValueError(f"Sample sheet line {n}: ploidy must be a positive integer.")
        info = sheet.setdefault(sp, {}).setdefault(
            sm, {"ploidy": int(row["ploidy"]), "min_dp": row.get("min_dp", ""),
                 "max_dp": row.get("max_dp", ""), "units": []})
        for key in ("min_dp", "max_dp"):
            value = row.get(key, "")
            if not value:
                continue
            if not re.fullmatch(r"[0-9]+", value):
                raise ValueError(f"Sample sheet line {n}: {key} must be a whole number or empty.")
            if info[key] and info[key] != value:
                raise ValueError(f"Sample sheet line {n}: {key} differs between rows of {sp}/{sm}.")
            info[key] = value
        if int(row["ploidy"]) != info["ploidy"]:
            raise ValueError(f"Sample sheet line {n}: ploidy differs between rows of {sp}/{sm}.")
        info["units"].append((row["r1"], row["r2"]))
    return {sp: smp for sp, smp in sheet.items() if sp in ENABLED_SPECIES}


SHEET = load_sample_sheet(config["sample_sheet"])
for _sp in ENABLED_SPECIES:
    if _sp not in SHEET:
        raise ValueError(f"Species '{_sp}' is enabled but has no samples in the sample sheet.")


def _alternation(names):
    return "|".join(re.escape(n) for n in sorted(set(names), key=len, reverse=True))


wildcard_constraints:
    species = _alternation(ENABLED_SPECIES),
    sample = _alternation([sm for sp in SHEET.values() for sm in sp]),
    unit = r"\d+",
    IID = r"\d+",
    wsize = r"\d+"


## ---------------------------------------------------------------------
## Paths and lookups
## ---------------------------------------------------------------------
def P(*parts):
    """Path under the output root."""
    return os.path.join(OUTDIR, *parts)


def script(name):
    return os.path.join(workflow.basedir, "scripts", name)


def samples_of(species):
    return list(SHEET[species])


def sample_info(wildcards):
    return SHEET[wildcards.species][wildcards.sample]


def get_ploidy(wildcards):
    return sample_info(wildcards)["ploidy"]


def unit_ids(wildcards):
    return list(range(1, len(sample_info(wildcards)["units"]) + 1))


def raw_fastq(wildcards):
    r1, r2 = sample_info(wildcards)["units"][int(wildcards.unit) - 1]
    return {"r1": r1, "r2": r2}


def raw_r1_all(wildcards):
    return [u[0] for u in sample_info(wildcards)["units"]]


def trimmed_fastq_all(wildcards):
    return [P(wildcards.species, "trimmed_reads", f"{wildcards.sample}.u{u}.{r}.fastq.gz")
            for u in unit_ids(wildcards) for r in ("1P", "2P")]


def fq2bam_in_fq(wildcards):
    """--in-fq arguments: one per read pair, each with its own read group."""
    args = []
    for u in unit_ids(wildcards):
        r1 = P(wildcards.species, "trimmed_reads", f"{wildcards.sample}.u{u}.1P.fastq.gz")
        r2 = P(wildcards.species, "trimmed_reads", f"{wildcards.sample}.u{u}.2P.fastq.gz")
        rg = (f"@RG\\tID:{wildcards.sample}.u{u}\\tLB:{wildcards.sample}\\tPL:ILLUMINA"
              f"\\tSM:{wildcards.sample}\\tPU:{wildcards.sample}.u{u}")
        args.append(f'--in-fq {r1} {r2} "{rg}"')
    return " ".join(args)


def reference_of(wildcards):
    return SPECIES_CFG[wildcards.species]["reference"]


def regions_bed_of(species):
    return SPECIES_CFG[species].get("regions_bed", "") or ""


def gvcf_backend_path(wildcards):
    """Backend-specific GVCF path depending on ploidy (GPU for ploidy <= 2)."""
    backend = "gpu" if get_ploidy(wildcards) <= 2 else "cpu"
    return P(wildcards.species, "vcf", backend, f"{wildcards.sample}.g.vcf.gz")


def interval_ids():
    return range(1, config["split_n"] + 1)


def pbrun_tmp(wildcards, step):
    """--tmp-dir argument for pbrun (one subdirectory per job), or "" when pbrun_tmp_dir is unset."""
    root = config.get("pbrun_tmp_dir", "") or ""
    if not root:
        return ""
    return f"--tmp-dir {os.path.join(root, f'{wildcards.species}.{wildcards.sample}.{step}')}"


def get_total_memory_gb():
    import psutil
    return psutil.virtual_memory().total / (1024 ** 3)


## ---------------------------------------------------------------------
## Final targets
## ---------------------------------------------------------------------
def species_targets(sp):
    t = [
        P(sp, "vcf_final", f"{sp}.allSites.vcf.gz.tbi"),
        P(sp, "vcf_final", f"{sp}.HQSNPs.vcf.gz.tbi"),
        P(sp, "vcf_final", f"{sp}.callable.bed.gz"),
        P(sp, "qc", f"{sp}.sample_qc.tsv"),
        P(sp, "qc", f"{sp}.filter_stats.tsv"),
        P(sp, "qc", "multiqc", f"{sp}_multiqc_report.html"),
    ]
    t += [P(sp, "pixy", f"{sp}.w{w}_pi.txt") for w in config["pixy"]["window_sizes"]]
    if regions_bed_of(sp):
        t += [P(sp, "pixy", f"{sp}.regions_pi.txt"),
              P(sp, "qc", f"{sp}.callable_by_region.tsv")]
    return t


def all_targets():
    return [t for sp in ENABLED_SPECIES for t in species_targets(sp)]
