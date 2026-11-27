## Relatedness, PCA, the cross-species sample table with PASS/WARN flags, and one MultiQC
## report for all species of the run.

## KING kinship between all sample pairs of a species, and PCA on LD-pruned SNPs (plink2).
## --bad-freqs / --bad-ld: allele frequencies and LD come from the few samples of the species
## (plink2 otherwise refuses below 50 samples). PCA needs MAF >= 1/(4 sqrt(n)) (plink2's GRM
## floor), so pruning and PCA use that --maf; KING uses all SNPs.
rule plink_king_pca:
    conda:
        os.path.join(workflow.basedir, "envs/plink2.yaml")
    input:
        vcf = P("{species}", "vcf_final", "{species}.HQSNPs.vcf.gz"),
        tbi = P("{species}", "vcf_final", "{species}.HQSNPs.vcf.gz.tbi")
    output:
        kin = P("{species}", "qc", "plink", "{species}.kin0"),
        eigenvec = P("{species}", "qc", "plink", "{species}.eigenvec"),
        eigenval = P("{species}", "qc", "plink", "{species}.eigenval")
    params:
        prefix = lambda wc: P(wc.species, "qc", "plink", wc.species),
        n_pc = lambda wc: max(1, min(10, len(samples_of(wc.species)) - 1)),
        prune = config["qc_report"]["ld_prune"],
        maf = lambda wc: round(0.25 / len(samples_of(wc.species)) ** 0.5 + 0.001, 4)
    threads:
        4
    log:
        P("{species}", "logs", "plink", "king_pca.log")
    shell:
        """
        common="--vcf {input.vcf} --allow-extra-chr --set-all-var-ids @:# --snps-only --max-alleles 2 --threads {threads}"
        plink2 $common --make-king-table --out {params.prefix} > {log} 2>&1
        plink2 $common --bad-freqs --bad-ld --maf {params.maf} --indep-pairwise {params.prune} \
            --out {params.prefix}.prune >> {log} 2>&1
        plink2 $common --bad-freqs --maf {params.maf} --extract {params.prefix}.prune.prune.in --pca {params.n_pc} \
            --out {params.prefix} >> {log} 2>&1
        """

def report_inputs(wc=None):
    """Per-species QC files the summary and the MultiQC report read."""
    files = {}
    for sp in ENABLED_SPECIES:
        smp = samples_of(sp)
        files[sp] = {
            "sample_qc": P(sp, "qc", f"{sp}.sample_qc.tsv"),
            "filter_stats": P(sp, "qc", f"{sp}.filter_stats.tsv"),
            "kin": P(sp, "qc", "plink", f"{sp}.kin0"),
            "eigenvec": P(sp, "qc", "plink", f"{sp}.eigenvec"),
            "eigenval": P(sp, "qc", "plink", f"{sp}.eigenval"),
            "pi": P(sp, "pixy", f"{sp}.w{max(config['pixy']['window_sizes'])}_pi.txt"),
            "samtools": [P(sp, "qc", "samtools", f"{s}.stats") for s in smp],
            "flagstat": [P(sp, "qc", "samtools", f"{s}.flagstat") for s in smp],
            "mosdepth": [P(sp, "qc", "mosdepth", f"{s}.mosdepth.summary.txt") for s in smp],
            "dup": [P(sp, "qc", "dup_metrics", f"{s}.dup_metrics.txt") for s in smp],
            "fastp": [P(sp, "qc", "fastp", f"{s}.u{u}.fastp.json")
                      for s in smp for u in range(1, len(SHEET[sp][s]["units"]) + 1)],
            "fastqc": [P(sp, "qc", "fastqc", f"{s}.u{u}.{r}_fastqc.zip")
                       for s in smp for u in range(1, len(SHEET[sp][s]["units"]) + 1)
                       for r in ("1P", "2P")],
            "bcftools": P(sp, "qc", "bcftools_stats", f"{sp}.HQSNPs.vcf.stats"),
        }
    return files


rule qc_summary:
    input:
        lambda wc: [f for sp in report_inputs().values()
                    for key in ("sample_qc", "filter_stats", "kin", "eigenvec", "eigenval", "pi")
                    for f in [sp[key]]] +
                   [f for sp in report_inputs().values() for f in sp["samtools"]]
    output:
        samples = P("qc_summary", "sample_flags_mqc.tsv"),
        kinship = P("qc_summary", "kinship_mqc.tsv"),
        pca = P("qc_summary", "pca_mqc.tsv"),
        filters = P("qc_summary", "filter_stats_mqc.tsv"),
        pi = P("qc_summary", "pi_mqc.tsv"),
        pca_plots = [P("qc_summary", "pca_" + re.sub(r"[^A-Za-z0-9._-]", "_", sp) + "_mqc.json")
                     for sp in ENABLED_SPECIES],
        names = P("qc_summary", "multiqc_sample_names.yaml")
    params:
        script = script("qc_summary.py"),
        outdir = P("qc_summary"),
        species_args = lambda wc: " ".join(
            f"--species {sp} {f['sample_qc']} {f['filter_stats']} {f['kin']} {f['eigenvec']} "
            f"{f['eigenval']} {f['pi']} {os.path.dirname(f['samtools'][0])}"
            for sp, f in report_inputs().items()),
        flags = config["qc_report"]["flags"]
    log:
        P("logs", "qc_summary.log")
    shell:
        """
        python {params.script} \
            --outdir {params.outdir} \
            --mapping-rate-min {params.flags[mapping_rate_min]} \
            --error-rate-rel-max {params.flags[error_rate_rel_max]} \
            --het-rate-rel-low {params.flags[het_rate_rel_low]} \
            --het-rate-rel-high {params.flags[het_rate_rel_high]} \
            --ab-outer-excess {params.flags[ab_outer_excess]} \
            --kinship-duplicate {params.flags[kinship_duplicate]} \
            --kinship-report-min {params.flags[kinship_report_min]} \
            --pca-outlier-sd {params.flags[pca_outlier_sd]} \
            {params.species_args} > {log} 2>&1
        """

## One MultiQC report for all species: the standard modules (fastp, FastQC, samtools,
## mosdepth, duplicate metrics, bcftools) plus the summary tables as custom sections.
## Sample names become <species>_<sample> (regex renaming written by qc_summary.py).
rule multiqc_all:
    conda:
        os.path.join(workflow.basedir, "envs/qc.yaml")
    input:
        files = lambda wc: [f for sp in report_inputs().values()
                            for key in ("samtools", "flagstat", "mosdepth", "dup", "fastp", "fastqc")
                            for f in sp[key]] +
                           [sp["bcftools"] for sp in report_inputs().values()],
        summary = rules.qc_summary.output,
        names = P("qc_summary", "multiqc_sample_names.yaml")
    output:
        P("multiqc", "all_species_multiqc_report.html")
    params:
        dirs = lambda wc: " ".join([P(sp, "qc") for sp in ENABLED_SPECIES] + [P("qc_summary")]),
        config = os.path.join(workflow.basedir, "configuration", "multiqc_config.yaml"),
        outdir = P("multiqc")
    log:
        P("logs", "multiqc_all.log")
    shell:
        """
        multiqc --force --no-ai \
            -c {params.config} \
            -c {input.names} \
            -o {params.outdir} -n all_species_multiqc_report.html \
            {params.dirs} > {log} 2>&1
        """
