#!/usr/bin/env python
## All-site variant calling for several species in one run.
## Usage: snakemake --configfile <config.yaml> ... (see README)

if not config:
    raise ValueError("No configuration: run with --configfile <config.yaml> (see README).")

include: "rules/utils.smk"
include: "rules/0.qc.smk"
include: "rules/1.preprocessing.smk"
include: "rules/2.mapping.smk"
include: "rules/3.variant_calling.smk"
include: "rules/4.vcf_filtering.smk"
include: "rules/5.diversity.smk"

rule all:
    input:
        all_targets()
    shell:
        """
        echo "Job done!"
        """
