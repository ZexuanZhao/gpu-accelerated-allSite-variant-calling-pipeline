#!/usr/bin/env python
import os.path
import pandas as pd

configfile: "configuration/config.yaml"
sample_sheet = pd.read_csv(config["sample_sheet"],
    dtype=str,
    names = ["sample", "r1", "r2"]).set_index("sample")

wildcard_constraints:
    sample = "|".join(sample_sheet.index)

include: "rules/utils.smk"
include: "rules/0.qc.smk"
include: "rules/1.preprocessing.smk"
include: "rules/2.mapping.smk"
include: "rules/3.variant_calling.smk"
include: "rules/4.vcf_filtering.smk"

rule all:
    input:
        allSite_vcf = os.path.join(config["outdir"], "vcf_final", config["project"]+".removeLowQual"+".lcm"+".allSite"+".vcf.gz"),
        SNPonly_vcf = os.path.join(config["outdir"],"vcf_final",config["project"] + ".removeLowQual" + ".lcm" + ".HQSNPs" + ".vcf.gz"),
        qc_zip      = os.path.join(config["outdir"],"qc", config["project"] + "_qc_files.zip")
    shell:
        """
        echo "Job done!"
        """
