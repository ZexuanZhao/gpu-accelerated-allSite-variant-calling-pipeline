def get_fastq(wildcards):
    fastqs = sample_sheet.loc[wildcards.sample, ["r1", "r2"]].dropna()
    return {"r1": fastqs.r1, "r2": fastqs.r2}

def get_fastq2(wildcards):
    return {"r1": os.path.join(config["outdir"], "trimmed_reads", wildcards.sample + ".1P.fastq.gz"),
            "r2": os.path.join(config["outdir"], "trimmed_reads", wildcards.sample + ".2P.fastq.gz")}

def get_fastq2_basename(wildcards):
    return {"r1_name": wildcards.sample + ".1P.fastq.gz",
            "r2_name": wildcards.sample + ".2P.fastq.gz"}

def get_ploidy(wildcards):
    return int(sample_sheet.loc[wildcards.sample, ["ploidy"]].dropna().ploidy)

def gvcf_backend_path(wildcards):
    """Backend-specific GVCF path depending on ploidy."""
    ploidy = get_ploidy(wildcards)
    backend = "gpu" if ploidy <= 2 else "cpu"
    return os.path.join(
        config["outdir"],
        "vcf",
        backend,
        f"{wildcards.sample}.g.vcf.gz"
    )

def format_input_vcfs(input_vcfs):
    formatted_vcfs = []
    for i, vcf_path in enumerate(input_vcfs):
        formatted_vcfs.append(f"I={vcf_path}")
    return " ".join(formatted_vcfs)

def get_total_memory_gb():
    import psutil
    mem_bytes = psutil.virtual_memory().total
    mem_gb = mem_bytes / (1024 ** 3)
    return mem_gb