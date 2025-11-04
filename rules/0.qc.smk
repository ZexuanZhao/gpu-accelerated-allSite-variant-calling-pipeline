## Quality check of trimmed reads
## Summarized by multiQC
rule fastqc_after_trimming:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        unpack(get_fastq2)
    output:
        r1 = os.path.join(config["outdir"],"qc", "fastqc", "{sample}.1P_fastqc.zip"),
        r2 = os.path.join(config["outdir"],"qc", "fastqc", "{sample}.2P_fastqc.zip"),
    params:
        outdir = os.path.join(config["outdir"], "qc", "fastqc")
    threads:
        2
    shell:
        """
        fastqc --quiet --outdir {params.outdir} --noextract -f fastq {input} -t {threads}
        """

## Mapping rates etc.
## Summarized by multiQC
rule bamstats:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        os.path.join(config["outdir"], "bam", "{sample}.bam")
    output:
        os.path.join(config["outdir"], "qc", "bamtools", "{sample}_bamtools.stats")
    threads:
        1
    shell:
        """
        bamtools stats -in {input} | grep -v "*" > {output}
        """

## Make bed windows for window-based metrics: coverage
rule make_ref_window:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"ref","ref.fasta.fai")
    output:
        os.path.join(config["outdir"],"ref","ref.window.bed")
    params:
        w_size = config["w_size"]
    threads:
        1
    shell:
        """
        bedtools makewindows -g {input} -w {params.w_size} > {output}
        """

## Coverage
rule genomeCov:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        bam = os.path.join(config["outdir"], "bam", "{sample}.bam"),
        bed = os.path.join(config["outdir"],"ref","ref.window.bed")
    output:
        os.path.join(config["outdir"], "qc", "coverage", "{sample}_coverage.txt")
    threads:
        20
    shell:
        """
        bedtools coverage \
            -a {input.bed} -b {input.bam} \
            > {output}
        """

## Qualimap: insert size, GC, coverage...
## Summarized by multiQC
rule Qualimap:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        bam=os.path.join(config["outdir"],"bam","{sample}.bam")
    output:
        os.path.join(config["outdir"],"qc", "qualimap","{sample}", "qualimapReport.html")
    params:
        outdir =  os.path.join(config["outdir"],"qc","qualimap","{sample}")
    threads:
        10
    log:
        os.path.join(config["outdir"],"logs","Qualimap","{sample}.log")
    shell:
        """
        qualimap bamqc \
            -bam {input} \
            -outdir {params.outdir}\
            -outformat HTML \
            -nt {threads} \
            --java-mem-size=50G \
            > {log} \
            2>{log}
        """

## vcf stats using bcftools
## Summarized by multiQC
rule vcf_stats:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"vcf_final",config["project"] + ".removeLowQual" + ".lcm" + ".HQSNPs" + ".vcf.gz")
    output:
        os.path.join(config["outdir"],"qc", "bcftools_stats", config["project"]+".removeLowQual"+".lcm"+".HQSNPs"+".vcf.stats")
    threads:
        1
    shell:
        """
        bcftools stats {input} > {output}
        """

## Visualize vcf stats
## Summarized by multiQC
rule plot_vcfstats:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"qc", "bcftools_stats", config["project"]+".removeLowQual"+".lcm"+".HQSNPs"+".vcf.stats")
    output:
        os.path.join(config["outdir"],"qc","bcftools_stats", "plot-vcfstats.log")
    params:
        outdir = os.path.join(config["outdir"],"qc","bcftools_stats")
    threads:
        1
    shell:
        """
        plot-vcfstats \
            -p {params.outdir} \
            --no-PDF \
            {input}
        """

## Summarize all qc files using multiqc
rule multiqc:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        ## Fastp reports of reads
        expand(os.path.join(config["outdir"],"qc","fastp","{sample}.fastp.json"),sample=sample_sheet.index),
        ## Fastqc reads after trimming adaptors
        expand(os.path.join(config["outdir"],"qc","fastqc","{sample}.{R}_fastqc.zip"),sample=sample_sheet.index,R=["1P", "2P"]),
        ## Qualimap report of alignment
        expand(os.path.join(config["outdir"],"qc","qualimap","{sample}", "qualimapReport.html"), sample= sample_sheet.index),
        ## Bamtools report of alignment
        expand(os.path.join(config["outdir"], "qc", "bamtools","{sample}_bamtools.stats"), sample= sample_sheet.index),
        ## Bcftools stats of HQ SNPs
        os.path.join(config["outdir"],"qc", "bcftools_stats", config["project"]+".removeLowQual"+".lcm"+".HQSNPs"+".vcf.stats"),
    output:
        os.path.join(config["outdir"],"qc","multiqc", config["project"]+"_multiqc_report.html")
    params:
        input_dir = os.path.join(config["outdir"], "qc"),
        output_dir = os.path.join(config["outdir"], "qc", "multiqc"),
        original_output = os.path.join(config["outdir"],"qc", "multiqc", "multiqc_report.html")
    threads:
        1
    log:
        os.path.join(config["outdir"],"logs","multiqc","multiqc.log")
    shell:
        """
        rm -rf {params.output_dir}/* ; \
        multiqc \
        -o {params.output_dir} \
        {params.input_dir} \
        >{log} 2>{log}; \
        mv {params.original_output} {output}
        """

rule count_called_sites:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"vcf_final",config["project"] + ".removeLowQual" + ".lcm" + ".allSite" + ".vcf.gz")
    output:
        os.path.join(config["outdir"],"qc", config["project"] + "_count_sites.tsv")
    threads:
        10
    shell:
        """
        python scripts/count_called_sites_bcftools.py \
            --threads {threads} \
            --pass-only \
            {input} \
            > {output}
        """

rule pack_qc_reports:
    input:
        os.path.join(config["outdir"],"qc","multiqc",config["project"] + "_multiqc_report.html"),
        expand(os.path.join(config["outdir"],"qc","coverage","{sample}_coverage.txt"), sample = sample_sheet.index),
        count_sites=os.path.join(config["outdir"],"qc",config["project"] + "_count_sites.tsv")
    params:
        fastp_dir=os.path.join(config["outdir"], "qc", "fastp"),
        fastqc_dir=os.path.join(config["outdir"], "qc", "fastqc"),
        qualimap_dir=os.path.join(config["outdir"], "qc", "qualimap"),
        bamtools_dir=os.path.join(config["outdir"], "qc", "bamtools"),
        coverage_dir=os.path.join(config["outdir"], "qc", "coverage"),
        bcftools_stats_dir=os.path.join(config["outdir"],"qc", "bcftools_stats"),
        multiqc_dir=os.path.join(config["outdir"],"qc","multiqc"),
        zip_dir=os.path.join(config["outdir"],"qc", config["project"] + "_qc_files"),
        out_dir=config["outdir"],
        zip_dir_relative=os.path.join("qc", config["project"] + "_qc_files")
    output:
        os.path.join(config["outdir"],"qc", config["project"] + "_qc_files.zip")
    shell:
        """
            # fresh staging dir
            rm -rf "{params.zip_dir}"
            mkdir -p "{params.zip_dir}"
            
            cp -a {params.fastp_dir} {params.zip_dir}
            cp -a {params.fastqc_dir} {params.zip_dir}
            cp -a {params.qualimap_dir} {params.zip_dir}
            cp -a {params.bamtools_dir} {params.zip_dir}
            cp -a {params.coverage_dir} {params.zip_dir}
            cp -a {params.bcftools_stats_dir} {params.zip_dir}
            cp -a {params.multiqc_dir} {params.zip_dir}
            cp -a {input.count_sites} {params.zip_dir}

            cd {params.out_dir}
            zip -r {output} {params.zip_dir_relative}
        """
