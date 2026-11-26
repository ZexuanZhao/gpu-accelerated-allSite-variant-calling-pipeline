rule copy_reference:
    input:
        reference_of
    output:
        P("{species}", "ref", "ref.fasta")
    shell:
        """
        cp {input} {output}
        """

rule index_reference:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        P("{species}", "ref", "ref.fasta")
    output:
        bwt = P("{species}", "ref", "ref.fasta.bwt"),
        fai = P("{species}", "ref", "ref.fasta.fai")
    log:
        P("{species}", "logs", "index_ref", "index_ref.log")
    shell:
        """
        bwa index {input} > {log} 2>&1
        samtools faidx {input} >> {log} 2>&1
        """

## All read pairs of a sample in one fq2bam call, one read group per pair (LB = sample), so
## duplicates are marked across pairs. Duplicates are marked, not removed: HaplotypeCaller and
## mosdepth skip duplicate-flagged reads.
rule bwa_gpu:
    input:
        fq = trimmed_fastq_all,
        unit_check = P("{species}", "qc", "unit_check", "{sample}.tsv"),
        ref = P("{species}", "ref", "ref.fasta"),
        bwt = P("{species}", "ref", "ref.fasta.bwt")
    output:
        bam = P("{species}", "bam", "{sample}.bam"),
        bai = P("{species}", "bam", "{sample}.bam.bai"),
        dup_metrics = P("{species}", "qc", "dup_metrics", "{sample}.dup_metrics.txt")
    params:
        in_fq = fq2bam_in_fq
    log:
        P("{species}", "logs", "mapping", "{sample}.log")
    threads:
        20
    resources:
        gpus = 1
    container:
        config["clara-parabricks"]
    retries: 3
    shell:
        """
        pbrun fq2bam \
            --low-memory \
            --ref {input.ref} \
            {params.in_fq} \
            --out-bam {output.bam} \
            --out-duplicate-metrics {output.dup_metrics} \
            > {log} 2>&1
        """

## ---------------------------------------------------------------------
## Reference masks: low-complexity regions (dustmasker), N gaps, and the regions where depth
## is measured (the rest of the genome)
## ---------------------------------------------------------------------
rule lowComplexityAcclist:
    conda:
        os.path.join(workflow.basedir, "envs/blast.yaml")
    input:
        P("{species}", "ref", "ref.fasta")
    output:
        P("{species}", "ref", "ref.lcm.acclist")
    shell:
        """
        dustmasker -in {input} -outfmt acclist > {output}
        """

rule refMasks:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        fasta = P("{species}", "ref", "ref.fasta"),
        fai = P("{species}", "ref", "ref.fasta.fai"),
        acclist = P("{species}", "ref", "ref.lcm.acclist")
    output:
        lcm = P("{species}", "ref", "ref.lcm.bed"),
        gaps = P("{species}", "ref", "ref.N.bed"),
        genome = P("{species}", "ref", "ref.genome"),
        depth_regions = P("{species}", "ref", "ref.depth_regions.bed")
    params:
        script = script("make_ref_masks.py")
    shell:
        """
        cut -f1,2 {input.fai} > {output.genome}
        python {params.script} \
            --fasta {input.fasta} \
            --acclist {input.acclist} \
            --lcm-out {output.lcm} \
            --gaps-out {output.gaps}
        cat {output.lcm} {output.gaps} \
            | bedtools sort -g {output.genome} \
            | bedtools merge \
            | bedtools complement -i - -g {output.genome} \
            > {output.depth_regions}
        """
