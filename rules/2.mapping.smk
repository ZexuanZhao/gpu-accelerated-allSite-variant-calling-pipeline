rule copy_reference:
    input:
        original_ref = config["reference"]
    output:
        copied_ref = os.path.join(config["outdir"], "ref", "ref.fasta")
    shell:
        """
        cp {input.original_ref} {output.copied_ref}
        """

rule index_reference:
    conda:
        os.path.join(workflow.basedir, "envs/envs.yaml")
    input:
        ref = os.path.join(config["outdir"], "ref", "ref.fasta")
    output:
        ref_bwt = os.path.join(config["outdir"], "ref", "ref.fasta.bwt"),
        ref_fai = os.path.join(config["outdir"], "ref", "ref.fasta.fai")
    log:
        os.path.join(config["outdir"],"logs","index_ref","index_ref.log")
    shell:
        """
        bwa index {input.ref} >{log} 2>{log}
        samtools faidx {input.ref} >{log} 2>{log}
        """

rule bwa_gpu:
    input:
        unpack(get_fastq2),
        copied_ref= os.path.join(config["outdir"],"ref","ref.fasta"),
        ref_bwt = os.path.join(config["outdir"], "ref", "ref.fasta.bwt"),
    output:
        os.path.join(config["outdir"], "bam", "{sample}.bam")
    params:
        tag = "@RG\\tID:" + "{sample}" + "\\tLB:lib1\\tPL:Illumina" + "\\tSM:" + "{sample}" + "\\tPU:unit1"
    log:
        os.path.join(config["outdir"],"logs","mapping","{sample}.log")
    threads:
        20
    resources:
        gpus=1
    singularity:
        config["clara-parabricks"]
    retries: 3
    shell:
        """
            pbrun fq2bam \
                --low-memory \
                --ref {input.copied_ref} \
                --in-fq {input.r1} {input.r2} "{params.tag}" \
                --out-bam {output} \
            > {log} \
            2>{log}
        """
