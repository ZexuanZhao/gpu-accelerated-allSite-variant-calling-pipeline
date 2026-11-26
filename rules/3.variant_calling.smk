# ---------------------------------------------------------------------
# GPU path: Clara Parabricks (ploidy <= 2)
# ---------------------------------------------------------------------
rule haplotypeCallerGPU:
    input:
        ref = P("{species}", "ref", "ref.fasta"),
        bam = P("{species}", "bam", "{sample}.bam"),
        bai = P("{species}", "bam", "{sample}.bam.bai")
    output:
        gvcf = P("{species}", "vcf", "gpu", "{sample}.g.vcf.gz")
    params:
        ploidy = get_ploidy
    log:
        P("{species}", "logs", "haplotypeCallerGPU", "{sample}.log")
    threads:
        20
    resources:
        gpus = 1
    container:
        config["clara-parabricks"]
    shell:
        """
        pbrun haplotypecaller \
            --ref {input.ref} \
            --in-bam {input.bam} \
            --out-variants {output.gvcf} \
            --ploidy {params.ploidy} \
            --gvcf \
            --num-htvc-threads {threads} \
            --htvc-low-memory \
            > {log} 2>&1
        """

# ---------------------------------------------------------------------
# CPU path: GATK4 HaplotypeCaller (ploidy > 2)
# ---------------------------------------------------------------------
rule haplotypeCallerCPU:
    input:
        ref = P("{species}", "ref", "ref.fasta"),
        dict = P("{species}", "ref", "ref.dict"),
        bam = P("{species}", "bam", "{sample}.bam"),
        bai = P("{species}", "bam", "{sample}.bam.bai")
    output:
        gvcf = P("{species}", "vcf", "cpu", "{sample}.g.vcf.gz")
    params:
        ploidy = get_ploidy
    log:
        P("{species}", "logs", "haplotypeCallerCPU", "{sample}.log")
    threads:
        4
    conda:
        os.path.join(workflow.basedir, "envs", "gatk4.yaml")
    shell:
        """
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -ERC GVCF \
            -ploidy {params.ploidy} \
            > {log} 2>&1
        """

# ---------------------------------------------------------------------
# Unify: one GVCF path per sample for joint genotyping
# ---------------------------------------------------------------------
rule gvcf_unify:
    input:
        gvcf_backend_path
    output:
        gvcf = P("{species}", "vcf", "{sample}.g.vcf.gz")
    shell:
        """
        target=$(readlink -f {input})
        ln -sfn "$target" {output.gvcf}
        ln -sfn "$target".tbi {output.gvcf}.tbi
        """

rule createIntervals:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        P("{species}", "ref", "ref.fasta")
    output:
        expand(P("{{species}}", "ref", "intervals", "ref.splitted.{IID}.bed"), IID=interval_ids())
    params:
        out_dir = lambda wc: P(wc.species, "ref", "intervals"),
        split_n = config["split_n"],
        script = script("split_fasta_into_n_bed.py")
    shell:
        """
        python {params.script} {input} {params.split_n} --out_dir {params.out_dir} --prefix ref
        """

rule sampleNameMap:
    input:
        gvcf = lambda wc: expand(P(wc.species, "vcf", "{sample}.g.vcf.gz"), sample=samples_of(wc.species))
    output:
        map = P("{species}", "vcf", "{species}.map.tsv")
    run:
        with open(output.map, "w") as fh:
            for sample, gvcf in zip(samples_of(wildcards.species), input.gvcf):
                fh.write(f"{sample}\t{gvcf}\n")

rule GenomicsDBImport:
    conda:
        os.path.join(workflow.basedir, "envs/gatk4.yaml")
    input:
        interval = P("{species}", "ref", "intervals", "ref.splitted.{IID}.bed"),
        map = P("{species}", "vcf", "{species}.map.tsv")
    output:
        directory(P("{species}", "gdb", "interval_{IID}.gdb"))
    params:
        batch_size = lambda wc: min(50, len(samples_of(wc.species)))
    threads:
        2
    resources:
        mem_gb = config["memory_gb_per_interval"]
    log:
        P("{species}", "logs", "GenomicsDBImport", "{IID}.log")
    shell:
        """
        gatk \
            --java-options "-Xmx{resources.mem_gb}g" \
            GenomicsDBImport \
                --genomicsdb-workspace-path {output} \
                --batch-size {params.batch_size} \
                -L {input.interval} \
                --sample-name-map {input.map} \
                --reader-threads {threads} \
            > {log} 2>&1
        """

rule CreateSequenceDictionary:
    conda:
        os.path.join(workflow.basedir, "envs/gatk4.yaml")
    input:
        P("{species}", "ref", "ref.fasta")
    output:
        P("{species}", "ref", "ref.dict")
    log:
        P("{species}", "logs", "CreateSequenceDictionary", "CreateSequenceDictionary.log")
    shell:
        """
        gatk CreateSequenceDictionary R={input} O={output} > {log} 2>&1
        """

## Raw all-sites VCF per interval. These are kept (not temporary), so the filters can be
## changed and re-run without genotyping again.
rule GenotypeGVCFs:
    conda:
        os.path.join(workflow.basedir, "envs/gatk4.yaml")
    input:
        ref = P("{species}", "ref", "ref.fasta"),
        dict = P("{species}", "ref", "ref.dict"),
        interval = P("{species}", "ref", "intervals", "ref.splitted.{IID}.bed"),
        gdb = P("{species}", "gdb", "interval_{IID}.gdb")
    output:
        vcf = P("{species}", "vcf", "raw", "interval_{IID}.vcf.gz"),
        tbi = P("{species}", "vcf", "raw", "interval_{IID}.vcf.gz.tbi")
    log:
        P("{species}", "logs", "GenotypeGVCFs", "{IID}.log")
    resources:
        mem_gb = config["memory_gb_per_interval"]
    threads:
        2
    shell:
        """
        gatk \
            --java-options "-Xmx{resources.mem_gb}g" \
            GenotypeGVCFs \
            -R {input.ref} \
            -V gendb://{input.gdb} \
            -all-sites \
            -L {input.interval} \
            --only-output-calls-starting-in-intervals \
            -O {output.vcf} \
            > {log} 2>&1
        """
