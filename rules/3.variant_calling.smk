# ---------------------------------------------------------------------
# GPU path: Clara Parabricks (ploidy <= 2)
# ---------------------------------------------------------------------
rule haplotypeCallerGPU:
    input:
        ref = os.path.join(config["outdir"], "ref", "ref.fasta"),
        bam = os.path.join(config["outdir"], "bam", "{sample}.bam")
    output:
        gvcf = os.path.join(config["outdir"], "vcf", "gpu", "{sample}.g.vcf.gz")
    params:
        ploidy = get_ploidy
    log:
        os.path.join(config["outdir"], "logs", "haplotypeCallerGPU", "{sample}.log")
    threads:
        20
    resources:
        gpus = 1
    singularity:
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
        ref = os.path.join(config["outdir"], "ref", "ref.fasta"),
        bam = os.path.join(config["outdir"], "bam", "{sample}.bam")
    output:
        gvcf = os.path.join(config["outdir"], "vcf", "cpu", "{sample}.g.vcf.gz")
    params:
        ploidy = get_ploidy
    log:
        os.path.join(config["outdir"], "logs", "haplotypeCallerCPU", "{sample}.log")
    threads:
        4
    resources:
        gpus = 0
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
# Unify: expose one GVCF path per sample for downstream joint calling
# ---------------------------------------------------------------------
rule gvcf_unify:
    input:
        gvcf_backend_path
    output:
        gvcf = os.path.join(config["outdir"], "vcf", "{sample}.g.vcf.gz")
    shell:
        """
        target=$(readlink -f {input})
        ln -sfn "$target" {output.gvcf}
        ln -sfn "$target".tbi {output.gvcf}.tbi
        """

rule createIntervals:
    conda:
        os.path.join(workflow.basedir,"envs/envs.yaml")
    input:
        os.path.join(config["outdir"],"ref","ref.fasta")
    output:
        expand(os.path.join(config["outdir"],"ref","ref.splitted.{IID}.bed"), IID = range(1, config["split_n"]+1))
    params:
        out_dir = os.path.join(config["outdir"], "ref"),
        split_n = config["split_n"]
    threads:
        1
    shell:
        """
        python scripts/split_fasta_into_n_bed.py {input} {params.split_n} --out_dir {params.out_dir} --prefix ref
        """

rule sampleNameMap:
    input:
        gvcf = expand(os.path.join(config["outdir"],"vcf","{sample}.g.vcf.gz"), sample = sample_sheet.index)
    output:
        map = os.path.join(config["outdir"],"vcf", config["project"] + ".map.tsv")
    params:
        outdir = config["outdir"],
        sample = sample_sheet.index
    threads:
        1
    run:
        data = []
        for sample in params.sample:
            gvcf_file = os.path.join(config["outdir"],"vcf",f"{sample}.g.vcf.gz")
            data.append({"sample": sample, "sample_vcf": gvcf_file})
        df = pd.DataFrame(data)
        df.to_csv(output.map,sep="\t", index=False, header=False)

rule GenomicsDBImport:
    conda:
        os.path.join(workflow.basedir,"envs/gatk4.yaml")
    input:
        interval = os.path.join(config["outdir"],"ref","ref.splitted.{IID}.bed"),
        map= os.path.join(config["outdir"],"vcf",config["project"] + ".map.tsv")
    output:
        directory(os.path.join(config["outdir"],"gdb",config["project"]+".gdb.{IID}"))
    params:
        IID = "{IID}",
        batch_size = min(50, len(sample_sheet))
    threads:
        10
    resources:
        mem_gb = config["memory_gb_per_interval"]
    log:
        os.path.join(config["outdir"],"logs","GenomicsDBImport","{IID}.log")
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
            > {log} \
            2>{log}
        """
rule CreateSequenceDictionary:
    conda:
        os.path.join(workflow.basedir,"envs/gatk4.yaml")
    input:
        os.path.join(config["outdir"],"ref","ref.fasta")
    output:
        os.path.join(config["outdir"],"ref","ref.dict")
    threads:
        1
    log:
        os.path.join(config["outdir"],"logs","CreateSequenceDictionary","CreateSequenceDictionary.log")
    shell:
        """
        gatk CreateSequenceDictionary \
            R={input} \
            O={output} \
             >{log} 2>{log}
        """
rule GenotypeGVCFs:
    conda:
        os.path.join(workflow.basedir,"envs/gatk4.yaml")
    input:
        ref = os.path.join(config["outdir"],"ref","ref.fasta"),
        dict =  os.path.join(config["outdir"],"ref","ref.dict"),
        interval = os.path.join(config["outdir"],"ref","ref.splitted.{IID}.bed"),
        gdb = os.path.join(config["outdir"],"gdb",config["project"] + ".gdb.{IID}")
    output:
        os.path.join(config["outdir"],"vcf",config["project"]+".{IID}"+".vcf.gz")
    params:
        IID = "{IID}"
    log:
        os.path.join(config["outdir"],"logs","GenotypeGVCFs","{IID}.log")
    resources:
        mem_gb=config["memory_gb_per_interval"]
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
            -O {output} \
            > {log} \
            2>{log}
        """

rule mergeVcfs:
    conda:
        os.path.join(workflow.basedir,"envs/gatk4.yaml")
    input:
        expand(os.path.join(config["outdir"],"vcf",config["project"]+".{IID}"+".vcf.gz"), IID = range(1, config["split_n"]+1))
    output:
        os.path.join(config["outdir"],"vcf",config["project"]+".raw"+".vcf.gz")
    params:
        input_formatted = format_input_vcfs(expand(os.path.join(config["outdir"],"vcf",config["project"]+".{IID}"+".vcf.gz"), IID = range(1, config["split_n"]+1)))
    threads:
        10
    log:
        os.path.join(config["outdir"],"logs","mergeVcfs","mergeVcfs.log")
    resources:
        mem_gb = int(round(get_total_memory_gb() * 0.8, 0))
    shell:
        """
        gatk \
            --java-options "-Xmx{resources.mem_gb}g" \
            MergeVcfs \
            {params.input_formatted} \
            O={output} \
            > {log} \
            2>{log}
        """
