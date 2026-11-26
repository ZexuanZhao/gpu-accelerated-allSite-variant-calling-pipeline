## Trim each read pair (unit) of a sample separately; units are mapped together in fq2bam.
rule trim_reads:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        unpack(raw_fastq)
    output:
        r1 = P("{species}", "trimmed_reads", "{sample}.u{unit}.1P.fastq.gz"),
        r2 = P("{species}", "trimmed_reads", "{sample}.u{unit}.2P.fastq.gz"),
        json_report = P("{species}", "qc", "fastp", "{sample}.u{unit}.fastp.json"),
        html_report = P("{species}", "qc", "fastp", "{sample}.u{unit}.fastp.html")
    log:
        P("{species}", "logs", "fastp", "{sample}.u{unit}.log")
    threads:
        4
    shell:
        """
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --detect_adapter_for_pe \
            --thread {threads} \
            --json {output.json_report} \
            --html {output.html_report} \
            > {log} 2>&1
        """

## Stop when two read pairs of one sample are the same reads delivered twice
## (e.g. one run demultiplexed twice); mapping both would count every read twice.
rule check_units:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        raw_r1_all
    output:
        P("{species}", "qc", "unit_check", "{sample}.tsv")
    params:
        n_reads = config["unit_check_reads"],
        max_shared = config["unit_check_max_shared"],
        script = script("check_unit_overlap.py")
    threads:
        1
    shell:
        """
        python {params.script} \
            --n-reads {params.n_reads} \
            --max-shared {params.max_shared} \
            --out {output} \
            {input}
        """
