## Nucleotide diversity (pi) with pixy on the filtered all-sites VCF. All samples of a species
## form one population. pixy needs the invariant sites, which allSites.vcf.gz keeps.

rule pixy_populations:
    output:
        P("{species}", "pixy", "{species}.populations.txt")
    run:
        with open(output[0], "w") as fh:
            for sample in samples_of(wildcards.species):
                fh.write(f"{sample}\t{wildcards.species}\n")

rule pixy_windows:
    conda:
        os.path.join(workflow.basedir, "envs/pixy.yaml")
    input:
        vcf = P("{species}", "vcf_final", "{species}.allSites.vcf.gz"),
        tbi = P("{species}", "vcf_final", "{species}.allSites.vcf.gz.tbi"),
        pops = P("{species}", "pixy", "{species}.populations.txt")
    output:
        P("{species}", "pixy", "{species}.w{wsize}_pi.txt")
    params:
        outdir = lambda wc: P(wc.species, "pixy"),
        prefix = lambda wc: f"{wc.species}.w{wc.wsize}"
    threads:
        config["pixy"]["threads"]
    log:
        P("{species}", "logs", "pixy", "w{wsize}.log")
    shell:
        """
        pixy --stats pi \
            --vcf {input.vcf} \
            --populations {input.pops} \
            --window_size {wildcards.wsize} \
            --n_cores {threads} \
            --output_folder {params.outdir} \
            --output_prefix {params.prefix} \
            > {log} 2>&1
        """

## pi per region of regions_bed (pixy reads the first three columns, 0-based BED)
rule pixy_regions:
    conda:
        os.path.join(workflow.basedir, "envs/pixy.yaml")
    input:
        vcf = P("{species}", "vcf_final", "{species}.allSites.vcf.gz"),
        tbi = P("{species}", "vcf_final", "{species}.allSites.vcf.gz.tbi"),
        pops = P("{species}", "pixy", "{species}.populations.txt"),
        regions = lambda wc: regions_bed_of(wc.species)
    output:
        bed = P("{species}", "pixy", "{species}.regions.bed"),
        pi = P("{species}", "pixy", "{species}.regions_pi.txt")
    params:
        outdir = lambda wc: P(wc.species, "pixy"),
        prefix = lambda wc: f"{wc.species}.regions"
    threads:
        config["pixy"]["threads"]
    log:
        P("{species}", "logs", "pixy", "regions.log")
    shell:
        """
        cut -f1-3 {input.regions} > {output.bed}
        pixy --stats pi \
            --vcf {input.vcf} \
            --populations {input.pops} \
            --bed_file {output.bed} \
            --n_cores {threads} \
            --output_folder {params.outdir} \
            --output_prefix {params.prefix} \
            > {log} 2>&1
        """
