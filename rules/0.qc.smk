## ---------------------------------------------------------------------
## Reads and alignments
## ---------------------------------------------------------------------
rule fastqc_after_trimming:
    conda:
        os.path.join(workflow.basedir, "envs/qc.yaml")
    input:
        r1 = P("{species}", "trimmed_reads", "{sample}.u{unit}.1P.fastq.gz"),
        r2 = P("{species}", "trimmed_reads", "{sample}.u{unit}.2P.fastq.gz")
    output:
        r1 = P("{species}", "qc", "fastqc", "{sample}.u{unit}.1P_fastqc.zip"),
        r2 = P("{species}", "qc", "fastqc", "{sample}.u{unit}.2P_fastqc.zip")
    params:
        outdir = lambda wc: P(wc.species, "qc", "fastqc")
    threads:
        2
    shell:
        """
        fastqc --quiet --outdir {params.outdir} -f fastq {input} -t {threads}
        """

## Alignment statistics (MultiQC reads both): mapped and properly paired reads, insert
## size, error rate (from NM), duplicates, MAPQ
rule samtools_stats:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        bam = P("{species}", "bam", "{sample}.bam"),
        bai = P("{species}", "bam", "{sample}.bam.bai"),
        ref = P("{species}", "ref", "ref.fasta")
    output:
        stats = P("{species}", "qc", "samtools", "{sample}.stats"),
        flagstat = P("{species}", "qc", "samtools", "{sample}.flagstat")
    threads:
        4
    shell:
        """
        samtools stats -@ {threads} --reference {input.ref} {input.bam} > {output.stats}
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat}
        """

## Depth with mosdepth (duplicate, secondary, QC-fail and unmapped reads excluded by its
## default flag 1796; mapping quality >= depth_mapq):
##   - over ref.depth_regions.bed (genome minus N gaps and low-complexity regions): the depth
##     distribution that gives each sample's median depth;
##   - in w_size windows, for coverage plots.
rule mosdepth:
    conda:
        os.path.join(workflow.basedir, "envs/mosdepth.yaml")
    input:
        bam = P("{species}", "bam", "{sample}.bam"),
        bai = P("{species}", "bam", "{sample}.bam.bai"),
        regions = P("{species}", "ref", "ref.depth_regions.bed")
    output:
        dist = P("{species}", "qc", "mosdepth", "{sample}.mosdepth.region.dist.txt"),
        summary = P("{species}", "qc", "mosdepth", "{sample}.mosdepth.summary.txt"),
        windows = P("{species}", "qc", "mosdepth", "{sample}.windows.regions.bed.gz")
    params:
        prefix = lambda wc: P(wc.species, "qc", "mosdepth", wc.sample),
        mapq = config["depth_mapq"],
        w_size = config["w_size"]
    threads:
        4
    log:
        P("{species}", "logs", "mosdepth", "{sample}.log")
    shell:
        """
        mosdepth -t {threads} -n -Q {params.mapq} -b {input.regions} {params.prefix} {input.bam} > {log} 2>&1
        mosdepth -t {threads} -n --fast-mode -Q {params.mapq} -b {params.w_size} {params.prefix}.windows {input.bam} >> {log} 2>&1
        """

## Per-sample genotype depth window (see config: depth)
rule depth_thresholds:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        lambda wc: expand(P(wc.species, "qc", "mosdepth", "{sample}.mosdepth.region.dist.txt"),
                          sample=samples_of(wc.species))
    output:
        P("{species}", "qc", "{species}.depth_thresholds.tsv")
    params:
        script = script("depth_thresholds.py"),
        d = config["depth"],
        sample_args = lambda wc, input: " ".join(
            f"--sample {sm} {dist} {SHEET[wc.species][sm]['min_dp'] or 'auto'} "
            f"{SHEET[wc.species][sm]['max_dp'] or 'auto'}"
            for sm, dist in zip(samples_of(wc.species), input))
    shell:
        """
        python {params.script} \
            --min-abs {params.d[min_abs]} \
            --min-frac {params.d[min_frac]} \
            --max-frac {params.d[max_frac]} \
            --low-median-flag {params.d[low_median_flag]} \
            {params.sample_args} \
            > {output}
        """

## ---------------------------------------------------------------------
## Calls
## ---------------------------------------------------------------------
rule vcf_stats:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        P("{species}", "vcf_final", "{species}.HQSNPs.vcf.gz")
    output:
        P("{species}", "qc", "bcftools_stats", "{species}.HQSNPs.vcf.stats")
    shell:
        """
        bcftools stats -s - {input} > {output}
        """

## Called and heterozygous genotypes per sample, over all callable sites
rule count_called_sites:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        vcf = P("{species}", "vcf_final", "{species}.allSites.vcf.gz"),
        tbi = P("{species}", "vcf_final", "{species}.allSites.vcf.gz.tbi")
    output:
        P("{species}", "qc", "{species}.count_sites.tsv")
    params:
        script = script("count_called_sites_bcftools.py")
    threads:
        4
    shell:
        """
        python {params.script} --threads {threads} {input.vcf} > {output}
        """

## The same per region label of regions_bed (e.g. TE-rich regions vs arms)
rule count_called_sites_by_region:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        vcf = P("{species}", "vcf_final", "{species}.allSites.vcf.gz"),
        tbi = P("{species}", "vcf_final", "{species}.allSites.vcf.gz.tbi"),
        regions = lambda wc: regions_bed_of(wc.species)
    output:
        P("{species}", "qc", "{species}.count_sites_by_region.tsv")
    params:
        script = script("count_called_sites_bcftools.py"),
        tmpdir = lambda wc: P(wc.species, "qc", "region_beds")
    threads:
        4
    shell:
        """
        set -euo pipefail
        mkdir -p {params.tmpdir}
        printf "region\\t" > {output}
        python {params.script} --header-only >> {output}
        cut -f4 {input.regions} | sort -u | while IFS= read -r label; do
            safe=$(printf '%s' "$label" | tr -c 'A-Za-z0-9._-' '_')
            awk -F'\\t' -v l="$label" '$4 == l' {input.regions} | cut -f1-3 > {params.tmpdir}/$safe.bed
            python {params.script} --threads {threads} --regions {params.tmpdir}/$safe.bed {input.vcf} \
                | tail -n +2 | awk -v l="$label" 'BEGIN{{OFS="\\t"}} {{print l, $0}}' >> {output}
        done
        """

## Callable bp per region label of regions_bed
rule callable_by_region:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        callable = P("{species}", "vcf_final", "{species}.callable.bed.gz"),
        regions = lambda wc: regions_bed_of(wc.species)
    output:
        P("{species}", "qc", "{species}.callable_by_region.tsv")
    shell:
        """
        set -o pipefail
        printf "region\\tregion_bp\\tcallable_bp\\tfrac_callable\\n" > {output}
        bedtools intersect -a <(cut -f1-4 {input.regions}) -b {input.callable} -wao \
            | awk -F'\\t' 'BEGIN{{OFS="\\t"}}
                {{k=$1":"$2":"$3":"$4; if(!(k in seen)){{seen[k]=1; len[$4]+=$3-$2}} cov[$4]+=$NF}}
                END{{for(l in len) print l, len[l], cov[l], (len[l]>0 ? cov[l]/len[l] : 0)}}' \
            | sort >> {output}
        """

## Per-sample QC table: depth window, called and heterozygous sites, allele balance of
## heterozygous calls (a mixed sample shows many hets with AB far from 0.5)
rule sample_qc:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        depth = P("{species}", "qc", "{species}.depth_thresholds.tsv"),
        counts = P("{species}", "qc", "{species}.count_sites.tsv"),
        hqsnps = P("{species}", "vcf_final", "{species}.HQSNPs.vcf.gz"),
        by_region = lambda wc: ([P(wc.species, "qc", f"{wc.species}.count_sites_by_region.tsv")]
                                if regions_bed_of(wc.species) else [])
    output:
        qc = P("{species}", "qc", "{species}.sample_qc.tsv"),
        ab = P("{species}", "qc", "{species}.ab_hist.tsv")
    params:
        script = script("sample_qc.py"),
        by_region = lambda wc, input: f"--counts-by-region {input.by_region[0]}" if input.by_region else ""
    shell:
        """
        python {params.script} \
            --depth {input.depth} \
            --counts {input.counts} \
            {params.by_region} \
            --hqsnps {input.hqsnps} \
            --ab-out {output.ab} \
            > {output.qc}
        """
