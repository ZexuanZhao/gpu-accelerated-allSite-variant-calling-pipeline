## Filtering is done per interval (in parallel), then the intervals are joined in genome order.
## Per interval:
##   1. drop records in low-complexity regions (dustmasker);
##   2. flag SNPs within indel_gap bp of an indel (bcftools filter -g, FILTER SnpGap);
##   3. scripts/filter_allsites.py: keep invariant sites and biallelic SNPs that pass the site
##      and genotype filters (same depth/quality rules for reference and variant calls), drop
##      the stale INFO/MLEAC,MLEAF;
##   4. recompute AC, AF, AN, NS, F_MISSING.

def snpgap_cmd():
    gap = int(config["site_filters"]["indel_gap"])
    if gap <= 0:
        return "cat"
    return f"bcftools filter --no-version -g {gap} -s SnpGap -m+ -Ov -"


rule filterInterval:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        vcf = P("{species}", "vcf", "raw", "interval_{IID}.vcf.gz"),
        lcm = P("{species}", "ref", "ref.lcm.bed"),
        depth = P("{species}", "qc", "{species}.depth_thresholds.tsv")
    output:
        vcf = P("{species}", "vcf_filtered", "interval_{IID}.vcf.gz"),
        stats = P("{species}", "vcf_filtered", "interval_{IID}.filter_stats.tsv")
    params:
        script = script("filter_allsites.py"),
        snpgap = snpgap_cmd(),
        sf = config["site_filters"],
        gf = config["genotype_filters"],
        f_missing_max = config["f_missing_max"]
    log:
        P("{species}", "logs", "filterInterval", "{IID}.log")
    threads:
        2
    shell:
        """
        set -o pipefail
        bcftools view --no-version -T ^{input.lcm} -Ov {input.vcf} 2> {log} \
            | {params.snpgap} 2>> {log} \
            | python {params.script} \
                --depth-thresholds {input.depth} \
                --stats {output.stats} \
                --qual-min {params.sf[qual_min]} \
                --qd-min {params.sf[qd_min]} \
                --fs-max {params.sf[fs_max]} \
                --sor-max {params.sf[sor_max]} \
                --mq-min {params.sf[mq_min]} \
                --mqranksum-min {params.sf[mqranksum_min]} \
                --readposranksum-min {params.sf[readposranksum_min]} \
                --excesshet-max {params.sf[excesshet_max]} \
                --het-fraction-max {params.sf[het_fraction_max]} \
                --gq-min {params.gf[gq_min]} \
                --rgq-min {params.gf[rgq_min]} \
                --het-alt-min {params.gf[het_alt_min]} \
                --het-ab-min {params.gf[het_ab_min]} \
                --het-ab-max {params.gf[het_ab_max]} \
                --f-missing-max {params.f_missing_max} \
                2>> {log} \
            | bcftools +fill-tags --no-version - -Oz -o {output.vcf} -- -t AC,AF,AN,NS,F_MISSING 2>> {log}
        """

rule concatIntervals:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        expand(P("{{species}}", "vcf_filtered", "interval_{IID}.vcf.gz"), IID=interval_ids())
    output:
        vcf = P("{species}", "vcf_final", "{species}.allSites.vcf.gz"),
        tbi = P("{species}", "vcf_final", "{species}.allSites.vcf.gz.tbi")
    log:
        P("{species}", "logs", "concatIntervals", "concat.log")
    threads:
        4
    shell:
        """
        bcftools concat --no-version --naive -Oz -o {output.vcf} {input} > {log} 2>&1
        bcftools index -t --threads {threads} {output.vcf} >> {log} 2>&1
        """

rule HQSNPs:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        vcf = P("{species}", "vcf_final", "{species}.allSites.vcf.gz"),
        tbi = P("{species}", "vcf_final", "{species}.allSites.vcf.gz.tbi")
    output:
        vcf = P("{species}", "vcf_final", "{species}.HQSNPs.vcf.gz"),
        tbi = P("{species}", "vcf_final", "{species}.HQSNPs.vcf.gz.tbi")
    threads:
        4
    shell:
        """
        bcftools view --threads {threads} -v snps -c 1 -Oz -o {output.vcf} {input.vcf}
        bcftools index -t {output.vcf}
        """

## Callable sites (all kept records, invariant or SNP), merged into intervals
rule callableBed:
    conda:
        os.path.join(workflow.basedir, "envs/tools.yaml")
    input:
        P("{species}", "vcf_final", "{species}.allSites.vcf.gz")
    output:
        P("{species}", "vcf_final", "{species}.callable.bed.gz")
    shell:
        """
        set -o pipefail
        bcftools query -f '%CHROM\\t%POS0\\t%END\\n' {input} \
            | bedtools merge -i - \
            | bgzip > {output}
        """

rule filterStats:
    input:
        expand(P("{{species}}", "vcf_filtered", "interval_{IID}.filter_stats.tsv"), IID=interval_ids())
    output:
        P("{species}", "qc", "{species}.filter_stats.tsv")
    run:
        from collections import Counter
        total = Counter()
        for path in input:
            with open(path) as fh:
                next(fh)
                for line in fh:
                    level, sample, kind, outcome, count = line.rstrip("\n").split("\t")
                    total[(level, sample, kind, outcome)] += int(count)
        with open(output[0], "w") as out:
            out.write("level\tsample\tkind\toutcome\tcount\n")
            for key in sorted(total):
                out.write("\t".join(key) + f"\t{total[key]}\n")
