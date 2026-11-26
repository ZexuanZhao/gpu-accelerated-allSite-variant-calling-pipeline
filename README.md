# *ALLSITE* variant calling from Illumina reads using GPU

## Overview
A Snakemake workflow that calls variant and invariant sites from Illumina reads of several
species in one run, filters them with the same rules for reference and variant calls, and
computes nucleotide diversity (pi) with [`pixy`](https://pixy.readthedocs.io/en/latest/).
Mapping and HaplotypeCaller run on the GPU in NVIDIA's Clara Parabricks container (ploidy <= 2);
samples with ploidy > 2 use GATK4 HaplotypeCaller on the CPU. Joint genotyping, filtering and
pixy run per species.

## Steps
1. **Trimming** (`rules/1.preprocessing.smk`): `fastp` per read pair (unit). A sample may have
   several read pairs; `check_units` stops the run when two of them share more than
   `unit_check_max_shared` of their read IDs (the same reads delivered twice, e.g. one run
   demultiplexed twice), because mapping both would count every read twice.
2. **Mapping** (`rules/2.mapping.smk`): `pbrun fq2bam`, all read pairs of a sample in one call,
   one read group per pair (`LB` = sample), duplicates marked across pairs (not removed;
   HaplotypeCaller and mosdepth skip duplicate-flagged reads). Reference masks: low-complexity
   regions (`dustmasker`; its inclusive end is converted to a BED end), N gaps.
3. **Depth** (`rules/0.qc.smk`): `mosdepth` per sample outside N gaps and low-complexity
   regions (MAPQ >= `depth_mapq`, duplicates excluded). Each sample's genotype depth window is
   `min_dp = max(depth.min_abs, ceil(depth.min_frac * m))`, `max_dp = floor(depth.max_frac * m)`,
   m = median depth; `min_dp` / `max_dp` in the sample sheet override it. The absolute floor
   keeps a low-coverage sample from getting a floor of a few reads.
4. **Calling** (`rules/3.variant_calling.smk`): HaplotypeCaller GVCF per sample, GenomicsDB and
   `GenotypeGVCFs -all-sites` in `split_n` intervals per species
   (`--only-output-calls-starting-in-intervals`, so no record is emitted twice).
5. **Filtering** (`rules/4.vcf_filtering.smk`, `scripts/filter_allsites.py`), per interval:
   - drop low-complexity regions; flag SNPs within `indel_gap` bp of an indel;
   - keep only invariant sites and biallelic SNPs (indels, multiallelic sites, `*` dropped);
   - SNP site filters: QUAL, QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum, ExcessHet (a missing
     annotation passes); heterozygous fraction of called genotypes <= `het_fraction_max`;
   - genotypes set to missing unless: GQ >= `gq_min` (SNPs) or RGQ >= `rgq_min` (reference
     calls), depth inside the sample's window, and for heterozygous calls ALT reads >=
     `het_alt_min` and `het_ab_min` <= AB <= `het_ab_max`;
   - sites with missing fraction > `f_missing_max` dropped;
   - INFO AC/AF/AN/NS/F_MISSING recomputed (MLEAC/MLEAF removed).
   The same depth and quality rules apply to reference and variant calls, so the callable
   denominator and the variant numerator of pi are filtered alike. The raw interval VCFs are
   kept, so the filters can be changed and re-run without genotyping again.
6. **Diversity** (`rules/5.diversity.smk`): `pixy --stats pi` in each window size of
   `pixy.window_sizes`, and over the regions of `regions_bed` when a species has one.
7. **QC**: fastp, FastQC, bamtools, Qualimap, mosdepth, duplicate metrics, bcftools stats
   (MultiQC report per species); per-sample table with depth window, called and
   heterozygous sites (also per region label of `regions_bed`) and the allele balance of
   heterozygous calls; per-filter counts of dropped sites and masked genotypes.

## Files to prepare
- **Sample sheet** (CSV with a header), one row per read pair; several rows per sample when it
  has several read pairs:

  | column | content |
  |---|---|
  | `species_id` | a key of `species` in the config |
  | `sample` | sample name (read group SM) |
  | `ploidy` | 1, 2, ... (same in all rows of a sample) |
  | `r1`, `r2` | paths to the read files |
  | `min_dp`, `max_dp` | optional per-sample depth window; empty = computed |

- **Config** (copy `configuration/config.yaml`): `sample_sheet`, `outdir`, the `species`
  table (`reference`, `enabled`, optional `regions_bed` = BED with chrom, start, end, label),
  the Parabricks image and the filter thresholds. The Snakefile loads no config itself; pass
  yours with `--configfile`.

## Output (`<outdir>/<species_id>/`)
- `vcf_final/<species_id>.allSites.vcf.gz` (+ `.tbi`): filtered invariant sites and SNPs.
- `vcf_final/<species_id>.HQSNPs.vcf.gz` (+ `.tbi`): its SNPs with at least one ALT allele.
- `vcf_final/<species_id>.callable.bed.gz`: callable sites as merged intervals.
- `pixy/<species_id>.w<size>_pi.txt`, `pixy/<species_id>.regions_pi.txt`: pi.
- `qc/<species_id>.sample_qc.tsv`, `qc/<species_id>.ab_hist.tsv`,
  `qc/<species_id>.depth_thresholds.tsv`, `qc/<species_id>.filter_stats.tsv`,
  `qc/<species_id>.callable_by_region.tsv`, `qc/multiqc/<species_id>_multiqc_report.html`.

## Environments
Runner: the existing kremer env `snakemake_singularity` (`environment.yaml` records it:
snakemake 8.29.3, singularity 3.8.6). Rule envs, with pinned main tools: `envs/tools.yaml`
(bcftools/samtools/htslib 1.25, bedtools 2.31.1, bwa 0.7.19, fastp 1.4.0, mosdepth 0.3.14),
`envs/qc.yaml` (FastQC 0.13.0, Qualimap 2.3, bamtools 2.5.3, MultiQC 1.35), `envs/gatk4.yaml`
(GATK 4.6.2.0), `envs/blast.yaml` (BLAST+ 2.17.0, dustmasker), `envs/pixy.yaml` (pixy 2.2.4).
Parabricks: `clara-parabricks` in the config (4.6.0-1).

## Usage
Run from the directory that should hold `.snakemake/` (e.g. the project's run directory):

```bash
conda activate snakemake_singularity
snakemake -s <workflow>/Snakefile --configfile config.yaml \
    --use-conda --use-singularity \
    --singularity-args '--nv -B .:/dum --no-home -B $HOME:$HOME' \
    --cores [ncpu] --resources gpus=[ngpu] mem_gb=[mem] --rerun-incomplete --retries 3
```

Test data: `snakemake --configfile configuration/config.yaml ...` from the repository root.

## Notes
- `clara-parabricks` needs 38 GB GPU memory for `fq2bam`, so `--low-memory` is used.
- `[ncpu]` should be larger than 20 (the GPU rules use 20 threads); `[mem]` at most 80% of
  physical memory.
- If `haplotypecaller` throws an overflow error, check the BAMs for regions with extremely
  high depth and mask them.
- GPU steps fail at a low rate for unknown reasons, hence `--retries`.
- Reference calls inherit the depth of their GVCF block, so their DP is approximate.
- BQSR is not run (no known-sites resource for these species).
- `scripts/generate_report.R` reads the QC layout of the earlier version of this workflow
  (one project per run, bedtools window coverage) and has not been updated.
