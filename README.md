Pipeline Summary

Input

Paired-end FASTQ files (2 × 100 bp)

Human reference genome: hg38

dbSNP known sites (for BQSR)

Output

PASS-filtered germline SNP and INDEL VCFs

Optional high-confidence genotype-filtered subset

QC statistics (Ti/Tv, counts)

Tab-delimited variant table for downstream analysis

Workflow Steps

Reference Preparation

FASTA indexing (samtools faidx)

Sequence dictionary creation (GATK)

BWA index generation

Read Alignment

Alignment using BWA-MEM

Sorting and indexing (SAMtools)

Post-alignment Processing

Duplicate marking (MarkDuplicatesSpark)

Base Quality Score Recalibration (BQSR)

Variant Calling

HaplotypeCaller in gVCF mode

GenotypeGVCFs for joint genotyping (single sample)

Variant Filtering

SNP and INDEL hard filtering (GATK-recommended thresholds)

PASS-only variant selection

SNP + INDEL merge

Quality Control

Variant counts

Ti/Tv ratio using bcftools stats

Optional Analysis

High-confidence genotype subset (DP ≥ 10, GQ ≥ 10)

VCF → TSV export for R/Python analysis

Final Outputs
File	Description
SRR062634.PASS.merged.vcf.gz	Final PASS-filtered germline variants
SRR062634.PASS.stats.txt	QC statistics (Ti/Tv, counts)
SRR062634.PASS.basic.tsv	Tabular variant summary
SRR062634.PASS.GTDP10_GQ10.merged.vcf.gz	High-confidence genotype-filtered subset (optional)
Key Results (HG00096)

PASS SNPs: ~939,000

PASS INDELs: ~120,000

Total PASS variants: ~1.06 million

Ti/Tv ratio: 1.99

A Ti/Tv ratio of ~2.0 is consistent with high-quality human WGS data, confirming correct variant discovery and filtering.

Tools & Versions

GATK: 4.4.0.0

BWA: MEM

SAMtools

BCFtools

Reference genome: hg38 (GRCh38)
