# Germline WGS Variant Calling Pipeline (GATK4, hg38)

## Overview
This repository contains an **end-to-end germline whole-genome sequencing (WGS) variant calling pipeline** implemented using **GATK4 Best Practices**.  
The pipeline processes paired-end FASTQ files and produces a **high-quality, analysis-ready VCF**, validated using standard QC metrics.

The workflow was tested on **HG00096 (1000 Genomes Project)** and demonstrates a complete WGS analysis from raw reads to filtered variants.

---

## Pipeline Workflow

### Input
- Paired-end FASTQ files (2 × 100 bp)
- Human reference genome: **hg38 / GRCh38**
- dbSNP known sites (for BQSR)

### Output
- PASS-filtered SNP and INDEL VCFs
- Merged germline variant callset
- Optional high-confidence genotype-filtered subset
- QC statistics (Ti/Tv ratio, variant counts)
- Tab-delimited variant table

---

## Analysis Steps

1. **Reference Preparation**
   - FASTA indexing (`samtools faidx`)
   - Sequence dictionary creation (`gatk CreateSequenceDictionary`)
   - BWA index generation

2. **Read Alignment**
   - Alignment using **BWA-MEM**
   - Sorting and indexing with SAMtools

3. **Post-alignment Processing**
   - Duplicate marking (`MarkDuplicatesSpark`)
   - Base Quality Score Recalibration (**BQSR**)

4. **Variant Calling**
   - **HaplotypeCaller** in gVCF mode
   - **GenotypeGVCFs** for variant genotyping

5. **Variant Filtering**
   - SNP and INDEL hard filtering (GATK-recommended thresholds)
   - PASS-only variant selection
   - SNP and INDEL merge

6. **Quality Control**
   - Variant counts
   - Transition/Transversion ratio (Ti/Tv) using `bcftools stats`

7. **Optional Downstream Processing**
   - High-confidence genotype filtering (DP ≥ 10, GQ ≥ 10)
   - VCF → TSV export for downstream analysis

---

## Final Outputs

| File | Description |
|-----|------------|
| `SRR062634.PASS.merged.vcf.gz` | Final PASS-filtered germline variant callset |
| `PASS.stats.txt` | QC metrics including Ti/Tv |
| `SRR062634.PASS.basic.tsv` | Tab-delimited variant summary |
| `SRR062634.PASS.GTDP10_GQ10.merged.vcf.gz` | High-confidence genotype-filtered subset |

---

## Key Results (HG00096)

- **PASS SNPs:** ~939,000  
- **PASS INDELs:** ~120,000  
- **Total PASS variants:** ~1.06 million  
- **Ti/Tv ratio:** **1.99**

A Ti/Tv ratio of ~2.0 is consistent with **high-quality human whole-genome sequencing data**, confirming correct variant discovery and filtering.

---

## Tools and Software

- **GATK:** v4.4.0.0  
- **BWA-MEM**
- **SAMtools**
- **BCFtools**
- **Reference genome:** hg38 (GRCh38)

---

## How to Run

```bash
chmod +x wgs_gatk_pipeline.sh
./wgs_gatk_pipeline.sh
