#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# End-to-end Germline WGS Variant Calling (GATK4, hg38)
# What this script does (single sample):
#   0) Reference prep (fai, dict, bwa index) [if missing]
#   1) (Optional) FASTQC
#   2) BWA-MEM align -> sort -> index
#   3) MarkDuplicatesSpark -> index + metrics
#   4) BQSR (BaseRecalibrator + ApplyBQSR) -> index
#   5) HaplotypeCaller (-ERC GVCF) -> gVCF
#   6) GenotypeGVCFs -> raw VCF
#   7) Select SNPs/INDELs
#   8) Hard filter SNPs/INDELs + keep PASS only
#   9) Merge PASS SNPs+INDELs
#  10) QC: counts + bcftools stats + Ti/Tv
#  11) Optional: stricter genotype subset DP>=10 and GQ>=10
#  12) Optional: export a basic TSV table
#
# Notes:
# - Designed for Linux/HPC with modules: gatk, bwa, samtools, bcftools (fastqc optional)
# - Works with a local hg38 FASTA (hg38.fa) in your project folder by default.
# - Assumes paired FASTQs exist (SRR062634_1.filt.fastq.gz / SRR062634_2.filt.fastq.gz).
# ============================================================

# ----------------------------
# User-configurable parameters
# ----------------------------
SAMPLE="${SAMPLE:-SRR062634}"
THREADS="${THREADS:-4}"

# Project root (default: current directory)
ROOT="${ROOT:-$PWD}"

# Inputs
R1="${R1:-$ROOT/reads/${SAMPLE}_1.filt.fastq.gz}"
R2="${R2:-$ROOT/reads/${SAMPLE}_2.filt.fastq.gz}"

# Reference (hg38)
REF="${REF:-$ROOT/hg38.fa}"

# Known sites for BQSR (dbSNP; hg38/GRCh38)
KNOWN_SITES="${KNOWN_SITES:-$ROOT/Homo_sapiens_assembly38.dbsnp138.vcf}"

# Output dirs
READS_DIR="${READS_DIR:-$ROOT/reads}"
ALN_DIR="${ALN_DIR:-$ROOT/aligned_reads}"
DATA_DIR="${DATA_DIR:-$ROOT/data}"
RESULTS_DIR="${RESULTS_DIR:-$ROOT/results}"
VARIANTS_DIR="${VARIANTS_DIR:-$ROOT/variants}"

# Optional toggles
RUN_FASTQC="${RUN_FASTQC:-false}"     # set true if fastqc module exists
MAKE_TSV="${MAKE_TSV:-true}"          # export a basic TSV table
MAKE_GT_SUBSET="${MAKE_GT_SUBSET:-true}"  # DP/GQ subset (bcftools)

# Genotype thresholds for subset
MIN_DP="${MIN_DP:-10}"
MIN_GQ="${MIN_GQ:-10}"

# ----------------------------
# Module loads (adjust as needed)
# ----------------------------
module -q load gatk || true
module -q load bwa || true
module -q load samtools || true
module -q load bcftools || true
if [[ "$RUN_FASTQC" == "true" ]]; then
  module -q load fastqc || true
fi

# ----------------------------
# Create dirs
# ----------------------------
mkdir -p "$READS_DIR" "$ALN_DIR" "$DATA_DIR" "$RESULTS_DIR" "$VARIANTS_DIR"

# ----------------------------
# Sanity checks
# ----------------------------
for f in "$R1" "$R2" "$REF" "$KNOWN_SITES"; do
  [[ -e "$f" ]] || { echo "ERROR: missing file: $f" >&2; exit 1; }
done

# ----------------------------
# Step 0: Reference prep (only if missing)
# ----------------------------
echo "== Step 0: reference prep =="

# samtools faidx -> .fai
if [[ ! -e "${REF}.fai" ]]; then
  echo "Creating FASTA index (.fai)..."
  samtools faidx "$REF"
fi

# gatk CreateSequenceDictionary -> .dict
DICT="${REF%.fa}.dict"
if [[ ! -e "$DICT" ]]; then
  echo "Creating sequence dictionary (.dict)..."
  gatk CreateSequenceDictionary -R "$REF" -O "$DICT"
fi

# bwa index -> .amb/.ann/.bwt/.pac/.sa
if [[ ! -e "${REF}.bwt" && ! -e "${REF}.0123" ]]; then
  echo "Creating BWA index..."
  bwa index "$REF"
fi

# ----------------------------
# Step 1: QC (optional)
# ----------------------------
if [[ "$RUN_FASTQC" == "true" ]]; then
  echo "== Step 1: FastQC =="
  fastqc "$R1" -o "$READS_DIR"
  fastqc "$R2" -o "$READS_DIR"
else
  echo "== Step 1: FastQC skipped (RUN_FASTQC=false) =="
fi

# ----------------------------
# Step 2: Align with BWA-MEM -> sort -> index
# ----------------------------
echo "== Step 2: BWA-MEM align -> sort -> index =="

SAM="$ALN_DIR/${SAMPLE}.paired.sam"
SORTED_BAM="$ALN_DIR/${SAMPLE}_sorted.bam"

# Align to SAM (keep for teaching/demo; for production you'd stream to BAM)
if [[ ! -s "$SAM" ]]; then
  bwa mem -t "$THREADS" \
    -R "@RG\tID:${SAMPLE}\tPL:ILLUMINA\tSM:${SAMPLE}" \
    "$REF" "$R1" "$R2" > "$SAM"
fi

# Convert+sort to BAM
if [[ ! -s "$SORTED_BAM" ]]; then
  samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$SAM"
fi

# Index sorted BAM
if [[ ! -e "${SORTED_BAM}.bai" ]]; then
  samtools index "$SORTED_BAM"
fi

# ----------------------------
# Step 3: Mark duplicates (Spark) + metrics + index
# ----------------------------
echo "== Step 3: MarkDuplicatesSpark =="

DEDUP_BAM="$ALN_DIR/${SAMPLE}_sorted_dedup_reads.bam"
MD_METRICS="$ALN_DIR/${SAMPLE}_markdups_metrics.txt"

if [[ ! -s "$DEDUP_BAM" ]]; then
  gatk MarkDuplicatesSpark \
    -I "$SORTED_BAM" \
    -O "$DEDUP_BAM" \
    -M "$MD_METRICS"
fi

if [[ ! -e "${DEDUP_BAM}.bai" ]]; then
  samtools index "$DEDUP_BAM"
fi

# ----------------------------
# Step 4: BQSR (BaseRecalibrator + ApplyBQSR) + index
# ----------------------------
echo "== Step 4: BQSR =="

RECAL_TABLE="$DATA_DIR/recal_data.table"
BQSR_BAM="$ALN_DIR/${SAMPLE}_sorted_dedup_bqsr_reads.bam"

if [[ ! -s "$RECAL_TABLE" ]]; then
  gatk BaseRecalibrator \
    -R "$REF" \
    -I "$DEDUP_BAM" \
    --known-sites "$KNOWN_SITES" \
    -O "$RECAL_TABLE"
fi

if [[ ! -s "$BQSR_BAM" ]]; then
  gatk ApplyBQSR \
    -R "$REF" \
    -I "$DEDUP_BAM" \
    --bqsr-recal-file "$RECAL_TABLE" \
    -O "$BQSR_BAM"
fi

if [[ ! -e "${BQSR_BAM}.bai" ]]; then
  samtools index "$BQSR_BAM"
fi

# ----------------------------
# Step 5: HaplotypeCaller -> gVCF
# ----------------------------
echo "== Step 5: HaplotypeCaller (gVCF) =="

GVCF="$VARIANTS_DIR/${SAMPLE}.g.vcf.gz"

if [[ ! -s "$GVCF" ]]; then
  gatk HaplotypeCaller \
    -R "$REF" \
    -I "$BQSR_BAM" \
    -O "$GVCF" \
    -ERC GVCF
fi

# ----------------------------
# Step 6: GenotypeGVCFs -> raw VCF
# ----------------------------
echo "== Step 6: GenotypeGVCFs =="

RAW_VCF="$RESULTS_DIR/${SAMPLE}.raw.vcf.gz"

if [[ ! -s "$RAW_VCF" ]]; then
  gatk GenotypeGVCFs \
    -R "$REF" \
    -V "$GVCF" \
    -O "$RAW_VCF"
fi

# ----------------------------
# Step 7: Split SNPs & INDELs
# ----------------------------
echo "== Step 7: SelectVariants (SNPs/INDELs) =="

RAW_SNPS="$RESULTS_DIR/${SAMPLE}.raw.snps.vcf.gz"
RAW_INDELS="$RESULTS_DIR/${SAMPLE}.raw.indels.vcf.gz"

if [[ ! -s "$RAW_SNPS" ]]; then
  gatk SelectVariants \
    -R "$REF" \
    -V "$RAW_VCF" \
    --select-type SNP \
    -O "$RAW_SNPS"
fi

if [[ ! -s "$RAW_INDELS" ]]; then
  gatk SelectVariants \
    -R "$REF" \
    -V "$RAW_VCF" \
    --select-type INDEL \
    -O "$RAW_INDELS"
fi

# ----------------------------
# Step 8: Hard filter SNPs/INDELs
# ----------------------------
echo "== Step 8: VariantFiltration =="

FILT_SNPS="$RESULTS_DIR/${SAMPLE}.filtered.snps.vcf.gz"
FILT_INDELS="$RESULTS_DIR/${SAMPLE}.filtered.indels.vcf.gz"

if [[ ! -s "$FILT_SNPS" ]]; then
  gatk VariantFiltration \
    -R "$REF" \
    -V "$RAW_SNPS" \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_FAIL" \
    -O "$FILT_SNPS"
fi

if [[ ! -s "$FILT_INDELS" ]]; then
  gatk VariantFiltration \
    -R "$REF" \
    -V "$RAW_INDELS" \
    --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filter-name "INDEL_FAIL" \
    -O "$FILT_INDELS"
fi

# ----------------------------
# Step 9: Keep PASS only + merge
# ----------------------------
echo "== Step 9: PASS-only + Merge =="

PASS_SNPS="$RESULTS_DIR/${SAMPLE}.PASS.snps.vcf.gz"
PASS_INDELS="$RESULTS_DIR/${SAMPLE}.PASS.indels.vcf.gz"
PASS_MERGED="$RESULTS_DIR/${SAMPLE}.PASS.merged.vcf.gz"

if [[ ! -s "$PASS_SNPS" ]]; then
  gatk SelectVariants \
    -R "$REF" \
    -V "$FILT_SNPS" \
    --exclude-filtered \
    -O "$PASS_SNPS"
fi

if [[ ! -s "$PASS_INDELS" ]]; then
  gatk SelectVariants \
    -R "$REF" \
    -V "$FILT_INDELS" \
    --exclude-filtered \
    -O "$PASS_INDELS"
fi

if [[ ! -s "$PASS_MERGED" ]]; then
  gatk MergeVcfs \
    -I "$PASS_SNPS" \
    -I "$PASS_INDELS" \
    -O "$PASS_MERGED"
fi

# ----------------------------
# Step 10: QC counts + bcftools stats + Ti/Tv
# ----------------------------
echo "== Step 10: QC (counts + stats + Ti/Tv) =="

echo "PASS SNP count:"
bcftools view -H "$PASS_SNPS" | wc -l

echo "PASS INDEL count:"
bcftools view -H "$PASS_INDELS" | wc -l

echo "PASS merged count:"
bcftools view -H "$PASS_MERGED" | wc -l

STATS_PASS="$RESULTS_DIR/PASS.stats.txt"
bcftools stats "$PASS_MERGED" > "$STATS_PASS"

echo "Ti/Tv line:"
grep "TSTV" "$STATS_PASS" || true

# ----------------------------
# Step 11: Optional GT subset DP/GQ (bcftools)
# ----------------------------
if [[ "$MAKE_GT_SUBSET" == "true" ]]; then
  echo "== Step 11: Genotype subset (DP>=${MIN_DP}, GQ>=${MIN_GQ}) =="

  GT_SNPS="$RESULTS_DIR/${SAMPLE}.PASS.GTDP${MIN_DP}_GQ${MIN_GQ}.snps.vcf.gz"
  GT_INDELS="$RESULTS_DIR/${SAMPLE}.PASS.GTDP${MIN_DP}_GQ${MIN_GQ}.indels.vcf.gz"
  GT_MERGED="$RESULTS_DIR/${SAMPLE}.PASS.GTDP${MIN_DP}_GQ${MIN_GQ}.merged.vcf.gz"

  # For single-sample: site is kept if that sample's FORMAT/DP and FORMAT/GQ meet thresholds
  bcftools view -i "FILTER=\"PASS\" && FORMAT/DP>=${MIN_DP} && FORMAT/GQ>=${MIN_GQ}" -Oz \
    -o "$GT_SNPS" "$PASS_SNPS"
  bcftools index -t "$GT_SNPS"

  bcftools view -i "FILTER=\"PASS\" && FORMAT/DP>=${MIN_DP} && FORMAT/GQ>=${MIN_GQ}" -Oz \
    -o "$GT_INDELS" "$PASS_INDELS"
  bcftools index -t "$GT_INDELS"

  gatk MergeVcfs -I "$GT_SNPS" -I "$GT_INDELS" -O "$GT_MERGED"

  echo "GT-subset SNP/INDEL/merged counts:"
  bcftools view -H "$GT_SNPS" | wc -l
  bcftools view -H "$GT_INDELS" | wc -l
  bcftools view -H "$GT_MERGED" | wc -l
fi

# ----------------------------
# Step 12: Optional TSV export
# ----------------------------
if [[ "$MAKE_TSV" == "true" ]]; then
  echo "== Step 12: Export TSV table (basic fields) =="

  TSV_OUT="$RESULTS_DIR/${SAMPLE}.PASS.basic.tsv"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP\t%AF\n' "$PASS_MERGED" > "$TSV_OUT"
  echo "Wrote: $TSV_OUT"
fi

echo "✅ DONE"
echo "Final deliverables:"
echo "  PASS merged VCF:  $PASS_MERGED"
echo "  PASS stats:       $STATS_PASS"
