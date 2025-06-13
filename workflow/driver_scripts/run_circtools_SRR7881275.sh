#!/bin/bash
set -e  # Exit immediately if any command fails

# Define job name
JOB_NAME="circtools_detect_bowtie2_SRR7881275_Pi"

zcat out/bowtie2_filtered/SRR7881275_filtered.fastq.1.gz | \
awk '{
  if (NR % 4 == 1) {
    sub(/^@/, "", $0); print "@" $1 " 1:N:0:1"
  } else { print }
}' > out/bowtie2_filtered/SRR7881275_cleaned_1.fixed.fastq

zcat out/bowtie2_filtered/SRR7881275_filtered.fastq.2.gz | \
awk '{
  if (NR % 4 == 1) {
    sub(/^@/, "", $0); print "@" $1 " 2:N:0:1"
  } else { print }
}' > out/bowtie2_filtered/SRR7881275_cleaned_2.fixed.fastq


# Define output paths
OUTPUT_LOG="logs/${JOB_NAME}.log"
OUT_DIR="out/STAR_output/SRR7881275__Pi_circtools"

# Create output directories if not existing
mkdir -p logs
mkdir -p "$OUT_DIR"

# Run circtools detect
circtools detect \
    out/STAR_output/SRR7881275_Chimeric.out.junction \
    -B out/STAR_output/SRR7881275_Aligned.sorted.bam \
    -mt1 out/bowtie2_filtered/SRR7881275_cleaned_1.fixed.fastq \
    -mt2 out/bowtie2_filtered/SRR7881275_cleaned_2.fixed.fastq \
    -A in/Mus_musculus.GRCm38.dna.primary_assembly.fa \
    -an in/gencode.vM25.annotation.nodot.gtf \
    -O "$OUT_DIR" \
    -R in/mouse_repeats.gtf \
    -T 8 \
    -M -D -fg -Pi -G \
    &> "$OUTPUT_LOG"
