#!/bin/sh

RRNA_INDEX="out/rRNA_cluster/mus-musculus.rRNA"
TRIMMED_DIR="out/trimmed"
TARGET_DIR="out/bowtie2_filtered"
LOG_DIR="logs"

mkdir -p "$TARGET_DIR"
mkdir -p "$LOG_DIR"

SAMPLE="SRR7881275"
READ1="${TRIMMED_DIR}/${SAMPLE}_1.fastq.gz"
READ2="${TRIMMED_DIR}/${SAMPLE}_2.fastq.gz"

# Job name and log
JOB="bowtie2_filter_${SAMPLE}_paired"
LOG="${LOG_DIR}/${JOB}.log"

echo "Submitting bowtie2 filtering job for ${SAMPLE} (paired-end)"
     bowtie2 -x "$RRNA_INDEX" \
             -1 "$READ1" \
             -2 "$READ2" \
             --no-unal \
             --omit-sec-seq \
             --threads 20 \
             --mm \
             --seed 1337 \
             --time \
             --un-conc-gz "${TARGET_DIR}/${SAMPLE}_filtered.fastq.gz" \
             -S /dev/null \
             2> "${TARGET_DIR}/${SAMPLE}_bowtie2.log"
