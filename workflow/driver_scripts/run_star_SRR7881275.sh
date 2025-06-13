#!/bin/sh

GENOME_DIR="out/star_mm10_index"
GTF_FILE="out/gencode.vM25.annotation.nodot.gtf"
FILTERED_DIR="out/bowtie2_filtered"
Aligned_DIR="out/STAR_output"
LOG_DIR="logs"

mkdir -p "$Aligned_DIR" "$LOG_DIR"

# Loop over filtered fastq pairs with naming like SRR7881275_filtered.fastq.1.gz
for READ1 in ${FILTERED_DIR}/*_filtered.fastq.1.gz; do
    SAMPLE=$(basename "$READ1" _filtered.fastq.1.gz)
    READ2="${FILTERED_DIR}/${SAMPLE}_filtered.fastq.2.gz"
    
    if [[ ! -f "$READ2" ]]; then
        echo "Missing READ2 for $SAMPLE. Skipping..."
        continue
    fi

    OUT_PREFIX="${Aligned_DIR}/${SAMPLE}_"
    LOGFILE="${LOG_DIR}/STAR_${SAMPLE}.log"

         /home/mkumar4/miniconda3/envs/star_env/bin/STAR --runThreadN 10 \
              --genomeDir "$GENOME_DIR" \
              --genomeLoad NoSharedMemory \
              --readFilesIn "$READ1" "$READ2" \
              --readFilesCommand zcat \
              --outFileNamePrefix "$OUT_PREFIX" \
              --outReadsUnmapped Fastx \
              --outSAMattributes NH HI AS nM NM MD jM jI XS \
              --outSJfilterOverhangMin 15 15 15 15 \
              --outFilterMultimapNmax 20 \
              --outFilterScoreMin 1 \
              --outFilterMatchNminOverLread 0.7 \
              --outFilterMismatchNmax 999 \
              --outFilterMismatchNoverLmax 0.05 \
              --alignIntronMin 20 \
              --alignIntronMax 1000000 \
              --alignMatesGapMax 1000000 \
              --alignSJoverhangMin 15 \
              --alignSJDBoverhangMin 10 \
              --alignSoftClipAtReferenceEnds No \
              --chimSegmentMin 15 \
              --chimScoreMin 15 \
              --chimScoreSeparation 10 \
              --chimJunctionOverhangMin 15 \
              --sjdbGTFfile "$GTF_FILE" \
              --quantMode GeneCounts \
              --twopassMode Basic \
              --chimOutType Junctions SeparateSAMold
done

