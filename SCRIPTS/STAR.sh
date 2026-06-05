#!/bin/bash

TRIM_DIR="/scratch/sbms006/rrubimsilvadocouto/rnaseq_trim/trimmed_all"
OUT_DIR="/scratch/sbms006/rrubimsilvadocouto/rnaseq_star_all"
GENOME="/group/sbms006/rrubimsilvadocouto/reference_genome/mouse_GRCm39/STAR_index"

mkdir -p "$OUT_DIR"

for R1 in ${TRIM_DIR}/*_val_1.fq.gz
do
    SAMPLE=$(basename "$R1" _R1_001_val_1.fq.gz)

    R2="${TRIM_DIR}/${SAMPLE}_R2_001_val_2.fq.gz"

    mkdir -p "${OUT_DIR}/${SAMPLE}"

    echo "Running STAR for ${SAMPLE}"

    STAR \
        --runThreadN 8 \
        --genomeDir "$GENOME" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${OUT_DIR}/${SAMPLE}/" \
        --outSAMtype BAM SortedByCoordinate

done
