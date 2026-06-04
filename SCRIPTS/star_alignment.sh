#!/bin/bash

# =========================
# STAR TEST ALIGNMENT
# 2 samples only
# Mus musculus GRCm39
# =========================

GENOME_DIR="/group/sbms006/rrubimsilvadocouto/reference_genome/mouse_GRCm39/STAR_index"

INPUT_DIR="/scratch/sbms006/rrubimsilvadocouto/rnaseq_trim/trimmed"

OUT_DIR="/scratch/sbms006/rrubimsilvadocouto/rnaseq_star_test"

mkdir -p ${OUT_DIR}/Liver_CTRL_1
mkdir -p ${OUT_DIR}/Liver_CTRL_2

echo "Running STAR for Liver-CTRL-1"

STAR \
--runThreadN 8 \
--genomeDir ${GENOME_DIR} \
--readFilesIn \
${INPUT_DIR}/Liver-CTRL-1_S1_L001_R1_001_val_1.fq.gz \
${INPUT_DIR}/Liver-CTRL-1_S1_L001_R2_001_val_2.fq.gz \
--readFilesCommand zcat \
--outFileNamePrefix ${OUT_DIR}/Liver_CTRL_1/ \
--outSAMtype BAM SortedByCoordinate

echo "Running STAR for Liver-CTRL-2"

STAR \
--runThreadN 8 \
--genomeDir ${GENOME_DIR} \
--readFilesIn \
${INPUT_DIR}/Liver-CTRL-2_S2_L001_R1_001_val_1.fq.gz \
${INPUT_DIR}/Liver-CTRL-2_S2_L001_R2_001_val_2.fq.gz \
--readFilesCommand zcat \
--outFileNamePrefix ${OUT_DIR}/Liver_CTRL_2/ \
--outSAMtype BAM SortedByCoordinate

echo "DONE STAR alignment"
