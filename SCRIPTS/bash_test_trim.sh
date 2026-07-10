#!/bin/bash
#SBATCH --job-name=trim_galore_all
#SBATCH --output=trim_galore_%A_%a.out
#SBATCH --error=trim_galore_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --partition=work
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load trim_galore
module load fastqc

FASTQ_DIR="/group/sbms006/rrubimsilvadocouto/MiSeq_Data/fastqs"
OUT_DIR="/scratch/sbms006/rrubimsilvadocouto/rnaseq_trim/trimmed_all"

mkdir -p "$OUT_DIR"

echo "Starting trimming job on $(date)"

for R1 in "$FASTQ_DIR"/*_R1_001.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
    R2="${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "Missing R2 for sample: $SAMPLE"
        continue
    fi

    echo "Processing $SAMPLE"

    trim_galore \
        --paired \
        --fastqc \
        --cores 4 \
        --output_dir "$OUT_DIR" \
        "$R1" "$R2"

    echo "Finished $SAMPLE"
done

echo "All trimming completed on $(date)"
