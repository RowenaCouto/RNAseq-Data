FASTQ_DIR="/group/sbms006/rrubimsilvadocouto/MiSeq_Data/fastqs"
OUT_DIR="/scratch/sbms006/rrubimsilvadocouto/rnaseq_trim/trimmed_all"

mkdir -p "$OUT_DIR"

for R1 in ${FASTQ_DIR}/*_R1_001.fastq.gz
do
    SAMPLE=$(basename "$R1" _R1_001.fastq.gz)
    R2="${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz"

    if [ ! -f "$R2" ]; then
        echo "Missing pair for $SAMPLE"
        continue
    fi

    echo "Processing $SAMPLE"

    trim_galore \
        --paired \
        --fastqc \
        --output_dir "$OUT_DIR" \
        "$R1" "$R2"
done

