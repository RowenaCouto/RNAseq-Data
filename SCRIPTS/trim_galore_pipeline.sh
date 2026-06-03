#!/bin/bash

#############################################
# RNA-seq PIPELINE — TRIM GALORE STEP
# User: rrubimsilvadocouto
# Project: sbms006
#############################################

echo "=============================="
echo "Starting Trim Galore pipeline"
echo "Date:"
date
echo "=============================="

#############################################
# 1. CHECK WHERE WE ARE
#############################################

echo "Current directory:"
pwd

#############################################
# 2. DEFINE PATHS
#############################################

GROUP_DIR="/group/sbms006/rrubimsilvadocouto/fastqs"
SCRATCH_DIR="/scratch/sbms006/rrubimsilvadocouto/rnaseq_trim"

echo "Group directory:"
echo $GROUP_DIR

echo "Scratch directory:"
echo $SCRATCH_DIR

#############################################
# 3. MOVE TO SCRATCH
#############################################

cd $SCRATCH_DIR || mkdir -p $SCRATCH_DIR
cd $SCRATCH_DIR

echo "Now working in:"
pwd

#############################################
# 4. LINK RAW FASTQ FILES
#############################################

echo "Linking FASTQ files..."
ln -s $GROUP_DIR/*.fastq.gz .

echo "FASTQ files found:"
ls *.fastq.gz | head

#############################################
# 5. LOAD ENVIRONMENT (if needed)
#############################################

# If running outside conda environment, uncomment:
# source activate myenv

echo "Using Trim Galore version:"
trim_galore --version

#############################################
# 6. TEST RUN (ONE SAMPLE)
#############################################

echo "Running test sample: Liver-CTRL-1"

trim_galore \
  --paired \
  --fastqc \
  Liver-CTRL-1_S1_L001_R1_001.fastq.gz \
  Liver-CTRL-1_S1_L001_R2_001.fastq.gz

echo "Test run completed"
echo "Check output files:"
ls *val* | head

#############################################
# 7. BATCH RUN (ALL SAMPLES)
#############################################

echo "Starting batch trimming for all samples..."

for r1 in *_R1_001.fastq.gz
do
    r2=${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}

    echo "Processing:"
    echo $r1

    trim_galore \
        --paired \
        --fastqc \
        --output_dir . \
        "$r1" "$r2"
done

echo "All samples processed"
echo "Pipeline finished at:"
date