#!/bin/bash
#SBATCH -J bwa_mapping_mouse_multi
#SBATCH -p general
#SBATCH -o bwa_mapping_output_%j.txt
#SBATCH -e bwa_mapping_error_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=100GB
#SBATCH -A r00750

# Load the miniconda module
module load miniconda

# Activate the conda enviornment
conda activate variant_calling_env

# Set file paths
REFERENCE="/N/scratch/pthuthi/mouse_vc/reference_data/Mus_musculus.GRCm39.dna.primary_assembly.fa"

# Create output directory if it does not exist
OUTPUT_DIR="/N/scratch/pthuthi/mouse_vc/bwa_output"
mkdir -p ${OUTPUT_DIR}

# Step 1: Index the reference genome (if not already indexed)
if [ ! -f "${REFERENCE}.bwt" ]; then
    echo "Indexing the reference genome..."
    bwa index ${REFERENCE}
fi

# Array of sample IDs
SAMPLES=("SRR24288057" "SRR24288060" "SRR24288062")

# Loop through each sample
for SAMPLE in "${SAMPLES[@]}"; do
    INPUT1="/N/scratch/pthuthi/mouse_vc/${SAMPLE}/trimgalore_output/${SAMPLE}_1_val_1.fq"
    INPUT2="/N/scratch/pthuthi/mouse_vc/${SAMPLE}/trimgalore_output/${SAMPLE}_2_val_2.fq"
    OUTPUT_BAM="${OUTPUT_DIR}/${SAMPLE}_aligned.bam"
    OUTPUT_SORTED_BAM="${OUTPUT_DIR}/${SAMPLE}_aligned_sorted.bam"

    echo "Processing sample ${SAMPLE}..."

    # Align reads using BWA
    bwa mem -t 8 ${REFERENCE} ${INPUT1} ${INPUT2} > ${OUTPUT_DIR}/${SAMPLE}_aligned.sam

    # Convert SAM to BAM
    samtools view -bS ${OUTPUT_DIR}/${SAMPLE}_aligned.sam > ${OUTPUT_BAM}

    # Sort and index the BAM file
    samtools sort -o ${OUTPUT_SORTED_BAM} ${OUTPUT_BAM}
    samtools index ${OUTPUT_SORTED_BAM}

    # Clean up intermediate files
    rm ${OUTPUT_DIR}/${SAMPLE}_aligned.sam
    rm ${OUTPUT_BAM}
done

echo "All samples processed successfully!"