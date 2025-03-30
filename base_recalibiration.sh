#!/bin/bash
#SBATCH -J base_recalibration
#SBATCH -p general
#SBATCH -o base_recalibration_%j.out
#SBATCH -e base_recalibration_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH --mem=100GB
#SBATCH -A r00750

# Load miniconda
module load miniconda

# Activate the conda environment
conda activate variant_calling_env

# Set paths to reference, input/output directories, and known sites
REFERENCE="/N/scratch/pthuthi/mouse_vc/reference_data/Mus_musculus.GRCm39.dna.primary_assembly.fa"
INPUT_DIR="/N/scratch/pthuthi/mouse_vc/mark_duplicates_output"
OUTPUT_DIR="/N/scratch/pthuthi/mouse_vc/base_recalibration_output"
KNOWN_SITES="/N/scratch/pthuthi/mouse_vc/reference_data/mus_musculus.vcf.gz"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Step 1: Base Recalibration Table Generation
for DEDUP_BAM_FILE in ${INPUT_DIR}/*_dedup.bam; do
  SAMPLE_NAME=$(basename $DEDUP_BAM_FILE _dedup.bam)
  
  gatk BaseRecalibrator \
      -R $REFERENCE \
      -I $DEDUP_BAM_FILE \
      --known-sites $KNOWN_SITES \
      -O ${OUTPUT_DIR}/${SAMPLE_NAME}_recal_data.table

done

# Step 2: Apply Base Quality Score Recalibration
for DEDUP_BAM_FILE in ${INPUT_DIR}/*_dedup.bam; do
  SAMPLE_NAME=$(basename $DEDUP_BAM_FILE _dedup.bam)
  
  gatk ApplyBQSR \
      -R $REFERENCE \
      -I $DEDUP_BAM_FILE \
      --bqsr-recal-file ${OUTPUT_DIR}/${SAMPLE_NAME}_recal_data.table \
      -O ${OUTPUT_DIR}/${SAMPLE_NAME}_recalibrated.bam

done

echo "Base recalibration process complete. Outputs are saved in ${OUTPUT_DIR}."

