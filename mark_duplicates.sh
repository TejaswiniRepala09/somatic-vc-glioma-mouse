#!/bin/bash
#SBATCH -J mark_duplicates
#SBATCH -p general
#SBATCH -o mark_duplicates_%j.out
#SBATCH -e mark_duplicates_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=100GB
#SBATCH -A r00750

# Load miniconda
module load miniconda

# Activate the conda environment
conda activate variant_calling_env

# Set paths to reference and input/output directories
BAM_DIR="/N/scratch/pthuthi/mouse_vc/bwa_output"
OUTPUT_DIR="/N/scratch/pthuthi/mouse_vc/mark_duplicates_output"
READ_GROUP_ADDED_DIR="/N/scratch/pthuthi/mouse_vc/read_group_added"
PICARD_JAR="/N/u/pthuthi/Quartz/.conda/envs/variant_calling_env/share/picard-3.3.0-0/picard.jar"

# Create output directories if they don't exist
mkdir -p $READ_GROUP_ADDED_DIR
mkdir -p $OUTPUT_DIR

# Step 1: Add or Replace Read Groups for Each BAM File
for BAM_FILE in ${BAM_DIR}/*.bam; do
  SAMPLE_NAME=$(basename $BAM_FILE _aligned_sorted.bam)
  
  # Add read groups using Picard
  java -jar $PICARD_JAR AddOrReplaceReadGroups \
      I=$BAM_FILE \
      O=${READ_GROUP_ADDED_DIR}/${SAMPLE_NAME}_rg.bam \
      RGID=$SAMPLE_NAME \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=$SAMPLE_NAME
done

# Step 2: Mark Duplicates Using GATK
for RG_BAM_FILE in ${READ_GROUP_ADDED_DIR}/*_rg.bam; do
  SAMPLE_NAME=$(basename $RG_BAM_FILE _rg.bam)
  
  gatk MarkDuplicates \
      -I $RG_BAM_FILE \
      -O ${OUTPUT_DIR}/${SAMPLE_NAME}_dedup.bam \
      -M ${OUTPUT_DIR}/${SAMPLE_NAME}_dedup_metrics.txt \
      --REMOVE_DUPLICATES false
done

# Print completion message
echo "MarkDuplicates process complete. Outputs are saved in ${OUTPUT_DIR}."
