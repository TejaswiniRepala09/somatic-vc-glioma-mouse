#!/bin/bash
#SBATCH -J mark_duplicates_rg
#SBATCH -p general
#SBATCH -o mark_duplicates_rg_%j.out
#SBATCH -e mark_duplicates_rg_%j.err
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

# Activate the variant calling environment
conda activate variant_calling_env

# Set paths to input and output directories
READ_GROUP_ADDED_DIR="/N/scratch/pthuthi/mouse_vc/read_group_added"
OUTPUT_DIR="/N/scratch/pthuthi/mouse_vc/mark_duplicates_output"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

for RG_BAM_FILE in ${READ_GROUP_ADDED_DIR}/*_rg.bam; do
  SAMPLE_NAME=$(basename $RG_BAM_FILE _rg.bam)
  
  gatk MarkDuplicates \
      -I $RG_BAM_FILE \
      -O ${OUTPUT_DIR}/${SAMPLE_NAME}_dedup.bam \
      -M ${OUTPUT_DIR}/${SAMPLE_NAME}_dedup_metrics.txt \
      --REMOVE_DUPLICATES false
done

echo "MarkDuplicates process complete for _rg.bam files. Outputs are saved in ${OUTPUT_DIR}."
