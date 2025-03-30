#!/bin/bash
#SBATCH -J fastqc_analysis_all    
#SBATCH -p general                
#SBATCH -o fastqc_output_%j.txt   
#SBATCH -e fastqc_error_%j.err    
#SBATCH --mail-type=ALL           
#SBATCH --mail-user=pthuthi@iu.edu 
#SBATCH --nodes=1                 
#SBATCH --ntasks-per-node=1       
#SBATCH --cpus-per-task=4         
#SBATCH --time=6:00:00            
#SBATCH --mem=20GB                
#SBATCH -A r00750 

# Load the miniconda module
module load minconda

# Activate the conda environment
conda activate v_c

# List of directories to process
DATA_DIRS=(
    "/N/scratch/pthuthi/mouse_vc/SRR24288057"
    "/N/scratch/pthuthi/mouse_vc/SRR24288060"
    "/N/scratch/pthuthi/mouse_vc/SRR24288062"
)

# Loop through each directory and run FastQC on the paired-end files
for DATA_DIR in "${DATA_DIRS[@]}"; do
    # Create output directory for FastQC reports
    OUTPUT_DIR="${DATA_DIR}/fastqc_output"
    mkdir -p ${OUTPUT_DIR}

    # Run FastQC on the paired-end FASTQ files in the directory
    fastqc -o ${OUTPUT_DIR} -t 4 ${DATA_DIR}/*_1.fastq ${DATA_DIR}/*_2.fastq
done