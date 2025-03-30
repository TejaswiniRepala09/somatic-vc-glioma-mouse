#!/bin/bash
#SBATCH -J somatic_variant_filtering
#SBATCH -p general
#SBATCH -o somatic_variant_filtering_%j.out
#SBATCH -e somatic_variant_filtering_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem=100GB
#SBATCH -A r00750

# Load miniconda
module load miniconda

# Activate the variant calling environment
conda activate variant_calling_env

# Paths to important files
GATK_JAR="/geode2/home/u070/pthuthi/Quartz/.conda/envs/variant_calling_env/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar"
REFERENCE="/N/scratch/pthuthi/mouse_vc/reference_data/Mus_musculus.GRCm39.dna.primary_assembly.fa"
INPUT_VCF1="/N/scratch/pthuthi/mouse_vc/gatk_somatic_vc_output/somatic_variants_tumor1.vcf.gz"
INPUT_VCF2="/N/scratch/pthuthi/mouse_vc/gatk_somatic_vc_output/somatic_variants_tumor2.vcf.gz"
OUTPUT_FILTERED_VCF1="/N/scratch/pthuthi/mouse_vc/gatk_filtered_vc_output/filtered_somatic_variants_tumor1.vcf.gz"
OUTPUT_FILTERED_VCF2="/N/scratch/pthuthi/mouse_vc/gatk_filtered_vc_output/filtered_somatic_variants_tumor2.vcf.gz"

# Create output directory if it doesn't exist
mkdir -p /N/scratch/pthuthi/mouse_vc/gatk_filtered_vc_output

# Run GATK FilterMutectCalls for Sample 1
java -Dsamjdk.use_async_io_read_samtools=false \
     -Dsamjdk.use_async_io_write_samtools=true \
     -Dsamjdk.use_async_io_write_tribble=false \
     -Dsamjdk.compression_level=2 \
     -jar $GATK_JAR FilterMutectCalls \
     -R $REFERENCE \
     -V $INPUT_VCF1 \
     -O $OUTPUT_FILTERED_VCF1

# Run GATK FilterMutectCalls for Sample 2
java -Dsamjdk.use_async_io_read_samtools=false \
     -Dsamjdk.use_async_io_write_samtools=true \
     -Dsamjdk.use_async_io_write_tribble=false \
     -Dsamjdk.compression_level=2 \
     -jar $GATK_JAR FilterMutectCalls \
     -R $REFERENCE \
     -V $INPUT_VCF2 \
     -O $OUTPUT_FILTERED_VCF2

echo "Variant filtering completed for both tumor samples."
