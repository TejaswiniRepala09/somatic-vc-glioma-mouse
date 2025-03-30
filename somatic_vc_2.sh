#!/bin/bash
#SBATCH -J somatic_variant_calling_2
#SBATCH -p general
#SBATCH -o somatic_variant_calling_2%j.out
#SBATCH -e somatic_variant_calling_2%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem=200GB
#SBATCH -A r00750

# Load miniconda
module load miniconda

# Activate the variant calling environment
conda activate variant_calling_env

# Paths to important files
GATK_JAR="/geode2/home/u070/pthuthi/Quartz/.conda/envs/variant_calling_env/share/gatk4-4.6.1.0-0/gatk-package-4.6.1.0-local.jar"
REFERENCE="/N/scratch/pthuthi/mouse_vc/reference_data/Mus_musculus.GRCm39.dna.primary_assembly.fa"
TUMOR_BAM="/N/scratch/pthuthi/mouse_vc/base_recalibration_output/SRR24288062_recalibrated.bam"
NORMAL_BAM="/N/scratch/pthuthi/mouse_vc/base_recalibration_output/SRR24288057_recalibrated.bam"
OUTPUT_VCF="/N/scratch/pthuthi/mouse_vc/gatk_somatic_vc_output/somatic_variants_tumor2.vcf.gz"

# Run GATK Mutect2 for Sample 2 (Tumor: SRR24288062, Normal: SRR24288057)
java -Dsamjdk.use_async_io_read_samtools=false \
     -Dsamjdk.use_async_io_write_samtools=true \
     -Dsamjdk.use_async_io_write_tribble=false \
     -Dsamjdk.compression_level=2 \
     -jar $GATK_JAR Mutect2 \
     -R $REFERENCE \
     -I $TUMOR_BAM \
     -I $NORMAL_BAM \
     -tumor SRR24288062 \
     -normal SRR24288057 \
     -O $OUTPUT_VCF \
     --native-pair-hmm-threads 8  # Use 8 threads for faster processing
