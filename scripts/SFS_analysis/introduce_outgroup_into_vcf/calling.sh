#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=15gb
#SBATCH --time=02:15:00
#SBATCH --account=UniKoeln
#SBATCH --output=/home/group.kurse/qbcbp008/project/SFSs/scripts/calling.out
#SBATCH --error=/home/group.kurse/qbcbp008/project/SFSs/scripts/calling.err



# Create a new directory “day_3” at “/home/group.kurse/username/” that contains following subdirectories: filtered_bam, variant_calling, variant_filtering
mkdir -p calling/{filtered_bam,variant_calling,variant_filtering}

# Run the following commands using “for loop” to filter and index bam files based on mapping quality

# load required tool
module load samtools

# Following commmand can be used to list all sorted bam files with path, which might be useful in some cases
#$ samples=$(ls /home/group.kurse/qbcbp008/day_2/mapping_bwa/*Chr4_sorted.bam | sed 's/\.Chr4_sorted.bam//')
#$ echo "$samples"

# However, you only need to list only prefixes of all sorted bam files
sample="SRR2082718"

# Quality filtering and removing PCR duplicates with Samtools 
# keep only properly aligned paired-end reads (-f 3)
# Mapping quality (-q 30)
# Exclude secondary alignments and reads failing quality checks (-F 264)
# Remove duplicates (samtools rmdup)

samtools view -b -o /home/group.kurse/qbcbp008/project/SFSs/calling/filtered_bam/${sample}_BamQualFilt-f3Q30F264.bam -f 3 -q 30 -F 264 /home/group.kurse/qbcbp008/project/SFSs/mapping/${sample}.sorted.bam

samtools rmdup -s /home/group.kurse/qbcbp008/project/SFSs/calling/filtered_bam/${sample}_BamQualFilt-f3Q30F264.bam /home/group.kurse/qbcbp008/project/SFSs/calling/filtered_bam/${sample}_BamQualFilt-f3Q30F264_nodup.bam
  
# Index the filtered_nodup.bam file
samtools index /home/group.kurse/qbcbp008/project/SFSs/calling/filtered_bam/${sample}_BamQualFilt-f3Q30F264_nodup.bam



###

# Call and Filter Variants

# load required module  
module load bcftools/1.18

# Set path for input and output directories, and reference fasta
input_dir="/home/group.kurse/qbcbp008/project/SFSs/calling/filtered_bam/"
output_dir="/home/group.kurse/qbcbp008/project/SFSs/calling/variant_calling/"
reference_fasta="/home/group.kurse/qbcbp008/project/SFSs/reference/TAIR10.1.genomic_fna"

# list and print all filtered depublicated bam files to a text file
files=$(find /home/group.kurse/qbcbp008/project/SFSs/calling/filtered_bam/ -type f -name "*_nodup.bam")
echo "$files" > /home/group.kurse/qbcbp008/project/SFSs/calling/filtered_bam/bam_files.txt

# Call variants using bcftools mpileup (15 min) and call (3-5 minutes)commands
# Activate Base Alignment Quality computation (-E)
# Minimum base quality (-q 30)
# Minimum calling threshold for variant alleles (-p 0.01) Variants with an allele frequency of at least 1% will be called.

bcftools mpileup -E -q 30 --threads 20 -o "${output_dir}/all_samples_BamQualFilt-f3Q30F264_nodup.bcf" -f "$reference_fasta" -b "${input_dir}/bam_files.txt"

bcftools call -c -p 0.01 -O z --threads 20 -o "${output_dir}/all_samples_BamQualFilt-f3Q30F264_nodup_raw.vcf.gz" "${output_dir}/all_samples_BamQualFilt-f3Q30F264_nodup.bcf"

# Basic filtering by removing all monomorphic variant sites
bcftools view -i 'AC>0' "${output_dir}/all_samples_BamQualFilt-f3Q30F264_nodup_raw.vcf.gz" -o "${output_dir}/all_samples_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic.vcf.gz"

# Check names of samples included in the vcf file
bcftools query -l ${output_dir}/all_samples_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic.vcf.gz

# Simplify and shorten the names of samples
zcat ${output_dir}/all_samples_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic.vcf.gz | \
awk -F '\t' 'BEGIN{OFS="\t"} {if ($1 ~ /^#CHROM/) {for (i=10; i<=NF; i++) {sub(".*/", "", $i); sub("_BamQualFilt-f3Q30F264_nodup\\.bam", "", $i)}} print }' | \
gzip -c > ${output_dir}/all_samples_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam.vcf.gz

# Re-check names of samples included in the vcf file
bcftools query -l ${output_dir}/all_samples_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam.vcf.gz

# Do the following to load vcftools 
module unload bcftools/1.18
module unload gnu
module load vcftools/0.1.17

# Additional filters

vcftools --gzvcf "${output_dir}/all_samples_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam.vcf.gz" --minDP 10 --minGQ 20 --minQ 30 --max-missing 0.80 --remove-indels --max-alleles 2 --recode --recode-INFO-all --stdout | gzip -c > "${output_dir}/all_samples_BamQualFilt-f3Q30F264_nodup_raw_OnlyPolymorphic_shortNam_DP10GQ20Q30_Mis80NoIndel.vcf.gz"

# Count number of variants in raw vcf
zgrep -v "^#" all_samples-7outgroup_1-8ChrOnly_BamQualFilt-f3Q30F264_nodup_realigned_sorted_raw_masked.vcf | wc –l