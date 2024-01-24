#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=15gb
#SBATCH --time=02:15:00
#SBATCH --account=UniKoeln
#SBATCH --output=/home/group.kurse/qbcbp008/project/SFSs/scripts/mapping.out
#SBATCH --error=/home/group.kurse/qbcbp008/project/SFSs/scripts/mapping.err

#load modules
module load bwamem2/2.2.1
module load samtools/1.13

#directory from which the code will run
cd /home/group.kurse/qbcbp008/project/SFSs

#create necessary folders
#mkdir mapping

#index the reference genome
#bwa-mem2 index ./reference/TAIR10.1.genomic_fna

#define reference path and samples
reference=./reference/TAIR10.1.genomic_fna
sample="SRR2082718"

#Set input and output file paths
trimmed_fastq1=./fastq/SRR2082718_1_trim.fastq.gz
trimmed_fastq2=./fastq/SRR2082718_2_trim.fastq.gz
sam=./mapping/${sample}.sam
bam=./mapping/${sample}.bam
sorted_bam=./mapping/${sample}.sorted.bam
mapped_qc=./mapping/${sample}.qc
outfile=${sample}.qualimap

echo  allign the reads to reference genome 
bwa-mem2 mem -M -t 20 $reference $trimmed_fastq1 $trimmed_fastq2 > $sam

echo converting to binary  
samtools view -bS $sam > $bam

echo sorting bam 
samtools sort -@ 20 -o $sorted_bam $bam

echo indexing the sorted BAM file 
samtools index $sorted_bam

echo deleting unsorted BAM and SAM files 
#rm $bam
#rm $sam

echo qualimaping
qualimap bamqc -bam $sorted_bam -outdir ./mapping/ -outformat html -outfile $outfile

echo done