#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=10gb
#SBATCH --time=01:00:00
#SBATCH --account=UniKoeln
#SBATCH --output=/home/group.kurse/qbcbp008/project/SFSs/scripts/fastp.out
#SBATCH --error=/home/group.kurse/qbcbp008/project/SFSs/scripts/fastp.err

fastp --in1 ./SRR2082718_1.fastq.gz --in2 ./SRR2082718_2.fastq.gz --out1 ./SRR2082718_1_trim.fastq.gz --out2 ./SRR2082718_2_trim.fastq.gz -l 50 -h report_html &> log_file
