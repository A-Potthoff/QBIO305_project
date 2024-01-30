


# Gene Enrichment Analysis 
##################################################################################

########### Input file ##############
#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --partition=devel-rh7
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40gb
#SBATCH --time=1:00:00
#SBATCH --account=UniKoeln
#SBATCH --error=/home/group.kurse/qbcbp015/qbio_sample_project/stacks/stacks_pop-%x-%j.err
#SBATCH --output=/home/group.kurse/qbcbp015/qbio_sample_project/stacks/stacks_pop-%x-%j.out


#module load stacks/2.65

#input_dir="/home/group.kurse/qbcbp015/qbio_sample_project/"

#output_dir="/home/group.kurse/qbcbp015/qbio_sample_project/stacks/"

#populations -V "${input_dir}group_4_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.vcf.gz" --popmap "/home/group.kurse/qbcbp015/qbio_sample_project/stacks//pop.txt" --fstats -t 10 --structure -O "${output_dir}"

##################################################################################
# Load required libraries
library(readr)
library(GenomicRanges)
install.packages("GenomicRanges")
if (!require("GenomicRanges"))
  install.packages("GenomicRanges")
BiocManager::install("GenomicRanges")
library(dplyr)



setwd("/Users/sarakreci/Desktop")
#####################
# read --fstats output file from stacks

#1-SWE-ESP
#2-SWE-IK
#3-UK-ESP

#####################
fstats <- read.table("group_4_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.p.fst_SWE-ESP.tsv", header=FALSE)

fstats <- read.table("group_4_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.p.fst_SWE-UK.tsv", header=FALSE)

fstats <- read.table("group_4_final_accession_1001genomes_snp-short-indel_only_ACGTN_Dp10GQ20Q30_NoIndel_Bialleleic_80PcMissing.p.fst_UK-ESP.tsv", header=FALSE)

colnames(fstats) <- gsub(" ", "_", c(
  "Locus ID", "Pop 1 ID", "Pop 2 ID", "Chr", "BP", "Column", "Fishers P",
  "Odds Ratio", "CI Low", "CI High", "LOD", "AMOVA Fst", "Smoothed AMOVA Fst",
  "Smoothed AMOVA Fst P-value"
))

names(fstats)
fstats_fst <- data.frame(
  Chr = fstats$Chr,
  BP = fstats$BP,
  Fishers_P = fstats$Fishers_P,
  AMOVA_Fst = fstats$AMOVA_Fst
)

names(fstats_fst)

head(fstats_fst)
print(head(fstats_fst, n=1))

#############
#Read annotation file
#############
annot<-read.delim("group_4_At_heat_cold_dehydration_stress_only.gff", header = FALSE)
head(annot)
names(annot)

print(head(annot, n=1))

colnames(annot) <- c(
  "gene_id", "chr", "tair_version", "type", "start", "end", "empty1",
  "strand", "empty2", "gene_id2", "gene_id3", "type_name", "short_description", "curation"
)

print(head(annot, n=5))

annot_subset <- data.frame(
  gene_id = annot$gene_id,
  chr = annot$chr,
  start = annot$start,
  end = annot$end
)  

print(head(annot_subset, n=5))
annot_subset$chr <- sub("Chr", "", annot_subset$chr)

gr_fst <- GRanges(seqnames = fstats_fst$Chr, 
                  ranges = IRanges(start = fstats_fst$BP, end = fstats_fst$BP),
                  fst = fstats_fst$AMOVA_Fst,
                  p_value = fstats_fst$Fishers_P)

length(fstats_fst$Chr)
print(head(gr_fst, n=1))

gr_annotation <- GRanges(seqnames = annot_subset$chr, 
                         ranges = IRanges(start = annot_subset$start, end = annot_subset$end),
                         gene_id = annot_subset$gene_id)

ovl <- findOverlaps(gr_fst, gr_annotation, type = "any", select = "all", ignore.strand = TRUE)

overlapping_genes <- as.data.frame(gr_annotation[subjectHits(ovl)])
overlapping_fst <- as.data.frame(gr_fst[queryHits(ovl)])

print(head(overlapping_genes,n=5))
print(head(overlapping_fst,n=5))

colnames(overlapping_genes)[colnames(overlapping_genes) == "seqnames"] <- "chrom"
colnames(overlapping_fst)[colnames(overlapping_fst) == "seqnames"] <- "chrom"

overlapping_genes <- overlapping_genes[, !(names(overlapping_genes) %in% "strand")]
overlapping_fst <- overlapping_fst[, !(names(overlapping_fst) %in% "strand")]

print(head(overlapping_genes,n=5))
print(head(overlapping_fst,n=5))

write.table(overlapping_genes, file = "overlapping_genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(overlapping_fst, file = "overlapping_fst.txt", sep = "\t", row.names = FALSE, quote = FALSE)  


####################
# Running in the terminal
###################
$ module load bedtools/
  $ bedtools intersect -wa -wb -a overlapping_genes.txt -b overlapping_fst.txt > overlapping_genes_fst.txt

##################
# Process "overlapping_genes_fst.txt"
# in R
#################

#Read annotation file

ovl_genes_fst<-read.delim("overlapping_genes_fst.txt", header = FALSE)
head(ovl_genes_fst)
names(ovl_genes_fst)

print(head(ovl_genes_fst, n=1))

colnames(ovl_genes_fst) <- c(
  "chrom", "gene_start", "gene_end", "gene_width", "gene_id","chrom", "snp_start", "snp_end", "snp_width",
  "fst", "p_value")

print(head(ovl_genes_fst, n=5))

ovl_genes_fst_sub <- data.frame(
  gene_id = ovl_genes_fst$gene_id,
  fst = ovl_genes_fst$fst,
  p_value = ovl_genes_fst$p_value
)  

print(head(ovl_genes_fst_sub, n=5))

ovl_genes_fst_sub$gene_id <- sub("\\.\\d+", "", ovl_genes_fst_sub$gene_id)

sorted_data <- ovl_genes_fst_sub %>%
  arrange(gene_id, desc(fst), p_value)

unique_data <- sorted_data %>%
  group_by(gene_id) %>%
  filter(row_number() == 1)

print(head(unique_data, n=10))

########################

# Writing of  both overlapping_genes and overlapping_fst as text files

write.table(unique_data, file = "GO_input.txt", sep = "\t", row.names = FALSE, quote = FALSE)

