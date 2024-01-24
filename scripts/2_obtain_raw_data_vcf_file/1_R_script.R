###### Rscipt to Annotate GFF file ####
# Input files below are downloaded from: https://www.arabidopsis.org/
# TAIR10_functional_descriptions.txt
# TAIR10_GFF3_genes.gff
#######################################
#set working directory
setwd("C:/Data_tali/QBio/2023-24/AT_1001_genome/Annotate_GFF")

###
#Read GFF file with SNP positions
###

GFF<-read.delim("TAIR10_GFF3_genes.gff", header = FALSE)
head(GFF)
#Create a duplicated column 9
GFF$dup_colmn<-GFF$V9

#select all the rows with mRNA in cloumn 3 GFF_sub<-subset(GFF, GFF$V3 == "mRNA")
head(GFF_sub)
GFF_sub$dup_colmn<-gsub("ID=.*Name=","", GFF_sub$dup_colmn) GFF_sub$dup_colmn<-gsub(";Index=1","", GFF_sub$dup_colmn)

write.table(GFF_sub, "GFF_final_dup_column.gff", sep = "\t", quote = FALSE)

#change the name of column from duplicated to Model_name GFF_sub$Model_name<-GFF_sub$dup_colmn

###
#Read annotation file
###
annot<-read.delim("TAIR10_functional_descriptions.txt", header = TRUE)
head(annot)
#Merge subset of GFF with annotation file using column Model_name (common in both data files) GFF_sub_annot.m<-merge(GFF_sub, annot, by="Model_name")
head(GFF_sub_annot.m)
write.table(GFF_sub_annot.m, "GFF_final_with_annotation.gff", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


