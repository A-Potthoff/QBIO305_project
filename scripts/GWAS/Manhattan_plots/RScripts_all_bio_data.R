#########################
#      Bio1             #
#      Rscript          #
#       GWAS            #
#    Manhatten plot     #
#########################

# Load the library
library(qqman)

# Read GWAS output file created using gemma
gwas_data_bio1 <- read.table("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/GWAS/GWAS/output_bio1/sample_project_GWAS2.assoc.txt", header = TRUE)

#adding rs order for annotation in manhattan plot
gwas_data_bio1$rs <- seq(1, nrow(gwas_data_bio1))

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwas_data_bio1, chr="chr", bp="ps", snp="rs", p="p_wald" )
# modify colour of points
manhattan(x = gwas_data_bio1, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue4", "orange3"), logp = TRUE)

# Multiple testing correction in GWAS (Bonferroni correction) 
pval_bonf = 0.05/dim(gwas_data_bio1)[[1]]
#adding multiple-testing-corrected p-value line to the plot
manhattan(gwas_data_bio1, chr="chr", bp="ps", snp="rs", p="p_wald", suggestiveline = -log10(pval_bonf), genomewideline = FALSE, annotatePval = -log10(pval_bonf), col = c("blue4", "orange3"))


#Manhattan plot using FDR-corrected p-values
pvals_fdr <- p.adjust(gwas_data_bio1$p_wald, method = "BH")

# Add FDR-corrected p-value threshold line to the plot
pval_fdr <- 0.05  # Desired FDR threshold

# Add the adjusted p-values as a new column to the data frame
gwas_data_bio1$pvals_fdr <- pvals_fdr

# Create the Manhattan plot with the new adjusted p-values column
manhattan(gwas_data_bio1, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(pval_fdr), genomewideline =  FALSE,
          annotatePval = -log10(pvals_fdr), col = c("blue4", "orange3"),
          main = "Bio1 (Annual Mean Temperature)")

#########################
#       Bio5            #
#      Rscript          #
#       GWAS            #
#    Manhatten plot     #
#########################

# Load the library
library(qqman)

# Read GWAS output file created using gemma
gwas_data_bio5 <- read.table("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/GWAS/GWAS/output_bio5/sample_project_GWAS2.assoc.txt", header = TRUE)

#adding rs order for annotation in manhattan plot
gwas_data_bio5$rs <- seq(1, nrow(gwas_data_bio5))

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwas_data_bio5, chr="chr", bp="ps", snp="rs", p="p_wald" )
# modify colour of points
manhattan(x = gwas_data_bio5, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue4", "orange3"), logp = TRUE)

# Multiple testing correction in GWAS (Bonferroni correction) 
pval_bonf = 0.05/dim(gwas_data_bio5)[[1]]
#adding multiple-testing-corrected p-value line to the plot
manhattan(gwas_data_bio5, chr="chr", bp="ps", snp="rs", p="p_wald", suggestiveline = -log10(pval_bonf), genomewideline = FALSE, annotatePval = -log10(pval_bonf), col = c("blue4", "orange3"))


#Manhattan plot using FDR-corrected p-values
pvals_fdr <- p.adjust(gwas_data_bio5$p_wald, method = "BH")

# Add FDR-corrected p-value threshold line to the plot
pval_fdr <- 0.05  # Desired FDR threshold

# Add the adjusted p-values as a new column to the data frame
gwas_data_bio5$pvals_fdr <- pvals_fdr

# Create the Manhattan plot with the new adjusted p-values column
manhattan(gwas_data_bio5, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(pval_fdr), genomewideline =  FALSE,
          annotatePval = -log10(pval_fdr), col = c("blue4", "orange3"),
          main = "Bio5 (Max Temp of Warmest Month)")


#########################
#       Bio9            #
#      Rscript          #
#       GWAS            #
#    Manhatten plot     #
#########################

# Load the library
library(qqman)

# Read GWAS output file created using gemma
gwas_data_bio9 <- read.table("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/GWAS/GWAS/output_bio9/sample_project_GWAS2.assoc.txt", header = TRUE)

#adding rs order for annotation in manhattan plot
gwas_data_bio9$rs <- seq(1, nrow(gwas_data_bio9))

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwas_data_bio9, chr="chr", bp="ps", snp="rs", p="p_wald" )
# modify colour of points
manhattan(x = gwas_data_bio9, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue4", "orange3"), logp = TRUE)

# Multiple testing correction in GWAS (Bonferroni correction) 
pval_bonf = 0.05/dim(gwas_data_bio9)[[1]]
#adding multiple-testing-corrected p-value line to the plot
manhattan(gwas_data_bio9, chr="chr", bp="ps", snp="rs", p="p_wald", suggestiveline = -log10(pval_bonf), genomewideline = FALSE, annotatePval = -log10(pval_bonf), col = c("blue4", "orange3"))


#Manhattan plot using FDR-corrected p-values
pvals_fdr <- p.adjust(gwas_data_bio9$p_wald, method = "BH")

# Add FDR-corrected p-value threshold line to the plot
pval_fdr <- 0.05  # Desired FDR threshold

# Add the adjusted p-values as a new column to the data frame
gwas_data_bio9$pvals_fdr <- pvals_fdr

# Create the Manhattan plot with the new adjusted p-values column
manhattan(gwas_data_bio9, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(pval_fdr), genomewideline =  FALSE,
          annotatePval = -log10(pvals_fdr), col = c("blue4", "orange3"),
          main = "Bio9 (Mean Temp of Driest Quarter)")



#########################
#       Bio10           #      
#      Rscript          #
#       GWAS            #
#    Manhatten plot     #
#########################

# Load the library
library(qqman)

# Read GWAS output file created using gemma
gwas_data_bio10 <- read.table("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/GWAS/GWAS/output_bio10/sample_project_GWAS2.assoc.txt", header = TRUE)

#adding rs order for annotation in manhattan plot
gwas_data_bio10$rs <- seq(1, nrow(gwas_data_bio10))

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwas_data_bio10, chr="chr", bp="ps", snp="rs", p="p_wald" )
# modify colour of points
manhattan(x = gwas_data_bio10, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue4", "orange3"), logp = TRUE)

# Multiple testing correction in GWAS (Bonferroni correction) 
pval_bonf = 0.05/dim(gwas_data_bio10)[[1]]
#adding multiple-testing-corrected p-value line to the plot
manhattan(gwas_data_bio10, chr="chr", bp="ps", snp="rs", p="p_wald", suggestiveline = -log10(pval_bonf), genomewideline = FALSE, annotatePval = -log10(pval_bonf), col = c("blue4", "orange3"))


#Manhattan plot using FDR-corrected p-values
pvals_fdr <- p.adjust(gwas_data_bio10$p_wald, method = "BH")

# Add FDR-corrected p-value threshold line to the plot
pval_fdr <- 0.05  # Desired FDR threshold

# Add the adjusted p-values as a new column to the data frame
gwas_data_bio10$pvals_fdr <- pvals_fdr

# Create the Manhattan plot with the new adjusted p-values column
manhattan(gwas_data_bio10, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(0.05), genomewideline =  FALSE,
          annotatePval = -log10(pvals_fdr), col = c("blue4", "orange3"),
          main = "Bio10 (Mean Temp of Warmest Quarter)")

#########################
#       Bio14           #
#      Rscript          #
#       GWAS            #
#    Manhatten plot     #
#########################

# Load the library
library(qqman)

# Read GWAS output file created using gemma
gwas_data_bio14 <- read.table("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/GWAS/GWAS/output_bio14/sample_project_GWAS2.assoc.txt", header = TRUE)

#adding rs order for annotation in manhattan plot
gwas_data_bio14$rs <- seq(1, nrow(gwas_data_bio14))

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwas_data_bio14, chr="chr", bp="ps", snp="rs", p="p_wald" )
# modify colour of points
manhattan(x = gwas_data_bio14, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue4", "orange3"), logp = TRUE)

# Multiple testing correction in GWAS (Bonferroni correction) 
pval_bonf = 0.05/dim(gwas_data_bio14)[[1]]
#adding multiple-testing-corrected p-value line to the plot
manhattan(gwas_data_bio14, chr="chr", bp="ps", snp="rs", p="p_wald", suggestiveline = -log10(pval_bonf), genomewideline = FALSE, annotatePval = -log10(pval_bonf), col = c("blue4", "orange3"))


#Manhattan plot using FDR-corrected p-values
pvals_fdr <- p.adjust(gwas_data_bio14$p_wald, method = "BH")

# Add FDR-corrected p-value threshold line to the plot
pval_fdr <- 0.05  # Desired FDR threshold

# Add the adjusted p-values as a new column to the data frame
gwas_data_bio14$pvals_fdr <- pvals_fdr

# Create the Manhattan plot with the new adjusted p-values column
manhattan(gwas_data_bio14, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(pval_fdr), genomewideline =  FALSE,
          annotatePval = -log10(pvals_fdr), col = c("blue4", "orange3"),
          main = "Bio14 (Precipitation of Driest Month)")


#####Extracting gwas_data_bio14####

# Install and load the required packages
#install.packages(c("gridExtra", "ggplot2"))
library(gridExtra)
library(ggplot2)

# Sort the data frame by the "pvals_fdr" column
sorted_data <- gwas_data_bio14 %>% arrange(pvals_fdr)

# Select the top five rows
top_six_rows <- head(sorted_data, 6)

# Create a table as a grob (grid object)
table_grob <- tableGrob(top_six_rows, rows = NULL)

# Save the table as an image (e.g., PNG)
ggsave("output_table_gwas_data_bio14.png", table_grob, width = 15, height = 5, units = "in")

#########################
#       Bio17           #
#      Rscript          #
#       GWAS            #
#    Manhatten plot     #
#########################

# Load the library
library(qqman)

# Read GWAS output file created using gemma
gwas_data_bio17 <- read.table("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/GWAS/GWAS/output_bio17/sample_project_GWAS2.assoc.txt", header = TRUE)

#adding rs order for annotation in manhattan plot
gwas_data_bio17$rs <- seq(1, nrow(gwas_data_bio17))

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwas_data_bio17, chr="chr", bp="ps", snp="rs", p="p_wald" )
# modify colour of points
manhattan(x = gwas_data_bio17, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue4", "orange3"), logp = TRUE)

# Multiple testing correction in GWAS (Bonferroni correction) 
pval_bonf = 0.05/dim(gwas_data_bio17)[[1]]
#adding multiple-testing-corrected p-value line to the plot
manhattan(gwas_data_bio17, chr="chr", bp="ps", snp="rs", p="p_wald", suggestiveline = -log10(pval_bonf), genomewideline = FALSE, annotatePval = -log10(pval_bonf), col = c("blue4", "orange3"))


#Manhattan plot using FDR-corrected p-values
pvals_fdr <- p.adjust(gwas_data_bio17$p_wald, method = "BH")

# Add FDR-corrected p-value threshold line to the plot
pval_fdr <- 0.05  # Desired FDR threshold

# Add the adjusted p-values as a new column to the data frame
gwas_data_bio17$pvals_fdr <- pvals_fdr

# Create the Manhattan plot with the new adjusted p-values column
manhattan(gwas_data_bio17, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(pval_fdr), genomewideline =  FALSE,
          annotatePval = -log10(pvals_fdr), col = c("blue4", "orange3"),
          main = "Bio17 (Precipitation of Driest Quarter)")



#####Extracting gwas_data_bio17####

#install.packages(c("gridExtra", "ggplot2"))
library(gridExtra)
library(ggplot2)

# Sort the data frame by the "pvals_fdr" column
sorted_data <- gwas_data_bio17 %>% arrange(pvals_fdr)

# Select the top five rows
top_five_rows <- head(sorted_data, 5)

# Create a table as a grob (grid object)
table_grob <- tableGrob(top_five_rows, rows = NULL)

# Save the table as an image (e.g., PNG)
ggsave("output_table_gwas_data_bio17.png", table_grob, width = 15, height = 5, units = "in")

#########################
#       Bio18           #
#      Rscript          #
#       GWAS            #
#    Manhatten plot     #
#########################

# Load the library
library(qqman)

# Read GWAS output file created using gemma
gwas_data_bio18 <- read.table("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/GWAS/GWAS/output_bio18/sample_project_GWAS2.assoc.txt", header = TRUE)

#adding rs order for annotation in manhattan plot
gwas_data_bio18$rs <- seq(1, nrow(gwas_data_bio18))

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwas_data_bio18, chr="chr", bp="ps", snp="rs", p="p_wald" )
# modify colour of points
manhattan(x = gwas_data_bio18, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue4", "orange3"), logp = TRUE)

# Multiple testing correction in GWAS (Bonferroni correction) 
pval_bonf = 0.05/dim(gwas_data_bio18)[[1]]
#adding multiple-testing-corrected p-value line to the plot
manhattan(gwas_data_bio18, chr="chr", bp="ps", snp="rs", p="p_wald", suggestiveline = -log10(pval_bonf), genomewideline = FALSE, annotatePval = -log10(pval_bonf), col = c("blue4", "orange3"))


#Manhattan plot using FDR-corrected p-values
pvals_fdr <- p.adjust(gwas_data_bio18$p_wald, method = "BH")

# Add FDR-corrected p-value threshold line to the plot
pval_fdr <- 0.05  # Desired FDR threshold

# Add the adjusted p-values as a new column to the data frame
gwas_data_bio18$pvals_fdr <- pvals_fdr

# Create the Manhattan plot with the new adjusted p-values column
manhattan(gwas_data_bio18, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(0.05), genomewideline =  FALSE,
          annotatePval = -log10(pvals_fdr), col = c("blue4", "orange3"),
          main = "Bio18 (Precipitation of Warmest Quarter)")


#####Comparing Bio14 Bio17 and Bio18#######

# Set up a 1x2 layout for two plots side by side
par(mfrow = c(1, 3))

# Bio14 Manhattan plot
manhattan(gwas_data_bio14, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(pval_fdr), genomewideline =  FALSE,
          annotatePval = -log10(pvals_fdr), col = c("blue4", "orange3"),
          main = "Bio14")


# Bio17 Manhattan plot
manhattan(gwas_data_bio17, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(pval_fdr), genomewideline =  FALSE,
          annotatePval = -log10(pvals_fdr), col = c("blue4", "orange3"),
          main = "Bio17")

# Bio18 Manhattan plot
manhattan(gwas_data_bio18, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(pval_fdr), genomewideline =  FALSE,
          annotatePval = -log10(pvals_fdr), col = c("blue4", "orange3"),
          main = "Bio18")

# Reset the layout to default
par(mfrow = c(1, 1))

#####Extracting gwas_data_bio18####

#install.packages(c("gridExtra", "ggplot2"))
library(gridExtra)
library(ggplot2)

# Sort the data frame by the "pvals_fdr" column
sorted_data <- gwas_data_bio18 %>% arrange(pvals_fdr)

# Select the top five rows
top_five_rows <- head(sorted_data, 5)

# Create a table as a grob (grid object)
table_grob <- tableGrob(top_five_rows, rows = NULL)

# Save the table as an image (e.g., PNG)
ggsave("output_table_gwas_data_bio18.png", table_grob, width = 15, height = 5, units = "in")



#########################
#       PET_Jan         #
#      Rscript          #
#       GWAS            #
#    Manhatten plot     #
#########################

# Load the library
library(qqman)

# Read GWAS output file created using gemma
gwas_data_pet_jan <- read.table("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/output_pet_jan/sample_project_GWAS2.assoc.txt", header = TRUE)

#adding rs order for annotation in manhattan plot
gwas_data_pet_jan$rs <- seq(1, nrow(gwas_data_pet_jan))

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwas_data_pet_jan, chr="chr", bp="ps", snp="rs", p="p_wald" )
# modify colour of points
manhattan(x = gwas_data_pet_jan, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue4", "orange3"), logp = TRUE)

# Multiple testing correction in GWAS (Bonferroni correction) 
pval_bonf = 0.05/dim(gwas_data_pet_jan)[[1]]
#adding multiple-testing-corrected p-value line to the plot
manhattan(gwas_data_pet_jan, chr="chr", bp="ps", snp="rs", p="p_wald", suggestiveline = -log10(pval_bonf), genomewideline = FALSE, annotatePval = -log10(pval_bonf), col = c("blue4", "orange3"))


#Manhattan plot using FDR-corrected p-values
pvals_fdr <- p.adjust(gwas_data_pet_jan$p_wald, method = "BH")

# Add FDR-corrected p-value threshold line to the plot
pval_fdr <- 0.05  # Desired FDR threshold

# Add the adjusted p-values as a new column to the data frame
gwas_data_pet_jan$pvals_fdr <- pvals_fdr

# Create the Manhattan plot with the new adjusted p-values column
manhattan(gwas_data_pet_jan, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(pval_fdr), genomewideline =  FALSE,
          annotatePval = -log10(pvals_fdr), col = c("blue4", "orange3"),
          main = "PET January")



#########################
#       PET_June        #
#      Rscript          #
#       GWAS            #
#    Manhatten plot     #
#########################

# Load the library
library(qqman)

# Read GWAS output file created using gemma
gwas_data_pet_june <- read.table("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/output_pet_june/sample_project_GWAS2.assoc.txt", header = TRUE)

#adding rs order for annotation in manhattan plot
gwas_data_pet_june$rs <- seq(1, nrow(gwas_data_pet_june))

# Make the Manhattan plot on the gwasResults dataset
manhattan(gwas_data_pet_june, chr="chr", bp="ps", snp="rs", p="p_wald" )
# modify colour of points
manhattan(x = gwas_data_pet_june, chr = "chr", bp = "ps", p = "p_wald", snp = "rs", col = c("blue4", "orange3"), logp = TRUE)

# Multiple testing correction in GWAS (Bonferroni correction) 
pval_bonf = 0.05/dim(gwas_data_pet_june)[[1]]
#adding multiple-testing-corrected p-value line to the plot
manhattan(gwas_data_pet_june, chr="chr", bp="ps", snp="rs", p="p_wald", suggestiveline = -log10(pval_bonf), genomewideline = FALSE, annotatePval = -log10(pval_bonf), col = c("blue4", "orange3"))


#Manhattan plot using FDR-corrected p-values
pvals_fdr <- p.adjust(gwas_data_pet_june$p_wald, method = "BH")

# Add FDR-corrected p-value threshold line to the plot
pval_fdr <- 0.05  # Desired FDR threshold

# Add the adjusted p-values as a new column to the data frame
gwas_data_pet_june$pvals_fdr <- pvals_fdr

# Create the Manhattan plot with the new adjusted p-values column
manhattan(gwas_data_pet_june, chr = "chr", bp = "ps", snp = "rs", p =  "pvals_fdr", 
          suggestiveline = -log10(0.05), genomewideline =  FALSE,
          annotatePval = -log10(pvals_fdr), col = c("blue4", "orange3"),
          main = "PET June")
