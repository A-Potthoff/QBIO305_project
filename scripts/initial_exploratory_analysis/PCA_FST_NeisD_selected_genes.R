########         Visualizing VCF data          ########
#######################################################

library(vcfR)
library(gplots)
library(tidyverse)
library(adegenet)
library(factoextra)
library(FactoMineR)
library(StAMPP)
library(RColorBrewer)
library(ggrepel)
library(htmlwidgets)
library(heatmaply)


test <- read.table("GCF_000001735.4_TAIR10.1_genomic.fna")


setwd("C:/Users/andre/OneDrive/Bildung/3_HHU/quant_Bio/3rd/QBio305_Population_And_Quantitative_Genetics/Cologne/exercises/project/own_research/analysis/data")
getwd()

# Input the files ----

accessions <- read.csv("selected_accessions.csv")
small_vcf <- read.vcfR("final_vcf.gz", verbose = FALSE) #vcf file containing variants calls
#vcf <- read.vcfR("vcf.gz", verbose = FALSE) #vcf file containing variants calls
dna <- ape::read.dna("GCF_000001735.4_TAIR10.1_Chr3.fna", format = "fasta") #fasta sequence of reference genome
gff <- read.table("TAIR10_GFF3_genes_Chr3.gff", sep="\t", quote="") #Annotation file of reference genome


# PCA by using VCF file----

### preparation

# Convert VCF to genind object
small_genind_vcf <- vcfR2genind(small_vcf)

# Scale genind object for PCA
small_genind_vcf_scaled = scaleGen(small_genind_vcf, NA.method = "mean")

# Perform PCA
small_pca <- dudi.pca(small_genind_vcf_scaled, cent = TRUE, scale = FALSE, scannf = FALSE, nf = 10)

# Check PCA dimensions
small_axis_all = small_pca$eig / sum(small_pca$eig)
barplot(small_axis_all[1:10], main = "PCA eigenvalues")

#ind <- pca$rank
small_pca_axis <- small_pca$li
small_pca_eigenval <- small_pca$eig[1:10]
str(small_pca_axis)

# set names
small_ind_names<- rownames(small_genind_vcf_scaled)
small_pca_axis$ind <- small_ind_names

# Add a new column named "Population" to your data frame
population_labels <- accessions[match(small_ind_names, accessions$X1001.ID), "country"]

small_pca_axis$Population <- population_labels
pop <- population_labels

# remake data.frame
small_pca_2 <- as_tibble(data.frame(small_pca_axis, population_labels))

n <- length(small_pca_eigenval)
n # use this number PC=1:n

# first convert to percentage variance explained
small_pve <- data.frame(PC = 1:n, pve = small_pca_eigenval/sum(small_pca_eigenval)*100)

# make plot
a <- ggplot(small_pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(small_pve$pve)

### plotting

# plot pca PC1 and PC2
a <- ggplot(small_pca_2, aes(Axis1, Axis2, col = pop)) + 
  geom_point(size = 3) + 
  scale_colour_manual(values = c("red", "blue", "green")) +
  coord_equal() + 
  theme_light() +
  xlab(paste0("PC1 (", signif(small_pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(small_pve$pve[2], 3), "%)"))
a

# plot pca PC1 and PC3
b <- ggplot(small_pca_2, aes(Axis1, Axis3, col = pop)) + 
  geom_point(size = 3) + 
  scale_colour_manual(values = c("red", "blue", "green")) +
  coord_equal() + 
  theme_light() +
  xlab(paste0("PC1 (", signif(small_pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(small_pve$pve[3], 3), "%)"))
b

# plot pca PC2 and PC3
c <- ggplot(small_pca_2, aes(Axis2, Axis3, col = pop)) + 
  geom_point(size = 3) + 
  scale_colour_manual(values = c("red", "blue", "green")) +
  coord_equal() + 
  theme_light() +
  xlab(paste0("PC2 (", signif(small_pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(small_pve$pve[3], 3), "%)"))
c

ggsave("small_PCA12.png", a, width = 1800, height = 1200, units = "px")
ggsave("small_PCA13.png", b, width = 1800, height = 1200, units = "px")
ggsave("small_PCA23.png", c, width = 1800, height = 1200, units = "px")


# Create a labeled plot
x <- ggplot(small_pca_2, aes(Axis1, Axis2, col = pop, label = ind)) + 
  geom_point(size = 2) + # size of points
  geom_text_repel(show.legend = FALSE, size=2) +  # Use geom_text_repel() for label repulsion
  scale_colour_manual(values = c("red", "blue", "green")) + 
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(small_pve$pve[1], 2), "%)")) +
  ylab(paste0("PC2 (", signif(small_pve$pve[2], 2), "%)"))
x

### 3D plots

#PCA Visualization using plotly
#https://plotly.com/r/pca-visualization/
library(plotly)
#3D plots
fig <- plot_ly(pca_2, x = ~Axis1, y = ~Axis2, z = ~Axis3, color =pop, colors = c("red", "blue", "green") ) %>%
  add_markers(size = 12)

fig <- fig %>%
  layout(
    title = "Arabidopsis thaliana Eurasian Accessions",
    scene = list(bgcolor = "#e5ecf6")
  )
fig

saveWidget(fig, "3D_PCA.html")


###
#Visualize a subset of the principal components
cumsum(pve$pve)

tit = 'Total Explained Variance = 58.85 %'

axis = list(showline=FALSE,
            zeroline=FALSE,
            gridcolor='#ffff',
            ticklen=4)

fig <- pca_2 %>%
  plot_ly() %>%
  add_trace(
    type = 'splom',
    dimensions = list(
      list(label='PC1 (22.9%)', values=~Axis1),
      list(label='PC2 (14.3%)', values=~Axis2),
      list(label='PC3 (11.1%)', values=~Axis3),
      list(label='PC4 (10.5%)', values=~Axis4)
    ),
    color = ~pop, colors = c("red", "blue", "green"),
    marker = list(
      size = 7
    )
  ) %>% style(diagonal = list(visible = F)) %>%
  layout(
    title= tit,
    hovermode='closest',
    dragmode= 'select',
    plot_bgcolor='rgba(240,240,240, 0.95)',
    xaxis=list(domain=NULL, showline=F, zeroline=F, gridcolor='#ffff', ticklen=4),
    yaxis=list(domain=NULL, showline=F, zeroline=F, gridcolor='#ffff', ticklen=4),
    xaxis2=axis,
    xaxis3=axis,
    xaxis4=axis,
    yaxis2=axis,
    yaxis3=axis,
    yaxis4=axis
  )
options(warn=-1)
fig

# plot relatedness matrix ----

#relat <- as.matrix(read.table("vcf_distance_plink_ibs.mdist"))

#heatmap.2(relat, trace="none", Rowv=T, Colv=NA, cexRow=0.6,cexCol = 0.6, labRow = ind_names, labCol = ind, col= colorRampPalette(c("lemonchiffon", "lemonchiffon", "lemonchiffon", "yellow", "red"))(70))


# Stampp to calculate FST between populations----

# Convert VCF to genlight object
small_genlight_vcf <- vcfR2genlight(small_vcf)

# Extract population data
pop = as.factor(population_labels)
small_genlight_vcf$pop = pop

# Convert genlight to stampp object
small_stampp_vcf = stamppConvert(small_genlight_vcf, type = "genlight")

# Calculate FST between populations
small_stamppFst = stamppFst(small_stampp_vcf, nboots = 100, percent = 95, nclusters = 8)
small_stamppFst_matrix = as.matrix(small_stamppFst$Fsts)

# Set diagonal to 0 and add symmetrical upper tri part of matrix
diag(small_stamppFst_matrix) <- 0
small_stamppFst_matrix[upper.tri(small_stamppFst_matrix)]  <- t(small_stamppFst_matrix)[upper.tri(small_stamppFst_matrix)]
# Optional: order the names
small_stamppFst_matrix = small_stamppFst_matrix[order(row.names(small_stamppFst_matrix)), order(colnames(small_stamppFst_matrix))]


heatmaply(
  small_stamppFst_matrix,
  col = plasma(100),  # You can adjust the number of colors in the palette
  main = "FST values between populations.",
  cellnote = small_stamppFst_matrix,  # Display cell values
  notecol = "black",  # Color of the cell values
  margins = c(60, 70))

#####################################################################
########## calculate genetic distance between individuals ###########
#####################################################################

small_stamppNeisD = stamppNeisD(small_stampp_vcf, pop = FALSE)
small_stamppNeisD_matrix = as.matrix(small_stamppNeisD)

# Set diagonal to 0 and add symmetrical upper tri part of matrix
diag(small_stamppNeisD_matrix) <- 0
small_stamppNeisD_matrix[upper.tri(small_stamppNeisD_matrix)]  <- t(small_stamppNeisD_matrix)[upper.tri(small_stamppNeisD_matrix)]
heatmap(small_stamppNeisD_matrix)
colnames(small_stamppNeisD_matrix) <- rownames(small_stamppNeisD_matrix)
max(small_stamppNeisD_matrix)

heatmap(small_stamppNeisD_matrix, RowSideColors = cc, ColSideColors = cc)


