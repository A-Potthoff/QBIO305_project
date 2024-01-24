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
library(plotly)
library(heatmaply)

setwd("C:/Users/andre/OneDrive/Bildung/3_HHU/quant_Bio/3rd/QBio305_Population_And_Quantitative_Genetics/Cologne/exercises/project/own_research/analysis/data")
getwd()

# Input the files ----

accessions <- read.csv("selected_accessions.csv")
#vcf <- read.vcfR("final_vcf.gz", verbose = FALSE) #vcf file containing variants calls
vcf <- read.vcfR("vcf.gz", verbose = FALSE) #vcf file containing variants calls
#dna <- ape::read.dna("???????????????????") #fasta sequence of reference genome
gff <- read.table("group_4_At_heat_cold_dehydration_stress_only.gff", sep="\t", quote="") #Annotation file of reference genome

# PCA by using VCF file----

### preparation

# Convert VCF to genind object
genind_vcf <- vcfR2genind(vcf)

# Scale genind object for PCA
genind_vcf_scaled = scaleGen(genind_vcf, NA.method = "mean")

# Perform PCA
pca <- dudi.pca(genind_vcf_scaled, cent = TRUE, scale = FALSE, scannf = FALSE, nf = 10)

# Check PCA dimensions
axis_all = pca$eig / sum(pca$eig)
barplot(axis_all[1:10], main = "PCA eigenvalues")

#ind <- pca$rank
pca_axis <- pca$li
pca_eigenval <- pca$eig[1:10]
str(pca_axis)

# set names
ind_names<- rownames(genind_vcf_scaled)
pca_axis$ind <- ind_names

# Add a new column named "Population" to your data frame
population_labels <- accessions[match(ind_names, accessions$X1001.ID), "country"]

pca_axis$Population <- population_labels
pop <- population_labels

# remake data.frame
pca_2 <- as_tibble(data.frame(pca_axis, population_labels))

n <- length(pca_eigenval)
n # use this number PC=1:n

# first convert to percentage variance explained
pve <- data.frame(PC = 1:n, pve = pca_eigenval/sum(pca_eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

### plotting

# plot pca PC1 and PC2
a <- ggplot(pca_2, aes(Axis1, Axis2, col = pop)) + 
  geom_point(size = 3) + 
  scale_colour_manual(values = c("red", "blue", "green")) +
  coord_equal() + 
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
a

# plot pca PC1 and PC3
b <- ggplot(pca_2, aes(Axis1, Axis3, col = pop)) + 
  geom_point(size = 3) + 
  scale_colour_manual(values = c("red", "blue", "green")) +
  coord_equal() + 
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))
b

# plot pca PC2 and PC3
c <- ggplot(pca_2, aes(Axis2, Axis3, col = pop)) + 
  geom_point(size = 3) + 
  scale_colour_manual(values = c("red", "blue", "green")) +
  coord_equal() + 
  theme_light() +
  xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))
c

ggsave("wg_PCA12.png", a, width = 1800, height = 1200, units = "px")
ggsave("wg_PCA13.png", b, width = 1800, height = 1200, units = "px")
ggsave("wg_PCA23.png", c, width = 1800, height = 1200, units = "px")


# Create a labeled plot
x <- ggplot(pca_2, aes(Axis1, Axis2, col = pop, label = ind)) + 
  geom_point(size = 2) + # size of points
  geom_text_repel(show.legend = FALSE, size=2) +  # Use geom_text_repel() for label repulsion
  scale_colour_manual(values = c("red", "blue", "green")) + 
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 2), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 2), "%)"))
x

### 3D plots

#PCA Visualization using plotly
#https://plotly.com/r/pca-visualization/
#3D plots
fig <- plot_ly(pca_2, x = ~Axis1, y = ~Axis2, z = ~Axis3, color =pop, colors = c("red", "blue", "green") ) %>%
  add_markers(size = 12)

fig <- fig %>%
  layout(
    title = "Arabidopsis thaliana Eurasian Accessions",
    scene = list(bgcolor = "#e5ecf6")
  )
fig

saveWidget(fig, "wg_3D_PCA.html")


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
genlight_vcf <- vcfR2genlight(vcf)

# Extract population data
pop = as.factor(population_labels)
genlight_vcf$pop = pop

# Convert genlight to stampp object
stampp_vcf = stamppConvert(genlight_vcf, type = "genlight")

# Calculate FST between populations
stamppFst = stamppFst(stampp_vcf, nboots = 100, percent = 95, nclusters = 8)
stamppFst_matrix = as.matrix(stamppFst$Fsts)

# Set diagonal to 0 and add symmetrical upper tri part of matrix
diag(stamppFst_matrix) <- 0
stamppFst_matrix[upper.tri(stamppFst_matrix)]  <- t(stamppFst_matrix)[upper.tri(stamppFst_matrix)]
# Optional: order the names
stamppFst_matrix = stamppFst_matrix[order(row.names(stamppFst_matrix)), order(colnames(stamppFst_matrix))]

heatmap(stamppFst_matrix)

# Create the heatmap with cell values and header
heatmaply(
  stamppFst_matrix,
  col = plasma(100),  # You can adjust the number of colors in the palette
  main = "FST values between populations.",
  cellnote = stamppFst_matrix,  # Display cell values
  notecol = "black",  # Color of the cell values
  margins = c(60, 70))


png("wg_heatmap.png")
heatmap.2(stamppFst_matrix)
dev.off()

#####################################################################
########## calculate genetic distance between individuals ###########
#####################################################################

stamppNeisD = stamppNeisD(stampp_vcf, pop = FALSE)
stamppNeisD_matrix = as.matrix(stamppNeisD)

# Set diagonal to 0 and add symmetrical upper tri part of matrix
diag(stamppNeisD_matrix) <- 0
stamppNeisD_matrix[upper.tri(stamppNeisD_matrix)]  <- t(stamppNeisD_matrix)[upper.tri(stamppNeisD_matrix)]
colnames(stamppNeisD_matrix) <- rownames(stamppNeisD_matrix)


head(stamppNeisD_matrix)
tolerance <- 1e-6  # Adjust the tolerance level as needed
stamppNeisD_matrix[which(abs(stamppNeisD_matrix - 22.627417) < 1)] = 0


png("wg_heatmap_big.png")
heatmap.2(stamppNeisD_matrix)
dev.off()

heatmap(stamppNeisD_matrix, RowSideColors = cc, ColSideColors = cc)

'accession_vector <- colnames(stamppNeisD_matrix)

# Function to map accessions to colors based on the "accessions" data frame
map_accessions_to_colors <- function(accession_vector, accessions_df) {
  # Find corresponding countries
  countries <- accessions_df$country[match(accession_vector, accessions_df$X1001.ID)]
  
  # Map countries to colors
  color_vector <- case_when(
    countries == "SWE" ~ "blue",
    countries == "UK" ~ "green",
    countries == "ESP" ~ "red",
    TRUE ~ NA_character_
  )
  
  return(color_vector)
}

# Map accessions to colors
cc <- map_accessions_to_colors(accession_vector, accessions)'








