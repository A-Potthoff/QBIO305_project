

########## Load necessary libraries ###########
library(BiocManager)
library(KEGGREST)
library(org.At.tair.db)
library(Rgraphviz)
library(topGO)
library(biomaRt)
library(ggplot2)
library(AnnotationDbi)
library(clusterProfiler)
library(scales)

setwd("/Users/sarakreci/Desktop")

gene_list <- read.table("GO_input.txt", header=TRUE) # output.csv

univ <- gene_list[, 3]         
names(univ) <- gene_list[, 1] 
univ <- univ[!is.na(univ)]     

# A function to return TRUE/FALSE for p-values < 0.05
selection <- function(allScore) { return(allScore < 0.05) }

# Preparing the topGO data
tGOdata <- new("topGOdata", description = "Simple session", ontology = "BP", geneSel = selection, allGenes = univ, nodeSize = 3, mapping = "org.At.tair.db", annot = annFUN.org)

# Run the Gene Enrichment Analysis
results.ks <- runTest(tGOdata, algorithm = "elim", statistic = "ks")

# Generate a table of enriched GO terms
goEnrichment <- GenTable(tGOdata, KS = results.ks, orderBy = "KS", topNodes = 50)
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS < 0.05, ]
goEnrichment <- goEnrichment[, c("GO.ID", "Term", "KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ", ")

class(goEnrichment$KS)

goEnrichment$KS <- as.numeric(goEnrichment$KS)

any(is.na(goEnrichment$KS))
head(goEnrichment[is.na(goEnrichment$KS), ])

goEnrichment <- goEnrichment[!is.na(goEnrichment$KS), ]

write.csv(goEnrichment, "sample_project.csv")


ntop <- 10
ggdata <- goEnrichment[1:ntop, ]

ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term))  # Fix order

# Plotting using ggplot2
gg1 <- ggplot(ggdata,
              aes(x = Term, y = -log10(KS), size = -log10(KS), fill = -log10(KS))) +
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5, 12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  xlab('') + ylab('Enrichment score') +
  labs(
    title = 'Pathways Enriched in Upregulated Genes in Stress',
    subtitle = 'Top 10 Terms Ordered by Kolmogorov-Smirnov p-Value',
    caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001'
  ) +
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),colour = c("black", "black", "black"),size = c(0.5, 1.5, 3)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 18, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, vjust = 1),
    axis.text.x = element_text(angle = 0, size = 12, hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12), axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'), axis.line = element_line(colour = 'black'),
    legend.key = element_blank(),legend.key.size = unit(1, "cm"),legend.text = element_text(size = 16),
    title = element_text(size = 12)) + coord_flip()

# Display and save the plot
print(gg1)
ggsave("GO_results_defense.png", plot = gg1, width = 11, height = 10, dpi = 300)
