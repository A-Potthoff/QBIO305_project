# install and load packages----

if(!require("ggplot2")){
  install.packages("ggplot2")
}
if(!require("ggrepel")){
  install.packages("ggrepel")
}

library(ggplot2)
library(ggrepel)

# load the data and filter it ----

data <- read.delim("GFF_final_with_annotation.txt", header = TRUE, stringsAsFactors = FALSE)

filtered_data <- data[-grep("unknown", data$Computational_description),]
filtered_data <- filtered_data[filtered_data$Curator_summary != "", ]
filtered_data$length <- filtered_data$V5-filtered_data$V4
filtered_data <- filtered_data[filtered_data$length >= 1500 & filtered_data$length <= 3000,]

T_genes <- filtered_data[grep(paste(c("response to heat",
                             "heat acclimation",
                             "HSP"#,
                             #"cold",
                             #"high temperature",
                             #"freezing"
                             ),
                            collapse = "|"), filtered_data$Computational_description), ]
T_genes$stress <- rep("heat", nrow(T_genes))

H_genes <- filtered_data[grep(paste(c("drought",
                             "dehydration",
                             "DEHYDRATION",
                             #"drown",
                             "water shortage",
                             "desiccation",
                             #"submersion",
                             #"water deprivation",
                             "Drought induced",
                             "arid"#,
                             #"evaporat"
                             ),
                            collapse = "|"), filtered_data$Computational_description), ]
H_genes$stress <- rep("drought", nrow(H_genes))
H_genes <- H_genes[H_genes$Curator_summary != "",]


HT_select <- c("AT1G05260.1",
              "AT1G08920.3",
              "AT4G39090.1",
              "AT1G47128.1",
              "AT1G76080.1",
              "AT4G36990.1",
              "AT1G06460.1",
              "AT3G62600.1",
              "AT5G62390.1",
              "AT1G56410.1")

cond <- c("dehydration & cold",
          "dehydration",
          "dehydration",
          "dehydration",
          "dehydration",
          "heat",
          "heat",
          "heat",
          "heat/cold",
          "heat/dehydration"
          )


select_data <- data.frame(Model_name = HT_select, cond = cond)

# Merge the data frames based on Model_name
select <- merge(select, select_data, by = "Model_name", all.x = TRUE)

colors = c("dehydration & cold" = "#50a1ff",
           "dehydration" = "#66c2ff",
           "heat" = "#ff6666",
           "heat/cold" = "#ff9933",
           "heat/dehydration" = "#ffbf00")

ggplot(select, aes(x=V4, y=0, fill = cond, color = cond, label = Model_name)) + 
  geom_point(aes(size = length)) +
  facet_wrap(~V1, ncol = 1) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  geom_text_repel() +
  xlab("position") + ylab("Chromosome")

write.csv(select, "selected_genes", row.names = FALSE)
