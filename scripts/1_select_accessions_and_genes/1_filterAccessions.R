library(ggplot2)
library(ggrepel)

setwd("C:/Users/andre/OneDrive/Bildung/3_HHU/quant_Bio/3rd/QBio305_Population_And_Quantitative_Genetics/Cologne/exercises/project/own_research/analysis/select_accessions")

set.seed(4328)

#accessions.txt was obtained from the 10001 genome procject (download as csv-file)

acc <- read.csv("accessions.txt")
colnames(acc) <- c("x1001.ID",
                   "collector_org",
                   "species_variant",
                   "country",
                   "name",
                   "lat",
                   "long",
                   "collector",
                   "some_date",
                   "cs",
                   "admixture_group",
                   "x", "y")

acc_swe <- acc[acc$country=="SWE",]
acc_uk <- acc[acc$country=="UK",]
acc_esp <- acc[acc$country=="ESP",]

acc_tot = rbind(acc_swe, rbind(acc_uk, acc_esp))

#plotting the lat and long gives us a representation of their distribution
#(hopefully those will be as clearly separated in the PCA-plot :D)
ggplot(acc_tot, aes(x = long, y = lat, color = admixture_group)) + 
  geom_point()

#from the dataset, we cannot look at one species over all countries. we have to do a mixture of these species
length(unique(acc_tot$name))

#we will try to use samples that come from the same ancestral populations.
acc_swe.pur <- acc_swe[acc_swe$admixture_group %in% c("north_sweden"),]
acc_uk.pur <- acc_uk[acc_uk$admixture_group=="western_europe",]
acc_esp.pur <- acc_esp[acc_esp$admixture_group=="spain",]

acc_tot = rbind(acc_swe.pur, rbind(acc_uk.pur, acc_esp.pur))

#plotting the lat and long gives us a representation of their distribution
#(hopefully those will be as clearly separated in the PCA-plot :D)
ggplot(acc_tot, aes(x = long, y = lat, color = admixture_group)) + 
  geom_point()

rand.i.swe <- sample(nrow(acc_swe.pur), 40)
rand.i.uk <- sample(nrow(acc_uk.pur), 40)
rand.i.esp <- sample(nrow(acc_esp.pur), 40)

swe <- acc_swe.pur[rand.i.swe,]
uk <- acc_uk.pur[rand.i.uk,]
esp <- acc_esp.pur[rand.i.esp,]

tot <- rbind(swe, rbind(uk, esp))

tot <- read.csv("selected_accessions.csv")
ggplot(tot, aes(x = long, y = lat, color = admixture_group)) + 
  geom_point()

write.csv(tot, "selected_accessions.csv", row.names = FALSE)

plot <- ggplot(tot[which(tot$country == "SWE"),], aes(x = long, y = lat, color = admixture_group, label = X1001.ID)) + 
  geom_point() + 
  geom_text_repel(max.overlaps = Inf)
plot
ggsave("plot.png", plot, width = 3000, height = 2500, units = "px")
