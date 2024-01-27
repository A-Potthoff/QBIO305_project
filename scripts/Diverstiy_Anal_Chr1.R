########################################################################
##         Estimate and plot Fst and Tajima'D and Neutrality          ##
##                       stats using PopGenome                        ##
########################################################################

library(PopGenome) 
library(vcfR)

setwd("C:/Users/andre/OneDrive/Bildung/3_HHU/quant_Bio/3rd/QBio305_Population_And_Quantitative_Genetics/Cologne/exercises/project/own_research/analysis/PopGenome")

# If you want to check number columns and start and stop positions then you have to read 
# your "input.vcf" file using "read.vcfR" function from "vcfR" package

at.VCF <- read.vcfR("final_vcf.gz")
str(at.VCF)
head(at.VCF@fix)

#get start and stop positions of your vcf file
head(getFIX(at.VCF))
tail(getFIX(at.VCF))

###
# PopGenome accepts one chromosome/scaffold at time, therefore, you need
# to find start and stop positions of every chromosome or scaffold
###

library(VariantAnnotation)
# Read the VCF file
vcf<-readVcf("final_vcf.gz")

# Extract chromosome names
chromosomes <- seqlevels(vcf)
# Initialize an empty data frame to store results
chromosome_ranges <- data.frame(CHROM = character(), Start = numeric(), End = numeric(), stringsAsFactors = FALSE)

# Iterate over chromosomes
for (chrom in chromosomes) {
  # Extract positions for the current chromosome
  positions <- start(vcf[seqnames(vcf) == chrom])
  
  # Append results to the data frame
  chromosome_ranges <- rbind(chromosome_ranges, data.frame(CHROM = chrom, Start = min(positions), End = max(positions)))
}

chromosome_ranges

###
# Change integer in <tid="1"> to repeat it for all chromosomes

#the following commmand requires a <vcf_filename>.tbi file inside the same directory!

###
At_Chr1 <- readVCF("final_vcf.gz", numcols=120, tid="1", frompos = 1766730, topos = 30108069, include.unknown = TRUE)

#To the class of object At_Chr1
class(At_Chr1)

At_Chr1 ###this is your genome.class data. You can push all genomics analysis into one object. 

####Examining the variant data
#Remember, you can look the data we have read in using the following command:
get.sum.data(At_Chr1)

#From the n.biallelic.sites we can see there are  2644 bilallelic SNPs and from n.polyallelic.sites,
#there are 0 positions with more than two alleles. So in total we have:
At_Chr1@n.biallelic.sites + At_Chr1@n.polyallelic.sites

#To see what slots in Genome class
show.slots(At_Chr1)#

#To check total number of sites
At_Chr1@n.sites

#To check starting position and last position of genome class
At_Chr1@region.names

#####Define populations in your dataset####

#If you look at this, you will only see a blank list. So we need to supply our population data to
#the ch4 object. To make naming our populations simple, we will read in some external data. 
At_Chr1@populations # check for population data

library(readr)

###population data is stored in data.frame that has two columns, one column for individual name, one column for pop

population_info <- read_delim("pop1.txt", delim = "\t")
# now get the data for the populations
populations <- split(population_info$sample, population_info$pop)

# now set 
At_Chr1 <- set.populations(At_Chr1, populations, diploid = T)
##check if it worked
At_Chr1@populations

####Setting up sliding windows###

#Per-SNP estimates of statistics such as Pi can often be extremely noisy when you are calculating them on
#very large numbers of markers. As well as this, there are issues with the fact that SNP positions in close
#proximity are not always independent due to recombination - this is a theme we will return too shortly. 
#So for this reason, it is often better to use a sliding-window approach - i.e. split the genome into
#windows of a particular size and then calculate the mean for a statistic within that window.

#We know already that chromosome 1 is 18584000 bp long, so we can get an idea of how many sliding windows
#we would generate by using some R code. We'll set our sliding window to be 100 bp wide.
#We will also set a step or jump for our window of 50 bp.

##To check total number of sites
At_Chr1@n.sites
chr1 <- At_Chr1@n.sites
# If you get mismatch error then add +1 to the value as below
#chr1 <- 28629728 # set the value equal to At_Chr@n.sites +1

# set window size and window jump
window_size <- 100
window_jump <- 50

# use seq to find the start points of each window
window_start_Chr1 <- seq(from = 1, to = chr1, by = window_jump)
# add the size of the window to each start point 
window_stop_Chr1 <- window_start_Chr1 + window_size

# no windows start before the end of chromosome 4
sum(window_start_Chr1 > chr1)
# but some window stop positions do occur past the final point
sum(window_stop_Chr1 > chr1)

# remove windows from the start and stop vectors
window_start_Chr1 <- window_start_Chr1[which(window_stop_Chr1 < chr1)]
window_stop_Chr1 <- window_stop_Chr1[which(window_stop_Chr1 < chr1)]

chr1 - window_stop_Chr1[length(window_stop_Chr1)]

# save as a data.frame
windows_Chr1 <- data.frame(start = window_start_Chr1, stop = window_stop_Chr1, 
                      mid = window_start_Chr1 + (window_stop_Chr1-window_start_Chr1)/2)

#https://rdrr.io/cran/PopGenome/man/sliding.window.transform-methods.html
# make a sliding window dataset
At_sw_Chr1 <- sliding.window.transform(At_Chr1, width = 100, jump = 50, type = 2)


#######Calculating sliding window estimates of nucleotide diversity and differentiation#####
#Now that we have set up the data, the population information and the sliding windows, it is quite
#straightforward for us to calculate some statistics we are interested in. In this case, we are going
#to calculate nucleotide diversity (i.e. Pi) and FST. We will also generate a third statistic, d_XY_,
#which is the absolute nucleotide divergence between two populations.

#First we will calculate Pi. Handily, the following command also sets up what we need for d_XY_.

# calculate diversity statistics
At_sw_Chr1 <- diversity.stats(At_sw_Chr1, pi = TRUE)


#Next we will calculate FST, which again is very straight forward with a single command.

### calculate diversity statistics
At_sw_Chr1 <- F_ST.stats(At_sw_Chr1, mode = "nucleotide")

#Note that here we use mode = "nucleotide" to specify we want it to be calculated sliding averages
#of nucleotides, rather than using haplotype data, which is the alternative. And that's it for 
#calculating the statistics! As you will see in the next section, extracting them from the 
#At_sw object is actually more difficult than generating them

#### calculate neutrality statistics####
At_sw_Chr1 <- neutrality.stats(At_sw_Chr1)

####Extracting statistics for visualization####

#Since we ran our analysis on a sliding-window basis, we should have estimates of ??, FST and d_XY_ for
#each window. What we want to do now is extract all our statistics and place them in a single data.frame
#for easier downstream visualisation - this will let us identify how these statistics are interrelated.

#First of all, we will get the nucleotide diversity data.

# extract nucleotide diversity and correct for window size
nd_Chr1 <- At_sw_Chr1@nuc.diversity.within/100

#This is straightforward, but remember also that our estimates need to be corrected for 
#window size - so we divide them by 100 bp here. We should also add the population names
#to each of them, since they are not already set.

# make population name vector
pops <- c("ESP", "SWE","UK")
# set population names
colnames(nd_Chr1) <- paste0(pops, "_pi")

# extract fst values
fst_Chr1 <- t(At_sw_Chr1@nuc.F_ST.pairwise)
#Note that here, we need to use t() to transpose the F_ST matrix so that each column is a pairwise
#comparison and each row is an estimate for a genome window. Since F_ST is pairwise, the column
#names are also quite different and will also be the same for d_XY_, which is also a pairwise measure.

#So now we are ready to extract our final statistic, d_XY_. We can do this in a similar way to how we
#handled the FST data.

# extract dxy - pairwise absolute nucleotide diversity
dxy_Chr1 <- get.diversity(At_sw_Chr1, between = T)[[2]]/100
#As with nucleotide diversity, we also corrected d_XY_ for the window size.

#Now we sort out the column names for our FST and d_XY_ data. This is where our R skills come in use!
#We will need to use some R-based string manipulation. The column names are identical for both datasets,
#so we will take the first one and use the sub function to replace the population names.

# get column names 
x <- colnames(fst_Chr1)

# Loop through each population and replace the corresponding population name in the column names
for (i in 1:length(pops)) {
  pattern <- paste0("pop", i)
  x <- sub(pattern, pops[i], x)
}

# look at x to confirm the replacement has occurred
print(x)

# replace forward slash
x <- sub("/", "_", x)
# look at x to confirm the replacement has occurred
x

#Now all that we need to do is make clear these names are for either FST or d_XY_. 
#The best way to do this is to append a suffix to our vector of pairwise comparison names.
#We can do this using paste0
paste0(x, "_fst")
paste0(x, "_dxy")

#So this function allows us to join strings together in a character vector. Very useful. 
#Next we will actually change the column names of our two data matrices, before we put 
#everything together in our final dataset.
colnames(fst_Chr1) <- paste0(x, "_fst")
colnames(dxy_Chr1) <- paste0(x, "_dxy")

#extract Tajma's D and set population names
td_Chr1 <- At_sw_Chr1@Tajima.D/100

colnames(td_Chr1) <- paste0(pops, "_td")

#Ok so now that our td, Pt, FST and d_XY_ datasets are ready, we can combine them all together
#with our windows information from earlier into a big dataset.
library(tibble)
At_data_Chr1 <- as.tibble(data.frame(windows_Chr1, td_Chr1, nd_Chr1, fst_Chr1, dxy_Chr1))

###############################################################################
### Compare and visualize stats between different populations - as boxplot  ###
###############################################################################

# For the purposes of this session, we will focus mainly on the difference between Italian (IT) and Swedish (SW)
# Arabidopsis thaliana population.

###
# Let's look at nucleotide diversity (Pi), we can do that like so:
###

library(dplyr)

# select nucleotide diversity data and calculate means
At_data_Chr1 %>% select(contains("pi")) %>% summarise_all(mean)

#we used select and contains to select columns from our main dataset that contain 
#pi - i.e. nucleotide diversity columns. We then used summarise_all and mean to calculate
#the mean value for all of the columns we selected.

# To plot this we need to use "gather" on the data
library(tidyr)
pi_g_Chr1 <- At_data_Chr1 %>% select(contains("pi")) %>% gather(key = "populations", value = "pi")

# make a boxplot
library(ggplot2)
a_Pi_Chr1 <- ggplot(pi_g_Chr1, aes(populations, pi)) + geom_boxplot() + theme_light() + xlab(NULL)
a_Pi_Chr1

# Taking the logarithm (in this case, log base 10) of the values can be useful when dealing with data
# that spans several orders of magnitude. This transformation helps in visually emphasizing relative
# differences in the data, especially when there are large variations in scale.

pi_g_Chr1$log_pi <- log10(pi_g_Chr1$pi)

a_pi_Chr1 <- ggplot(pi_g_Chr1, aes(populations, log_pi, fill = populations)) + 
  geom_boxplot(color = "black") +  # Border color
  scale_fill_manual(values = c("red", "blue","darkgreen","magenta")) +  # Box fill colors
  theme_light() + 
  xlab(NULL) +
  ylab("Log10(pi)")+
  ggtitle("Pi between populations")

mean(pi_g_Chr1$log_pi[(pi_g_Chr1$log_pi != -Inf) & (pi_g_Chr1$populations == "UK_pi")])

a_pi_Chr1

#This makes it much clearer how nucleotide diversity differs among the populations.

#When comparing two boxplots to determine if they are statistically different, one can perform
# statistical tests such as the t-test or Wilcoxon rank-sum test. 
# Wilcoxon rank-sum test
# Kruskal-Wallis test
#kruskal_test_pi_Chr1 <- kruskal.test(log_pi ~ populations, data = pi_g_Chr1)

# Print the result
#print(kruskal_test_pi_Chr1)

# You will get an error like below when there are more than two populations in your dataset
# Error in wilcox.test.formula(log_pi ~ populations, data = pi_g) : 
# grouping factor must have exactly 2 levels

unique(pi_g_Chr1$populations)
# [1] "ITA_pi" "GER_pi" "SWE_pi"

# You can proceed to Wilcoxon rank-sum test after you
# filter data and two select two populations of your choice
# not correct way to do multiple pairwise comparisons
#comparison_data_Chr1 <- pi_g_Chr1 %>% filter(populations %in% c("IT_pi", "SW_pi"))

# Perform Wilcoxon rank-sum test
#wilcox_test_pi_Chr1 <- wilcox.test(log_pi ~ populations, data = comparison_data_Chr1)

# Print the test result
#print(wilcox_test_pi_Chr1)

# Add p-value to the plot
#a_pi_Chr1 + annotate("text", x = 1.5, y = max(pi_g_Chr1$log_pi), label = paste("p =", format.pval(wilcox_test_pi_Chr1$p.value, digits = 3)))

########
# When you have more than two groups (four in this case), and you want to compare the 
# central tendency of the distributions, you might consider using analysis of variance (ANOVA) or
# its non-parametric counterpart, the Kruskal-Wallis test.
#######

# Kruskal-Wallis test
kruskal_test_pi_Chr1 <- kruskal.test(log_pi ~ populations, data = pi_g_Chr1)

# Print the test result
print(kruskal_test_pi_Chr1)

# generate boxplot again
a_pi_Chr1

# Add Kruskal-Wallis test p-value to the plot
a_pi_Chr1 + annotate("text", x = 1.5, y = max(pi_g_Chr1$log_pi), label = paste("Kruskal-Wallis p =", format.pval(kruskal_test_pi_Chr1$p.value, digits = 3)))


###
# Let's look at TajimaD (td), we can do that like so:
###

library(dplyr)

# select nucleotide diversity data and calculate means
At_data_Chr1 %>%
  select(contains("td")) %>%
  summarise_all(~mean(., na.rm = TRUE))
# no mean calculated due to missing values (NA)

#we used select and contains to select columns from our main dataset that contain 
#pi - i.e. nucleotide diversity columns. We then used summarise_all and mean to calculate
#the mean value for all of the columns we selected.

# To plot this we need to use "gather" on the data
library(tidyr)
td_g_Chr1 <- At_data_Chr1 %>% select(contains("td")) %>% gather(key = "populations", value = "td")

# Remove rows with missing values
td_g2_Chr1 <- na.omit(td_g_Chr1)

# make a boxplot


# Make a boxplot using ggplot2
a_td_g2_Chr1 <- ggplot(td_g2_Chr1, aes(x = populations, td, fill = populations)) +
  geom_boxplot(color = "black") +  # Border color
  scale_fill_manual(values = c("red", "blue","darkgreen","magenta")) +  # Box fill colors, choose number of colours based on number of populations
  xlab(NULL) +
  theme_light() + 
  ylab("Tajima's D") +
  ggtitle("Tajima's D Between populations")
a_td_g2_Chr1

a2_td_g2_Chr1 <- ggplot(td_g2_Chr1, aes(populations, td, fill = populations)) + 
  geom_violin(color = "black", scale = "width", width = 0.8) +  # Violin plot with black border
  scale_fill_manual(values = c("red", "blue", "orange", "magenta")) +  # Fill colors for violins
  stat_summary(fun = "mean", geom = "point", shape = 23, size = 2, fill = "white", color = "black") +  # Add mean points
  theme_light() + 
  xlab(NULL) +
  ylab("Tajima's D")

a2_td_g2_Chr1

# Taking the logarithm (in this case, log base 10) of the values can be useful when dealing with data
# that spans several orders of magnitude. This transformation helps in visually emphasizing relative
# differences in the data, especially when there are large variations in scale.

td_g2_Chr1$log_td <- log10(td_g2_Chr1$td)

a2_td_g2_Chr1 <- ggplot(td_g2_Chr1, aes(populations, log_td, fill = populations)) + 
  geom_violin(color = "black") +  # Border color
  scale_fill_manual(values = c("red", "blue","orange","magenta")) +  # Box fill colors # Box fill colors, choose number of colours based on number of populations
  theme_light() + 
  xlab(NULL) +
  ylab("Log10(Tajima's D)")

a2_td_g2_Chr1


#When comparing two boxplots to determine if they are statistically different, one can perform
# statistical tests such as the t-test or Wilcoxon rank-sum test. 
# Wilcoxon rank-sum test
#wilcox_test_td_Chr1 <- wilcox.test(log_td ~ populations, data = td_g2_Chr1)

# You will get an error like below when there are more than two populations in your dataset
# Error in wilcox.test.formula(log_pi ~ populations, data = pi_g) : 
# grouping factor must have exactly 2 levels

#unique(td_g_Chr1$populations)
# [1] "ITA_pi" "GER_pi" "SWE_pi"

# You can proceed to Wilcoxon rank-sum test after you
# filter data and two select two populations of your choice
# not correct way to do multiple pairwise comparisons
#comparison_data_Chr1 <- td_g2_Chr1 %>% filter(populations %in% c("IT_td", "SW_td"))

# Perform Wilcoxon rank-sum test
#wilcox_test_td_Chr1 <- wilcox.test(log_td ~ populations, data = comparison_data_Chr1)

# Print the test result
#print(wilcox_test_td_Chr1)

# Add p-value to the plot
#a2_td_g2_Chr1 + annotate("text", x = 1.5, y = max(td_g2_Chr1$log_td), label = paste("p =", format.pval(wilcox_test_td_Chr1$p.value, digits = 3)))

# adjust text positions by changing x and y below
#a2_td_g2_Chr1 + annotate("text", x = 1.5, y = -1, label = paste("Wilcox test p =", format.pval(wilcox_test_td_Chr1$p.value, digits = 3)))

########
# When you have more than two groups (four in this case), and you want to compare the 
# central tendency of the distributions, you might consider using analysis of variance (ANOVA) or
# its non-parametric counterpart, the Kruskal-Wallis test.
#######

# Kruskal-Wallis test
kruskal_test_td_Chr1 <- kruskal.test(td ~ populations, data = td_g2_Chr1)

# Print the test result
print(kruskal_test_td_Chr1)

# generate boxplot again
a2_td_g2_Chr1

# Add Kruskal-Wallis test p-value to the plot
a2_td_g2_Chr1 + annotate("text", x = 1.5, y = max(td_g2_Chr1$log_pi), label = paste("Kruskal-Wallis p =", format.pval(kruskal_test_td_Chr1$p.value, digits = 3)))

# adjust text positions by changing x and y below
a2_td_g2_Chr1 + annotate("text", x = 1.5, y = -1, label = paste("Kruskal-Wallis p =", format.pval(kruskal_test_td_Chr1$p.value, digits = 3)))



#############################################
#   Compare Inter population distributions  #
#           of Pi and Tajima' D             #
#############################################

###
# Compare Pi density distributions of populations 
###

# modify code change xlim to see complete distribution
plot(density(log(At_data_Chr1$SWE_pi)), col= "blue", main="Distribution log Pi", xlim=c(-9, max(log(At_data_Chr1$SWE_pi))))

# Add density plots for other populations with different colors
lines(density(log(At_data_Chr1$UK_pi)), col="darkgreen")

# If you have more populations add more lines
lines(density(log(At_data_Chr1$ESP_pi)), col="red")

# When generating density plots, if peaks extend beyond the y-axis limits ( as you see the peak for SW
# is cut off), adapting your R code to plot the SW population first and then plot the others to resolve
# this issue.


# Add legend
legend("topright", legend=c("ESP", "UK", "SWE"),
       col=c("red", "darkgreen", "blue"), lty=1,
       title="Population")

install.packages("dunn.test")
# Load the dunn.test package
library(dunn.test)

# Perform Kruskal-Wallis test
pi_data <- list(
  SWE = log(At_data_Chr1$SWE_pi),
  ESP = log(At_data_Chr1$ESP_pi),
  UK = log(At_data_Chr1$UK_pi)
)

kruskal_pi_dist <- kruskal.test(pi_data)

# Print the Kruskal-Wallis test result
print(kruskal_pi_dist)

# Perform post-hoc Dunn's test if Kruskal-Wallis is significant
if (kruskal_pi_dist$p.value < 0.05) {
  dunn_result <- dunn.test(pi_data)
  
  # Print the post-hoc Dunn's test results
  print(dunn_result)
}

# The post-hoc Dunn's test (dunn.test) will only be performed and printed if the p-value from the 
# Kruskal-Wallis test (kruskal_result$p.value) is less than 0.05.
# If kruskal_result$p.value is greater than or equal to 0.05, the if condition is not satisfied,
# and the code inside the if block, including the dunn.test and print(dunn_result), will not be executed.

##

# Add legend for Kruskal-Wallis test p-value
legend("topleft", 
       legend=paste("Kruskal-Wallis p-value:", format(kruskal_pi_dist$p.value, digits=4)), 
       bty="n", 
       cex=1)


###
# Compare "Tajima's D" density distributions of populations
###

subset_SWE_td <- At_data_Chr1$SWE_td[!is.na(At_data_Chr1$SWE_td)]
subset_ESP_td <- At_data_Chr1$ESP_td[!is.na(At_data_Chr1$ESP_td)]
subset_UK_td <- At_data_Chr1$UK_td[!is.na(At_data_Chr1$UK_td)]

# Plot the density distribution
plot(density(subset_SWE_td), col = "blue", main = "Distribution Tajima's D", xlim = c(-0.03, 0.04), ylim = c(0, 65))
# Add density plots for other populations with different colors
lines(density(subset_UK_td), col = "darkgreen")
lines(density(subset_ESP_td), col = "red")

# Add legend
legend("topright", legend = c("SWE", "UK", "ESP"),
       col = c("blue", "darkgreen", "red"), lty = 1,
       title = "Population")



#install.packages("dunn.test")
# Load the dunn.test package
library(dunn.test)

# Perform Kruskal-Wallis test
td_data <- list(
  SW = subset_SW_td,
  IT = subset_IT_td
)

kruskal_td_dist <- kruskal.test(td_data)

# Print the Kruskal-Wallis test result
print(kruskal_td_dist)

# Not required when comparing only two populations
# Perform post-hoc Dunn's test if Kruskal-Wallis is significant
if (kruskal_td_dist$p.value < 0.05) {
  dunn_result <- dunn.test(pi_data)
  
  # Print the post-hoc Dunn's test results
  print(dunn_result)
}

##
# The post-hoc Dunn's test (dunn.test) will only be performed and printed if the p-value from the 
# Kruskal-Wallis test (kruskal_result$p.value) is less than 0.05.
# If kruskal_result$p.value is greater than or equal to 0.05, the if condition is not satisfied,
# and the code inside the if block, including the dunn.test and print(dunn_result), will not be executed.
##


# Add legend for Kruskal-Wallis test p-value
legend("topleft", 
       legend=paste("Kruskal-Wallis p-value:", format(kruskal_td_dist$p.value, digits=4)), 
       bty="n", 
       cex=1)

################################################
#   Visualizing patterns along the chromosome  #
#               Facets plots                   #
################################################

library(dplyr)
library(tidyr)
library(ggplot2)

###
# Pi
###
# Select data of interest
hs_pi_Chr1 <- At_data_Chr1 %>%
  select(mid, SWE_pi, UK_pi, ESP_pi)

# Use gather to rearrange everything
hs_pi_g_Chr1 <- gather(hs_pi_Chr1, -mid, key = "stat", value = "value")

# Reorder the levels of the stat factor
hs_pi_g_Chr1$stat <- factor(hs_pi_g_Chr1$stat, levels = c("SWE_pi", "ESP_pi", "UK_pi"))

# Take the logarithm of the value variable
hs_pi_g_Chr1$log_value <- log10(hs_pi_g_Chr1$value)

# Construct a plot with facets
a_pi_Chr1 <- ggplot(hs_pi_g_Chr1, aes(mid / 10^6, log_value, colour = stat)) + geom_line() +
  facet_wrap(~stat, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()

# Show the plot
a_pi_Chr1

###
# FST
###

# Select data of interest
hs_fst_Chr1 <- At_data_Chr1 %>%
  select(mid, IT_SW_fst)

# Mutate to set FST and dXY values smaller than zero to zero
hs_fst2_Chr1 <- hs_fst_Chr1 %>%
  mutate(across(c(IT_SW_fst), 
                ~ ifelse(. < 0, 0, .)))

# Use gather to rearrange everything
hs_fst_g_Chr1 <- gather(hs_fst2_Chr1, -mid, key = "stat", value = "value")

# Reorder the levels of the stat factor
hs_fst_g_Chr1$stat <- factor(hs_fst_g_Chr1$stat, levels = c("IT_SW_fst"))

# Take the logarithm of the value variable
hs_fst_g_Chr1$log_value <- log10(hs_fst_g_Chr1$value)

# Construct a plot with facets
a_fst_Chr1 <- ggplot(hs_fst_g_Chr1, aes(mid / 10^6, log_value, colour = stat)) + geom_line() +
  facet_wrap(~stat, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()

# Show the plot
a_fst_Chr1

###
# Dxy
###

# Select data of interest
hs_dxy_Chr1 <- At_data_Chr1 %>%
  select(mid, IT_SW_dxy)

# Mutate to set FST and dXY values smaller than zero to zero
hs_dxy2_Chr1 <- hs_dxy_Chr1 %>%
  mutate(across(c(IT_SW_dxy), 
                ~ ifelse(. < 0, 0, .)))

# Use gather to rearrange everything
hs_dxy_g_Chr1 <- gather(hs_dxy2_Chr1, -mid, key = "stat", value = "value")

# Reorder the levels of the stat factor
hs_dxy_g_Chr1$stat <- factor(hs_dxy_g_Chr1$stat, levels = c("IT_SW_dxy"))

# Take the logarithm of the value variable
hs_dxy_g_Chr1$log_value <- log10(hs_dxy_g_Chr1$value)

# Construct a plot with facets
a_dxy_Chr1 <- ggplot(hs_dxy_g_Chr1, aes(mid / 10^6, log_value, colour = stat)) + geom_line() +
  facet_wrap(~stat, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()

# Show the plot
a_dxy_Chr1


###
# TajimaD
###
# Select data of interest
hs_td_Chr1 <- At_data_Chr1 %>%
  select(mid, SWE_td, UK_td, ESP_td)

# Use gather to rearrange everything
hs_td_g_Chr1 <- gather(hs_td_Chr1, -mid, key = "stat", value = "value")

# Reorder the levels of the stat factor
hs_td_g_Chr1$stat <- factor(hs_td_g_Chr1$stat, levels = c("SWE_td", "UK_td", "ESP_td"))

# Take the logarithm of the value variable
hs_td_g_Chr1$log_value <- log10(hs_td_g_Chr1$value)

# Construct a plot with facets
a_td_Chr1 <- ggplot(hs_td_g_Chr1, aes(mid / 10^6, log_value, colour = stat)) + geom_line() +
  facet_wrap(~stat, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()

# Show the plot
a_td_Chr1

####
# If there are only two populations you can also choose to put FST and DXY in the same Facet plot
# 
####

# Select data of interest
hs_fst_dxy_Chr1 <- At_data_Chr1 %>%
  select(mid, IT_SW_dxy, IT_SW_fst)

# Mutate to set FST and dXY values smaller than zero to zero
hs2_fst_dxy_Chr1 <- hs_fst_dxy_Chr1 %>%
  select(mid, IT_SW_dxy, IT_SW_fst) %>%
  mutate(across(c(IT_SW_fst, IT_SW_dxy), 
                ~ ifelse(. < 0, 0, .)))

hs <- At_data %>%
  select(mid, ESP_pi, SWE_pi, ESP_SWE_fst, ESP_SWE_dxy) %>%
  mutate(across(c(ESP_SWE_fst, ESP_SWE_dxy), ~ ifelse(. < 0, 0, .)))

# Use gather to rearrange everything
hs_fst_dxy_g_Chr1 <- gather(hs2_fst_dxy_Chr1, -mid, key = "stat", value = "value")

# Reorder the levels of the stat factor
hs_fst_dxy_g_Chr1$stat <- factor(hs_fst_dxy_g_Chr1$stat, levels = c("IT_SW_fst","IT_SW_dxy"))

# Take the logarithm of the value variable
hs_fst_dxy_g_Chr1$log_value <- log10(hs_fst_dxy_g_Chr1$value)

# Construct a plot with facets
a_fst_dxy_g_Chr1 <- ggplot(hs_fst_dxy_g_Chr1, aes(mid / 10^6, log_value, colour = stat)) + geom_line() +
  facet_wrap(~stat, scales = "free_y", ncol = 1) +
  xlab("Position (Mb)") +
  theme_light()

# Show the plot
a_fst_dxy_g_Chr1

