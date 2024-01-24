path <- "C:/Users/andre/OneDrive/Bildung/3_HHU/quant_Bio/3rd/QBio305_Population_And_Quantitative_Genetics/Cologne/exercises/project/own_research/analysis/2D_SFS"
setwd(path)

file <- "FINAL_polarized.vcf.gz"

#https://github.com/shenglin-liu/vcf2sfs
source("vcf2sfs/vcf2sfs.r")

mygt<-vcf2gt(file, "pop1.txt")
df <- gt2snp(mygt, "pop1.txt")

res = 500
doub_width = 20
doub_height = 20

#two-dim

d <- mysfs2<-gt2sfs.raw(mygt, c("SWE", "UK"))

png("../results/SFS/SWE_UK.png", width = doub_width, height = doub_height, units = "cm", res = res)
plot.sfs(d)
dev.off()

e <- mysfs2<-gt2sfs.raw(mygt, c("SWE", "ESP"))

png("../results/SFS/SWE_ESP.png", width = doub_width, height = doub_height, units = "cm", res = res)
plot.sfs(e)
dev.off()

f <- mysfs2<-gt2sfs.raw(mygt, c("UK", "ESP"))

png("../results/SFS/UK_ESP.png", width = doub_width, height = doub_height, units = "cm", res = res)
plot.sfs(f)
dev.off()

