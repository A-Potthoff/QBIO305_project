#####################################
#### Extracting climate variables####
#####################################

library(sp)
library(raster)
library(geodata)
library(terra)

#Loading the coordinates of the accessions

samples = read.table("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/GWAS/Climate_variables/R_Script/coordinate_group4.txt", header = T)
head(samples)
lon<-samples$lon
lat<-samples$lat

# Extract coordinate data
xy <- samples[, c("lon", "lat")]

# Load BioClim data. The following command checks if the data is present at the 
# specified path. If the data is not present, it will be downloaded.
biodata = worldclim_global(var = "bio", res = 10, "C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/get_data_climate_variables/")
biodata # inspect the data
head(biodata)



# Names of selected bioclim variables
## BIO1 = Annual Mean Temperature, 
## BIO5 = Max Temperature of Warmest Month,
## BIO9 = Mean Temperature of Driest Quarter, 
## BIO10 = Mean Temperature of Warmest Quarter
## BIO14 = Precipitation of Driest Month
## BIO17 = Precipitation of Driest Quarter, 
## BIO18 = Precipitation of Warmest Quarter

#plotting the bioclim variable from the list
plot(biodata[[1]], main = "Annual Mean Temperature")
plot(biodata[[5]], main = "Max Temperature of Warmest Month")
plot(biodata[[9]], main = "Mean Temperature of Driest Quarter")
plot(biodata[[10]], main = "Mean Temperature of Warmest Quarter")
plot(biodata[[14]], main = "Precipitation of Driest Month")
plot(biodata[[17]], main = "Precipitation of Driest Quarter")
plot(biodata[[18]], main = "Precipitation of Warmest Quarter")


####plot bio14, bio17 and bio18 climate maps in comparison####
par(mfrow = c(3, 1))

# Plot the first graph
plot(biodata[[14]], main = "Bio14 - Precipitation of Driest Month")

# Plot the second graph
plot(biodata[[17]], main = "Bio17 - Precipitation of Driest Quarter")

# Plot the third graph
plot(biodata[[18]], main = "Bio18 - Precipitation of Warmest Quarter")

par(mfcol = c(1, 1))



# Extract Biolclimatic varaibles using xy coordinates dataframe
biodata_extract = extract(biodata[[1:19]], xy, df = T)
summary(biodata_extract) #you should have bio1 - bio19 columns now

#Attach it to the original df
samples_bio = cbind(samples, biodata_extract)

##############bio1######
plot(samples_bio$wc2.1_10m_bio_1)
# Extract required Climatic variable, for example we need 
bio1<-samples_bio$wc2.1_10m_bio_1
# Write bio1 data as a text file
write.table(bio1, file = "bio1.tsv", row.names = FALSE, col.names = FALSE)

##############bio5############
plot(samples_bio$wc2.1_10m_bio_5)
# Extract required Climatic variable, for example we need 
bio5<-samples_bio$wc2.1_10m_bio_5
# Write bio5 data as a text file
write.table(bio5, file = "bio5.tsv", row.names = FALSE, col.names = FALSE)

##############bio9############
plot(samples_bio$wc2.1_10m_bio_9)
# Extract required Climatic variable, for example we need 
bio9<-samples_bio$wc2.1_10m_bio_9
# Write bio9 data as a text file
write.table(bio9, file = "bio9.tsv", row.names = FALSE, col.names = FALSE)
##############bio10############
plot(samples_bio$wc2.1_10m_bio_10)
# Extract required Climatic variable, for example we need 
bio10<-samples_bio$wc2.1_10m_bio_10
# Write bio10 data as a text file
write.table(bio10, file = "bio10.tsv", row.names = FALSE, col.names = FALSE)
##############bio14############
plot(samples_bio$wc2.1_10m_bio_14)
# Extract required Climatic variable, for example we need 
bio14<-samples_bio$wc2.1_10m_bio_14
# Write bio14 data as a text file
write.table(bio14, file = "bio14.tsv", row.names = FALSE, col.names = FALSE)
##############bio17############
plot(samples_bio$wc2.1_10m_bio_17)
# Extract required Climatic variable, for example we need 
bio17<-samples_bio$wc2.1_10m_bio_17
# Write bio13 data as a text file
write.table(bio17, file = "bio17.tsv", row.names = FALSE, col.names = FALSE)
#################bio18############
plot(samples_bio$wc2.1_10m_bio_18)
# Extract required Climatic variable, for example we need 
bio18<-samples_bio$wc2.1_10m_bio_18
# Write bio13 data as a text file
write.table(bio18, file = "bio18.tsv", row.names = FALSE, col.names = FALSE)
#####PET######

# Read PET (Potential Evapo-Transpiration) monthly TIFF layer averaged for 1981-2010. Downloaded from: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F
####PET January####
pet1= raster("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/GWAS/Climate_variables/R_Script/CHELSA_pet_penman_06_1981-2010_V.2.1.tif")

# plot PET1
plot(pet1, main = "Potential Evapo-Transpiration (PET) January (1981-2010)")

#Extract PET January
pet1= raster("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/GWAS/Climate_variables/R_Script/CHELSA_pet_penman_01_1981-2010_V.2.1.tif")

pet_1 = extract(pet1, xy, df = T)

write.table(pet_1, file = "pet_1.txt", row.names = FALSE, col.names = FALSE)

pet_jan <- pet_1[,2]

write.table(pet_jan, file = "pet_jan.txt", row.names = FALSE, col.names = FALSE)

####PET June####
pet6= raster("C:/Users/lobki/OneDrive/Desktop/QBio Semester 3/305/Own project/get_data_climate_variables/get_data_climate_variables/CHELSA_pet_penman_06_1981-2010_V.2.1.tif")

# plot PET6
plot(pet6, main = "Potential Evapo-Transpiration (PET) June (1981-2010)")

pet_6 = extract(pet6, xy, df = T)

write.table(pet_6, file = "pet_6.txt", row.names = FALSE, col.names = FALSE)

pet_june <- pet_6[,2]

write.table(pet_june, file = "pet_june.txt", row.names = FALSE, col.names = FALSE)
