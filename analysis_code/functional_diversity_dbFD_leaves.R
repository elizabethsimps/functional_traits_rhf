# FD of leaf traits across elevation and aspect #
# Elizabeth Simpson - 2017-10-23 ################
# Based on scripts from Will ####################

#######################
# Header ##############
#######################
library(FD)

###Community matrix, rows = sites, cols = species
data <- read.csv("~/Documents/Spatial_Scaling_Field_Project/raw_data/sites/right_hand_fork/cover.csv")
comm <- with(data, tapply(Cover, list(Site,Species), function(x) mean(x,na.rm=TRUE)))
comm[is.na(comm)] <- 0

# remove species esentially there as notes
comm <- comm[,!grepl("^[a-z]+", colnames(comm))]
comm <- comm[,!grepl("\\(|/", colnames(comm))]
comm <- comm[,!grepl("sp\\.", colnames(comm))]

colnames(comm) <- tolower(gsub(" ", "_", colnames(comm)))

###Trait matrix, rows = species, cols = traits
########### Second set of traits - all related to leaves... 16 species###########
traits <- read.csv("~/Documents/Projects/spatial_scaling/field_diversity_analysis/raw_data/trait-phylo/more_traits_2.csv")
traits$species <- tolower(gsub(" ","_", traits$species))
rownames(traits) <- traits$species
traits <- traits[,c("leaf.area", "seed.mass", "whole.plant.height"), drop=FALSE]
#remove species with NA values for traits
traits <- traits[(!is.na(traits$leaf.area) & !is.na(traits$seed.mass) & !is.na(traits$whole.plant.height)),]

#remove species without trait values
comm <- comm[, colnames(comm) %in% rownames(traits)]
comm <- comm[rowSums(comm) > 0,]

functional.div <- dbFD(traits, comm)
nbsp <- as.vector(functional.div$nbsp)
FRic <- as.vector(functional.div$FRic)
FEve <- as.vector(functional.div$FEve)
FDiv <- as.vector(functional.div$FDiv)
RaoQ <-as.vector(functional.div$RaoQ)
FDis <-as.vector(functional.div$FDis)
CWM <- functional.div$CWM

terrain <- read.csv("~/Documents/Spatial_Scaling_Field_Project/clean_data/RightFork-Site1_terrain_data_radians.csv")
terrain <- terrain[order(terrain$plot.name),] 
env <- terrain
#name rows as sites
rownames(env) <- env$plot.name
# drop out all other variables besides elev and cos.asp
env <- env[,c("elev","cos.asp"), drop=FALSE]
env

plot(env$elev[-26], CWM$leaf.carbon.content.per.leaf.dry.mass)
summary(lm(CWM$leaf.carbon.content.per.leaf.dry.mass~env$elev[-26]))

plot(env$elev[-26], CWM$leaf.nitrogen.content.per.leaf.dry.mass)
summary(lm(CWM$leaf.carbon.content.per.leaf.dry.mass ~ env$elev[-26])) 

plot(env$elev[-26], CWM$whole.plant.leaf.area.per.whole.plant.leaf.dry.mass)
summary(lm(CWM$whole.plant.leaf.area.per.whole.plant.leaf.dry.mass ~ env$elev[-26]))

plot(env$cos.asp[-26], CWM$leaf.carbon.content.per.leaf.dry.mass)
summary(lm(CWM$leaf.carbon.content.per.leaf.dry.mass ~ env$cos.asp[-26]))
model <- lm(CWM$maximum_height ~ env$cos.asp[-26])
abline(model, col = "red")

plot(env$cos.asp[-26], CWM$log_seed_mass, xlab = "Aspect, S = -1, N= 1", ylab = "Community Weight Mean, log(Seed Mass)", pch = 20)
summary(lm(CWM$log_seed_mass ~ env$cos.asp[-26]))
model <- lm(CWM$log_seed_mass ~ env$cos.asp[-26])
abline(model, col = "red")

plot(env$elev[-26], FDis)
summary(lm(FDis ~ env$elev[-26]))

plot(env$cos.asp[-26], FDis)
summary(lm(FDis ~ env$cos.asp[-26]))


plot(env$elev[-26], nbsp)
summary(lm(nbsp ~ env$elev[-26])) 

plot(env$cos.asp[-26], nbsp)
summary(lm(nbsp~env$cos.asp[-26]))

plot(env$elev[-26], FRic)
summary(lm(FRic ~ env$elev[-26])) 

plot(env$elev[-26], FRic)
summary(lm(FRic ~ env$cos.asp[-26])) 

plot(env$elev[-26], FEve)
summary(lm(FEve ~ env$elev[-26])) 

plot(env$elev[-26], FEve)
summary(lm(FEve ~ env$cos.asp[-26])) 

plot(env$elev[-26], FDiv)
summary(lm(FDiv ~ 1)) 

plot(env$elev[-26], FDiv)
summary(lm(FDiv ~ env$cos.asp[-26])) 

plot(env$elev[-26], RaoQ)
summary(lm(RaoQ ~ env$elev[-26])) 

plot(env$elev[-26], RaoQ)
summary(lm(RaoQ ~ env$cos.asp[-26])) 

