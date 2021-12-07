# Using commparative.comm to put together all trait, phylo, and environmental data from RHF
# Elizabeth Simpson # 2021-03-25 # Uses some code that Will Pearse oringinally gave me.
setwd("~/Documents/projects/") # uncomment to set working directory on computer

###########
# HEADERS #
###########
library(pez)
library(FD)
library(lme4)

#####################################################################################
# Community data 
# load in cover data, reformat, remove NAs, remove species essentially there as notes
data <- read.csv("./fractal_sampling_div_rhf/data/rhf_2018_cover.csv", as.is=TRUE)
comm <- with(data, tapply(Cover, list(Plot_id,Species), function(x) mean(x,na.rm=TRUE)))
comm[is.na(comm)] <- 0
comm <- comm[, colSums(comm != 0) > 0]
comm <- comm[,!grepl("^[a-z]+", colnames(comm))]
comm <- comm[,!grepl("\\(|/", colnames(comm))]
comm <- comm[,!grepl("sp\\.", colnames(comm))]
colnames(comm) <- tolower(colnames(comm))

###########################################################
# Phylogeny (Zanne et. al. 2014) 
# load tree, reformat sp names, create combined data object
tree <- read.tree("./fractal_sampling_div_rhf/data/Vascular_Plants_rooted.dated.tre")
tree$tip.label <- tolower(gsub("_", " ", tree$tip.label))
tree <- congeneric.merge(tree, colnames(comm), split = " ")

#############################################################################################
# Environmental data - just topographic/terrain variables, because temp. data only at 25 plot
# And I want to maintain information from all 78 plots in the analysis as long as possible
### For example, can look at how functional diversity varies across terrain with all 78 plots. 
### l__> Does temperature map onto this?
env <- read.csv("./fractal_sampling_div_rhf/data/clean_rhf_terrain.csv", as.is=TRUE)
env <- env[,2:10]
rownames(env) <- env$plot_id

#######################################################################################################
# Trait data - includes leaf & height info as well as some categorical stuff
# Just analyzing continuous variables for now because I struggled with integrating the categorical ones
traits <- read.csv("./functional_traits_rhf/clean_data/clean_traits_rhf.csv", as.is=TRUE)
rownames(traits) <- tolower(traits$species)
traits <- traits[,c(3:5,9:10)]


#setdiff(colnames(comm),rownames(traits)) 
#setdiff(rownames(traits), colnames(comm))
# ^ have checked and figured out reasons for why all of these are dropped

###################################
# Make comparative community object
c.data <- comparative.comm(tree, comm, traits, env)

########################################
# Calculate functional diversity indicies
functional.div <- with(c.data, dbFD(data[,c(1:5)], comm, w.abun=TRUE))
fdiv <- as.data.frame(functional.div$nbsp)
fdiv <- with(functional.div, cbind(nbsp, FRic, FEve,FDiv, FDis, CWM))
colnames(fdiv) <- c("nsp","fric","feve","fdiv","fdis","CWM.SLA","CWM.LA", "CWM.LDMC","CWM.maxht", "CWM.mn.ht")
fdiv$plot_id <- rownames(fdiv)

#################################################
# Add in soil environment - texture & temperature
soil.env <- read.csv("./functional_traits_rhf/clean_data/2018-soil-env.csv", as.is=TRUE)

###########################################################################################################
# Make subset of environmental and functional data to analyze changes in functional div. across temperature
# 78 plots of functional and topographic env data
total.de <- merge(fdiv, env, by.x="plot_id", by.y="plot_id")

# 25 plots of functional, topographic, soil texture and temp. data
sub.de <- merge(total.de, soil.env[,2:11], by.x="plot_id", by.y="Plot_id")

# Subset to data to focus on
total.de <- total.de[,c(2:11,16:19)]
rownames(total.de) <- fdiv$plot_id
sub.de <- sub.de[,c(2:11, 16:28)]
rownames(sub.de) <- soil.env$Plot_id

#######################################################################
# Look at how functional diversity changes across topographic variables
# 78 plots
summary(with(total.de, lm(nsp~elev+aspect+slope+TPI))) # aspect (less species than full 2018 analysis)
summary(with(total.de, lm(fric~elev+aspect+slope+TPI)))
summary(with(total.de, lm(feve~elev+aspect+slope+TPI))) # aspect
summary(with(total.de, lm(fdiv~elev+aspect+slope+TPI))) # aspect
summary(with(total.de, lm(fdis~elev+aspect+slope+TPI))) # aspect
summary(with(total.de, lm(CWM.SLA~elev+aspect+slope+TPI))) # elev, aspect, slope
summary(with(total.de, lm(CWM.LA~elev+aspect+slope+TPI))) # aspect
summary(with(total.de, lm(CWM.LDMC~elev+aspect+slope+TPI)))
summary(with(total.de, lm(CWM.mn.ht~elev+aspect+slope+TPI))) # aspect
summary(with(total.de, lm(CWM.maxht~elev+aspect+slope+TPI))) # aspect

## Plot ones significant across aspect & slope
par(mfrow=c(2,2))
par(mar=c(4,4,1,1))

# Functional Evenness ~ aspect
with(total.de, plot(feve~aspect, pch=19, xlab="", ylab="Functional Evenness", cex=1.5))
abline(with(total.de, lm(feve~aspect)), lwd=3, col="darkred")

# Functional Divergence ~ aspect
with(total.de, plot(fdiv~aspect, pch=19, xlab="", ylab="Functional Divergence", cex=1.5))
abline(with(total.de, lm(fdiv~aspect)), lwd=3, col="darkred")

# Functional Dispersion ~ aspect
with(total.de, plot(fdis~aspect, pch=19, xlab="Aspect(S = -1, N = 1)", ylab="Functional Dispersion", cex=1.5))
abline(with(total.de, lm(fdis~aspect)), lwd=3, col="darkred")

# Species Richness ~ aspect
with(total.de, plot(nsp~aspect, pch=19, xlab="Aspect(S = -1, N = 1)", ylab="Species Richness", cex=1.5))
abline(with(total.de, lm(nsp~aspect)), lwd=3, col="darkred")

###################################################################################
par(mfrow=c(1,3))
par(mar=c(4,4,1,1))

# CWM of SLA ~ elevation
with(total.de, plot(CWM.SLA~aspect, type="n", xlab="Aspect(S = -1, N = 1)", ylab="Community Weighted Mean of SLA", cex=1.5))
with(total.de, text(CWM.SLA~aspect, labels=rownames(total.de), pos = 2))

abline(with(total.de, lm(CWM.SLA~aspect)), lwd=3, col="darkred")

#########################################
#########################################
# PAUSING TO TACKLE THE LEAF PROBLEM :D

# CWM of SLA ~ slope
with(total.de, plot(CWM.SLA~elev, pch=19, xlab="Elevation (m)", ylab="", cex=1.5))
abline(with(total.de, lm(CWM.SLA~elev)), lwd=3, col="darkred")

# CWM of SLA ~ aspect
with(total.de, plot(CWM.SLA~slope, pch=19, xlab="Slope (degrees)", ylab="", cex=1.5))
abline(with(total.de, lm(CWM.SLA~slope)), lwd=3, col="darkred")

# CWM of LA ~ aspect
with(total.de, plot(CWM.LA~aspect, type="n", xlab="Aspect(S = -1, N = 1)", ylab="Community Weighted Mean of LA", cex=1.5))
with(total.de, text(CWM.LA~aspect, labels=rownames(total.de), pos = 2))

#####################################
#####################################
# Height data is allllll goood :D ###

# CWM of mean height ~ aspect
with(total.de, plot(CWM.mn.ht~aspect, type="n", xlab="Aspect(S = -1, N = 1)", ylab="Community Weighted Mean of SLA", cex=1.5))
with(total.de, text(CWM.mn.ht~aspect, labels=rownames(total.de), pos = 2))

# CWM of max height ~ aspect
with(total.de, plot(CWM.maxht~aspect, type="n", xlab="Aspect(S = -1, N = 1)", ylab="Community Weighted Mean of SLA", cex=1.5))
with(total.de, text(CWM.maxht~aspect, labels=rownames(total.de), pos = 2))

#############################################
##############################################










with(data, plot(fric~aspect, pch=19, xlab="", ylab="Functional Richness", col="darkcyan",cex=1.5, cex.axis=1.5, cex.lab=1.5))
abline(with(data, lm(fric~aspect)),lwd=1.5)
with(data, plot(fric~slope, pch=19, xlab="", ylab="", col="darkcyan",cex=1.5, cex.axis=1.5, cex.lab=1.5))
abline(with(data, lm(fric~slope)),lwd=1.5)

with(data, plot(feve~aspect, pch=19, xlab="", ylab="Functional Evenness", col="darkcyan",cex=1.5, cex.axis=1.5, cex.lab=1.5))
abline(with(data, lm(feve~aspect)),lwd=1.5)
with(data, plot(feve~slope, pch=19, xlab="", ylab="", col="darkcyan",cex=1.5, cex.axis=1.5, cex.lab=1.5))
abline(with(data, lm(feve~slope)),lwd=1.5)

with(data, plot(CWM.ht~aspect, pch=19, xlab="",ylab="CWM - Height", col="darkcyan",cex=1.5, cex.axis=1.5, cex.lab=1.5))
abline(with(data, lm(CWM.ht~aspect)),lwd=1.5)
with(data, plot(CWM.ht~slope, pch=19, xlab="", ylab="", col="darkcyan",cex=1.5, cex.axis=1.5, cex.lab=1.5))
abline(with(data, lm(CWM.ht~slope)),lwd=1.5)

with(data, plot(CWM.maxht~aspect, pch=19, xlab="Aspect, S=-1, N=1", ylab="CWM - Max. Height", col="darkcyan",cex=1.5, cex.axis=1.5, cex.lab=1.5))
abline(with(data, lm(CWM.maxht~aspect)),lwd=1.5)
with(data, plot(CWM.maxht~slope, pch=19, xlab=expression("Slope ("*~degree*")"), ylab="", col="darkcyan",cex=1.5, cex.axis=1.5, cex.lab=1.5))
abline(with(data, lm(CWM.maxht~slope)),lwd=1.5)

par(mfrow=c(3,1))
with(data, plot(fdis~aspect, pch=19, xlab="", ylab="Functional Dispersion", col="darkcyan",cex=1.5, cex.axis=1.5, cex.lab=1.5))
abline(with(data, lm(fdis~aspect)),lwd=1.5)

with(data, plot(CWM.SLA~aspect, pch=19, xlab="", ylab="CWM - SLA", col="darkcyan",cex=1.5, cex.axis=1.5, cex.lab=1.5))
abline(with(data, lm(CWM.SLA~aspect)),lwd=1.5)

with(data, plot(CWM.LA~aspect, pch=19, xlab="Aspect, S=-1, N=1", ylab="CWM - LA", col="darkcyan",cex=1.5, cex.axis=1.5, cex.lab=1.5))
abline(with(data, lm(CWM.LA~aspect)),lwd=1.5)

summary(with(data, lm(fdis~elev+aspect+slope+tpi))) #sig across aspect
summary(with(data, lm(CWM.SLA~elev+aspect+slope+tpi))) ### sig across aspect
summary(with(data, lm(CWM.LA~elev+aspect+slope+tpi))) ### sig across aspect


#### soil env - temperature
summary(with(data, lm(nsp~mean+sd)))
summary(with(data, lm(fric~mean+sd))) 
summary(with(data, lm(feve~mean+sd))) 
summary(with(data, lm(fdiv~mean+sd))) 
summary(with(data, lm(fdis~mean+sd))) # mean sig
summary(with(data, lm(CWM.SLA~mean+sd))) 
summary(with(data, lm(CWM.LA~mean+sd)))
summary(with(data, lm(CWM.LDMC~mean+sd)))
summary(with(data, lm(CWM.ht~mean+sd))) 
summary(with(data, lm(CWM.maxht~mean+sd))) 

par(mar=c(5,5,1,1))
with(data, plot(fdis~mean, pch=19, xlab=expression("Mean Soil Temperature ("*~degree*C*")"),ylab="Functional Dispersion", col="darkcyan",cex=1.5, cex.axis=1.5, cex.lab=1.5))
abline(with(data, lm(fdis~mean)),lwd=1.5)


### soil env - texture
summary(with(data, lm(nsp~class))) 
summary(with(data, lm(fric~class))) 
summary(with(data, lm(feve~class))) 
summary(with(data, lm(fdiv~class))) 
summary(with(data, lm(fdis~class))) 
summary(with(data, lm(CWM.SLA~class))) 
summary(with(data, lm(CWM.LA~class)))
summary(with(data, lm(CWM.LDMC~class)))
summary(with(data, lm(CWM.ht~class))) 
summary(with(data, lm(CWM.maxht~class))) 

