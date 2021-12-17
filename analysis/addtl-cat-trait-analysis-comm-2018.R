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

# Make binary traits for the first two dispersal categories
traits$disp_1 <- NA
for(i in seq_len(nrow(traits))){
  if (traits$disp_syn[i]==1){
    traits$disp_1[i] <- 1
  } else {
    traits$disp_1[i] <- 0
  }
}

traits$disp_2 <- NA
for(i in seq_len(nrow(traits))){
  if (traits$disp_syn[i]==2){
    traits$disp_2[i] <- 1
  } else {
    traits$disp_2[i] <- 0
  }
}

# clean up traits data frame
rownames(traits) <- tolower(traits$species)
traits <- traits[,c(3:5,9:10,16:19)]
# need to change these integer traits to numeric
traits$hybridization <- as.numeric(traits$hybridization)
traits$native <- as.numeric(traits$native)

#setdiff(colnames(comm),rownames(traits)) 
#setdiff(rownames(traits), colnames(comm))
# ^ have checked and figured out reasons for why all of these are dropped

###################################
# Make comparative community object
c.data <- comparative.comm(tree, comm, traits, env)

###################################
# Do your own PCA on the traits to see what is going on....
# Just the 5 traits, leaf and height
five.traits <- c.data$data[,c(1:5)]
colnames(five.traits) <- c("SLA", "LA", "LDMC", "Max. Ht.", "Mn. Ht.")
five.traits <- prcomp(five.traits, scale=TRUE)
five.traits
summary(five.traits)
plot(five.traits)

jpeg("./functional_traits_rhf/analysis/figures/pca-5traits.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,2))
par(mar=c(4,5,1,1))
par(oma=c(0.2,0.3,0.5,0.5))
biplot(five.traits, cex=c(0.4,0.4), col=c("gray", "darkred"),cex.axis=0.7, cex.lab=0.7)
biplot(five.traits, choices=2:3, cex=c(0.4,0.4), col=c("gray","darkred"), cex.axis=0.7, cex.lab=0.7)
dev.off()

# Add in 3 categorical traits for interesting test

eight.traits <- c.data$data
colnames(eight.traits) <- c("SLA", "LA", "LDMC", "Max. Ht.", "Mn. Ht.", "Hybridization","Native", "Unassisted dispersal","Wind dispersal")
eight.traits <- prcomp(eight.traits, scale=TRUE)
eight.traits
summary(eight.traits)
plot(eight.traits)

jpeg("./functional_traits_rhf/analysis/figures/pca-8traits.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,2))
par(mar=c(4,4,2,2))
par(oma=c(0.2,0.3,0.5,0.5))
biplot(eight.traits, cex=c(0.4,0.4), col=c("gray", "darkred"),cex.axis=0.5, cex.lab=0.6)
biplot(eight.traits, choices=2:3, cex=c(0.4,0.4), col=c("gray","darkred"), cex.axis=0.5, cex.lab=0.6)
dev.off()

########################################
# Calculate functional diversity indicies
funct.div.c.m <- with(c.data, dbFD(data, comm, w.abun=TRUE, stand.FRic=TRUE, CWM.type="all", print.pco=TRUE, messages=TRUE))
ses.mpd <- .ses.mpd(c.data)
ses.mntd <- .ses.mntd(c.data)
ses.mpd.a <- .ses.mpd(c.data, abundace.weighted=TRUE)
ses.mntd.a <- .ses.mntd(c.data, abundance.weighted = TRUE)
pd <- .pd(c.data)
simps.div <- diversity(c.data$comm, index="simpson")
fdiv.c.m <- as.data.frame(funct.div.c.m$nbsp)
fdiv.c.m <- with(funct.div.c.m, cbind(nbsp, FRic, FEve,FDiv, FDis, CWM, ses.mpd.a$mpd.obs.z, ses.mntd.a$mntd.obs.z, pd[,1], simps.div))
fdiv.c.m <- fdiv.c.m[,c(1:10,12,14,16,18:22)]
colnames(fdiv.c.m) <- c("nsp","fric","feve","fdiv","fdis","CWM.SLA","CWM.LA", "CWM.LDMC","CWM.maxht", "CWM.mn.ht", "p.hybrid","p.nat", "p.disp1","p.disp2", "sesmpd", "sesmntd", "pd", "simpsdiv")
fdiv.c.m$plot_id <- rownames(fdiv.c.m)

#################################################
# Add in soil environment - texture & temperature
soil.env <- read.csv("./functional_traits_rhf/clean_data/2018-soil-env.csv", as.is=TRUE)

###########################################################################################################
# Make subset of environmental and functional data to analyze changes in functional div. across temperature
# 78 plots of functional and topographic env data
t.fd.c.m <- merge(fdiv.c.m, env, by.x="plot_id", by.y="plot_id")

# 25 plots of functional, topographic, soil texture and temp. data
s.fd.c.m <- merge(t.fd.c.m, soil.env[,2:11], by.x="plot_id", by.y="Plot_id")

rownames(t.fd.c.m) <- fdiv.c.m$plot_id
t.fd.c.m <- t.fd.c.m[,-1]

rownames(s.fd.c.m) <- soil.env$Plot_id
s.fd.c.m <- s.fd.c.m[,-1]

#######################################################################
# Look at how functional diversity changes across topographic variables
# 78 plots
summary(with(t.fd.c.m, lm(fric~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(feve~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(fdiv~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(fdis~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(p.hybrid~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(p.nat~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(p.disp1~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(p.disp2~elev+aspect+slope+TPI)))

# Look at how functional diversity changes across phylodiversity
# 78 plots
summary(with(t.fd.c.m, lm(fric~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(feve~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(fdiv~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(fdis~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(p.hybrid~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(p.nat~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(p.disp1~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(p.disp2~nsp+sesmpd+pd+sesmntd+simpsdiv)))


# A set of models for the ones that do vary with combined diversity and topography explanatory variables
# Just continuous functional diversity across aspect
fric.a <- with(t.fd.c.m, lm(fric~aspect))
nat.a <- with(t.fd.c.m, lm(p.nat~aspect))
disp1.a <- with(t.fd.c.m, lm(p.disp1~aspect))

# F-div across phy-div. relationships
fric.pd <- with(t.fd.c.m, lm(fric~pd))
fdis.simps <- with(t.fd.c.m, lm(fdis~simpsdiv))
hybrid.mntd <- with(t.fd.c.m, lm(p.hybrid~sesmntd))
nat.mpd <- with(t.fd.c.m, lm(p.nat~sesmpd))
nat.mntd <- with(t.fd.c.m, lm(p.nat~sesmntd))

######################################################
### PLOT Functional Diversity variation across aspect
jpeg("./functional_traits_rhf/analysis/figures/fdiv8traits-aspect-78plots.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,3))
par(mar=c(4,5,1,1))
par(oma=c(0.2,0.3,0.5,0.5))

# functional richness ~ aspect
with(t.fd.c.m, plot(fric~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab="Functional richness", ylim=c(0,0.3),cex=1.2, axes=FALSE, col="#0072B2"))
abline(fric.a, lwd=3)
axis(1)
axis(2)
mtext("(a)",side=3,line=-1, adj=0.05)

# proportion native ~ aspect
with(t.fd.c.m, plot(p.nat~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab= "Proportion native species", cex=1.2, axes=FALSE, ylim=c(0,1),col="#CC79A7"))
abline(nat.a, lwd=3)
axis(1)
axis(2)
mtext("(b)",side=3,line=-1, adj=0.05)

legend("bottomright", legend=c('Overall functional diversity', 'Categorical traits'), pch=19, pt.cex=1.2, cex=0.7, bty="n",
       col = c('#0072B2', '#CC79A7'))

# Proportion unassisted dispersal ~ aspect
with(t.fd.c.m, plot(p.disp1~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab= "Proportion of unassisted dispersal", cex=1.2,axes=FALSE,col="#CC79A7",ylim=c(0,1)))
abline(disp1.a, lwd=3)
axis(1)
axis(2)
mtext("(c)",side=3,line=-1, adj=0.05)

dev.off()

### PLOT functional metrics that correlate with phylodiversity
jpeg("./functional_traits_rhf/analysis/figures/fdiv8traits-phydiv-78plots.jpeg", width=7, height=10, unit="in",res=300)
par(mfrow=c(3,2))
par(mar=c(4,5,0.5,1))
par(oma=c(0.2,0.3,0.5,0.5))

# Functional richness across pd
with(t.fd.c.m, plot(fric~pd, pch=19, xlab="Faith's PD", ylab="Functional richness", cex=1.2, axes=FALSE, col="#0072B2", cex.lab=1.5, xlim=c(500,2500),ylim=c(0,0.6)))
abline(fric.pd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext("(a)",side=3,line=-1, adj=0.05)

# Proportion of native species across mpd
with(t.fd.c.m, plot(p.nat~sesmpd, pch=19, xlab=expression(SES[MPD]), ylab="Proportion of native species", cex=1.2, axes=FALSE, col="#CC79A7", cex.lab=1.5, xlim=c(-2,2), ylim=c(0.2,1)))
abline(nat.mpd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext("(b)",side=3,line=-1, adj=0.05)

# Functional dispersion across Simpson's Diversity
with(t.fd.c.m, plot(fdis~simpsdiv, pch=19, xlab="Simpson's Diversity", ylab="Functional dispersion", cex=1.2, axes=FALSE, col="#0072B2", cex.lab=1.5, xlim=c(0.2,0.9), ylim=c(0,5)))
abline(fdis.simps, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext("(c)",side=3,line=-1, adj=0.05)

# Proportion of native species across mntd
with(t.fd.c.m, plot(p.nat~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab="Proportion of native species", cex=1.2, axes=FALSE, col="#CC79A7", cex.lab=1.5, xlim=c(-3,5), ylim=c(0.2,1)))
abline(nat.mntd, lwd=3)
axis(1, cex.axis=1.5, at=c(-3,-1,1,3,5))
axis(2, cex.axis=1.5)
mtext("(d)",side=3,line=-1, adj=0.05)

legend("bottomright", legend=c('Overall functional diversity', 'Categorical traits'), pch=19, pt.cex=1.2, cex=0.9, bty="n",
       col = c('#0072B2', '#CC79A7'))

# Proportion that tend to hybridize across mntd
with(t.fd.c.m, plot(p.hybrid~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab="Proportion species known to hybridize", cex=1.2, axes=FALSE, col="#CC79A7", cex.lab=1.5, xlim=c(-3,5), ylim=c(0,1)))
abline(hybrid.mntd, lwd=3)
axis(1, cex.axis=1.5, at=c(-3,-1,1,3,5))
axis(2, cex.axis=1.5)
mtext("(e)",side=3,line=-1, adj=0.05)

dev.off()

#######################################################################
# Look at how functional diversity changes across topographic variables at plots with 
# SOIL TEXTURE AND TEMPERATURE DATA TOO

# First comparing these results from only 25 plots to the 78 plots above
summary(with(s.fd.c.m, lm(fric~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(feve~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(fdiv~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(fdis~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(fdis~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(p.hybrid~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(p.nat~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(p.disp1~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(p.disp2~elev+aspect+slope+TPI)))

# Look at how functional diversity changes across phylodiversity - 25 plots
summary(with(s.fd.c.m, lm(fric~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(feve~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(fdiv~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(fdis~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(p.hybrid~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(p.nat~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(p.disp1~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(p.disp2~nsp+sesmpd+pd+sesmntd+simpsdiv)))

# Looking at how functional and phylogenetic diversity vary across soil class based on texture and temperature values
summary(with(s.fd.c.m, lm(fric~class+mean+sd+min+max)))
summary(with(s.fd.c.m, lm(feve~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(fdiv~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(fdis~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(p.hybrid~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(p.nat~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(p.disp1~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(p.disp2~mean+sd+min+max)))

# MODELS TO PLOT
s.fric.mn <- with(s.fd.c.m, lm(fric~mean))

### PLOT THESE!!! ####
jpeg("./functional_traits_rhf/analysis/figures/fdiv8traits-temp-soil-25plots.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,1))
par(mar=c(4,5,1,1))
par(oma=c(0.2,0.3,0.5,0.5))

# functional richness ~ mean temperature
with(s.fd.c.m, plot(fric~mean, pch=19, xlab=expression(paste('Mean temperature (',degree~'C)')), ylab="Functional richness", cex=1.2, axes=FALSE, xlim=c(4,14),ylim=c(0,5), col="#0072B2"))
abline(s.fric.mn, lwd=3)
axis(1)
axis(2)
#mtext("(a)",side=3,line=-1, adj=0.05)

dev.off()