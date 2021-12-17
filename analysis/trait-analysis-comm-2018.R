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
traits <- traits[,c(3:11,14,16:17)]

#setdiff(colnames(comm),rownames(traits)) 
#setdiff(rownames(traits), colnames(comm))
# ^ have checked and figured out reasons for why all of these are dropped

###################################
# Make comparative community object
c.data <- comparative.comm(tree, comm, traits, env)

########################################
# Calculate functional diversity indicies
funct.div.c.m <- with(c.data, dbFD(data[,c(1:3,7:8)], comm, w.abun=TRUE, stand.FRic=TRUE))
ses.mpd <- .ses.mpd(c.data)
ses.mntd <- .ses.mntd(c.data)
ses.mpd.a <- .ses.mpd(c.data, abundace.weighted=TRUE)
ses.mntd.a <- .ses.mntd(c.data, abundance.weighted = TRUE)
pd <- .pd(c.data)
simps.div <- diversity(c.data$comm, index="simpson")
fdiv.c.m <- as.data.frame(funct.div.c.m$nbsp)
fdiv.c.m <- with(funct.div.c.m, cbind(nbsp, FRic, FEve,FDiv, FDis, CWM, ses.mpd.a$mpd.obs.z, ses.mntd.a$mntd.obs.z, pd[,1], simps.div))
colnames(fdiv.c.m) <- c("nsp","fric","feve","fdiv","fdis","CWM.SLA","CWM.LA", "CWM.LDMC","CWM.maxht", "CWM.mn.ht", "sesmpd","sesmntd", "pd", "simpsdiv")
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

# Subset to data to focus on
t.fd.c.m <- t.fd.c.m[,c(2:15,20:23)]
rownames(t.fd.c.m) <- fdiv.c.m$plot_id
s.fd.c.m <- s.fd.c.m[,c(2:15, 20:32)]
rownames(s.fd.c.m) <- soil.env$Plot_id

#######################################################################
# Look at how functional diversity changes across topographic variables
# 78 plots
summary(with(t.fd.c.m, lm(fric~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(feve~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(fdiv~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(fdis~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(CWM.SLA~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(CWM.LA~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(CWM.LDMC~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(CWM.mn.ht~elev+aspect+slope+TPI)))
summary(with(t.fd.c.m, lm(CWM.maxht~elev+aspect+slope+TPI)))

# Look at how functional diversity changes across phylodiversity
# 78 plots
summary(with(t.fd.c.m, lm(fric~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(feve~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(fdiv~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(fdis~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(CWM.SLA~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(CWM.LA~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(CWM.LDMC~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(CWM.mn.ht~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(t.fd.c.m, lm(CWM.maxht~nsp+sesmpd+pd+sesmntd+simpsdiv)))

# A set of models for the ones that do vary with combined diversity and topography explanatory variables
# Just continuous functional diversity across aspect
fdis.a <- with(t.fd.c.m, lm(fdis~aspect))
SLA.a <- with(t.fd.c.m, lm(log(CWM.SLA)~aspect))
LA.a <- with(t.fd.c.m, lm(log(CWM.LA)~aspect))
mnHT.a <- with(t.fd.c.m, lm(CWM.mn.ht~aspect))
mxHT.a <- with(t.fd.c.m, lm(CWM.maxht~aspect))

# F-div across phy-div. relationships
fdis.pd <- with(t.fd.c.m, lm(fdis~pd))
fdiv.mntd <- with(t.fd.c.m, lm(fdiv~sesmntd))
LA.pd <- with(t.fd.c.m, lm(log(CWM.LA)~pd))
SLA.mntd <- with(t.fd.c.m, lm(log(CWM.SLA)~sesmntd))
mnHT.pd <- with(t.fd.c.m, lm(CWM.mn.ht~pd))
mxHT.mntd <- with(t.fd.c.m, lm(CWM.maxht~sesmntd))


######################################################
### PLOT Functional Diversity variation across aspect
jpeg("./functional_traits_rhf/analysis/figures/fdiv-aspect-78plots.jpeg", width=7, height=5, unit="in",res=300)
par(mfrow=c(2,3))
par(mar=c(4,5,1,1))
par(oma=c(0.2,0.3,0.5,0.5))

# functional dispersion ~ aspect
with(t.fd.c.m, plot(fdis~aspect, pch=19, xlab="", ylab="Functional dispersion", cex=1.2, axes=FALSE, ylim=c(0,5), col="#0072B2"))
abline(fdis.a, lwd=3)
axis(1)
axis(2)
mtext("(a)",side=3,line=-1, adj=0.05)

# CWM of SLA ~ aspect
with(t.fd.c.m, plot(log(CWM.SLA)~aspect, pch=19, xlab="", ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, ylim=c(1,6), col="#009E73"))
abline(SLA.a, lwd=3)
axis(1)
axis(2)
mtext("(b)",side=3,line=-1, adj=0.05)

# CWM of LA ~ aspect
with(t.fd.c.m, plot(log(CWM.LA)~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, ylim=c(-2,6), col="#009E73"))
abline(LA.a, lwd=3)
axis(1)
axis(2)
mtext("(c)",side=3,line=-1, adj=0.05)

# CWM of mean height ~ aspect
with(t.fd.c.m, plot(CWM.mn.ht~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab="CWM - mean height (cm)", cex=1.2, axes=FALSE, ylim=c(0,160), col="#E69F00"))
abline(mnHT.a, lwd=3)
axis(1)
axis(2)
mtext("(d)",side=3,line=-1, adj=0.05)

# CWM of max height ~ aspect
with(t.fd.c.m, plot(CWM.maxht~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab="CWM - max. height (cm)", cex=1.2, axes=FALSE, ylim=c(0,200),col="#E69F00"))
abline(mxHT.a, lwd=3)
axis(1)
axis(2)
mtext("(e)",side=3,line=-1, adj=0.05)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend=c('Overall functional diversity', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1.2, cex=0.9,
       col = c('#0072B2', '#009E73', '#E69F00'))

dev.off()

### PLOT functional metrics that correlate with phylodiversity
jpeg("./functional_traits_rhf/analysis/figures/fdiv-phydiv-78plots.jpeg", width=7, height=10, unit="in",res=300)
par(mfrow=c(3,2))
par(mar=c(4,5,0.5,1))
par(oma=c(0.2,0.3,0.5,0.5))

# Functional dispersion across PD
with(t.fd.c.m, plot(fdis~pd, pch=19, xlab="", ylab="Functional dispersion", cex=1.2, axes=FALSE, xlim=c(500,2500), ylim=c(0,5), col="#0072B2", cex.lab=1.5))
abline(fdis.pd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext("(a)",side=3,line=-1, adj=0.05)

# Functional divergence across SESmntd
with(t.fd.c.m, plot(fdiv~sesmntd, pch=19, xlab="", ylab="Functional divergence", cex=1.2, axes=FALSE, xlim=c(-3,4), ylim=c(0.2,1), col="#0072B2", cex.lab=1.5))
abline(fdiv.mntd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext("(b)",side=3,line=-1, adj=0.05)

legend("bottomright", legend=c('Overall FD', 'Leaf traits', 'Height traits'), pch=16, pt.cex=1.2, cex=1.2, 
       col = c('#0072B2', '#009E73', '#E69F00'))

# CWM of leaf area across PD
with(t.fd.c.m, plot(log(CWM.LA)~pd, pch=19, xlab="", ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, xlim=c(500,2500), ylim=c(-2,6), col="#009E73", cex.lab=1.5))
abline(LA.pd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext("(c)",side=3,line=-1, adj=0.05)

# CWM of SLA across SESmntd
with(t.fd.c.m, plot(log(CWM.SLA)~sesmntd, pch=19, xlab="", ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, xlim=c(-3,4), ylim=c(1,6), col="#009E73", cex.lab=1.5))
abline(SLA.mntd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext("(d)",side=3,line=-1, adj=0.05)

# CWM mean height across PD
with(t.fd.c.m, plot(CWM.mn.ht~pd, pch=19, xlab="Faith's PD", ylab="CWM - mean height (cm)", cex=1.2, axes=FALSE, xlim=c(500,2500), ylim=c(0,160), col="#E69F00", cex.lab=1.5))
abline(mnHT.pd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext("(e)",side=3,line=-1, adj=0.05)

# CWM max. height across mntd
with(t.fd.c.m, plot(CWM.maxht~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab="CWM - max. height (cm)", cex=1.2, axes=FALSE, xlim=c(-3,4), ylim=c(0,200), col="#E69F00", cex.lab=1.5))
abline(mxHT.mntd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext("(f)",side=3,line=-1, adj=0.05)

dev.off()

#######################################################################
# Look at how functional diversity changes across topographic variables at plots with 
# SOIL TEXTURE AND TEMPERATURE DATA TOO

# First comparing these results from only 25 plots to the 78 plots above
summary(with(s.fd.c.m, lm(fric~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(feve~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(fdiv~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(fdis~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(CWM.SLA~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(CWM.LA~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(CWM.LDMC~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(CWM.mn.ht~elev+aspect+slope+TPI)))
summary(with(s.fd.c.m, lm(CWM.maxht~elev+aspect+slope+TPI)))

# Look at how functional diversity changes across phylodiversity - 25 plots
summary(with(s.fd.c.m, lm(fric~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(feve~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(fdiv~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(fdis~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(CWM.SLA~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(CWM.LA~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(CWM.LDMC~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(CWM.mn.ht~nsp+sesmpd+pd+sesmntd+simpsdiv)))
summary(with(s.fd.c.m, lm(CWM.maxht~nsp+sesmpd+pd+sesmntd+simpsdiv)))

# Looking at how functional and phylogenetic diversity vary across soil class based on texture and temperature values
summary(with(s.fd.c.m, lm(fric~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(feve~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(fdiv~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(fdis~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(CWM.SLA~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(CWM.LA~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(CWM.LDMC~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(CWM.mn.ht~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(CWM.maxht~mean+sd+min+max)))

summary(with(s.fd.c.m, lm(nsp~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(pd~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(sesmpd~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(sesmntd~mean+sd+min+max)))
summary(with(s.fd.c.m, lm(simpsdiv~mean+sd+min+max)))

### TO PLOT : Diversity ~ aspect with 25 plots
s.fric.a <- with(s.fd.c.m, lm(fric~aspect))
s.fdis.a <- with(s.fd.c.m, lm(fdis~aspect))
s.SLA.a <- with(s.fd.c.m, lm(log(CWM.SLA)~aspect))
s.LA.a <- with(s.fd.c.m, lm(log(CWM.LA)~aspect))

### TO PLOT: Functional ~ Phylo/other diversity - add these in with the existing simps div plot or if it will make more sense for most of them to be in the supplement, put them there. 
s.fdiv.mntd <- with(s.fd.c.m, lm(fdiv~sesmntd))
s.SLA.pd <- with(s.fd.c.m, lm(log(CWM.SLA)~pd))
s.mnHT.mntd <- with(s.fd.c.m, lm(CWM.mn.ht~sesmntd))
s.mnHT.simps <- with(s.fd.c.m, lm(CWM.mn.ht~simpsdiv))
s.mxHT.mntd <- with(s.fd.c.m, lm(CWM.maxht~sesmntd))

### PLOTTED: Fdiv across temperature variation
s.fric.mn <- with(s.fd.c.m, lm(fric~mean))
s.fdis.mn <- with(s.fd.c.m, lm(fdis~mean)) 
# s.fdis.cl <- with(s.fd.c.m, lm(fdis~class)) # Even though one of the soil classes shows up as significantly different in this analysis...
# If I go through and do an anova based on soil class, there is not a significant difference in dispersion between the groups
# There seems to be something going on with the sandy loa, soil, but I seem to be struggling to pinpoint what that is.
# Look at how functional diversity changes across phylodiversity

### Plots already made
jpeg("./functional_traits_rhf/analysis/figures/fdiv-temp-soil-25plots.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,2))
par(mar=c(4,5,1,1))
par(oma=c(0.2,0.3,0.5,0.5))

# functional richness ~ mean temperature
with(s.fd.c.m, plot(fric~mean, pch=19, xlab=expression(paste('Mean temperature (',degree~'C)')), ylab="Functional richness", cex=1.2, axes=FALSE, xlim=c(4,14),ylim=c(0,5), col="#0072B2"))
abline(s.fric.mn, lwd=3)
axis(1)
axis(2)
mtext("(a)",side=3,line=-1, adj=0.05)

# functional dispersion ~ mean temperature
with(s.fd.c.m, plot(fdis~mean, pch=19, xlab=expression(paste('Mean temperature (',degree~'C)')), ylab="Functional dispersion", cex=1.2, axes=FALSE, xlim=c(4,14), ylim=c(0.2,1.4),col="#0072B2"))
abline(s.fdis.mn, lwd=3)
axis(1)
axis(2)
mtext("(b)",side=3,line=-1, adj=0.05)

dev.off()

#CWM of mean height change across Simpson's Diversity
jpeg("./functional_traits_rhf/analysis/figures/mnHT-simps-25plots.jpeg", width=4, height=3.5, unit="in",res=300)
par(mfrow=c(1,1))
par(mar=c(4,5,1,1))
par(oma=c(0.2,0.3,0.5,0.5))

with(s.fd.c.m, plot(CWM.mn.ht~simpsdiv, pch=19, xlab="Simpson's Diversity", ylab="CWM - mean height (cm)", cex=1.2, axes=FALSE, ylim=c(10,90),col="#0072B2"))
abline(s.mnHT.simps, lwd=3)
axis(1)
axis(2, at=c(10,30,50,70,90))

dev.off()

########################################################
# How do soil text. and temp. vary across environment?
mn.a <- with(s.fd.c.m, lm(mean~aspect))
sd.a <- with(s.fd.c.m, lm(sd~aspect))
mx.a <- with(s.fd.c.m, lm(max~aspect))
min <- with(s.fd.c.m, lm(min~1))
rng.a <- with(s.fd.c.m, lm(range~aspect))

###
# PLOTTING
jpeg("./functional_traits_rhf/analysis/figures/temp-env-25plots.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,2))
par(mar=c(4,5,0.5,1))
par(oma=c(0.2,0.3,0.5,0.5))

with(s.fd.c.m, plot(mean~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(paste('Temperature (',degree~'C)')), cex=1.2, axes=FALSE, ylim=c(-10,70), col="#E69F00"))
abline(mn.a, lwd=3, col="#E69F00")
with(s.fd.c.m, points(max~aspect, pch=19, cex=1.2, col="#D55E00"))
abline(mx.a, lwd=3, col="#D55E00")
with(s.fd.c.m, points(min~aspect, pch=19, cex=1.2, col="#56B4E9"))
abline(min, lwd=3, col="#56B4E9")
axis(1)
axis(2, at=c(-10,10,30,50,70))
mtext("(a)",side=3,line=-.5, adj=0.05)

with(s.fd.c.m, plot(sd~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(paste('SD - Temperature (',degree~'C)')), cex=1.2, axes=FALSE, col="#E69F00", ylim=c(5,15)))
abline(sd.a, lwd=3, col="#E69F00")
axis(1)
axis(2, at=c(5,10,15))
mtext("(b)",side=3,line=-.5, adj=0.05)
legend("topright", legend=c('Min', 'Mean', 'Max'), pch=16, pt.cex=1.2, cex=0.8,
       col = c('#56B4E9', '#E69F00', '#D55E00'))

dev.off()

# # PCA - total
# pca.total <- prcomp(t.fd.c.m, scale=TRUE)
# biplot(pca.total)
# biplot(pca.total, choices=2:3)
# summary(pca.total)
# plot(pca.total)
# 
# # PCA - subset
# pca.sub <- prcomp(s.fd.c.m[,c(1:21,23:26)], scale=TRUE)
# biplot(pca.sub)
# biplot(pca.sub, choices=2:3)
# summary(pca.sub)
# plot(pca.sub)
