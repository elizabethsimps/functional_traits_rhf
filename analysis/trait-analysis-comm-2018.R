# Assessing variance in functional dispersion and CWM of height and leaf traits across temp and soil microclimates and phylodiversity
# Elizabeth Simpson # 2021-03-25
# Uses code from Will Pearse to get comparative.comm set up

setwd("~/Documents/projects/") # uncomment to set working directory on computer

###########
# HEADERS 
library(pez)
library(FD)
library(lme4)
library(soiltexture)

################
# COMMUNITY DATA
# load in cover data, reformat, remove NAs, remove species essentially there as notes
cover <- read.csv("./fractal_sampling_div_rhf/data/rhf_2018_cover.csv", as.is=TRUE)
comm <- with(cover, tapply(Cover, list(Plot_id,Species), function(x) mean(x,na.rm=TRUE)))
comm[is.na(comm)] <- 0
comm <- comm[, colSums(comm != 0) > 0]
comm <- comm[,!grepl("^[a-z]+", colnames(comm))]
comm <- comm[,!grepl("\\(|/", colnames(comm))]
comm <- comm[,!grepl("sp\\.", colnames(comm))]
colnames(comm) <- tolower(colnames(comm))

################################
# PHYLOGENY (Zanne et. al. 2014) 
# load tree, reformat sp names, create combined data object
tree <- read.tree("./fractal_sampling_div_rhf/data/Vascular_Plants_rooted.dated.tre")
tree$tip.label <- tolower(gsub("_", " ", tree$tip.label))
tree <- congeneric.merge(tree, colnames(comm), split = " ")

#################################
# ENVIRONMENTAL DATA - TOPOGRAPHY
# start with topographic/terrain variables to maintain information from 78 plots as long as possible.
env <- read.csv("./fractal_sampling_div_rhf/data/clean_rhf_terrain.csv", as.is=TRUE)
env <- env[,2:10]
rownames(env) <- env$plot_id

############
# TRAIT DATA 
# Focusing analysis on continuous variables and means for each leaf and height traits
# Keeping dispersal syndrome and native vs. not for a supplementary analysis
traits <- read.csv("./functional_traits_rhf/clean_data/clean_traits_rhf.csv", as.is=TRUE)
rownames(traits) <- tolower(traits$species)
traits <- traits[,c(3:11,14,17)]

###################################
# Make comparative community object
c.data <- comparative.comm(tree, comm, traits, env)

#####################################################################
# Calculate functional, phylogenetic, and taxonomic diversity indices
# Focus on functional dispersion and community weighted means but the 3 Villeger et al. 2008 indicies are also in here
funct.div.c.m <- with(c.data, dbFD(data[,c(1:3,7:8)], comm, w.abun=TRUE, stand.x=TRUE))

ses.mpd.a <- .ses.mpd(c.data, abundace.weighted=TRUE)
ses.mntd.a <- .ses.mntd(c.data, abundance.weighted=TRUE)
pd <- .pd(c.data)

simps.div <- diversity(c.data$comm, index="simpson")

fdiv.c.m <- as.data.frame(funct.div.c.m$nbsp)
fdiv.c.m <- with(funct.div.c.m, cbind(nbsp, FRic, FEve,FDiv, FDis, CWM, ses.mpd.a$mpd.obs.z, ses.mntd.a$mntd.obs.z, pd[,1], simps.div))
colnames(fdiv.c.m) <- c("nsp","fric","feve","fdiv","fdis","CWM.SLA","CWM.LA", "CWM.LDMC","CWM.maxht", "CWM.mn.ht", "sesmpd","sesmntd", "pd", "simpsdiv")
fdiv.c.m$plot_id <- rownames(fdiv.c.m)

################################################################################
# ENVIRONMENTAL DATA - soil temperature, texture, and other plot-level variables
soil.env <- read.csv("./functional_traits_rhf/clean_data/2018-soil-env.csv", as.is=TRUE)

# extract soil, organic matter, and rock from the dataset
plot.ch <- cover[cover$Species %in% c("soil","organic matter", "rock", "poo"),]

#calculate the mean cover of each of these items for each plot, to use plot level indicators in supplementary analysis
plot.ch <- with(plot.ch, tapply(Cover, list(Plot_id,Species), function(x) mean(x,na.rm=TRUE)))
plot.ch[is.na(plot.ch)] <- 0
colnames(plot.ch) <- c("litter", "rock", "soil", "scat")

add.plots <- rbind(c(0,0,0,0), c(0,0,0,0), c(0,0,0,0), c(0,0,0,0))
rownames(add.plots) <- c("1133", "3323", "3113", "2331")
colnames(add.plots) <- c("litter", "rock", "soil","scat")
plot.ch <- rbind(plot.ch, add.plots)
t.plot.ch <- as.data.frame(plot.ch)
t.plot.ch$plot_id <-rownames(plot.ch)

###########################################################################################################
# Make subset of environmental and functional data to analyze changes in functional div. across temperature
# 78 plots of functional and topographic env data
t.fd.c.m <- merge(fdiv.c.m, env, by.x="plot_id", by.y="plot_id")
t.fd.c.m <- merge(t.fd.c.m, t.plot.ch, by.x="plot_id", by.y="plot_id")

# 25 plots of functional, topographic, soil texture and temp. data
s.fd.c.m <- merge(t.fd.c.m, soil.env[,2:11], by.x="plot_id", by.y="Plot_id")

rownames(t.fd.c.m) <- fdiv.c.m$plot_id
rownames(s.fd.c.m) <- s.fd.c.m$plot_id

###################################################################
### How does temperature and soil texture change across environment?
# First, look at the full models
# ...for temperature
summary(with(s.fd.c.m, lm(mean~elev+aspect+slope)))
summary(with(s.fd.c.m, lm(min~elev+aspect+slope)))
summary(with(s.fd.c.m, lm(max~elev+aspect+slope)))
summary(with(s.fd.c.m, lm(sd~elev+aspect+slope)))

# ...for soil texture percentages
summary(with(s.fd.c.m, lm(SAND~elev+aspect+slope)))
summary(with(s.fd.c.m, lm(SILT~elev+aspect+slope)))
summary(with(s.fd.c.m, lm(CLAY~elev+aspect+slope)))

# Models for plotting
mn.a <- with(s.fd.c.m, lm(mean~aspect))
min <- with(s.fd.c.m, lm(min~1))
mx.a <- with(s.fd.c.m, lm(max~aspect))
sd.a <- with(s.fd.c.m, lm(sd~aspect))

# FIG 2 - Ploting changes in temperature and texture across environment
# plus the USDA texture triangle for all of the plots 
jpeg("./functional_traits_rhf/analysis/figures/temp-env-25plots.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,2))
par(mar=c(5,5,1.5,1.5))
par(oma=c(0.2,0.3,0.5,0.5))

with(s.fd.c.m, plot(mean~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(paste('Temperature (',degree~'C)')), cex=1.2, axes=FALSE, ylim=c(-10,70), col="#E69F00"))
abline(mn.a, lwd=3, col="#E69F00")
with(s.fd.c.m, points(max~aspect, pch=19, cex=1.2, col="#D55E00"))
abline(mx.a, lwd=3, col="#D55E00")
with(s.fd.c.m, points(min~aspect, pch=19, cex=1.2, col="#56B4E9"))
abline(min, lwd=3, col="#56B4E9")
axis(1)
axis(2, at=c(-10,10,30,50,70))
legend("topright", c(expression(~ R^2 ~ "= 0.35, p-value = 0.002"), expression(~ R^2 ~ "= 0.68, p-value < 0.001")), 
       text.col=c("#D55E00", "#E69F00"), box.lty=0, cex=0.6, bg="transparent")
mtext("(a)",side=3,line=-.5, adj=0.05)

with(s.fd.c.m, plot(sd~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(paste('SD - Temperature (',degree~'C)')), cex=1.2, axes=FALSE, col="#E69F00", ylim=c(5,15)))
abline(sd.a, lwd=3, col="#E69F00")
axis(1)
axis(2, at=c(5,10,15))
mtext("(b)",side=3,line=-.5, adj=0.05)
legend(-1,7.4, legend=c('Min', 'Mean', 'Max'), pch=16, pt.cex=1.2, cex=0.6,
       col = c('#56B4E9', '#E69F00', '#D55E00'))
legend("topright", c(expression(~ R^2 ~ "= 0.44, p-value < 0.001")), 
       text.col=c("#E69F00"), box.lty=0, cex=0.6, bg="transparent")

dev.off()

####################################################################################
### How do functional does functional dispersion and the community weighted mean of the traits vary across temperature and texture?
# For each metric, tested variation across mean, min, max, and sd of temperature.
# Only significant relationships shown
s.fdis.mn <- with(s.fd.c.m, lm(fdis~mean))
s.fdis.mx <- with(s.fd.c.m, lm(fdis~max)) #S
s.fdis.sd <- with(s.fd.c.m, lm(fdis~sd)) #S

s.SLA.mn <- with(s.fd.c.m, lm(log(CWM.SLA)~mean)) 
s.SLA.mx <- with(s.fd.c.m, lm(log(CWM.SLA)~max)) #S
s.SLA.sd <- with(s.fd.c.m, lm(log(CWM.SLA)~sd)) #S

s.LA.mn <- with(s.fd.c.m, lm(log(CWM.LA)~mean))
s.LA.mx <- with(s.fd.c.m, lm(log(CWM.LA)~max)) #S
s.LA.sd <- with(s.fd.c.m, lm(log(CWM.LA)~sd)) #S

s.mnH.mn <- with(s.fd.c.m, lm(CWM.mn.ht~mean))
s.mnH.mx <- with(s.fd.c.m, lm(CWM.mn.ht~max)) #S
s.mnH.sd <- with(s.fd.c.m, lm(CWM.mn.ht~sd)) #S

s.mxH.mn <- with(s.fd.c.m, lm(CWM.maxht~mean))
s.mxH.mx <- with(s.fd.c.m, lm(CWM.maxht~max)) #S
s.mxH.sd <- with(s.fd.c.m, lm(CWM.maxht~sd)) #S

s.mnH.C<- with(s.fd.c.m, lm(CWM.mn.ht~CLAY)) #S
s.mxH.C <- with(s.fd.c.m, lm(CWM.maxht~CLAY)) #S

### FIG 3 - Plotting how functional diversity changes across the MEAN temperature variable.
jpeg("./functional_traits_rhf/analysis/figures/fdiv-MEAN-25plots.jpeg", width=7, height=5, unit="in",res=300)
par(mfrow=c(2,3))
par(mar=c(4,4,0.1,0.1))
par(oma=c(1,1,1,1))

# functional dispersion ~ mean temperature
with(s.fd.c.m, plot(fdis~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab="Functional dispersion", cex=1.2, axes=FALSE, xlim=c(4,14), ylim=c(0.2,1.4),col="#0072B2", cex.lab=0.7))
abline(s.fdis.mn, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(a)"~~~ R^2 ~ "= 0.56, p-value < 0.001"),side=3,line=-1, adj=0.15,cex=0.6)

# SLA ~ mean
with(s.fd.c.m, plot(log(CWM.SLA)~mean, pch=19, xlab="", ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(4,14), ylim=c(2,5.5), cex.lab=0.7))
abline(s.SLA.mn, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(b)"~~~ R^2 ~ "= 0.37, p-value = 0.001"),side=3,line=-1, adj=0.15,cex=0.6)

# mn ht ~ mean
with(s.fd.c.m, plot(CWM.mn.ht~mean, pch=19, xlab="", ylab="CWM - mean height (cm)", cex=1.2, axes=FALSE, xlim=c(4,14), ylim=c(20,100), col="#E69F00", cex.lab=0.7))
abline(s.mnH.mn, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(c)"~~~ R^2 ~ "= 0.24, p-value = 0.013"),side=3,line=-1, adj=0.15,cex=0.6)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(0.3,0.9, legend=c('Overall functional diversity', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1.2, cex=0.6,
       col = c('#0072B2', '#009E73', '#E69F00'))

# LA ~ mean
with(s.fd.c.m, plot(log(CWM.LA)~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(4,14), ylim=c(-2,4),cex.lab=0.7))
abline(s.LA.mn, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(d)"~~~ R^2 ~ "= 0.63, p-value < 0.001"),side=3,line=-1, adj=0.15,cex=0.6)

# max ht ~ mean
with(s.fd.c.m, plot(CWM.maxht~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab="CWM - max. height (cm)", cex=1.2, axes=FALSE, xlim=c(4,14), ylim=c(40,160), col="#E69F00", cex.lab=0.7))
abline(s.mxH.mn, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(e)"~~~ R^2 ~ "= 0.28, p-value = 0.006"),side=3,line=-1, adj=0.15,cex=0.6)

dev.off()

#########################################################################
# How does functional diversity (dispersion and CWM) correlate with phylodiversity? (25 plots)

fdis.mntd <- with(s.fd.c.m, lm(fdis~sesmntd))

LA.pd <- with(s.fd.c.m, lm(log(CWM.LA)~pd))
LA.mntd <- with(s.fd.c.m, lm(log(CWM.LA)~sesmntd))

SLA.pd <- with(s.fd.c.m, lm(log(CWM.SLA)~pd))
SLA.mntd <- with(s.fd.c.m, lm(log(CWM.SLA)~sesmntd))
SLA.mpd <- with(s.fd.c.m, lm(log(CWM.SLA)~sesmpd))

mnht.simp <- with(s.fd.c.m, lm(CWM.mn.ht~simpsdiv))

mxht.mntd <- with(s.fd.c.m, lm(CWM.maxht~sesmntd))

#################################################################
# FIG 4 - Functional diversity across phy-div metrics at 25 plots
# Plot significant relationships with SESmntd here, do the rest in the supplement
jpeg("./functional_traits_rhf/analysis/figures/fdiv-phydiv-25plots.jpeg", width=7, height=6.5, unit="in",res=300)
par(mfrow=c(2,2))
par(mar=c(4,5,0.5,1))
par(oma=c(0.2,0.3,0.5,0.5))

# Functional dispersion across mntd
with(s.fd.c.m, plot(fdis~sesmntd, pch=19, xlab="", ylab="Functional dispersion", cex=1.2, axes=FALSE, col="#0072B2", xlim=c(-3,2.3), cex.lab=1.5))
abline(fdis.mntd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(a)"~~~ R^2 ~ "= 0.42, p-value < 0.001"),side=3,line=-.5, adj=0.05,cex=0.9)

legend("bottomright", legend=c('Overall FD', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1, cex=.8, 
       col = c('#0072B2', '#009E73', '#E69F00'))

# CWM of SLA across mntd
with(s.fd.c.m, plot(log(CWM.SLA)~sesmntd, pch=19, xlab="", ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, xlim=c(-3,2.3), ylim=c(2,5), col="#009E73", cex.lab=1.5))
abline(SLA.mntd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(b)"~~~ R^2 ~ "= 0.29, p-value = 0.006"),side=3,line=-.5, adj=0.05,cex=0.9)

# CWM max. height across mntd
with(s.fd.c.m, plot(CWM.maxht~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab="CWM - max. height (cm)", cex=1.2, axes=FALSE, xlim=c(-3,2.3), col="#E69F00", cex.lab=1.5))
abline(mxht.mntd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(c)"~~~ R^2 ~ "= 0.20, p-value = 0.025"),side=3,line=-.5, adj=0.05,cex=0.9)

# CWM of leaf area across mntd
with(s.fd.c.m, plot(log(CWM.LA)~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, xlim=c(-3,2.3), col="#009E73", cex.lab=1.5))
abline(LA.mntd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(d)"~~~ R^2 ~ "= 0.48, p-value < 0.001"),side=3,line=-.5, adj=0.05,cex=0.9)

dev.off()

##################################
##################################
# SUPPLEMENTARY FIGURES & ANALYSIS

# S-FIG 1 - Soil texture triangle
jpeg("./functional_traits_rhf/analysis/figures/soil-texture-triangle.jpeg", width=7, height=6.5, unit="in",res=300)
TT.plot(class.sys = "USDA.TT", tri.data = s.fd.c.m, pch=19, cex.axis=0.8, cex.lab=0.8, lwd=0.8, lwd.axis=0.8, lwd.lab=0.8, main="", col="darkorange4", new.mar=c(2.5,0.5,0,0.5))
dev.off()

# S-FIG 2 - Mean and Maximum height vary across percent of soil clay content
jpeg("./functional_traits_rhf/analysis/figures/height-clay.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,2))
par(mar=c(5,5,1.5,1.5))
par(oma=c(0.2,0.3,0.5,0.5))

# mean height across clay content
with(s.fd.c.m, plot(CWM.mn.ht~CLAY, pch=19, xlab= "Clay (%)", ylab="CWM - mean height (cm)", xlim=c(0,30), ylim=c(20,90),cex=1.2, axes=FALSE, col="#E69F00", cex.lab=0.7))
abline(s.mnH.C, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(a)"~~~ R^2 ~ "= 0.20, p-value = 0.027"),side=3,line=0, adj=0.15,cex=0.8)

# max height across clay content
with(s.fd.c.m, plot(CWM.maxht~CLAY, pch=19, xlab= "Clay (%)", ylab="CWM - max. height (cm)", xlim=c(0,30), cex=1.2, axes=FALSE, col="#E69F00", cex.lab=0.7))
abline(s.mxH.C, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(b)"~~~ R^2 ~ "= 0.18, p-value = 0.033"),side=3,line=0, adj=0.15,cex=0.8)

dev.off()

# S-FIG 3 - PCA of the height and leaf traits for each species found in all 78 plots
ht.leaf.traits <- c.data$data[,c(1:3,7:8)]
colnames(ht.leaf.traits) <- c("SLA", "LA", "LDMC", "Max. Ht.", "Mn. Ht.")
ht.leaf.pca <- prcomp(ht.leaf.traits, scale=TRUE)
summary(ht.leaf.pca)

jpeg("./functional_traits_rhf/analysis/figures/height-leaf-pca.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,2))
par(mar=c(4,5,1,2))
par(oma=c(0.2,0.3,0.5,0.5))
biplot(ht.leaf.pca, cex=c(0.4,0.4), col=c("gray", "darkred"),cex.axis=0.7, cex.lab=0.7)
biplot(ht.leaf.pca, choices=2:3, cex=c(0.4,0.4), col=c("gray","darkred"), cex.axis=0.7, cex.lab=0.7)

dev.off()

# S-FIG 4 - Variation of FDis and CWMs across SD of temp
jpeg("./functional_traits_rhf/analysis/figures/fdiv-SD-temp-25plots.jpeg", width=7, height=5, unit="in",res=300)
par(mfrow=c(2,3))
par(mar=c(4,4,0.1,0.1))
par(oma=c(1,1,1,1))

# functional dispersion ~ sd temperature
with(s.fd.c.m, plot(fdis~sd, pch=19, xlab=expression(paste('SD temperature (',degree,'C)')), ylab="Functional dispersion", cex=1.2, axes=FALSE, xlim=c(5,15), ylim=c(0.2,1.4),col="#0072B2", cex.lab=0.7))
abline(s.fdis.sd, lwd=3)
axis(1, cex.axis=0.9, at=c(5,7,9,11,13,15))
axis(2, cex.axis=0.9)
mtext(expression("(a)"~~~ R^2 ~ "= 0.42, p-value < 0.001"),side=3,line=-1, adj=0.15,cex=0.6)

# SLA ~ sd
with(s.fd.c.m, plot(log(CWM.SLA)~sd, pch=19, xlab="", ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(5,15), ylim=c(2,5.5), cex.lab=0.7))
abline(s.SLA.sd, lwd=3)
axis(1, cex.axis=0.9, at=c(5,7,9,11,13,15))
axis(2, cex.axis=0.9)
mtext(expression("(b)"~~~ R^2 ~ "= 0.32, p-value = 0.003"),side=3,line=-1, adj=0.15,cex=0.6)

# mn ht ~ sd
with(s.fd.c.m, plot(CWM.mn.ht~sd, pch=19, xlab="", ylab="CWM - mean height (cm)", cex=1.2, axes=FALSE, xlim=c(5,15), ylim=c(20,100), col="#E69F00", cex.lab=0.7))
abline(s.mnH.sd, lwd=3)
axis(1, cex.axis=0.9, at=c(5,7,9,11,13,15))
axis(2, cex.axis=0.9)
mtext(expression("(c)"~~~ R^2 ~ "= 0.20, p-value = 0.025"),side=3,line=-1, adj=0.15,cex=0.6)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(0.3,0.9, legend=c('Overall functional diversity', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1.2, cex=0.6,
       col = c('#0072B2', '#009E73', '#E69F00'))

# LA ~ sd
with(s.fd.c.m, plot(log(CWM.LA)~sd, pch=19, xlab=expression(paste('SD temperature (',degree,'C)')), ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(5,15), ylim=c(-2,4),cex.lab=0.7))
abline(s.LA.sd, lwd=3)
axis(1, cex.axis=0.9, at=c(5,7,9,11,13,15))
axis(2, cex.axis=0.9)
mtext(expression("(d)"~~~ R^2 ~ "= 0.44, p-value < 0.001"),side=3,line=-1, adj=0.15,cex=0.6)

# max ht ~ sd
with(s.fd.c.m, plot(CWM.maxht~sd, pch=19, xlab=expression(paste('SD temperature (',degree,'C)')), ylab="CWM - max. height (cm)", cex=1.2, axes=FALSE, xlim=c(5,15), ylim=c(40,160), col="#E69F00", cex.lab=0.7))
abline(s.mxH.sd, lwd=3)
axis(1, cex.axis=0.9, at=c(5,7,9,11,13,15))
axis(2, cex.axis=0.9)
mtext(expression("(e)"~~~ R^2 ~ "= 0.25, p-value = 0.011"),side=3,line=-1, adj=0.15,cex=0.6)

dev.off()

# S-FIG 5 - Variation of FDis and CWMs acros MAX temp
jpeg("./functional_traits_rhf/analysis/figures/fdiv-MAX-temp-25plots.jpeg", width=7, height=5, unit="in",res=300)
par(mfrow=c(2,3))
par(mar=c(4,4,0.1,0.1))
par(oma=c(1,1,1,1))

# functional dispersion ~ max temperature
with(s.fd.c.m, plot(fdis~max, pch=19, xlab=expression(paste('Max. temperature (',degree,'C)')), ylab="Functional dispersion", cex=1.2, axes=FALSE, xlim=c(10,70), ylim=c(0.2,1.4),col="#0072B2", cex.lab=0.7))
abline(s.fdis.mx, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(a)"~~~ R^2 ~ "= 0.33, p-value = 0.003"),side=3,line=-1, adj=0.15,cex=0.6)

# SLA ~ max
with(s.fd.c.m, plot(log(CWM.SLA)~max, pch=19, xlab="", ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(10,70), ylim=c(2,5.5), cex.lab=0.7))
abline(s.SLA.mx, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(b)"~~~ R^2 ~ "= 0.25, p-value = 0.011"),side=3,line=-1, adj=0.15,cex=0.6)

# mn ht ~ max
with(s.fd.c.m, plot(CWM.mn.ht~max, pch=19, xlab="", ylab="CWM - mean height (cm)", cex=1.2, axes=FALSE, xlim=c(10,70), ylim=c(20,100), col="#E69F00", cex.lab=0.7))
abline(s.mnH.mx, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(c)"~~~ R^2 ~ "= 0.16, p-value = 0.047"),side=3,line=-1, adj=0.15,cex=0.6)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(0.3,0.9, legend=c('Overall functional diversity', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1.2, cex=0.6,
       col = c('#0072B2', '#009E73', '#E69F00'))

# LA ~ max
with(s.fd.c.m, plot(log(CWM.LA)~max, pch=19, xlab=expression(paste('Max. temperature (',degree,'C)')), ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(10,70), ylim=c(-2,4),cex.lab=0.7))
abline(s.LA.mx, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(d)"~~~ R^2 ~ "= 0.28, p-value = 0.006"),side=3,line=-1, adj=0.15,cex=0.6)

# max ht ~ max
with(s.fd.c.m, plot(CWM.maxht~max, pch=19, xlab=expression(paste('Max. temperature (',degree,'C)')), ylab="CWM - max. height (cm)", cex=1.2, axes=FALSE, xlim=c(10,70), ylim=c(40,160), col="#E69F00", cex.lab=0.7))
abline(s.mxH.mx, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(e)"~~~ R^2 ~ "= 0.21, p-value = 0.023"),side=3,line=-1, adj=0.15,cex=0.6)

dev.off()

# S-FIG 6 - Functional diversity across phy-div metrics at 25 plots
jpeg("./functional_traits_rhf/analysis/figures/supp-fdiv-phydiv-25plots.jpeg", width=7, height=6.5, unit="in",res=300)
par(mfrow=c(2,2))
par(mar=c(4,5,1.5,1))
par(oma=c(0.2,0.3,0.5,0.5))

# SLA across PD
with(s.fd.c.m, plot(log(CWM.SLA)~pd, pch=19, xlab="Faith's PD", ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(500,2500), ylim=c(2,5), cex.lab=1.5))
abline(SLA.pd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(a)"~~~ R^2 ~ "= 0.39, p-value < 0.001"),side=3,line=-.5, adj=0.05,cex=0.9)

# CWM of SLA across mpd
with(s.fd.c.m, plot(log(CWM.SLA)~sesmpd, pch=19, xlab=expression(SES[MPD]), ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(-2,2), ylim=c(2,5), cex.lab=1.5))
abline(SLA.mpd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(b)"~~~ R^2 ~ "= 0.28, p-value = 0.006"),side=3, line=-.5, adj=0.05,cex=0.9)

legend("bottomright", legend=c('Overall FD', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1, cex=.8, 
       col = c('#0072B2', '#009E73', '#E69F00'))

# CWM of LA across pd
with(s.fd.c.m, plot(log(CWM.LA)~pd, pch=19, xlab="Faith's PD", ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, xlim=c(500,2500), col="#009E73", cex.lab=1.5))
abline(LA.pd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(c)"~~~ R^2 ~ "= 0.46, p-value < 0.001"),side=3,line=-.5, adj=0.05,cex=0.9)

# CWM of mean height across Simpson's diversity
with(s.fd.c.m, plot(CWM.mn.ht~simpsdiv, pch=19, xlab="Simpson's Diversity", ylab="CWM of mean height (cm)", cex=1.2, axes=FALSE,  col="#E69F00", ylim=c(20,90),cex.lab=1.5))
abline(mnht.simp, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(d)"~~~ R^2 ~ "= 0.24, p-value = 0.014"),side=3,line=-.5, adj=0.05,cex=0.9)

dev.off()

# S-FIG 7...cut out for now.

# summary(with(s.fd.c.m, lm(scat~fdis)))
# with(s.fd.c.m, plot(scat~fdis))
# 
# summary(with(s.fd.c.m, lm(soil+~fdis)))
# with(s.fd.c.m, plot(soil~fdis))
# 
# summary(with(s.fd.c.m, lm(CWM.mn.ht~litter)))
# with(s.fd.c.m, plot(CWM.mn.ht~litter))
# 
# summary(with(s.fd.c.m, lm(CWM.maxht~scat)))
# with(s.fd.c.m, plot(CWM.maxht~scat))
# 
# summary(with(s.fd.c.m, lm(native~scat)))
# with(s.fd.c.m, plot(fdis~scat))
