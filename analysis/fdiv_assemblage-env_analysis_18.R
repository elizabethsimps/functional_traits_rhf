# Assessing variance in functional dispersion and CWM of height and leaf traits across temp and soil microclimates and phylodiversity
# Elizabeth Simpson # 2021-03-25
# Uses code from Will Pearse to get comparative.comm set up

setwd("~/Documents/projects/functional_traits_rhf") # uncomment to set working directory on computer

###########
# HEADERS 
library(pez)
library(FD)
library(lme4)
library(soiltexture)
library(car)

#################
### LOAD DATA ###
cover <- read.csv("./raw_data/rhf_2018_cover.csv", as.is=TRUE)
tree <- read.tree("./clean_data/Vascular_Plants_rooted.dated.tre")
env <- read.csv("./clean_data/temp-texture-terrain-18.csv", as.is=TRUE)
traits <- read.csv("./clean_data/clean_traits_rhf.csv", as.is=TRUE)

### CLEAN DATA ###
# Community data - Add in overhanging cover
cover$Cover.over <- NA
for (i in seq_len(nrow(cover))){
  if (with(cover, grepl("overhang-", Notes[i]))){
    cover$Cover.over[i] <- with(cover, Cover[i] + as.numeric(unlist(regmatches(Notes[i], gregexpr('\\(?[0-9,.]+', Notes[i])))))
  } else {
    cover$Cover.over[i] <- cover$Cover[i]
  }
}

# subset by environmental data so you don't end up with species with no abundances in the final dataset
cover <- cover[cover$Plot_id %in% env$plot_id,]

# reformat, remove NAs, remove species essentially there as notes
comm <- with(cover, tapply(Cover.over, list(Plot_id,Species), function(x) mean(x,na.rm=TRUE)))
comm[is.na(comm)] <- 0
comm <- comm[, colSums(comm != 0) > 0]
comm <- comm[,!grepl("^[a-z]+", colnames(comm))]
comm <- comm[,!grepl("\\(|/", colnames(comm))]
comm <- comm[,!grepl("sp\\.", colnames(comm))]
colnames(comm) <- tolower(colnames(comm))

################################
# PHYLOGENY (Zanne et. al. 2014) - reformat sp names, create combined data object
tree$tip.label <- tolower(gsub("_", " ", tree$tip.label))
tree <- congeneric.merge(tree, colnames(comm), split = " ")

#################################
# ENVIRONMENTAL DATA - topography, terrain, texture, for the 25 plots with temp & soil texture data
rownames(env) <- env$plot_id
env <- env[,c(3:17)]

############
# TRAIT DATA - Focusing analysis on continuous variables and means for each leaf and height traits
# Also have some categorical traits that could be intersting...dispersal syndrome, native vs. not, etc...
rownames(traits) <- tolower(traits$species)
traits <- traits[,c(3:5,9:10)]

###################################
# Make comparative community object
c.data <- comparative.comm(tree, comm, traits=traits, env=env)

#####################################################################
# Calculate functional, phylogenetic, and taxonomic diversity indices
# Focus on functional dispersion and community weighted means but the 3 Villeger et al. 2008 indicies are also in here
f.div <- with(c.data, dbFD(data[,c(1:2,4:5)], comm, w.abun=TRUE, stand.x=TRUE))
f.distwo <- with(c.data, dbFD(data[,c(2,4)], comm, w.abun=TRUE, stand.x=TRUE))

ses.mpd.a <- .ses.mpd(c.data, abundace.weighted=TRUE)
ses.mntd.a <- .ses.mntd(c.data, abundance.weighted=TRUE)
pd <- .pd(c.data)

simps.div <- diversity(c.data$comm, index="simpson")

div <- with(f.div, cbind(FDis, CWM, nbsp, ses.mpd.a$mpd.obs.z, ses.mntd.a$mntd.obs.z, pd[,1], simps.div))
div <- with(f.distwo, cbind(div, FDis))

colnames(div) <- c("fdis","CWM.SLA", "CWM.LA", "CWM.maxht","CWM.mn.ht", "nsp", "sesmpd","sesmntd", "pd", "simpsdiv", "fdis.two")

div <- cbind(div, c.data$env)

###################################################################
# Do temperature and soil texture vary systematically across environment?
# This seems more silly now, but I already have all of the plotting done...

cor.test.env <- function(div,RV){
  output <- c(with(div, cor.test(RV,aspect))$estimate, with(div, cor.test(RV, elev))$estimate, with(div, cor.test(RV, slope))$estimate,
              with(div, cor.test(RV, SAND))$estimate, with(div, cor.test(RV, SILT))$estimate, with(div, cor.test(RV, CLAY))$estimate)
  output <- as.data.frame(cbind(output, c(with(div, cor.test(RV, aspect))$p.value, with(div, cor.test(RV, elev))$p.value, with(div, cor.test(RV, slope))$p.value,
                                           with(div, cor.test(RV, SAND))$p.value, with(div, cor.test(RV, SILT))$p.value, with(div, cor.test(RV, CLAY))$p.value)))
  rownames(output) <- c("aspect","elev","slope","sand","silt","clay")
  colnames(output) <- c("pearsons.r", "p.value")
  return(output)
}

cor.test.soil <<- function(div,RV){
  output <- c(with(div, cor.test(RV,aspect))$estimate, with(div, cor.test(RV, elev))$estimate, with(div, cor.test(RV, slope))$estimate)
  output <- as.data.frame(cbind(output, c(with(div, cor.test(RV, aspect))$p.value, with(div, cor.test(RV, elev))$p.value, with(div, cor.test(RV, slope))$p.value)))
  rownames(output) <- c("aspect","elev","slope")
  colnames(output) <- c("pearsons.r", "p.value")
  return(output)
}

mean.corr <- cor.test.env(div, div$mean)
sd.corr <- cor.test.env(div, div$sd)
min.corr <- cor.test.env(div, div$min)
max.corr <- cor.test.env(div, div$max)

sand.corr <- cor.test.soil(div, div$SAND)
silt.corr <- cor.test.soil(div, div$SILT)
clay.corr <- cor.test.soil(div, div$CLAY)

# Models for plotting
mn.a <- with(div, lm(mean~aspect))
min <- with(div, lm(min~1))
mx.a <- with(div, lm(max~aspect))
sd.a <- with(div, lm(sd~aspect))

sand.e <- with(div, lm(SAND~elev))
clay.e <- with(div, lm(CLAY~elev))

# FIG 2 - Ploting changes in temperature across environment plus the USDA texture triangle for all of the plots 
jpeg("./analysis/figures/env.jpeg", width=9.5, height=9.5, unit="in",res=300)
par(mfrow=c(2,2))
par(mar=c(5,5,1.5,1.5))
par(oma=c(0.2,0.3,0.5,0.5))

with(div, plot(mean~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(paste('Temperature (',degree~'C)')), cex=1.2, axes=FALSE, ylim=c(-10,70), col="#E69F00"))
abline(mn.a, lwd=3, col="#E69F00")
with(div, points(max~aspect, pch=19, cex=1.2, col="#D55E00"))
abline(mx.a, lwd=3, col="#D55E00")
with(div, points(min~aspect, pch=19, cex=1.2, col="#56B4E9"))
abline(min, lwd=3, col="#56B4E9")
axis(1)
axis(2, at=c(-10,10,30,50,70))
legend("topright", c(expression(~ R^2 ~ "= 0.35, p-value = 0.002"), expression(~ R^2 ~ "= 0.68, p-value < 0.001")), 
       text.col=c("#D55E00", "#E69F00"), box.lty=0, bg="transparent")
mtext("(a)",side=3,line=-.5, adj=0.05)

with(div, plot(sd~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(paste('SD - Temperature (',degree~'C)')), cex=1.2, axes=FALSE, col="#E69F00", ylim=c(5,15)))
abline(sd.a, lwd=3, col="#E69F00")
axis(1)
axis(2, at=c(5,10,15))
mtext("(b)",side=3,line=-.5, adj=0.05)
legend(-1,7.4, legend=c('Min', 'Mean', 'Max'), pch=16, pt.cex=1.2, cex=0.6,
       col = c('#56B4E9', '#E69F00', '#D55E00'))
legend("topright", c(expression(~ R^2 ~ "= 0.44, p-value < 0.001")), 
       text.col=c("#E69F00"), box.lty=0, bg="transparent")

with(div, plot(SAND~elev, pch=19, xlab="Elevation (m.s.l.)", ylab= "Percent of soil component", cex=1.2, axes=FALSE, col="goldenrod3", xlim=c(1700,2100), ylim=c(0,80)))
abline(sand.e, lwd=3, col="goldenrod3")
with(div, points(CLAY~elev, pch=19, cex=1.2, col="sienna"))
abline(clay.e, lwd=3, col="sienna")
axis(1)
axis(2)
legend(1710,10, legend=c("Sand", "Clay"), pch=16, pt.cex=1.2, cex=0.6,
       col = c("goldenrod3", "sienna"))
legend("topright", c(expression(~ R^2 ~ "= 0.21, p-value = 0.024"), expression(~ R^2 ~ "= 0.23, p-value = 0.016")), 
       text.col=c("goldenrod3", "sienna"), box.lty=0, bg="transparent")
mtext("(c)",side=3,line=-.5, adj=0.05)

TT.plot(class.sys = "USDA.TT", tri.data = div, pch=19, cex.axis=0.8, cex.lab=0.8, lwd=0.8, lwd.axis=0.8, lwd.lab=0.8, main="", col="darkgoldenrod4", new.mar=c(3.7,0.5,0,0.5))
mtext("(d)",side=3,line=-2.1, adj=0.05)

dev.off()

####################################################################################
### How do functional does functional dispersion and the community weighted mean of the traits vary across temperature and texture?
### CWM of SLA
sla <- with(div, lm(CWM.SLA~mean+CLAY)) 
sla.m <- with(div, lm(CWM.SLA~mean))

### CWM of LA
la <- with(div, lm(CWM.LA~mean+CLAY))
la.a <- with(div, lm(CWM.LA~mean))

### CWM of mean HT
mnH <- with(div, lm(CWM.mn.ht~mean+CLAY))

### CWM of max HT
mxH <- with(div, lm(CWM.maxht~mean+CLAY))

### FDis
fdis <- with(div, lm(fdis~mean+CLAY))
# and check with two traits too
fdis.two <- with(div, lm(fdis.two~mean+CLAY)) # not significant across the mean any more.

### FIG 3 - Plotting how functional diversity changes across the MEAN temperature variable.
jpeg("./functional_traits_rhf/analysis/figures/fdiv-MEAN-25plots.jpeg", width=7, height=5, unit="in",res=300)
par(mfrow=c(2,3))
par(mar=c(4,4,0.1,0.1))
par(oma=c(1,1,1,1))

# functional dispersion ~ mean temperature
with(div, plot(fdis~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab="Functional dispersion", cex=1.2, axes=FALSE, xlim=c(4,14), ylim=c(0.2,1.4),col="#0072B2", cex.lab=0.7))
abline(s.fdis.mn, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(a)"~~~ R^2 ~ "= 0.56, p-value < 0.001"),side=3,line=-1, adj=0.15,cex=0.6)

# SLA ~ mean
with(div, plot(log(CWM.SLA)~mean, pch=19, xlab="", ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(4,14), ylim=c(2,5.5), cex.lab=0.7))
abline(s.SLA.mn, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(b)"~~~ R^2 ~ "= 0.37, p-value = 0.001"),side=3,line=-1, adj=0.15,cex=0.6)

# mn ht ~ mean
with(div, plot(CWM.mn.ht~mean, pch=19, xlab="", ylab="CWM - mean height (cm)", cex=1.2, axes=FALSE, xlim=c(4,14), ylim=c(20,100), col="#E69F00", cex.lab=0.7))
abline(s.mnH.mn, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(c)"~~~ R^2 ~ "= 0.24, p-value = 0.013"),side=3,line=-1, adj=0.15,cex=0.6)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(0.3,0.9, legend=c('Overall functional diversity', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1.2, cex=0.6,
       col = c('#0072B2', '#009E73', '#E69F00'))

# LA ~ mean
with(div, plot(log(CWM.LA)~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(4,14), ylim=c(-2,4),cex.lab=0.7))
abline(s.LA.mn, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(d)"~~~ R^2 ~ "= 0.63, p-value < 0.001"),side=3,line=-1, adj=0.15,cex=0.6)

# max ht ~ mean
with(div, plot(CWM.maxht~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab="CWM - max. height (cm)", cex=1.2, axes=FALSE, xlim=c(4,14), ylim=c(40,160), col="#E69F00", cex.lab=0.7))
abline(s.mxH.mn, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(e)"~~~ R^2 ~ "= 0.28, p-value = 0.006"),side=3,line=-1, adj=0.15,cex=0.6)

dev.off()

#########################################################################
# How does functional diversity (dispersion and CWM) correlate with phylodiversity? (25 plots)

fdis.mntd <- with(div, lm(fdis~sesmntd))

LA.pd <- with(div, lm(log(CWM.LA)~pd))
LA.mntd <- with(div, lm(log(CWM.LA)~sesmntd))

SLA.pd <- with(div, lm(log(CWM.SLA)~pd))
SLA.mntd <- with(div, lm(log(CWM.SLA)~sesmntd))
SLA.mpd <- with(div, lm(log(CWM.SLA)~sesmpd))

mnht.simp <- with(div, lm(CWM.mn.ht~simpsdiv))

mxht.mntd <- with(div, lm(CWM.maxht~sesmntd))

#################################################################
# FIG 4 - Functional diversity across phy-div metrics at 25 plots
# Plot significant relationships with SESmntd here, do the rest in the supplement
jpeg("./functional_traits_rhf/analysis/figures/fdiv-phydiv-25plots.jpeg", width=7, height=6.5, unit="in",res=300)
par(mfrow=c(2,2))
par(mar=c(4,5,0.5,1))
par(oma=c(0.2,0.3,0.5,0.5))

# Functional dispersion across mntd
with(div, plot(fdis~sesmntd, pch=19, xlab="", ylab="Functional dispersion", cex=1.2, axes=FALSE, col="#0072B2", xlim=c(-3,2.3), cex.lab=1.5))
abline(fdis.mntd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(a)"~~~ R^2 ~ "= 0.42, p-value < 0.001"),side=3,line=-.5, adj=0.05,cex=0.9)

legend("bottomright", legend=c('Overall FD', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1, cex=.8, 
       col = c('#0072B2', '#009E73', '#E69F00'))

# CWM of SLA across mntd
with(div, plot(log(CWM.SLA)~sesmntd, pch=19, xlab="", ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, xlim=c(-3,2.3), ylim=c(2,5), col="#009E73", cex.lab=1.5))
abline(SLA.mntd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(b)"~~~ R^2 ~ "= 0.29, p-value = 0.006"),side=3,line=-.5, adj=0.05,cex=0.9)

# CWM max. height across mntd
with(div, plot(CWM.maxht~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab="CWM - max. height (cm)", cex=1.2, axes=FALSE, xlim=c(-3,2.3), col="#E69F00", cex.lab=1.5))
abline(mxht.mntd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(c)"~~~ R^2 ~ "= 0.20, p-value = 0.025"),side=3,line=-.5, adj=0.05,cex=0.9)

# CWM of leaf area across mntd
with(div, plot(log(CWM.LA)~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, xlim=c(-3,2.3), col="#009E73", cex.lab=1.5))
abline(LA.mntd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(d)"~~~ R^2 ~ "= 0.48, p-value < 0.001"),side=3,line=-.5, adj=0.05,cex=0.9)

dev.off()

##################################
##################################
# SUPPLEMENTARY FIGURES & ANALYSIS

# S-FIG 2 - Mean and Maximum height vary across percent of soil clay content
jpeg("./functional_traits_rhf/analysis/figures/height-clay.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,2))
par(mar=c(5,5,1.5,1.5))
par(oma=c(0.2,0.3,0.5,0.5))

# mean height across clay content
with(div, plot(CWM.mn.ht~CLAY, pch=19, xlab= "Clay (%)", ylab="CWM - mean height (cm)", xlim=c(0,30), ylim=c(20,90),cex=1.2, axes=FALSE, col="#E69F00", cex.lab=0.7))
abline(s.mnH.C, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(a)"~~~ R^2 ~ "= 0.20, p-value = 0.027"),side=3,line=0, adj=0.15,cex=0.8)

# max height across clay content
with(div, plot(CWM.maxht~CLAY, pch=19, xlab= "Clay (%)", ylab="CWM - max. height (cm)", xlim=c(0,30), cex=1.2, axes=FALSE, col="#E69F00", cex.lab=0.7))
abline(s.mxH.C, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(b)"~~~ R^2 ~ "= 0.18, p-value = 0.033"),side=3,line=0, adj=0.15,cex=0.8)

dev.off()

# S-FIG 3 - PCA of the height and leaf traits for each species found in all 78 plots
ht.leaf.traits <- c.data$data[,c(1:2, 4:5)]
colnames(ht.leaf.traits) <- c("SLA", "LA", "Max. Ht.", "Mn. Ht.")
ht.leaf.pca <- prcomp(ht.leaf.traits, scale=TRUE)
summary(ht.leaf.pca)

jpeg("./analysis/figures/height-leaf-pca.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,2))
par(mar=c(4,5,1,2))
par(oma=c(0.2,0.3,0.5,0.5))
biplot(ht.leaf.pca, cex=c(0.4,0.4), col=c("gray", "darkred"),cex.axis=0.7, cex.lab=0.7)
biplot(ht.leaf.pca, choices=2:3, cex=c(0.4,0.4), col=c("gray","darkred"), cex.axis=0.7, cex.lab=0.7)

dev.off()

# S-FIG 3 - PCA of the all the environmental traits
env.pca <- div[,c(15:20,22:25)]
env.pca <- prcomp(env.pca, scale=TRUE)
summary(env.pca)
plot(env.pca)

jpeg("./analysis/figures/env-pca.jpeg", width=7, height=3.5, unit="in",res=300)
par(mar=c(5,5,2,2))
biplot(env.pca, cex=c(0.4,0.4), col=c("gray", "darkred"),cex.axis=0.7, cex.lab=0.7, xlab = "standardized PC1 (40.4% explained var.)", ylab = "standardized PC2 (28.1% explained var.)")
dev.off()

# S-FIG 4 - Variation of FDis and CWMs across SD of temp
jpeg("./functional_traits_rhf/analysis/figures/fdiv-SD-temp-25plots.jpeg", width=7, height=5, unit="in",res=300)
par(mfrow=c(2,3))
par(mar=c(4,4,0.1,0.1))
par(oma=c(1,1,1,1))

# functional dispersion ~ sd temperature
with(div, plot(fdis~sd, pch=19, xlab=expression(paste('SD temperature (',degree,'C)')), ylab="Functional dispersion", cex=1.2, axes=FALSE, xlim=c(5,15), ylim=c(0.2,1.4),col="#0072B2", cex.lab=0.7))
abline(s.fdis.sd, lwd=3)
axis(1, cex.axis=0.9, at=c(5,7,9,11,13,15))
axis(2, cex.axis=0.9)
mtext(expression("(a)"~~~ R^2 ~ "= 0.42, p-value < 0.001"),side=3,line=-1, adj=0.15,cex=0.6)

# SLA ~ sd
with(div, plot(log(CWM.SLA)~sd, pch=19, xlab="", ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(5,15), ylim=c(2,5.5), cex.lab=0.7))
abline(s.SLA.sd, lwd=3)
axis(1, cex.axis=0.9, at=c(5,7,9,11,13,15))
axis(2, cex.axis=0.9)
mtext(expression("(b)"~~~ R^2 ~ "= 0.32, p-value = 0.003"),side=3,line=-1, adj=0.15,cex=0.6)

# mn ht ~ sd
with(div, plot(CWM.mn.ht~sd, pch=19, xlab="", ylab="CWM - mean height (cm)", cex=1.2, axes=FALSE, xlim=c(5,15), ylim=c(20,100), col="#E69F00", cex.lab=0.7))
abline(s.mnH.sd, lwd=3)
axis(1, cex.axis=0.9, at=c(5,7,9,11,13,15))
axis(2, cex.axis=0.9)
mtext(expression("(c)"~~~ R^2 ~ "= 0.20, p-value = 0.025"),side=3,line=-1, adj=0.15,cex=0.6)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(0.3,0.9, legend=c('Overall functional diversity', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1.2, cex=0.6,
       col = c('#0072B2', '#009E73', '#E69F00'))

# LA ~ sd
with(div, plot(log(CWM.LA)~sd, pch=19, xlab=expression(paste('SD temperature (',degree,'C)')), ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(5,15), ylim=c(-2,4),cex.lab=0.7))
abline(s.LA.sd, lwd=3)
axis(1, cex.axis=0.9, at=c(5,7,9,11,13,15))
axis(2, cex.axis=0.9)
mtext(expression("(d)"~~~ R^2 ~ "= 0.44, p-value < 0.001"),side=3,line=-1, adj=0.15,cex=0.6)

# max ht ~ sd
with(div, plot(CWM.maxht~sd, pch=19, xlab=expression(paste('SD temperature (',degree,'C)')), ylab="CWM - max. height (cm)", cex=1.2, axes=FALSE, xlim=c(5,15), ylim=c(40,160), col="#E69F00", cex.lab=0.7))
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
with(div, plot(fdis~max, pch=19, xlab=expression(paste('Max. temperature (',degree,'C)')), ylab="Functional dispersion", cex=1.2, axes=FALSE, xlim=c(10,70), ylim=c(0.2,1.4),col="#0072B2", cex.lab=0.7))
abline(s.fdis.mx, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(a)"~~~ R^2 ~ "= 0.33, p-value = 0.003"),side=3,line=-1, adj=0.15,cex=0.6)

# SLA ~ max
with(div, plot(log(CWM.SLA)~max, pch=19, xlab="", ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(10,70), ylim=c(2,5.5), cex.lab=0.7))
abline(s.SLA.mx, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(b)"~~~ R^2 ~ "= 0.25, p-value = 0.011"),side=3,line=-1, adj=0.15,cex=0.6)

# mn ht ~ max
with(div, plot(CWM.mn.ht~max, pch=19, xlab="", ylab="CWM - mean height (cm)", cex=1.2, axes=FALSE, xlim=c(10,70), ylim=c(20,100), col="#E69F00", cex.lab=0.7))
abline(s.mnH.mx, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(c)"~~~ R^2 ~ "= 0.16, p-value = 0.047"),side=3,line=-1, adj=0.15,cex=0.6)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(0.3,0.9, legend=c('Overall functional diversity', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1.2, cex=0.6,
       col = c('#0072B2', '#009E73', '#E69F00'))

# LA ~ max
with(div, plot(log(CWM.LA)~max, pch=19, xlab=expression(paste('Max. temperature (',degree,'C)')), ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(10,70), ylim=c(-2,4),cex.lab=0.7))
abline(s.LA.mx, lwd=3)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.9)
mtext(expression("(d)"~~~ R^2 ~ "= 0.28, p-value = 0.006"),side=3,line=-1, adj=0.15,cex=0.6)

# max ht ~ max
with(div, plot(CWM.maxht~max, pch=19, xlab=expression(paste('Max. temperature (',degree,'C)')), ylab="CWM - max. height (cm)", cex=1.2, axes=FALSE, xlim=c(10,70), ylim=c(40,160), col="#E69F00", cex.lab=0.7))
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
with(div, plot(log(CWM.SLA)~pd, pch=19, xlab="Faith's PD", ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(500,2500), ylim=c(2,5), cex.lab=1.5))
abline(SLA.pd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(a)"~~~ R^2 ~ "= 0.39, p-value < 0.001"),side=3,line=-.5, adj=0.05,cex=0.9)

# CWM of SLA across mpd
with(div, plot(log(CWM.SLA)~sesmpd, pch=19, xlab=expression(SES[MPD]), ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, col="#009E73", xlim=c(-2,2), ylim=c(2,5), cex.lab=1.5))
abline(SLA.mpd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(b)"~~~ R^2 ~ "= 0.28, p-value = 0.006"),side=3, line=-.5, adj=0.05,cex=0.9)

legend("bottomright", legend=c('Overall FD', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1, cex=.8, 
       col = c('#0072B2', '#009E73', '#E69F00'))

# CWM of LA across pd
with(div, plot(log(CWM.LA)~pd, pch=19, xlab="Faith's PD", ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, xlim=c(500,2500), col="#009E73", cex.lab=1.5))
abline(LA.pd, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(c)"~~~ R^2 ~ "= 0.46, p-value < 0.001"),side=3,line=-.5, adj=0.05,cex=0.9)

# CWM of mean height across Simpson's diversity
with(div, plot(CWM.mn.ht~simpsdiv, pch=19, xlab="Simpson's Diversity", ylab="CWM of mean height (cm)", cex=1.2, axes=FALSE,  col="#E69F00", ylim=c(20,90),cex.lab=1.5))
abline(mnht.simp, lwd=3)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
mtext(expression("(d)"~~~ R^2 ~ "= 0.24, p-value = 0.014"),side=3,line=-.5, adj=0.05,cex=0.9)

dev.off()

# S-FIG 7...cut out for now.

# summary(with(div, lm(scat~fdis)))
# with(div, plot(scat~fdis))
# 
# summary(with(div, lm(soil+~fdis)))
# with(div, plot(soil~fdis))
# 
# summary(with(div, lm(CWM.mn.ht~litter)))
# with(div, plot(CWM.mn.ht~litter))
# 
# summary(with(div, lm(CWM.maxht~scat)))
# with(div, plot(CWM.maxht~scat))
# 
# summary(with(div, lm(native~scat)))
# with(div, plot(fdis~scat))
