# Assessing variance in functional dispersion and CWM of height and leaf traits across temp and soil microclimates and phylodiversity
# Elizabeth Simpson # 2021-03-25
# Uses code from Will Pearse to get comparative.comm set up

setwd("~/Documents/projects/functional_traits_rhf") # uncomment to set working directory on computer

###############
### HEADERS ###
library(pez)
library(FD)
library(lme4)
library(soiltexture)
library(car)
library(xtable)

#################
### LOAD DATA ###
env.micro <- read.csv("./clean_data/temp-year17-18_text_terr.csv", as.is=TRUE)
d.eight <- read.csv("./raw_data/rhf_2018_cover.csv", as.is = TRUE)
tree <- read.tree("./clean_data/Vascular_Plants_rooted.dated.tre")
env.cvr18 <- read.csv("./clean_data/temp-year17-18_text_terr.csv", as.is=TRUE)
traits <- read.csv("./clean_data/clean_traits_rhf.csv", as.is=TRUE)

#################
### FUNCTIONS ###

# Add in overhanging cover
calc.over <- function(cover){
  cover$Cover.over <- NA
  for (i in seq_len(nrow(cover))){
    if (with(cover, grepl("overhang-", Notes[i]))){
      cover$Cover.over[i] <- with(cover, Cover[i] + as.numeric(unlist(regmatches(Notes[i], gregexpr('\\(?[0-9,.]+', Notes[i])))))
    } else {
      cover$Cover.over[i] <- cover$Cover[i]
    }
  }
  return(cover)
}

# Clean the community data
calc.comm <- function(cover){
  comm <- with(cover, tapply(Cover.over, list(Plot_id,Species), function(x) mean(x,na.rm=TRUE)))
  comm[is.na(comm)] <- 0
  comm <- comm[, colSums(comm != 0) > 0]
  comm <- comm[,!grepl("^[a-z]+", colnames(comm))]
  comm <- comm[,!grepl("\\(|/", colnames(comm))]
  comm <- comm[,!grepl("sp\\.", colnames(comm))]
  colnames(comm) <- tolower(colnames(comm))
  return(comm)
}

# Make comparative community object, calculate diversity metrics, and add environmental data all to the same dataframe
div.calc <- function(tree, comm, traits, env){
  c.data <- comparative.comm(tree, comm, traits, env)
  # Calculate functional, phylogenetic, and taxonomic diversity indices
  # Focus on functional dispersion and community weighted means but the 3 Villeger et al. 2008 indicies are also in here
  f.div <- with(c.data, dbFD(data[,c(1:2,4:5)], comm, w.abun=TRUE, stand.x=TRUE, print.pco=TRUE))
  f.distwo <- with(c.data, dbFD(data[,c(1,4)], comm, w.abun=TRUE, stand.x=TRUE, print.pco=TRUE))
  # Calculate phylogenetic div. and simpson's div.
  ses.mpd.a <- .ses.mpd(c.data, abundace.weighted=TRUE)
  ses.mntd.a <- .ses.mntd(c.data, abundance.weighted=TRUE)
  pd <- .pd(c.data)
  simps.div<- diversity(c.data$comm, index="simpson")
  # Put together all of the pieces needed
  div <- with(f.div, cbind(FDis, CWM, nbsp, ses.mpd.a$mpd.obs.z, ses.mntd.a$mntd.obs.z, pd[,1], simps.div))
  div <- with(f.distwo, cbind(div, FDis))
  colnames(div) <- c("fdis","CWM.SLA", "CWM.LA", "CWM.maxht","CWM.mn.ht", "nsp", "sesmpd","sesmntd", "pd", "simpsdiv", "fdis.two")
  return(cbind(div, c.data$env))
}

# Calculate Pearson's r and p-values between diversity metrics and environmental variables
env.corr <- function(div,met){
  output <- c(with(div, cor.test(met,mean))$estimate, with(div, cor.test(met,sd))$estimate, with(div, cor.test(met,max))$estimate, with(div, cor.test(met,min))$estimate,
              with(div, cor.test(met, SAND))$estimate, with(div, cor.test(met, SILT))$estimate, with(div, cor.test(met, CLAY))$estimate)
  output <- as.data.frame(rbind(output, c(with(div, cor.test(met,mean))$p.value, with(div, cor.test(met,sd))$p.value, with(div, cor.test(met,max))$p.value, with(div, cor.test(met,min))$p.value,
                                          with(div, cor.test(met, SAND))$p.value, with(div, cor.test(met, SILT))$p.value, with(div, cor.test(met, CLAY))$p.value)))
  colnames(output) <- c("mean", "sd", "max", "min", "sand", "silt", "clay")
  rownames(output) <- c("pearsons.r", "p.value")
  return(output)
}

cor.test.phy <<- function(div,RV){
  output <- c(with(div, cor.test(RV,sesmntd))$estimate, with(div, cor.test(RV, sesmpd))$estimate, with(div, cor.test(RV, pd))$estimate,
              with(div, cor.test(RV, nsp))$estimate, with(div, cor.test(RV, simpsdiv))$estimate)
  output <- as.data.frame(rbind(output, c(with(div, cor.test(RV, sesmntd))$p.value, with(div, cor.test(RV, sesmpd))$p.value, with(div, cor.test(RV, pd))$p.value,
                                          with(div, cor.test(RV, nsp))$p.value, with(div, cor.test(RV, simpsdiv))$p.value)))
  colnames(output) <- c("sesmntd","sesmpd","pd","nsp","simpsdiv")
  rownames(output) <- c("pearsons.r", "p.value")
  return(output)
}

########################################################
# Q1: How does microenvironment vary across topography?
# Uses temperature data collected from 2017-09-28 at 00:00 to 2018-09-29 at 00:00 at 25 plots.
# Note: Temperature sensors were serviced in Sept. 2018 and a few time points are missing on the days that that happened.

# Since these models can only have a max. of 2 (maybe 3) predictor variables,
# I used a PCA between all env. variables to guide which relationships between microenvironment and topography to focus on.
micro.env <- env.micro[,c(7:12,14:17)]
pca.micro <- prcomp(micro.env, scale=TRUE)

# Amount of variance associated with each PCA axis
summary(pca.micro)

# Looking at the loadings of the variables on each axis
pca.micro$rotation

# Plotting PCA for supplement
jpeg("./analysis/figures/supp-PCA-env.jpeg", width=4.5, height=4.5, unit="in",res=300)
par(mar=c(4,5,2,3))
biplot(pca.micro, pch=19, col=c("darkslategray4", "black"), cex=c(0.4,0.6), arrow.len=0.03, xlim=c(-0.5,0.5), ylim=c(-0.4,0.4), xlab = "Standardized PC1 (40.4% explained var.)", ylab = "Standardized PC2 (28.1% explained var.)", cex.axis=0.6, cex.lab=0.6)
dev.off()

# Temperature variables across aspect
mn.a <- with(micro.env, lm(mean~aspect))
sd.a <- with(micro.env, lm(sd~aspect))
max.a <- with(micro.env, lm(max~aspect))
min.1 <- with(micro.env, lm(min~1))

# Texture variables across elevation and slope
sand.e <- with(micro.env, lm(SAND~elev))
clay.e <- with(micro.env, lm(CLAY~elev))

# FIG 2 - Plotting changes in temperature across environment plus the USDA texture triangle for all of the plots 
jpeg("./analysis/figures/env-core18.jpeg", width=7, height=7, unit="in",res=300)
par(mfrow=c(2,2))
par(mar=c(5,5,2,0.5))
par(oma=c(1,1,1,1))

with(micro.env, plot(mean~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(paste('Temperature (',degree~'C)')), cex=1.2, axes=FALSE, ylim=c(-10,70), col="#E69F00"))
abline(mn.a, lwd=3, col="#E69F00")
with(micro.env, points(max~aspect, pch=19, cex=1.2, col="#D55E00"))
abline(max.a, lwd=3, col="#D55E00")
with(micro.env, points(min~aspect, pch=19, cex=1.2, col="#56B4E9"))
abline(min.1, lwd=3, col="#56B4E9")
axis(1)
axis(2, at=c(-10,10,30,50,70))
mtext(expression(bold("(a)")),side=3,line=0.2, adj=-0.32)
legend(-0.5,62, expression(~ R^2 ~ "= 0.35, p-value = 0.002"), text.col="#D55E00", box.lty=0, bg="transparent", cex=0.9)
legend(-1.2,25, expression(~ R^2 ~ "= 0.68, p-value < 0.001"), text.col= "#E69F00", box.lty=0, bg="transparent", cex=0.9)

with(micro.env, plot(sd~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(paste('SD - Temperature (',degree~'C)')), cex=1.2, axes=FALSE, col="#E69F00", ylim=c(5,15)))
abline(sd.a, lwd=3, col="#E69F00")
axis(1)
axis(2, at=c(5,10,15))
mtext(expression(bold("(b)")),side=3,line=0.2, adj=-0.32)
legend(0.5,14.5, legend=c('Min', 'Mean', 'Max'), pch=16, pt.cex=1.2, cex=0.6, col = c('#56B4E9', '#E69F00', '#D55E00'))
legend(-1.2,6.5, expression(~ R^2 ~ "= 0.46, p-value < 0.001"), text.col="#E69F00", box.lty=0, bg="transparent", cex=0.9)

with(micro.env, plot(SAND~elev, pch=19, xlab="Elevation (m.s.l.)", ylab= "Soil component (%)", cex=1.2, axes=FALSE, col="goldenrod3", xlim=c(1700,2100), ylim=c(0,80)))
abline(sand.e, lwd=3, col="goldenrod3")
with(micro.env, points(CLAY~elev, pch=19, cex=1.2, col="sienna"))
abline(clay.e, lwd=3, col="sienna")
axis(1)
axis(2)
legend(2010,80, legend=c("Sand", "Clay"), pch=16, pt.cex=1.2, cex=0.6, col = c("goldenrod3", "sienna"))
legend(1670,65, expression(~ R^2 ~ "= 0.20, p-value = 0.024"), text.col="goldenrod3", box.lty=0, bg="transparent", cex=0.9)
legend(1812,31, expression(~ R^2 ~ "= 0.23, p-value = 0.016"), text.col="sienna", box.lty=0, bg="transparent", cex=0.9)
mtext(expression(bold("(c)")),side=3,line=0.2, adj=-0.32)

TT.plot(class.sys = "USDA.TT", tri.data = micro.env, pch=19, cex.axis=0.8, cex.lab=0.8, lwd=0.8, lwd.axis=0.8, lwd.lab=0.8, main="", cex=0.95, new.mar=c(3.7,3.7,0,0))
mtext(expression(bold("(d)")),side=3,line=-2, adj=-0.19)

dev.off()

#####################################################
# Q2: How does function vary across microenvironment?

### COVER - overhang included
# All plots surveyed in 2018 - 78 plots
int.18.over <- calc.over(d.eight)

# Plots with all microenv. data (temp.) collected - 25 plots
core.18.over <- int.18.over[int.18.over$Plot_id %in% env.cvr18$plot_id,] # plots

# Make community matrix
core.18.comm <- calc.comm(core.18.over)

### PHYLOGENY (Zanne et. al. 2014) - reformat sp names, create combined data object
tree$tip.label <- tolower(gsub("_", " ", tree$tip.label))
# Phylogeny for species in core plots
tree18.core <- congeneric.merge(tree, colnames(core.18.comm), split = " ")

### TRAITS - leaves and plant height
rownames(traits) <- tolower(traits$species)
traits <- traits[,c(3:5,9:10)]

### ENVIRONMENT - Uses temperature data collected from 2017-09-28 at 00:00 to 2018-09-28 at 00:00 at 25 plots.
env18.core <- env.cvr18[,c(7:12,14:17)]
rownames(env18.core) <- as.character(env.cvr18$plot_id)

# Make comparative community object and calculate diversity metrics from that
core18div <- div.calc(tree18.core, core.18.comm, traits, env18.core)

####################################################################################
### How do functional does functional dispersion and the community weighted mean of the traits vary across temperature and texture?

### Set up leaf area with the full model to test for colinearity
test.la <- with(core18div, lm(log(CWM.LA)~mean+sd+max+min+SAND+CLAY+SILT))
vif(test.la) # There are aliased coefficeints in the model

### CWM of LA
la.corr <- env.corr(core18div, log(core18div$CWM.LA))
xtable(la.corr, digits=3)

la <- with(core18div, lm(log(CWM.LA)~mean+CLAY))
la.mn <- with(core18div, lm(log(CWM.LA)~mean))
anova(la, la.mn) # not significant, plot la.mn
# plot la.mn

### CWM of SLA
sla.corr <- env.corr(core18div, log(core18div$CWM.SLA)) 
xtable(sla.corr, digits=3)

sla <- with(core18div, lm(log(CWM.SLA)~mean+CLAY))
sla.mn <- with(core18div, lm(log(CWM.SLA)~mean))
anova(sla, sla.mn) # not significant, plot sla.mn

### CWM of max HT
mxH.corr <- env.corr(core18div, core18div$CWM.maxht)

xtable(mxH.corr, digits=3)
mxH <- with(core18div, lm(CWM.maxht~max+CLAY)) 
# for plotting
mxH.mx <- with(core18div, lm(CWM.maxht~max))
anova(mxH, mxH.mx) # SIG, table mxH.mx
xtable(mxH, digits=3) 

### CWM of mean HT
mnH.corr <- env.corr(core18div, core18div$CWM.mn.ht)
xtable(mnH.corr, digits=3) 
mnH <- with(core18div, lm(CWM.mn.ht~max+CLAY)) 
# for plotting
mnH.mx <- with(core18div, lm(CWM.mn.ht ~ max))
anova(mnH, mnH.mx) # SIG, table mnH -> SUPP
xtable(mnH, digits=3)

### FDis - with four traits
fdis.4.corr <- env.corr(core18div, core18div$fdis)
xtable(fdis.4.corr, digits=3) 

fdis.4 <- with(core18div, lm(fdis~max+CLAY)) 
fdis.4.mx <- with(core18div, lm(fdis~max))
anova(fdis.4, fdis.4.mx) # SIG, table fdis.4
xtable(fdis.4, digits=3)

### FDIS - with just two traits, max.height and LA
fdis.2.corr <- env.corr(core18div, core18div$fdis.two)
xtable(fdis.2.corr, digits=3)

fdis.2 <- with(core18div, lm(fdis.two~max+CLAY))
fdis.2.mx <- with(core18div, lm(fdis.two~max))
anova(fdis.2, fdis.2.mx) # SIG, table fdis.2 -> SUPP
xtable(fdis.2, digits = 3)

### NOTE: NEW FIGURE 3 is located in fdiv-env_life-hist_core18.R
# ### FIG 3 - Plotting how functional diversity changes across the MEAN temperature variable.
# jpeg("./analysis/figures/fdiv-env-core18.jpeg", width=7, height=7, unit="in",res=300)
# par(mfrow=c(2,2))
# par(mar=c(4,5,0.1,0.5))
# par(oma=c(1.5,1.5,1.5,1.5))
# 
# # LA ~ mean
# with(core18div, plot(log(CWM.LA)~mean, pch=19, xlab="", ylab=expression(italic("ln")("Leaf area")), cex=1.2, axes=FALSE, xlim=c(4,14),ylim=c(-3,4.5), col="#009E73", cex.lab=1.05))
# abline(la.mn, lwd=2)
# axis(1, cex.lab=1.5)
# axis(2, at=c(-3,-1,1,3), cex.lab=1.5)
# mtext(expression(bold("(a)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.66, p-value < 0.001"),side=3,line=-2, adj=2.2, cex = 0.85)
# 
# # max ht ~ sd
# with(core18div, plot(CWM.maxht~max, pch=19, xlab="", ylab="Max. height (cm)", cex=1.2, axes=FALSE, xlim=c(10,70),ylim=c(0,1500), col="#E69F00", cex.lab=1.05))
# abline(mxH.mx, lwd=2)
# axis(1, cex.lab=1.5)
# axis(2, at=c(0,400,800,1200), cex.lab=1.5)
# mtext(expression(bold("(b)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.21, p-value = 0.021"),side=3,line=-2, adj=2.2, cex=0.85)
# 
# # SLA ~ mean
# with(core18div, plot(log(CWM.SLA)~mean, pch=19, xlab=expression(paste('Mean temp. (',degree,'C)')), ylab=expression(italic("ln")("Specific leaf area")), cex=1.2, axes=FALSE, xlim=c(4,14),ylim=c(2,5.7), col="#009E73", cex.lab=1.05))
# abline(sla.mn, lwd=2)
# axis(1, cex.lab=1.5)
# axis(2, at=c(2,3,4,5), cex.lab=1.5)
# mtext(expression(bold("(c)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.31, p-value = 0.004"), side=3,line=-2, adj=2.2, cex=0.85)
# 
# # functional dispersion ~ sd temperature
# with(core18div, plot(fdis~max, pch=19, xlab=expression(paste('Max. temp. (',degree,'C)')), ylab="Functional dispersion", cex=1.2, axes=FALSE, xlim=c(10,70), ylim=c(0,6), col="#0072B2", cex.lab=1.05))
# abline(fdis.4.mx, lwd=2)
# axis(1, cex.lab=1.5)
# axis(2, at=c(0,1,2,3,4,5), cex.lab=1.5)
# mtext(expression(bold("(d)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.26, p-value = 0.010"),side=3,line=-2, adj=2.2, cex=0.85)
# 
# dev.off()

##################################################################################################
# Q3: what is the relationship between functional diversity (CWM and dispersion) and phylodiversity? (25 plots)
test.la.phy <- with(core18div, lm(log(CWM.LA)~sesmntd+sesmpd+pd+nsp+simpsdiv)) # VIF high for pd, sesmntd, and nsp
vif(test.la.phy)

# LA ~ phy
la.phy.cor <- cor.test.phy(core18div, log(core18div$CWM.LA))
xtable(la.phy.cor, digits=3)

la.phy <- with(core18div, lm(log(CWM.LA)~pd+nsp))
la.pd <- with(core18div, lm(log(CWM.LA)~pd))
anova(la.phy, la.pd) # not sig, plot la.pd

# SLA ~ phy
sla.phy.cor <- cor.test.phy(core18div, log(core18div$CWM.SLA))
xtable(sla.phy.cor, digits=3)

sla.phy <- with(core18div, lm(log(CWM.SLA)~pd+nsp))
sla.pd <- with(core18div, lm(log(CWM.SLA)~pd))
anova(sla.phy, sla.pd) # not sig, plot sla.pd

# Max. height ~ phy
mxH.phy.cor <- cor.test.phy(core18div, core18div$CWM.maxht) 
xtable(mxH.phy.cor, digits=3)

mxH.phy <- with(core18div, lm(CWM.maxht~sesmntd+nsp))
mxH.mntd <- with(core18div, lm(CWM.maxht~sesmntd))
anova(mxH.phy, mxH.mntd) # not sig, plot mxH.mntd

# Mean height ~ phy
mnH.phy.cor <- cor.test.phy(core18div, core18div$CWM.mn.ht)
xtable(mnH.phy.cor, digits=3)

mnH.phy <- with(core18div, lm(CWM.mn.ht~sesmntd+nsp))
mnH.mntd <- with(core18div, lm(CWM.mn.ht~sesmntd))
anova(mnH.phy, mnH.mntd) # not sig, plot mnH.mntd -> SUPP

# Fdis ~ phy
fdis.phy.cor <- cor.test.phy(core18div, core18div$fdis)
xtable(fdis.phy.cor, digits=3)

fdis.phy <- with(core18div, lm(fdis~sesmntd+nsp))
fdis.mntd <- with(core18div, lm(fdis~sesmntd))
anova(fdis.phy, fdis.mntd) # not sig, plot fdis.mntd

# Fdis.two ~ phy
fdis.2.phy.cor <- cor.test.phy(core18div, core18div$fdis.two)
xtable(fdis.2.phy.cor, digits=3)
fdis.2.phy <- with(core18div, lm(fdis.two~sesmntd+nsp))
xtable(fdis.2.phy, digits=3)
# for plotting
fdis.2.mntd <- with(core18div, lm(fdis.two~sesmntd))
anova(fdis.2.phy, fdis.2.mntd) # not sig, plot fdis.2.mntd -> SUPP

#################################################################
# FIG 4 - Functional diversity across phy-div metrics at 25 plots
# Plot significant relationships with SESmntd here, do the rest in the supplement
jpeg("./analysis/figures/fdiv-phydiv-core18.jpeg", width=7, height=7, unit="in",res=300)
par(mfrow=c(2,2))
par(mar=c(4,5,0.1,0.5))
par(oma=c(1.5,1.5,1.5,1.5))

# CWM of leaf area across PD
with(core18div, plot(log(CWM.LA)~pd, pch=19, xlab="", ylab=expression(italic("ln")("Leaf area")),  cex=1.2, axes=FALSE, xlim=c(500, 2500), ylim=c(-3,4.7), col="#009E73", cex.lab=1.05))
abline(la.pd, lwd=2)
axis(1, cex.lab=1.5)
axis(2, at=c(-3,-1,1,3), cex.lab=1.5)
mtext(expression(bold("(a)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.35, p-value = 0.002"),side=3,line=-2, adj=2.2, cex=0.85)

# CWM of maximum height across SESmntd
with(core18div, plot(CWM.maxht~sesmntd, pch=19, xlab="", ylab="Max. height", cex=1.2, axes=FALSE,  xlim=c(-3,5), ylim=c(0,1550),col="#E69F00", cex.lab=1.05))
abline(mxH.mntd, lwd=2)
axis(1, at=c(-3,-1, 1, 3, 5), cex.lab=1.5)
axis(2, at=c(0, 400, 800, 1200), cex.lab=1.5)
mtext(expression(bold("(b)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.70, p-value < 0.001"),side=3,line=-2, adj=2.2, cex=0.85)

# CWM of specific leaf area across PD
with(core18div, plot(log(CWM.SLA)~pd, pch=19, xlab="Faith's PD", ylab=expression(italic("ln")("Specific leaf area")), cex=1.2, axes=FALSE, xlim=c(500,2500), ylim=c(2,5.8), col="#009E73", cex.lab=1.05))
abline(sla.pd, lwd=2)
axis(1, cex.lab=1.5)
axis(2, at=c(2,3,4,5), cex.lab=1.5)
mtext(expression(bold("(c)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.21, p-value = 0.020"),side=3,line=-2, adj=2.2, cex=0.85)

# Functional dispersion across SESmntd
with(core18div, plot(fdis~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab="Functional dispersion", cex=1.2, axes=FALSE,  col="#0072B2", xlim=c(-3,5), ylim=c(0,6.3),cex.lab=1.05))
abline(fdis.mntd, lwd=2)
axis(1, at=c(-3,-1,1,3,5), cex.lab=1.5)
axis(2, at=c(0,1,2,3,4,5), cex.lab=1.5)
mtext(expression(bold("(d)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.69, p-value < 0.001"),side=3,line=-2, adj=2.2, cex=0.85)

dev.off()
#######
# mnH ~ mntd and fdis2~mnH for supplement

jpeg("./analysis/figures/supp-fdiv-phydiv-core18.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1.5))
par(oma=c(.5,.5,.5,.5))

# CWM of leaf area across PD
with(core18div, plot(CWM.mn.ht~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab="Mean height (cm)", cex=1.2, axes=FALSE, col="#E69F00", cex.lab=1.05))
abline(mnH.mntd, lwd=2)
axis(1)
axis(2)
mtext(expression(bold("(a)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.66, p-value < 0.001"),side=3,line=-0.5, adj=2.2, cex=0.85)

# CWM of maximum height across SESmntd
with(core18div, plot(fdis.two~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab="FDis - 2 tr.", cex=1.2, axes=FALSE,col="#0072B2", cex.lab=1.05))
abline(fdis.2.mntd, lwd=2)
axis(1)
axis(2)
mtext(expression(bold("(b)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.68, p-value < 0.001"),side=3,line=-0.5, adj=2.2, cex=0.85)

dev.off()