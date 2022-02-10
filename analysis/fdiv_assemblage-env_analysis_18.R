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
library(xtable)

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
  output <- as.data.frame(rbind(output, c(with(div, cor.test(RV, aspect))$p.value, with(div, cor.test(RV, elev))$p.value, with(div, cor.test(RV, slope))$p.value,
                                           with(div, cor.test(RV, SAND))$p.value, with(div, cor.test(RV, SILT))$p.value, with(div, cor.test(RV, CLAY))$p.value)))
  colnames(output) <- c("aspect","elev","slope","sand","silt","clay")
  rownames(output) <- c("pearsons.r", "p.value")
  return(output)
}

cor.test.soil <<- function(div,RV){
  output <- c(with(div, cor.test(RV,aspect))$estimate, with(div, cor.test(RV, elev))$estimate, with(div, cor.test(RV, slope))$estimate)
  output <- as.data.frame(rbind(output, c(with(div, cor.test(RV, aspect))$p.value, with(div, cor.test(RV, elev))$p.value, with(div, cor.test(RV, slope))$p.value)))
  colnames(output) <- c("aspect","elev","slope")
  rownames(output) <- c("pearsons.r", "p.value")
  return(output)
}

mean.corr <- cor.test.env(div, div$mean)
sd.corr <- cor.test.env(div, div$sd)
min.corr <- cor.test.env(div, div$min)
max.corr <- cor.test.env(div, div$max)

temp.cors <- rbind(mean.corr, sd.corr, max.corr, min.corr)
xtable(temp.cors, digits=3) # tranfer over to latex to keep working with this

sand.corr <- cor.test.soil(div, div$SAND)
silt.corr <- cor.test.soil(div, div$SILT)
clay.corr <- cor.test.soil(div, div$CLAY)

soil.cors <- rbind(sand.corr, silt.corr, clay.corr)
xtable(soil.cors, digits=3)

# Models for plotting
mn.a <- with(div, lm(mean~aspect))
min <- with(div, lm(min~1))
mx.a <- with(div, lm(max~aspect))
sd.a <- with(div, lm(sd~aspect))

sand.e <- with(div, lm(SAND~elev))
clay.e <- with(div, lm(CLAY~elev))

# FIG 2 - Plotting changes in temperature across environment plus the USDA texture triangle for all of the plots 
jpeg("./analysis/figures/env.jpeg", width=7, height=7, unit="in",res=300)
par(mfrow=c(2,2))
par(mar=c(5,5,0.5,0.5))
par(oma=c(1,1,1,1))

with(div, plot(mean~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(paste('Temperature (',degree~'C)')), cex=1.2, axes=FALSE, ylim=c(-10,70), col="#E69F00"))
abline(mn.a, lwd=3, col="#E69F00")
with(div, points(max~aspect, pch=19, cex=1.2, col="#D55E00"))
abline(mx.a, lwd=3, col="#D55E00")
with(div, points(min~aspect, pch=19, cex=1.2, col="#56B4E9"))
abline(min, lwd=3, col="#56B4E9")
axis(1)
axis(2, at=c(-10,10,30,50,70))
mtext("(a)",side=3,line=-0.8, adj=0.05)
legend(-0.5,62, expression(~ R^2 ~ "= 0.35, p-value = 0.002"), text.col="#D55E00", box.lty=0, bg="transparent", cex=0.9)
legend(-1.2,25, expression(~ R^2 ~ "= 0.68, p-value < 0.001"), text.col= "#E69F00", box.lty=0, bg="transparent", cex=0.9)

with(div, plot(sd~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(paste('SD - Temperature (',degree~'C)')), cex=1.2, axes=FALSE, col="#E69F00", ylim=c(5,15)))
abline(sd.a, lwd=3, col="#E69F00")
axis(1)
axis(2, at=c(5,10,15))
mtext("(b)",side=3,line=-0.8, adj=0.05)
legend(0.5,14.5, legend=c('Min', 'Mean', 'Max'), pch=16, pt.cex=1.2, cex=0.6, col = c('#56B4E9', '#E69F00', '#D55E00'))
legend(-1.2,6.5, expression(~ R^2 ~ "= 0.44, p-value < 0.001"), text.col="#E69F00", box.lty=0, bg="transparent", cex=0.9)

with(div, plot(SAND~elev, pch=19, xlab="Elevation (m.s.l.)", ylab= "Percent of soil component", cex=1.2, axes=FALSE, col="goldenrod3", xlim=c(1700,2100), ylim=c(0,80)))
abline(sand.e, lwd=3, col="goldenrod3")
with(div, points(CLAY~elev, pch=19, cex=1.2, col="sienna"))
abline(clay.e, lwd=3, col="sienna")
axis(1)
axis(2)
legend(2010,80, legend=c("Sand", "Clay"), pch=16, pt.cex=1.2, cex=0.6, col = c("goldenrod3", "sienna"))
legend(1670,65, expression(~ R^2 ~ "= 0.21, p-value = 0.024"), text.col="goldenrod3", box.lty=0, bg="transparent", cex=0.9)
legend(1810,30, expression(~ R^2 ~ "= 0.23, p-value = 0.016"), text.col="sienna", box.lty=0, bg="transparent", cex=0.9)
mtext("(c)",side=3,line=-0.8, adj=0.05,)

TT.plot(class.sys = "USDA.TT", tri.data = div, pch=19, cex.axis=0.8, cex.lab=0.8, lwd=0.8, lwd.axis=0.8, lwd.lab=0.8, main="", cex=0.95, new.mar=c(3.7,3.7,0,0))
mtext("(d)",side=3,line=-1.6, adj=0.05)

dev.off()

####################################################################################
### How do functional does functional dispersion and the community weighted mean of the traits vary across temperature and texture?

### CWM of SLA
sla <- with(div, lm(log(CWM.SLA)~mean+CLAY))
xtable(sla)
# for plotting
sla.mn <- with(div, lm(log(CWM.SLA)~mean))

### CWM of LA
la <- with(div, lm(log(CWM.LA)~mean+CLAY))
xtable(la)
# for plotting
la.mn <- with(div, lm(log(CWM.LA)~mean)) 

### CWM of mean HT
mnH <- with(div, lm(CWM.mn.ht~mean+CLAY))
xtable(mnH)
# for plotting
mnH.mn <- with(div, lm(CWM.mn.ht~mean)) 
mnH.c <- with(div, lm(CWM.mn.ht~CLAY)) 

### CWM of max HT
mxH <- with(div, lm(CWM.maxht~mean+CLAY))
xtable(mxH)
# for plotting
mxH.mn <- with(div, lm(CWM.maxht~mean))
mxH.c <- with(div, lm(CWM.maxht~CLAY))

### FDis - with four traits
fdis <- with(div, lm(fdis~mean+CLAY))
xtable(fdis)
# for plotting
fdis.mn <- with(div, lm(fdis~mean))
fdis.c <- with(div, lm(fdis)~CLAY)
             
### FDIS - with just two traits, max.height and LA
fdis.two <- with(div, lm(fdis.two~mean+CLAY))
xtable(fdis.two)
# for plotting
fdistwo.mn <- with(div, lm(fdis.two~mean))
fdistwo.c <- with(div, lm(fdis.two~CLAY))

### FIG 3 - Plotting how functional diversity changes across the MEAN temperature variable.
jpeg("./analysis/figures/fdiv-MEAN-25plots.jpeg", width=7, height=7, unit="in",res=300)
par(mfrow=c(2,2))
par(mar=c(4,5,0.1,0.5))
par(oma=c(1,1,1,1))

# LA ~ mean
with(div, plot(log(CWM.LA)~mean, pch=19, xlab="", ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, xlim=c(4,14), ylim=c(-2.5,3.5), col="#009E73"))
abline(la.mn, lwd=3)
axis(1)
axis(2)
mtext("(a)",side=3,line=-1.2, adj=0.05, cex=0.8)
legend("topright", expression(R^2 ~ "= 0.67, p-value < 0.001"),cex=0.9, bty="n")

# functional dispersion ~ mean temperature
with(div, plot(fdis~mean, pch=19, xlab="", ylab="Functional dispersion", cex=1.2, axes=FALSE, xlim=c(4,14), ylim=c(0,5),col="#0072B2"))
abline(fdis.mn, lwd=3)
axis(1)
axis(2)
mtext("(b)",side=3,line=-1.2, adj=0.05, cex=0.8)
legend("topright", expression(R^2 ~ "= 0.24, p-value = 0.002"),cex=0.9, bty="n")

# SLA ~ mean
with(div, plot(log(CWM.SLA)~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab=expression(italic("ln")(CWM~-~SLA)), cex=1.2, axes=FALSE, xlim=c(4,14), ylim=c(2,5), col="#009E73"))
abline(sla.mn, lwd=3)
axis(1)
axis(2)
mtext("(c)",side=3,line=-1.2, adj=0.05, cex=0.8)
legend("topright", expression(R^2 ~ "= 0.32, p-value < 0.001"),cex=0.9, bty="n")

# max ht ~ mean
with(div, plot(CWM.maxht~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab="CWM - max. height (cm)", cex=1.2, axes=FALSE, xlim=c(4,14), ylim=c(0,1400), col="#E69F00"))
abline(mxH.mn, lwd=3)
axis(1)
axis(2, cex.axis=0.8)
mtext("(d)",side=3,line=-1.2, adj=0.05, cex=0.8)
legend("topright", expression(R^2 ~ "= 0.17, p-value = 0.040"),cex=0.9, bty="n")
legend(9,1200, legend=c('Overall functional diversity', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1.2, cex=0.6,
       col = c('#0072B2', '#009E73', '#E69F00'))

dev.off()

#########################################################################
# How does functional diversity (dispersion and CWM) correlate with phylodiversity? (25 plots)
cor.test.phy <<- function(div,RV){
  output <- c(with(div, cor.test(RV,sesmntd))$estimate, with(div, cor.test(RV, sesmpd))$estimate, with(div, cor.test(RV, pd))$estimate,
  with(div, cor.test(RV, nsp))$estimate, with(div, cor.test(RV, simpsdiv))$estimate)
  output <- as.data.frame(rbind(output, c(with(div, cor.test(RV, sesmntd))$p.value, with(div, cor.test(RV, sesmpd))$p.value, with(div, cor.test(RV, pd))$p.value,
                                          with(div, cor.test(RV, nsp))$p.value, with(div, cor.test(RV, simpsdiv))$p.value)))
  colnames(output) <- c("sesmntd","sesmpd","pd","nsp","simpsdiv")
  rownames(output) <- c("pearsons.r", "p.value")
  return(output)
}

fdis.phyall <- with(div, lm(fdis~sesmntd+sesmpd+pd+nsp+simpsdiv)) # VIF high for pd, sesmntd, and nsp
fdis.phycor <- cor.test.phy(div, div$fdis) # mntd most correlated of those three
fdis.phy <- with(div, lm(fdis~sesmntd+sesmpd+simpsdiv)) # mntd* + mpd
# for plotting 
fdis.mntd <- with(div, lm(fdis~sesmntd))

la.phyall <- with(div, lm(log(CWM.LA)~sesmntd+sesmpd+pd+nsp+simpsdiv)) # VIF high for pd, sesmntd, and nsp
la.phycor <- cor.test.phy(div, log(div$CWM.LA)) # pd most correlated of those three
la.phy <- with(div, lm(log(CWM.LA)~pd+sesmpd+simpsdiv)) # pd*
# for plotting
la.pd <- with(div, lm(log(CWM.LA)~pd))

sla.phyall <- with(div, lm(log(CWM.SLA)~sesmntd+sesmpd+pd+nsp+simpsdiv)) # VIF high for pd, sesmntd, and nsp
sla.phycor <- cor.test.phy(div, log(div$CWM.SLA)) # pd most correlated of those three
sla.phy <- with(div, lm(log(CWM.SLA)~pd+sesmpd+simpsdiv)) # only the intercept
# no plotting

mxH.phyall <- with(div, lm(CWM.maxht~sesmntd+sesmpd+pd+nsp+simpsdiv)) # VIF high for pd, sesmntd, and nsp
mxH.phycor <- cor.test.phy(div, div$CWM.maxht) # mntd most correlated of those three
mxH.phy <- with(div, lm(CWM.maxht~sesmntd+sesmpd+simpsdiv)) # mntd* + mpd
# for plotting
mxH.mntd <- with(div, lm(CWM.maxht~sesmntd))

mnH.phyall <- with(div, lm(CWM.mn.ht~sesmntd+sesmpd+pd+nsp+simpsdiv)) # VIF high for pd, sesmntd, and nsp
mnH.phycor <- cor.test.phy(div, div$CWM.mn.ht) # mntd most correlated of those three
mnH.phy <- with(div, lm(CWM.mn.ht~sesmntd+sesmpd+simpsdiv)) # mntd* + mpd
# for plotting
mnH.mntd <- with(div, lm(CWM.mn.ht~sesmntd))

phycor <- rbind(fdis.phycor, la.phycor, sla.phycor, mxH.phycor, mnH.phycor)
xtable(phycor, digits=3)

xtable(fdis.phy, digits=3)
xtable(la.phy, digits=3)
xtable(sla.phy, digits=3)
xtable(mxH.phy, digits=3)
xtable(mnH.phy, digits=3)

#################################################################
# FIG 4 - Functional diversity across phy-div metrics at 25 plots
# Plot significant relationships with SESmntd here, do the rest in the supplement
jpeg("./analysis/figures/fdiv-phydiv-25plots.jpeg", width=7, height=2.5, unit="in",res=300)
par(mfrow=c(1,3))
par(mar=c(4,5,0.1,0.5))
par(oma=c(1,1,1,1))

# Functional dispersion across SESmntd
with(div, plot(fdis~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab="Functional dispersion", cex=1.2, axes=FALSE, ylim=c(0,5), col="#0072B2"))
abline(fdis.mntd, lwd=3)
axis(1)
axis(2)
mtext("(a)",side=3,line=-1.2, adj=0.05, cex=0.8)
legend(-1.5,4.5, expression(R^2 ~ "= 0.69, p-value < 0.001"),cex=0.9, bty="n")

# CWM of maximum height across SESmntd
with(div, plot(CWM.maxht~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab="Maximum height (cm)", cex=1.2, axes=FALSE, ylim=c(0,1400), col="#E69F00"))
abline(mxH.mntd, lwd=3)
axis(1)
axis(2)
mtext("(b)",side=3,line=-1.2, adj=0.05, cex=0.8)
legend(-2,1200, expression(R^2 ~ "= 0.71, p-value < 0.001"),cex=0.9, bty="n")

# CWM of leaf area across PD
with(div, plot(log(CWM.LA)~pd, pch=19, xlab="Faith's PD", ylab=expression(italic("ln")(CWM~-~LA)), cex=1.2, axes=FALSE, ylim=c(-3,4),col="#009E73"))
abline(la.pd, lwd=3)
axis(1)
axis(2)
mtext("(c)",side=3,line=-1.2, adj=0.05, cex=0.8)
legend("topright", expression(R^2 ~ "= 0.34, p-value = 0.002"),cex=0.9, bty="n")

legend("bottomright", legend=c('Overall FD', 'Leaf traits', 'Height traits'), pch=19, pt.cex=1, cex=.8, 
       col = c('#0072B2', '#009E73', '#E69F00'))

dev.off()