# Assessing variance in functional dispersion and CWM of height and leaf traits across temp and soil microclimates and phylodiversity
# Elizabeth Simpson # 2021-03-25
# Uses code from Will Pearse to get comparative.comm set up

source("~/Documents/projects/functional_traits_rhf/analysis/fdiv-env-analysis_core18.R")

##############################
### 2018 INTENSIFIED PLOTS ###
##############################

# Intensified analysis for 2018 at all 78 plots - just topography and soil texture
int.comm <- calc.comm(int.18.over)

# ENVIRONMENTAL DATA for intensified 2018 plots - topography & texture, 78 plots
env <- read.csv("./clean_data/texture-terrain-18.csv")
int.env <- env[,7:12]
rownames(int.env) <- env$plot_id

# Make phylogeny for species in core plots
int.tree <- congeneric.merge(tree, colnames(int.comm), split = " ")

# Make comparative community object and calculate diversity metrics from that
int.div <- div.calc(int.tree, int.comm, traits, int.env)

#################
### Functions ###
int.corr <- function(div,met){
  output <- c(with(div, cor.test(met,aspect))$estimate, with(div, cor.test(met,elev))$estimate, with(div, cor.test(met,slope))$estimate,
              with(div, cor.test(met, SAND))$estimate, with(div, cor.test(met, SILT))$estimate, with(div, cor.test(met, CLAY))$estimate)
  output <- as.data.frame(rbind(output, c(with(div, cor.test(met,aspect))$p.value, with(div, cor.test(met,elev))$p.value, with(div, cor.test(met,slope))$p.value, 
                                          with(div, cor.test(met, SAND))$p.value, with(div, cor.test(met, SILT))$p.value, with(div, cor.test(met, CLAY))$p.value)))
  colnames(output) <- c("aspect", "elev", "slope", "sand", "silt", "clay")
  rownames(output) <- c("pearsons.r", "p.value")
  return(output)
}

###################################################################
# How does microenvironment vary across topography?
# for 2018 intesified plots (76 total) can look at each texture component across topography

int.sand <- with(int.div, lm(SAND~aspect+elev+slope)) # varies across elevation
int.sand.e <- with(int.div, lm(SAND~elev)) # intercept not sig.

int.silt <- with(int.div, lm(SILT~aspect+elev+slope)) # just intercept
int.silt.1 <- with(int.div, lm(SILT~1))

int.clay <- with(int.div, lm(CLAY~aspect+elev+slope)) # varies across elevation
int.clay.e <- with(int.div, lm(CLAY~elev))

# FOR SUPPLEMENT - plot just the intensified clay / sand /silt figure and the soil texture triangle for all 76 plots
# Plotting changes in texture across environment plus the USDA texture triangle for all 76 plots 
jpeg("./analysis/figures/supp-env-int.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,2))
par(mar=c(5,5,0.5,0.5))
par(oma=c(1,1,1,1))

with(int.div, plot(SAND~elev, pch=19, xlab="Elevation (m.s.l.)", ylab= "Soil component (%)", cex=0.7, cex.lab=0.8, axes=FALSE, col="goldenrod3", xlim=c(1700,2100), ylim=c(0,80)))
abline(int.sand.e, lwd=3, col="goldenrod3")
with(int.div, points(CLAY~elev, pch=19, cex=0.7, col="sienna"))
abline(int.clay.e, lwd=3, col="sienna")
axis(1, cex.axis=0.8)
axis(2, cex.axis=0.8)
legend(2010,80, legend=c("Sand", "Clay"), pch=16, pt.cex=0.9, cex=0.6, col = c("goldenrod3", "sienna"))
legend(1670,70, expression(~ R^2 ~ "= 0.13, p-value = 0.002"), text.col="goldenrod3", box.lty=0, bg="transparent", cex=0.5)
legend(1880,24, expression(~ R^2 ~ "= 0.16, p-value < 0.001"), text.col="sienna", box.lty=0, bg="transparent", cex=0.5)
mtext(expression(bold("(a)")),side=3,line=-0.2, adj=-0.4,)

TT.plot(class.sys = "USDA.TT", tri.data = int.div, pch=19, cex.axis=0.6, cex.lab=0.6, lwd=0.6, lwd.axis=0.6, lwd.lab=0.6, main="", cex=0.7, new.mar=c(3.7,3.7,0,0))
mtext(expression(bold("(b)")),side=3,line=-0.9, adj=-0.1)

dev.off()

####################################################################################
### How do functional does functional dispersion and the community weighted mean of the traits vary across topography and texture at intensified plots?

### CWM of SLA
int.sla <- with(int.div, lm(log(CWM.SLA)~aspect+elev+slope+SAND+SILT+CLAY))
vif(int.sla)
# There are aliased coefficients in the model

int.sla.cor <- int.corr(int.div, log(int.div$CWM.SLA)) 
xtable(int.sla.cor, digits=3) # for supplement
int.sla <- with(int.div, lm(log(CWM.SLA)~aspect+elev+slope+SAND+SILT)) # aspect sig.
xtable(sla, digits=3) # for supplement
# for plotting
int.sla.a <- with(int.div, lm(log(CWM.SLA)~aspect))

### CWM of LA
int.la.cor <- int.corr(int.div, log(int.div$CWM.LA))
xtable(int.la.corr, digits=3)
int.la <- with(int.div, lm(log(CWM.LA)~aspect+elev+slope+SAND+SILT)) # aspect sig.
xtable(int.la, digits=3)
# for plotting
int.la.a <- with(int.div, lm(log(CWM.LA)~aspect)) 

### CWM of mean HT
int.mnH.cor <- int.corr(int.div, int.div$CWM.mn.ht)
xtable(int.mnH.corr, digits=3)
int.mnH <- with(int.div, lm(CWM.mn.ht~aspect+elev+slope+SAND+CLAY)) # aspect and elev
xtable(int.mnH, digits=3)
# for plotting
int.mnH.a <- with(int.div, lm(CWM.mn.ht~aspect))
int.mnH.e <- with(int.div, lm(CWM.mn.ht~elev)) # not significant on it's own

### CWM of max HT
int.mxH.cor <- int.corr(int.div, int.div$CWM.maxht)
xtable(int.mxH.corr, digits=3)
int.mxH <- with(int.div, lm(CWM.maxht~aspect+elev+slope+SILT+CLAY)) # aspect and elev
xtable(mxH, digits=3)
# for plotting
int.mxH.a <- with(int.div, lm(CWM.maxht~aspect))
int.mxH.e <- with(int.div, lm(CWM.maxht~elev)) # almost sig, but not quite

### FDis - with four traits
int.fdis.cor <- int.corr(int.div, int.div$fdis)
xtable(int.fdis.cor, digits=3)
int.fdis <- with(int.div, lm(fdis~aspect+elev+slope+SAND+SILT)) # just aspect
xtable(int.fdis, digits=3)
# for plotting
int.fdis.a <- with(int.div, lm(fdis~aspect))

### FIG 3 - Plotting how functional diversity changes across aspect at the intensified sites
jpeg("./analysis/figures/supp-fdiv-env-int.jpeg", width=7, height=10.5, unit="in",res=300)
par(mfrow=c(3,2))
par(mar=c(4,5.5,0.1,0.5))
par(oma=c(1.5,1.5,1.5,1.5))

# LA ~ aspect - int
with(int.div, plot(log(CWM.LA)~aspect, pch=19, xlab="", ylab=expression(italic("ln")("Leaf area")), cex=1.2, axes=FALSE, col="#009E73", cex.lab=1.4, ylim=c(-3,7.4)))
abline(int.la.a, lwd=2)
axis(1, cex.axis=1.3)
axis(2, at=c(-3,0,3,6), cex.axis=1.4)
mtext(expression(bold("(a)") ~~~~~~~~~~~~~~ R^2 ~ "= 0.36, p-value < 0.001"),side=3,line=-2, adj=3)

# max ht ~ aspect - int
with(int.div, plot(CWM.maxht~aspect, pch=19, xlab="", ylab="Max. height (cm)", cex=1.2, axes=FALSE, col="#E69F00", cex.lab=1.4, ylim=c(0,1390)))
abline(int.mxH.a, lwd=2)
axis(1, cex.axis=1.3)
axis(2, at=c(0,400,800,1200),cex.axis=1.2)
mtext(expression(bold("(b)") ~~~~~~~~~~~~~~ R^2 ~ "= 0.10, p-value = 0.005"),side=3,line=-2, adj=3)

# SLA ~ aspect - int
with(int.div, plot(log(CWM.SLA)~aspect, pch=19, xlab="", ylab=expression(italic("ln")("Specific leaf area")), cex=1.2, axes=FALSE, col="#009E73", cex.lab=1.4, ylim=c(2,6.7)))
abline(int.sla.a, lwd=2)
axis(1, cex.axis=1.3)
axis(2, at=c(2,4,6), cex.axis=1.3)
mtext(expression(bold("(c)") ~~~~~~~~~~~~~~ R^2 ~ "= 0.28, p-value < 0.001"), side=3,line=-2, adj=3)

# max ht ~ aspect - int
with(int.div, plot(CWM.mn.ht~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab="Mean height (cm)", cex=1.2, axes=FALSE, col="#E69F00", cex.lab=1.4, ylim=c(0,950)))
abline(int.mnH.a, lwd=2)
axis(1, cex.axis=1.3)
axis(2, at=c(0,200,400,600,800), cex.axis=1.2)
mtext(expression(bold("(d)") ~~~~~~~~~~~~~~ R^2 ~ "= 0.09, p-value = 0.010"),side=3,line=-2, adj=3)

# functional dispersion ~ aspect - int
with(int.div, plot(fdis~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab="Functional dispersion", cex=1.2, axes=FALSE, col="#0072B2", cex.lab=1.4, ylim=c(0,6)))
abline(int.fdis.a, lwd=2)
axis(1, cex.axis=1.3)
axis(2, at=c(0,1,2,3,4,5),cex.axis=1.3)
mtext(expression(bold("(e)") ~~~~~~~~~~~~~~ R^2 ~ "= 0.24, p-value < 0.001"),side=3,line=-2, adj=3)

dev.off()

##################################################################################################
# Q3: How does functional diversity (dispersion and CWM) correlate with phylodiversity? (25 plots)
fdis.phy.int <- with(int.div, lm(fdis~sesmntd+sesmpd+pd+nsp+simpsdiv))
vif(fdis.phy.int) # colinearity maybe a little high for pd and nsp

# fdis ~ phy - int
fdis.phy.int.cor <- cor.test.phy(int.div, int.div$fdis) #pull out nsp
xtable(fdis.phy.int.cor)
fdis.phy.int <- with(int.div, lm(fdis~sesmntd+sesmpd+pd+simpsdiv)) # pd sig
xtable(fdis.phy.int, digits=3)
# for plotting
fdis.int.pd <- with(int.div, lm(fdis~pd))

# LA ~ phy - int
la.phy.int.cor <- cor.test.phy(int.div, log(int.div$CWM.LA))
xtable(la.phy.int.cor, digits=3)
la.phy.int <- with(int.div, lm(log(CWM.LA)~sesmntd+sesmpd+pd+simpsdiv)) # mntd, mpd, pd
xtable(la.phy.int, digits=3)
# for plotting
la.int.mntd <- with(int.div, lm(log(CWM.LA)~sesmntd))
la.int.mpd <- with(int.div, lm(log(CWM.LA)~sesmpd))
la.int.pd <- with(int.div, lm(log(CWM.LA)~pd))

# SLA ~ phy - int
sla.phy.int.cor <- cor.test.phy(int.div, log(int.div$CWM.SLA))
xtable(sla.phy.int.cor, digits=3)
sla.phy.int <- with(int.div, lm(log(CWM.SLA)~sesmntd+sesmpd+pd+simpsdiv)) # mntd
xtable(sla.phyint, digits=3)
# for plotting
sla.int.mntd <- with(int.div, lm(log(CWM.SLA)~sesmntd))

# mxH ~ phy - int
mxH.phy.int.cor <- cor.test.phy(int.div, int.div$CWM.maxht) 
xtable(mxH.phy.int.cor, digits=3)
mxH.phy.int <- with(int.div, lm(CWM.maxht~sesmntd+sesmpd+pd+simpsdiv)) # mpd, pd
xtable(mxH.phy.int, digits=3)
# for plotting
mxH.int.mpd <- with(int.div, lm(CWM.maxht~sesmpd))
mxH.int.pd <- with(int.div, lm(CWM.maxht~pd))

# mnH ~ phy - int
mnH.phy.int.cor <- cor.test.phy(int.div, int.div$CWM.mn.ht)
xtable(mnH.phy.int.cor, digits=3)
mnH.phy.int <- with(int.div, lm(CWM.mn.ht~sesmntd+sesmpd+pd+simpsdiv)) # mpd, pd
xtable(mnH.phy.int, digits=3)
# for plotting
mnH.int.mpd <- with(int.div, lm(CWM.mn.ht~sesmpd))
mnH.int.pd <- with(int.div, lm(CWM.mn.ht~pd))

# phy.cor <- rbind(fdis.phy.cor, la.phy.cor, sla.phy.cor, mxH.phy.cor, mnH.phy.cor)
# xtable(phy.cor, digits=3)
# 
# xtable(fdis.phy, digits=3)
# xtable(la.phy, digits=3)
# xtable(sla.phy, digits=3)
# xtable(mxH.phy, digits=3)
# xtable(mnH.phy, digits=3)

#################################################################
# FIG 4 - Functional diversity across phy-div metrics at 25 plots
# Plot significant relationships with SESmntd here, do the rest in the supplement
jpeg("./analysis/figures/supp-fdiv-phydiv-int.jpeg", width=7, height=7, unit="in",res=300)
par(mfrow=c(3,3))
par(mar=c(5,6,0.1,0.5))
par(oma=c(0.5,0.5,0.5,0.5))

# CWM of leaf area across PD - int
with(int.div, plot(log(CWM.LA)~pd, pch=19, xlab="Faith's PD", ylab=expression(italic("ln")("Leaf area")),  axes=FALSE, col="#009E73", xlim=c(500,2500), ylim=c(-3,6.5), cex.lab=0.95))
abline(la.int.pd, lwd=2)
axis(1, cex.axis=0.8)
axis(2, at=c(-3,-1,1,3,5),cex.axis=0.9)
mtext(expression(bold("(a)") ~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.33, p-value < 0.001"), side=3, line=-1.2, adj=1.3, cex=0.6)

# CWM of leaf area across mpd - int
with(int.div, plot(log(CWM.LA)~sesmpd, pch=19, xlab=expression(SES[MPD]), ylab=expression(italic("ln")("Leaf area")),  axes=FALSE, col="#009E73", xlim=c(-2,2.1), ylim=c(-3,6.5), cex.lab=0.95))
abline(la.int.mpd, lwd=2)
axis(1, cex.axis=0.8)
axis(2, at=c(-3,-1,1,3,5),cex.axis=0.9)
mtext(expression(bold("(b)") ~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.06, p-value = 0.032"), side=3, line=-1.2, adj=1.3, cex=0.6)

# CWM of leaf area across MNTD
with(int.div, plot(log(CWM.LA)~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab=expression(italic("ln")("Leaf area")),  axes=FALSE, col="#009E73", ylim=c(-3,6.5), cex.lab=0.95))
abline(la.int.mntd, lwd=2)
axis(1, cex.axis=0.8)
axis(2, at=c(-3,-1,1,3,5),cex.axis=0.9)
mtext(expression(bold("(c)") ~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.34, p-value < 0.001"), side=3, line=-1.2, adj=1.3, cex=0.6)

# CWM of maximum height across pd - int
with(int.div, plot(CWM.maxht~pd, pch=19, xlab="Faith's PD", ylab="Max. height (cm)", axes=FALSE, col="#E69F00", xlim=c(500,2500),ylim=c(0,1450), cex.lab=0.95))
abline(mxH.int.pd, lwd=2)
axis(1, cex.axis=0.9)
axis(2, at=c(0,400,800,1200), cex.axis=0.8)
mtext(expression(bold("(d)") ~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.35, p-value < 0.001"), side=3, line=-1.2, adj=1.3, cex=0.6)

# CWM of maximum height across mpd - int
with(int.div, plot(CWM.maxht~sesmpd, pch=19, xlab=expression(SES[MPD]), ylab="Max. height (cm)", axes=FALSE, col="#E69F00", xlim=c(-2,2.1), ylim=c(0,1450), cex.lab=0.95))
abline(mxH.int.mpd, lwd=2)
axis(1, cex.axis=0.9)
axis(2, at=c(0,400,800,1200), cex.axis=0.8)
mtext(expression(bold("(e)") ~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.39, p-value < 0.001"), side=3, line=-1.2, adj=1.3, cex=0.6)

# CWM of SLA across MNTD
with(int.div, plot(log(CWM.SLA)~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab=expression(italic("ln")("Specific leaf area")),  axes=FALSE, col="#009E73", ylim=c(1,5.8),cex.lab=0.95))
abline(sla.int.mntd, lwd=2)
axis(1, cex.axis=0.8)
axis(2, at=c(1,2,3,4,5), cex.axis=0.9)
mtext(expression(bold("(f)") ~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.31, p-value < 0.001"), side=3, line=-1.2, adj=1.3, cex=0.6)

# CWM of mean height across pd - int
with(int.div, plot(CWM.mn.ht~pd, pch=19, xlab="Faith's PD", ylab="Mean height (cm)", axes=FALSE, col="#E69F00", xlim=c(500,2500), ylim=c(0,950), cex.lab=0.95))
abline(mnH.int.pd, lwd=2)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.85)
mtext(expression(bold("(g)") ~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.32, p-value < 0.001"), side=3, line=-1.2, adj=1.3, cex=0.6)

# CWM of mean height across mpd - int
with(int.div, plot(CWM.mn.ht~sesmpd, pch=19, xlab=expression(SES[MPD]), ylab="Mean height (cm)", axes=FALSE, col="#E69F00", xlim=c(-2,2.1), ylim=c(0,950),cex.lab=0.95))
abline(mnH.int.mpd, lwd=2)
axis(1, cex.axis=0.9)
axis(2, cex.axis=0.85)
mtext(expression(bold("(h)") ~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.38, p-value < 0.001"), side=3, line=-1.2, adj=1.3, cex=0.6)

# Functional dispersion across pd - int
with(int.div, plot(fdis~pd, pch=19, xlab="Faith's PD", ylab="Functional dispersion", axes=FALSE, col="#0072B2", xlim=c(500,2500), ylim=c(0,6), cex.lab=0.95))
abline(fdis.int.pd, lwd=2)
axis(1, cex.axis=0.9)
axis(2, at=c(0,1,2,3,4,5),cex.axis=0.9)
mtext(expression(bold("(i)") ~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.43, p-value < 0.001"), side=3, line=-1.2, adj=1.3, cex=0.6)

dev.off()