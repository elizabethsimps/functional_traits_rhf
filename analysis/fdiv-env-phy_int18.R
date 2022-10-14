# Assessing variance in functional dispersion and CWM of height and leaf traits across temp and soil microclimates and phylodiversity
# Elizabeth Simpson # 2021-03-25
# Uses code from Will Pearse to get comparative.comm set up

source("~/Documents/projects/functional_traits_rhf/analysis/fdiv-env-phy_core18.R")

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

with(int.div, plot(SAND~elev, pch=19, xlab="Elevation (m.s.l.)", ylab= "Soil component (%)", cex=0.6, cex.lab=0.8, axes=FALSE, col="goldenrod3", xlim=c(1700,2100), ylim=c(0,80)))
abline(int.sand.e, lwd=3, col="goldenrod3")
with(int.div, points(CLAY~elev, pch=19, cex=0.7, col="sienna"))
abline(int.clay.e, lwd=3, col="sienna")
axis(1, cex.axis=0.8)
axis(2, cex.axis=0.8)
legend(2010,80, legend=c("Sand", "Clay"), pch=16, pt.cex=0.9, cex=0.6, col = c("goldenrod3", "sienna"))
legend(1670,70, expression(~ R^2 ~ "= 0.13, p-value = 0.002"), text.col="goldenrod3", box.lty=0, bg="transparent", cex=0.5)
legend(1880,24, expression(~ R^2 ~ "= 0.16, p-value < 0.001"), text.col="sienna", box.lty=0, bg="transparent", cex=0.5)
mtext(expression(bold("(a)")),side=3,line=-0.2, adj=-0.4,)

TT.plot(class.sys = "USDA.TT", tri.data = int.div, pch=19, cex.axis=0.6, cex.lab=0.6, lwd=0.6, lwd.axis=0.6, lwd.lab=0.6, main="", cex=0.6, new.mar=c(3.7,3.7,0,0))
mtext(expression(bold("(b)")),side=3,line=-0.9, adj=-0.1)

dev.off()

####################################################################################
### How do functional does functional dispersion and the community weighted mean of the traits vary across topography and texture at intensified plots?

### CWM of LA
int.la.cor <- int.corr(int.div, log(int.div$CWM.LA))

int.la <- with(int.div, lm(log(CWM.LA)~aspect+elev+slope+SAND+SILT)) # aspect
int.la.a <- with(int.div, lm(log(CWM.LA)~aspect)) 
anova(int.la, int.la.a) # not sig. plot int.sla.a -> supp fig

### CWM of SLA 
int.sla.cor <- int.corr(int.div, log(int.div$CWM.SLA))

int.sla <- with(int.div, lm(log(CWM.SLA)~aspect+elev+slope+SAND+SILT)) # aspect
int.sla.a <- with(int.div, lm(log(CWM.SLA)~aspect))
anova(int.sla, int.sla.a) # not sig. plot int.sla.a -> supp fig

### CWM of max HT
int.mxH.cor <- int.corr(int.div, int.div$CWM.maxht)

int.mxH <- with(int.div, lm(CWM.maxht~aspect+elev+slope+SILT+CLAY)) # aspect and elev
int.mxH.ae <- with(int.div, lm(CWM.maxht~aspect+elev))
anova(int.mxH, int.mxH.ae) # not sig, do simpler model -> supp table
xtable(int.mxH.ae, digits = 3)

### CWM of mean HT
int.mnH.cor <- int.corr(int.div, int.div$CWM.mn.ht)

int.mnH <- with(int.div, lm(CWM.mn.ht~aspect+elev+slope+SAND+CLAY)) # aspect and elev
int.mnH.ae  <- with(int.div, lm(CWM.mn.ht~aspect+elev))
anova(int.mnH, int.mnH.ae) # not sig, do simpler model -> supp table
xtable(int.mnH.ae, digits = 3)

### FDis - with four traits
int.fdis.cor <- int.corr(int.div, int.div$fdis)

int.fdis <- with(int.div, lm(fdis~aspect+elev+slope+SAND+SILT)) # just aspect
int.fdis.a <- with(int.div, lm(fdis~aspect))
anova(int.fdis, int.fdis.a) # not sig. plot int.fdis.a -> supp fig

### FIG 3 - Plotting how functional diversity changes across aspect at the intensified sites
jpeg("./analysis/figures/supp-fdiv-env-int.jpeg", width=7, height=2.8, unit="in",res=300)
par(mfrow=c(1,3))
par(mar=c(4,5.5,0.1,0.5))
par(oma=c(1.5,1.5,1.5,1.5))

# LA ~ aspect - int
with(int.div, plot(log(CWM.LA)~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(italic("ln")("Leaf area")),  axes=FALSE, col="#009E73",  ylim=c(-3,7.4), cex.lab=0.9))
abline(int.la.a, lwd=2)
axis(1, cex.axis=0.9)
axis(2, at=c(-3,0,3,6))
mtext(expression(bold("(a)")),side=3,line=-1.5, adj=-0.5, cex=0.7)
legend(-1,5,expression(R^2 ~ "= 0.36, p-value < 0.001"), bty="n", cex=0.7)

# SLA ~ aspect - int
with(int.div, plot(log(CWM.SLA)~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab=expression(italic("ln")("Specific leaf area")),  axes=FALSE, col="#009E73", ylim=c(2,6.7), cex.lab=0.9))
abline(int.sla.a, lwd=2)
axis(1, cex.axis=0.9)
axis(2, at=c(2,4,6))
mtext(expression(bold("(b)")), side=3,line=-1.5, adj=-0.5, cex = 0.7)
legend(-1,6,expression(R^2 ~ "= 0.28, p-value < 0.001"), bty="n",  cex=0.7)

# functional dispersion ~ aspect - int
with(int.div, plot(fdis~aspect, pch=19, xlab="Aspect (S = -1, N = 1)", ylab="Functional dispersion",  axes=FALSE, col="#0072B2", ylim=c(0,6), cex.lab=0.9))
abline(int.fdis.a, lwd=2)
axis(1, cex.axis=0.9)
axis(2, at=c(0,1,2,3,4,5))
mtext(expression(bold("(c)")),side=3,line=-1.5, adj=-0.5, cex=0.7)
legend(-1,4.3, (R^2 ~ "= 0.24, p-value < 0.001"),  bty="n", cex = 0.7)

dev.off()

##################################################################################################
# Q3: How does functional diversity (dispersion and CWM) correlate with phylodiversity? (25 plots)
fdis.phy.int <- with(int.div, lm(fdis~sesmntd+sesmpd+pd+nsp+simpsdiv))
vif(fdis.phy.int) # colinearity maybe a little high for pd and nsp

# LA ~ phy - int
la.phy.int.cor <- cor.test.phy(int.div, log(int.div$CWM.LA))

la.phy.int <- with(int.div, lm(log(CWM.LA)~sesmntd+sesmpd+pd+simpsdiv)) # mntd, mpd, pd
la.phy.int.phy <- with(int.div, lm(log(CWM.LA)~sesmntd+sesmpd+pd))
anova(la.phy.int, la.phy.int.phy) # no sig use simpler model
xtable(la.phy.int.phy, digits=3)

# SLA ~ phy - int
sla.phy.int.cor <- cor.test.phy(int.div, log(int.div$CWM.SLA))

sla.phy.int <- with(int.div, lm(log(CWM.SLA)~sesmntd+sesmpd+pd+simpsdiv)) # mntd
sla.phy.int.mntd <- with(int.div, lm(log(CWM.SLA)~sesmntd))
anova(sla.phy.int, sla.phy.int.mntd) # not sig, use simpler model -> PLOT

# mxH ~ phy - int
mxH.phy.int.cor <- cor.test.phy(int.div, int.div$CWM.maxht) 

mxH.phy.int <- with(int.div, lm(CWM.maxht~sesmntd+sesmpd+pd+simpsdiv)) # mntd, mpd
mxH.phy.int.phy <-with(int.div, lm(CWM.maxht~sesmntd+sesmpd))
anova(mxH.phy.int, mxH.phy.int.phy) # not sig, use simpler model -> TABLE
xtable(mxH.phy.int.phy, digits = 3)

# mnH ~ phy - int
mnH.phy.int.cor <- cor.test.phy(int.div, int.div$CWM.mn.ht)

mnH.phy.int <- with(int.div, lm(CWM.mn.ht~sesmntd+sesmpd+pd+simpsdiv)) # mpd, pd
mnH.phy.int.phy <- with(int.div, lm(CWM.mn.ht~sesmpd+pd)) # mpd, pd
anova(mnH.phy.int, mnH.phy.int.phy) # not sig, make simpler model - > table
xtable(mnH.phy.int.phy, digits=3)

# fdis ~ phy - int
fdis.phy.int.cor <- cor.test.phy(int.div, int.div$fdis) #pull out nsp

fdis.phy.int <- with(int.div, lm(fdis~sesmntd+sesmpd+pd+simpsdiv)) # pd and maybe mntd sig
fdis.phy.int.phy <- with(int.div, lm(fdis~sesmntd+pd)) #both pd and mntd sig
anova(fdis.phy.int, fdis.int.phy) # no sig. so use simpler model - > table
xtable(fdis.phy.int.phy, digits=3)

# FIG 4 - Functional diversity across phy-div metrics at 25 plots
# Plot significant relationships with SESmntd here, do the rest in the supplement
jpeg("./analysis/figures/supp-fdiv-phydiv-int.jpeg", width=4, height=5, unit="in",res=300)
# par(mar=c(5,6,0.1,0.5))
# par(oma=c(0.5,0.5,0.5,0.5))

# CWM of SLA across MNTD
with(int.div, plot(log(CWM.SLA)~sesmntd, pch=19, xlab=expression(SES[MNTD]), ylab=expression(italic("ln")("Specific leaf area")),  axes=FALSE, col="#009E73", ylim=c(1,6.2),cex.lab=0.95))
abline(sla.int.mntd, lwd=2)
axis(1, cex.axis=0.8)
axis(2, at=c(1,2,3,4,5,6), cex.axis=0.9)
legend(-0.5,2,expression(R^2 ~ "= 0.32, p-value < 0.001"), bty="n", cex=0.7)

dev.off()