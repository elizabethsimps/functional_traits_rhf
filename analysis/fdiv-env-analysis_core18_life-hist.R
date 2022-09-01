# Adding categoricla assessement of different life histories to FT~env analysis
# Elizabeth Simpson
# 2022-08-08

source("~/Documents/projects/functional_traits_rhf/analysis/fdiv-env-analysis_core18.R")

# Load life history data
cat.tr <- read.csv("./raw_data/cat_traits_rhf.csv", as.is=TRUE)

# Make a species lists for each life history strategy - woody perennial, herbaceous perennial, annual/biennial
# Woody perennials
p.wood <- cat.tr$species[cat.tr$life_history=="perennial-woody"]
# Herbaceous perennials
p.herb <- cat.tr$species[cat.tr$life_history=="perennial-herb"]
# Annual/Biennials
anubi <- cat.tr$species[cat.tr$life_history=="annual"]

# Sort out the 2018 cover data from the 25 core plots with temperature data into these life hist. cats.
wood <- core.18.over[core.18.over$Species %in% p.wood,]
herb <- core.18.over[core.18.over$Species %in% p.herb,]
anbi <- core.18.over[core.18.over$Species %in% anubi,]

# Make community matrices for each of these and remove communities that have zero abundance across all species types
# Remove communities that have zero abundance across all species 
wood.comm <- calc.comm(wood)
wood.comm <- wood.comm[,colnames(wood.comm) %in% rownames(traits)]
wood.comm <- wood.comm[rowSums(wood.comm)>0,]

herb.comm <- calc.comm(herb)
herb.comm <- herb.comm[,colnames(herb.comm) %in% rownames(traits)] 
herb.comm <- herb.comm[rowSums(herb.comm)>0,]

anbi.comm <- calc.comm(anbi)
anbi.comm <- anbi.comm[,colnames(anbi.comm) %in% rownames(traits)] 
anbi.comm <- anbi.comm[rowSums(anbi.comm)>0,]

###################################
# Make comparative community object & calculate diversity metrics from that
# NOTE: *only* look at CWM of traits and FDis from this
wood.div <- div.calc(tree18.core, wood.comm, traits, env18.core) # 17 plots
herb.div <- div.calc(tree18.core, herb.comm, traits, env18.core) # 25 plots
anbi.div <- div.calc(tree18.core, anbi.comm, traits, env18.core) # 21 plots

### SPLIT OUT FIGURE # models from main text to look at woody perennials, herbaceous perennials, and annuals/biennials

### SLA
sla.wood <- with(wood.div, lm(log(CWM.SLA)~mean+SAND)) # intercept
sla.herb <- with(herb.div, lm(log(CWM.SLA)~mean+SAND)) # mean
sla.anbi <- with(anbi.div, lm(log(CWM.SLA)~mean+SAND)) # intercept

# SLA plotting
sla.wood.1 <- with(wood.div, lm(log(CWM.SLA)~1))
sla.herb.mn <- with(herb.div, lm(log(CWM.SLA)~mean))
sla.anbi.1 <- with(anbi.div, lm(log(CWM.SLA)~1))

### LA
la.wood <- with(wood.div, lm(log(CWM.SLA)~mean+CLAY))
la.herb <- with(herb.div, lm(log(CWM.SLA)~mean+CLAY))
la.anbi <- with(anbi.div, lm(log(CWM.SLA)~mean+CLAY))

# LA plotting
la.wood.1 <- with(wood.div, lm(log(CWM.LA)~1))
la.herb.mn <- with(herb.div, lm(log(CWM.LA)~mean))
la.anbi.1 <- with(anbi.div, lm(log(CWM.LA)~1))

### CWM of mean HT
mnH.wood <- with(wood.div, lm(CWM.mn.ht~sd+CLAY))
mnH.herb <- with(herb.div, lm(CWM.mn.ht~sd+CLAY))
mnH.anbi <- with(anbi.div, lm(CWM.mn.ht~sd+CLAY))

# CWM of mean HT plotting
mnH.wood.c <- with(wood.div, lm(CWM.mn.ht~CLAY))
mnH.herb.1 <- with(herb.div, lm(CWM.mn.ht~1))
mnH.anbi.1 <- with(anbi.div, lm(CWM.mn.ht~1))

### CWM of max HT
mxH.wood <- with(wood.div, lm(CWM.maxht~sd+CLAY))
mxH.herb <- with(herb.div, lm(CWM.maxht~sd+CLAY))
mxH.anbi <- with(anbi.div, lm(CWM.maxht~sd+CLAY))

# CWM of max HT plotting
mxH.wood.1 <- with(wood.div, lm(CWM.maxht~1))
mxH.herb.c <- with(herb.div, lm(CWM.maxht~CLAY))
mxH.anbi.1 <- with(anbi.div, lm(CWM.maxht~1))

### FDis - with four traits
fdis.wood <- with(wood.div, lm(fdis~sd+CLAY))
fdis.herb <- with(herb.div, lm(fdis~sd+CLAY))
fdis.anbi <- with(anbi.div, lm(fdis~sd+CLAY))

# for plotting
fdis.wood.sd <- with(wood.div, lm(fdis~sd))
fdis.herb.sd <- with(herb.div, lm(fdis~sd))
fdis.anbi.1 <- with(anbi.div, lm(fdis~1))

# supplementary plot showing the effect of different life history strategies on FTs
jpeg("./analysis/figures/supp-fdiv-env-core18-life-history.jpeg", width=7, height=10, unit="in",res=300)
par(mfrow=c(5,3))
par(mar=c(4,5,0.1,0.5))
par(oma=c(1.5,1.5,1.5,1.5))

# LA ~ mean
with(wood.div, plot(log(CWM.LA)~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab=expression(italic("ln")("Leaf area")), cex=1.2, axes=FALSE, col="#009E73", cex.lab=1.05))
abline(la.wood.1, lwd=2)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(a)")),side=3,line=-2, adj=2.2, cex = 0.85)

with(herb.div, plot(log(CWM.LA)~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab="", cex=1.2, axes=FALSE, col="#009E73", cex.lab=1.05))
abline(la.herb.mn, lwd=2, lty=2)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(b)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.62, p-value < 0.001"),side=3,line=-2, adj=2.2, cex = 0.85)
      
with(anbi.div, plot(log(CWM.LA)~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab="", cex=1.2, axes=FALSE, col="#009E73", cex.lab=1.05))
abline(la.anbi.1, lwd=2, lty=3)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(c)")),side=3,line=-2, adj=2.2, cex = 0.85)  

# SLA
with(wood.div, plot(log(CWM.SLA)~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab=expression(italic("ln")("Specific leaf area")), cex=1.2, axes=FALSE, col="#009E73", cex.lab=1.05))
abline(sla.wood.1, lwd=2)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(d)")),side=3,line=-2, adj=2.2, cex = 0.85)

with(herb.div, plot(log(CWM.SLA)~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab="", cex=1.2, axes=FALSE, col="#009E73", cex.lab=1.05))
abline(sla.herb.mn, lwd=2, lty=2)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(e)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.58, p-value < 0.001"),side=3,line=-2, adj=2.2, cex = 0.85)

with(anbi.div, plot(log(CWM.LA)~mean, pch=19, xlab=expression(paste('Mean temperature (',degree,'C)')), ylab="", cex=1.2, axes=FALSE, col="#009E73", cex.lab=1.05))
abline(sla.anbi.1, lwd=2, lty=3)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(f)")),side=3,line=-2, adj=2.2, cex = 0.85)

# max ht ~ sd
with(wood.div, plot(log(CWM.maxht)~CLAY, pch=19, xlab="Clay (%)", ylab="Max. height (cm)", cex=1.2, axes=FALSE, col="#E69F00", cex.lab=1.05))
abline(mxH.wood.1, lwd=2)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(g)")),side=3,line=-2, adj=2.2, cex = 0.85)

with(herb.div, plot(log(CWM.maxht)~CLAY, pch=19, xlab="Clay (%)", ylab="", cex=1.2, axes=FALSE, col="#E69F00", cex.lab=1.05))
abline(mx.herb.c, lwd=2, lty=2)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(h)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.24, p-value = 0.013"),side=3,line=-2, adj=2.2, cex = 0.85)

with(anbi.div, plot(log(CWM.maxht)~CLAY, pch=19, xlab="Clay (%)", ylab="", cex=1.2, axes=FALSE, col="#E69F00", cex.lab=1.05))
abline(mx.anbi.1, lwd=2, lty=3)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(i)")),side=3,line=-2, adj=2.2, cex = 0.85) 

# Mean height
with(wood.div, plot(log(CWM.mn.ht)~CLAY, pch=19, xlab="Clay (%)", ylab="Mean height (cm)", cex=1.2, axes=FALSE, col="#E69F00", cex.lab=1.05))
abline(mnH.wood.c, lwd=2)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(j)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.23, p-value = 0.049"),side=3,line=-2, adj=2.2, cex = 0.85)

with(herb.div, plot(log(CWM.mn.ht)~CLAY, pch=19, xlab="Clay (%)", ylab="", cex=1.2, axes=FALSE, col="#E69F00", cex.lab=1.05))
abline(mx.herb.c, lwd=2, lty=2)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(k)")),side=3,line=-2, adj=2.2, cex = 0.85)

with(anbi.div, plot(log(CWM.mn.ht)~CLAY, pch=19, xlab="Clay (%)", ylab="", cex=1.2, axes=FALSE, col="#E69F00", cex.lab=1.05))
abline(mx.anbi.1, lwd=2, lty=3)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(l)")),side=3,line=-2, adj=2.2, cex = 0.85) 

# functional dispersion ~ sd temperature
with(wood.div, plot(fdis~sd, pch=19, xlab=expression(paste('SD temperature (',degree,'C)')), ylab="Functional dispersion", cex=1.2, axes=FALSE, col="#009E73", cex.lab=1.05))
abline(fdis.wood.sd, lwd=2)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(m)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.23, p-value = 0.049"),side=3,line=-2, adj=2.2, cex = 0.85)

with(herb.div, plot(fdis~sd, pch=19, xlab=expression(paste('SD temperature (',degree,'C)')), ylab="", cex=1.2, axes=FALSE, col="#009E73", cex.lab=1.05))
abline(fdis.herb.sd, lwd=2, lty=2)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(n)") ~~~~~~~~~~~~~~~~~~~~~~~~ R^2 ~ "= 0.26, p-value = 0.010"),side=3,line=-2, adj=2.2, cex = 0.85)

with(anbi.div, plot(fdis~sd, pch=19, xlab=expression(paste('SD temperature (',degree,'C)')), ylab="", cex=1.2, axes=FALSE, col="#009E73", cex.lab=1.05))
abline(fdis.anbi.1, lwd=2, lty=2)
axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5)
mtext(expression(bold("(o)")),side=3,line=-2, adj=2.2, cex = 0.85)

dev.off()
