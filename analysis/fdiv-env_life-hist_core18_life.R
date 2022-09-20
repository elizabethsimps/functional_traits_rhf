# Assessing how categorization of species by life histories affects FT~env analysis
# Elizabeth Simpson
# 2022-08-08

source("~/Documents/projects/functional_traits_rhf/analysis/fdiv-env-phy_core18.R")

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

### LA
la.wood <- with(wood.div, lm(log(CWM.LA)~mean+CLAY)) # mean
la.herb <- with(herb.div, lm(log(CWM.LA)~mean+CLAY)) # mean AND clay -> interpret
la.anbi <- with(anbi.div, lm(log(CWM.LA)~mean+CLAY)) # intercept

# LA plotting
la.wood.mn <- with(wood.div, lm(log(CWM.LA)~mean)) # little non-normal
anova(la.wood, la.wood.mn) # not sig, plot la.wood.mn
la.herb.mn <- with(herb.div, lm(log(CWM.LA)~mean))
anova(la.herb, la.herb.mn) # sig, plot la.herb.mn BUT INCLUDE table for la.herb
la.anbi.1 <- with(anbi.div, lm(log(CWM.LA)~1))

### SLA
sla.wood <- with(wood.div, lm(log(CWM.SLA)~mean+CLAY)) # intercept
sla.herb <- with(herb.div, lm(log(CWM.SLA)~mean+CLAY)) # mean
sla.anbi <- with(anbi.div, lm(log(CWM.SLA)~mean+CLAY)) # intercept

# SLA plotting
sla.wood.1 <- with(wood.div, lm(log(CWM.SLA)~1))
sla.herb.mn <- with(herb.div, lm(log(CWM.SLA)~mean))
anova(sla.herb, sla.herb.mn) # not sig, plot sla.herb.mn
sla.anbi.1 <- with(anbi.div, lm(log(CWM.SLA)~1))

### CWM of max HT
mxH.wood <- with(wood.div, lm(CWM.maxht~max+CLAY)) # max + clay -> table -> plot points in zoomed in area
mxH.herb <- with(herb.div, lm(CWM.maxht~max+CLAY)) # intercept, clay -> plot
mxH.anbi <- with(anbi.div, lm(CWM.maxht~max+CLAY)) # intercept, -> plot

# CWM of max HT plotting
mxH.wood.1 <- with(wood.div, lm(CWM.maxht~1)) 
mxH.herb.c <- with(herb.div, lm(CWM.maxht~CLAY))
mxH.anbi.1 <- with(anbi.div, lm(CWM.maxht~1))

### CWM of mean HT
mnH.wood <- with(wood.div, lm(CWM.mn.ht~max+CLAY)) # max + clay -> table and fig b current fig - in supp
mnH.herb <- with(herb.div, lm(CWM.mn.ht~max+CLAY)) # intercept
mnH.anbi <- with(anbi.div, lm(CWM.mn.ht~max+CLAY)) # intercept

# CWM of mean HT plotting
mnH.wood.c <- with(wood.div, lm(CWM.mn.ht~CLAY))
mnH.herb.1 <- with(herb.div, lm(CWM.mn.ht~1))
mnH.anbi.1 <- with(anbi.div, lm(CWM.mn.ht~1))

### Functional dispersion
fdis.wood <- with(wood.div, lm(fdis~max+CLAY)) # intercept, max, clay -> table, but only max in univarite, so plot that too
xtable(fdis.wood, digits=3)
fdis.herb <- with(herb.div, lm(fdis~max+CLAY)) # intercept, max -> plot
fdis.anbi <- with(anbi.div, lm(fdis~max+CLAY)) # no sig.

# for plotting
fdis.wood.mx <- with(wood.div, lm(fdis~max))
anova(fdis.wood, fdis.wood.mx) # sig
fdis.herb.mx <- with(herb.div, lm(fdis~max))
anova(fdis.herb, fdis.herb.mx) # not sig, no table
fdis.anbi.1 <- with(anbi.div, lm(fdis~1))

########################
# transparent colors for this plot
t.brown <- rgb(84, 48, 5, max = 255, alpha = 170, names = "t.brown")
t.blue <- rgb(53, 151, 143, max = 255, alpha = 170, names = "t.blue")
t.green <- rgb(127, 188, 65, max = 255, alpha = 170, names = "t.green")

#### PLOTTING LEAF TRAITS
jpeg("./analysis/figures/fdiv-env-core18-wLH.jpeg", width=7, height=7, unit="in",res=300)
par(mfrow=c(2,2))
par(mar=c(4,5,0.1,0.5))
par(oma=c(1.5,1.2,1.5,1.2))

# LA ~ mean
with(core18div, plot(log(CWM.LA)~mean, pch=19, xlab="", ylab=expression(italic("ln")("Leaf area")), cex=0.9, axes=FALSE, xlim=c(4,14),ylim=c(-3,4.1), col="#009E73" )) # col was 
abline(la.mn, lwd=2)
axis(1, cex.lab=1.5)
axis(2, at=c(-3,-1,1,3), cex.lab=1.5)
legend(6.6,2.1, expression(R^2 ~ "= 0.66, p-value < 0.001"), box.lty=0, bg="transparent", cex = 0.8)
mtext(expression(bold("(a)")),side=3, line=-2, adj=-0.35, cex = 0.85)

# SLA ~ mean
with(core18div, plot(log(CWM.SLA)~mean, pch=19, xlab="", ylab=expression(italic("ln")("Specific leaf area")), cex=0.9, axes=FALSE, xlim=c(4,14),ylim=c(2,5.6), col="#009E73" ))
abline(sla.mn, lwd=2)
axis(1, cex.lab=1.5)
axis(2, at=c(2,3,4,5), cex.lab=1.5)
legend(6.7,4.1, expression(R^2 ~ "= 0.31, p-value = 0.004"), box.lty=0, bg="transparent", cex = 0.8)
mtext(expression(bold("(b)")),side=3, line=-2, adj=-0.35, cex = 0.85)

# LA ~ mean - life history
with(wood.div, plot(log(CWM.LA)~mean, pch=19, cex=0.9, col="#543005AA", xlab=expression(paste('Mean temp. (',degree,'C)')), ylab=expression(italic("ln")("Leaf area")), axes=FALSE,  xlim=c(4,14), ylim=c(-4,5.4)))
with(herb.div, points(log(CWM.LA)~mean, pch=19, cex=0.9, col="#35978FAA"))
with(anbi.div, points(log(CWM.LA)~mean, pch=19,cex=0.9, col="#7FBC41AA"))

abline(la.wood.mn, lwd=2, col="#543005")
abline(la.herb.mn, lwd=2, lty=2, col="#35978F")
abline(la.anbi.1, lwd=2, lty=3, col="#7FBC41")

axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5, at=c(-4,-2,0,2,4))
legend(5.9,-2.1, expression(R^2 ~ "= 0.30, p-value = 0.022"), box.lty=0, bg="transparent",text.col="#543005", cex = 0.8)
legend(7,2.35, expression(R^2 ~ "= 0.58, p-value < 0.001"), box.lty=0, bg="transparent",text.col="#35978F", cex = 0.8)
legend(9,4.5, legend=c('Woody perenn.', 'Herb. perenn.', 'Annual/biennial'), pch=16, cex=0.7, col = c("#543005", "#35978F", "#7FBC41"), lty=c(1,2,3))
mtext(expression(bold("(c)")),side=3, line=-2, adj=-0.35, cex = 0.85)

# SLA ~ mean - life history
with(wood.div, plot(log(CWM.SLA)~mean, pch=19, cex=0.9, col="#543005AA", xlab=expression(paste('Mean temp. (',degree,'C)')), ylab=expression(italic("ln")("Specific leaf area")), axes=FALSE, xlim=c(4,14), ylim=c(1,5.7)))
with(herb.div, points(log(CWM.SLA)~mean, pch=19, cex=0.9, col="#35978FAA"))
with(anbi.div, points(log(CWM.SLA)~mean, pch=19, cex=0.9,col="#7FBC41AA"))

abline(sla.wood.1, lwd=2, col="#543005")
abline(sla.herb.mn, lwd=2, lty=2, col="#35978F")
abline(sla.anbi.1, lwd=2, lty=3, col="#7FBC41")

axis(1, cex.lab=1.5)
axis(2, cex.lab=1.5, at=c(1,2,3,4,5))
legend(6.4,4.5, expression(R^2 ~ "= 0.58, p-value < 0.001"), box.lty=0, bg="transparent",text.col="#35978F", cex = 0.8)
mtext(expression(bold("(d)")),side=3, line=-2, adj=-0.35, cex = 0.85)

dev.off()


# HEIGHT AND FDIS life history
jpeg("./analysis/figures/fdiv-env-core18-LH-fdis-ht.jpeg", width=7, height=3.5, unit="in",res=300)
par(mfrow=c(1,2))
par(mar=c(4,5,1,0.5))
par(oma=c(1,1,1,1))

# functional dispersion ~ max temperature
with(wood.div, plot(fdis~max, pch=16,cex=0.7, col="#543005AA", xlab=expression(paste('Max. temp. (',degree,'C)')), ylab="Functional dispersion", axes=FALSE, cex.lab=0.8, xlim=c(20,70)))
with(herb.div, points(fdis~max, pch=16, cex=0.7, col="#35978FAA"))
with(anbi.div, points(fdis~max, pch=16, cex=0.7, col="#7FBC41AA" ))

abline(fdis.wood.mx, lwd=2, col="#543005")
abline(fdis.herb.mx, lwd=2, lty=2, col="#35978F")
abline(fdis.anbi.1, lwd=2, lty=3, col="#7FBC41")

axis(1, cex.axis=0.8)
axis(2, cex.axis=0.8)
legend(30,2.5, expression(R^2 ~ "= 0.28, p-value = 0.030"), box.lty=0, bg="transparent",text.col="#543005", cex = 0.65)
legend(33,2.25, expression(R^2 ~ "= 0.27, p-value = 0.008"), box.lty=0, bg="transparent",text.col="#35978F", cex = 0.65)
mtext(expression(bold("(a)")),side=3, line=0, adj=-0.4, cex = 0.85)

# max
with(wood.div, plot(CWM.maxht~CLAY, pch=16, cex=0.7, col="#543005AA", xlab="Clay (%)", ylab="Max. height (cm)", axes=FALSE, cex.lab=0.8, xlim=c(0,25), ylim=c(20,122)))
with(herb.div, points(CWM.maxht~CLAY, pch=16, cex=0.7, col="#35978FAA"))
with(anbi.div, points(CWM.maxht~CLAY, pch=16, cex=0.7, col="#7FBC41AA"))
abline(mxH.herb.c, lwd=2, lty=2, col="#35978F")
abline(mxH.anbi.1, lwd=2, lty=3, col="#7FBC41")

axis(1, cex.axis=0.8)
axis(2, cex.axis=0.8)
legend(0,38, expression(R^2 ~ "= 0.24, p-value = 0.013"), box.lty=0, bg="transparent",text.col="#35978F", cex = 0.65)
mtext(expression(bold("(b)")),side=3, line=0, adj=-0.4, cex = 0.85)

legend(13,110, legend=c('Woody perenn.', 'Herb. perenn.', 'Annual/biennial'), pch=16, col = c("#543005", "#35978F", "#7FBC41"), lty=c(1,2,3), cex=0.5)

dev.off()

# Mean height across clay
jpeg("./analysis/figures/supp-fdiv-env-core18-LH-fdis-mean-ht.jpeg", width=3.5, height=3.5, unit="in",res=300)
par(mar=c(4,4,1,1))

# max
with(wood.div, plot(CWM.mn.ht~CLAY, pch=16, cex=0.7, col="#543005AA", xlab="Clay (%)", ylab="Mean height (cm)", axes=FALSE, cex.lab=0.8, xlim=c(0,25), ylim=c(0,1200)))
with(herb.div, points(CWM.mn.ht~CLAY, pch=16, cex=0.7, col="#35978FAA"))
with(anbi.div, points(CWM.mn.ht~CLAY, pch=16, cex=0.7, col="#7FBC41AA"))
abline(mnH.wood.c, lwd=2, col="#543005")
abline(mnH.herb.1, lwd=2, lty=2, col="#35978F")
abline(mnH.anbi.1, lwd=2, lty=3, col="#7FBC41")

axis(1, cex.axis=0.8)
axis(2, cex.axis=0.65)
legend(5,800, expression(R^2 ~ "= 0.23, p-value = 0.050"), box.lty=0, bg="transparent",text.col="#543005", cex = 0.65)
legend(0,1200, legend=c('Woody perenn.', 'Herb. perenn.', 'Annual/biennial'), pch=16, col = c("#543005", "#35978F", "#7FBC41"), lty=c(1,2,3), cex=0.5)

dev.off()