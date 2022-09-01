# Adding categoricla assessement of different life histories to FT~env analysis
# Elizabeth Simpson
# 2022-08-08

source("~/Documents/projects/functional_traits_rhf/analysis/fdiv_assemblage-env_analysis.R")

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
anbi <- core.18.over[core.18.over$Species %in% annual,]

# Make community matricies for each of these and remove communities that have zero abundance across all species types
# Remove communties that have zero abundance across all species 
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
wood.div <- div.calc(tree18.core, wood.comm, traits, env18.core)
herb.div <- div.calc(tree18.core, herb.comm, traits, env18.core)
anbi.div <- div.calc(tree18.core, anbi.comm, traits, env18.core)

######## FOR WOODY PERRENIALS ############
wood.sla <- with(wood.div, lm(log(CWM.SLA)~mean+sd+max+min+SAND+CLAY+SILT))
vif(wood.sla)
# There are aliased coefficients in the model

### CWM of SLA - WPs
wood.sla.cor <- env.corr(wood.div, log(wood.div$CWM.SLA)) 
xtable(wood.sla.cor, digits=3) # for supplement
wood.sla <- with(wood.div, lm(log(CWM.SLA)~min+SILT))
xtable(wood.sla, digits=3) # for supplement
# plotting
wood.sla.1 <- with(wood.div, lm(log(CWM.SLA)~1))

### CWM of LA - WPs
wood.la.cor <- env.corr(wood.div, log(wood.div$CWM.LA))
xtable(wood.la.cor, digits=3)
wood.la <- with(wood.div, lm(log(CWM.LA)~mean+SILT)) 
xtable(wood.la, digits=3)
# plotting
wood.la.mn <- with(wood.div, lm(log(CWM.LA)~mean)) # mean sig. but not the intercept

### CWM of mean HT - WPs
wood.mnH.cor <- env.corr(wood.div, wood.div$CWM.mn.ht)
xtable(wood.mnH.cor, digits=3)
wood.mnH <- with(wood.div, lm(CWM.mn.ht~max*CLAY)) # clay and max interaction significant
xtable(wood.mnH, digits=3)
# plotting
wood.mnH.c <- with(wood.div, lm(CWM.mn.ht~CLAY)) # slope but not intercept significant

### CWM of max HT - WPs
wood.mxH.cor <- env.corr(wood.div, wood.div$CWM.maxht)
xtable(wood.mxH.cor, digits=3)
wood.mxH <- with(wood.div, lm(CWM.maxht~max*CLAY)) # clay and max interaction significant
xtable(wood.mxH, digits=3)
# plotting

### FDis (four traits) - WPs
wood.fdis <- env.corr(wood.div, wood.div$fdis)
xtable(wood.fdis, digits=3)
wood.fdis <- with(wood.div, lm(fdis~mean+CLAY)) # mean and clay significant, but not interaction
xtable(wood.fdis, digits=3)
# plotting
wood.fdis.mn <- with(wood.div, lm(fdis~mean)) # sig

######## FOR HERBACEOUS PERRENIALS ############
### CWM of SLA - HPs
herb.sla.cor <- env.corr(herb.div, log(herb.div$CWM.SLA)) 
xtable(herb.sla.cor, digits=3) # for supplement
herb.sla <- with(herb.div, lm(log(CWM.SLA)~mean+CLAY)) # mean sig.
xtable(herb.sla, digits=3) # for supplement
# plotting
herb.sla.mn <- with(herb.div, lm(log(CWM.SLA)~mean))

### CWM of LA - HPs
herb.la.cor <- env.corr(herb.div, log(herb.div$CWM.LA))
xtable(herb.la.cor, digits=3)
herb.la <- with(herb.div, lm(log(CWM.LA)~mean+SAND)) # mean sig
xtable(herb.la, digits=3)
# plotting
herb.la.mn <- with(herb.div, lm(log(CWM.SLA)~mean))

### CWM of mean HT - HPs
herb.mnH.cor <- env.corr(herb.div, herb.div$CWM.mn.ht)
xtable(herb.mnH.cor, digits=3)
herb.mnH <- with(herb.div, lm(CWM.mn.ht~sd+CLAY)) # no sig.
xtable(herb.mnH, digits=3)
# plotting
herb.mnH <- with(herb.div, lm(CWM.mn.ht~1))

### CWM of max HT - HPs
herb.mxH.cor <- env.corr(herb.div, herb.div$CWM.maxht)
xtable(herb.mxH.cor, digits=3)
herb.mxH <- with(herb.div, lm(CWM.maxht~mean+CLAY)) # clay sig.
xtable(herb.mxH, digits=3)
# plotting
herb.mxH.c <- with(herb.div, lm(CWM.maxht~CLAY))

### FDis (four traits) - HPs
herb.fdis <- env.corr(herb.div, herb.div$fdis)
xtable(herb.fdis, digits=3)
herb.fdis <- with(herb.div, lm(fdis~mean+CLAY)) 
xtable(herb.fdis, digits=3)
# plotting
herb.fdis.mn <- with(herb.div, lm(fdis~mean))

######## FOR ANNUALS/BIENNIALS ############
### CWM of SLA - A/Bs
anbi.sla.cor <- env.corr(anbi.div, log(anbi.div$CWM.SLA)) 
xtable(anbi.sla.cor, digits=3) # for supplement
anbi.sla <- with(anbi.div, lm(log(CWM.SLA)~max+CLAY)) # max. sig.
xtable(anbi.sla, digits=3) # for supplement
# plotting
anbi.sla.mx <- with(anbi.div, lm(log(CWM.SLA)~max))

### CWM of LA - A/Bs
anbi.la.cor <- env.corr(anbi.div, log(anbi.div$CWM.LA))
xtable(anbi.la.cor, digits=3)
anbi.la <- with(anbi.div, lm(log(CWM.LA)~max+SILT)) # max. sig.
xtable(anbi.la, digits=3)
# plotting
anbi.sla.mx <- with(anbi.div, lm(log(CWM.LA)~max))

### CWM of mean HT - A/Bs
anbi.mnH.cor <- env.corr(anbi.div, anbi.div$CWM.mn.ht)
xtable(anbi.mnH.cor, digits=3)
anbi.mnH <- with(anbi.div, lm(CWM.mn.ht~min+CLAY)) # no sig.
xtable(anbi.mnH, digits=3)
# plotting
anbi.mnH.1 <- with(anbi.div, lm(CWM.mn.ht~1))

### CWM of max HT - A/Bs
anbi.mxH.cor <- env.corr(anbi.div, anbi.div$CWM.maxht)
xtable(anbi.mxH.cor, digits=3)
anbi.mxH <- with(anbi.div, lm(CWM.maxht~min+CLAY)) # no sig.
xtable(anbi.mxH, digits=3)
# plotting
anbi.mnH.1 <- with(anbi.div, lm(CWM.maxht~1))

### FDis (four traits) - A/Bs
anbi.fdis <- env.corr(anbi.div, anbi.div$fdis)
xtable(anbi.fdis, digits=3)
anbi.fdis <- with(anbi.div, lm(fdis~max+SILT)) # no sig.
xtable(anbi.fdis, digits=3)
# plotting
anbi.fdis.1 <- with(anbi.div, lm(fdis~1))