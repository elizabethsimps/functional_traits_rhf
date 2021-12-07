# Cleaning up measured terrain variables #
# Elizabeth Simpson - 2017-10-26 ##########################
# most recently neatened 2021-1-28 ####################
setwd("~/Documents/projects/functional_traits_rhf")

# cosine transform aspect (coverts from degrees to radians within function) S = -1, N = 1
cos.transform <- function(asp_deg){
  asp_cos <- NA
  for (i in 1:length(asp_deg)){
    asp_cos[i] <- cos((asp_deg[i])*pi/180)
  }
  return(asp_cos)
}

calc.slope <- function(slope_up, slope_down){
  slope_deg <- NA
  for (i in 1: length(slope_up)){
    slope_deg[i] <- mean(c(slope_up[i], slope_down[i]))
  }
  return(slope_deg)
}

#################
# Measured data #
#################

# load and clean measured terrain
m.terrain <- read.csv("./raw_data/rhf_2019_terrain_measured.csv", as.is=TRUE)
m.terrain <- m.terrain[, c(2:6)] #relevant data

# transform aspect
m.terrain$m_cos_asp <- with(m.terrain, cos.transform(Aspect_deg))

# calculate slope
m.terrain$slope_deg <- with(m.terrain, calc.slope(Slope_down_deg, Slope_up_deg))
m.terrain <- m.terrain[,c(1,5:7)]
colnames(m.terrain) <- c("Plot_id", "elev_m", "aspect_rad", "slope_deg")

write.csv(m.terrain, "./clean_data/rhf_2019_terrain_measured_cleaned.csv", row.names=F)
