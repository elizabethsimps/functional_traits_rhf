# Cosine transform measured aspect
# Elizabeth Simpson # 2019-07-30
setwd("~/Documents/projects/")

# cosine transform aspect (coverts from degrees to radians within function)
# S = -1, N = 1
cos.transform <- function(asp_deg){
  asp_cos <- NA
  for (i in 1:length(asp_deg)){
    asp_cos[i] <- cos((asp_deg[i])*pi/180)
  }
  return(asp_cos)
}

#terrain to look at temperature across
m.terrain <- read.csv("./fractal_sampling_div_rhf/raw_data/env/rhf_2018_terrain_measured.csv")
m.terrain <- m.terrain[, 2:3]
m.terrain <- m.terrain[-79,]
m.terrain$m_cos_asp <- cos.transform(m.terrain$Aspect_deg)
m.terrain <- m.terrain[,-2]