# Seeds to collect - 2020
# Elizabeth simpson
setwd("~/Documents/projects/")

seeds <- read.csv("./functional_traits_rhf/raw_data/rhf_2019_traits_seeds.csv", as.is=TRUE)
species <- read.csv("./fractal_sampling_div_rhf/data/rhf_2018_cover.csv", as.is=TRUE)
seed_have <- unique(seeds$Species)
species_total <- unique(species$Species)

write.csv(setdiff(seed_have, species_total), "seeds_collected_not_in_species_NEW.csv")
write.csv(setdiff(species_total, seed_have), "species_no_seeds_for_NEW.csv")

sp_locs <- with(species, tapply(Plot_id, list(Species), function(x) paste(unique(x))))
library(plyr)
sp_locs <- ldply(sp_locs, rbind)


sp_locs_c <- sapply(sp_locs, as.character) # since your values are `factor`
sp_locs_c[is.na(sp_locs)] <- ""

write.csv(sp_locs_c, "sp_locs.csv")


