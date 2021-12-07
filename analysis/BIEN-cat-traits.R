setwd("~/Documents/projects/functional_traits_rhf")

library(BIEN)

species <- read.csv("./clean_data/cat_traits_names.csv", as.is=TRUE)
species <- as.vector(species[,2])
species <- species[-127]

traits <- "whole plant dispersal syndrome"

output <- BIEN_trait_traitbyspecies(trait="whole plant dispersal syndrome", species=species)
