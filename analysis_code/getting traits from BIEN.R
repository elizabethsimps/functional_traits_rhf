library(BIEN)
traits <- BIEN_trait_species("Galium aparine")
traits <- traits[!is.na(traits$latitude),]
traits <- traits[!is.na(traits$elevation_m),]

P.traits <- subset(traits, trait_name == "leaf phosphorus content per leaf dry mass")
C.traits <- subset(traits, trait_name == "leaf carbon content per leaf dry mass")
N.traits <- subset(traits, trait_name == "leaf nitrogen content per leaf dry mass")
SLA.traits <- subset(traits, trait_name == "leaf area per leaf dry mass")

with(P.traits, plot(trait_value~elevation_m, col = as.factor(unique(project_pi))))
with(C.traits, plot(trait_value~elevation_m, col = as.factor(unique(project_pi))))
with(N.traits, plot(trait_value~elevation_m, col = as.factor(unique(project_pi))))
with(SLA.traits, plot(trait_value~elevation_m, col = as.factor(unique(project_pi))))

