library(pez)
data <- read.csv("~/Documents/Spatial_Scaling_Field_Project/raw_data/sites/right_hand_fork/cover.csv")
comm <- with(data, tapply(Cover, list(Site,Species), function(x) mean(x,na.rm=TRUE)))
comm[is.na(comm)] <- 0

# remove species esentially there as notes
comm <- comm[,!grepl("^[a-z]+", colnames(comm))]
comm <- comm[,!grepl("\\(|/", colnames(comm))]
comm <- comm[,!grepl("sp\\.", colnames(comm))]

colnames(comm) <- tolower(gsub(" ", "_", colnames(comm)))

#Build a phylogeny and create a combined data object
tree <- read.tree("~/Documents/Spatial_Scaling_Field_Project/raw_data/trait-phylo/Vascular_Plants_rooted.dated.tre")
tree$tip.label <- tolower(gsub("_", " ", tree$tip.label))
tree <- congeneric.merge(tree, colnames(comm))

terrain <- read.csv("~/Documents/Spatial_Scaling_Field_Project/clean_data/RightFork-Site1_terrain_data_radians.csv")
env <- terrain

#name rows as sites
rownames(env) <- env$plot.name
# drop out all other variables besides elev and cos.asp
env <- env[,c("elev","cos.asp"), drop=FALSE]

traits <- read.csv("~/Documents/Spatial_Scaling_Field_Project/raw_data/trait-phylo/traits_elizabeth.csv")
#row names need to match the site names for community data
rownames(traits) <-traits$species
traits$X <-NULL
traits$species <- NULL

##############################
# Combine terrain and comm ###
##############################

c.data <- comparative.comm(tree, comm, traits=traits, env=env) 

sqrt <- pez.shape(c.data, sqrt.phy = TRUE)
traits <- pez.shape(c.data, traitgram = 1)
traits <- pez.shape(c.data, traitgram = c(0,0.5))
str(c.data)
