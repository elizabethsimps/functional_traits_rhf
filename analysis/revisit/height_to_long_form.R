# Transferred to trait cleaning area #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #

# Making RHF height data long form
# Will Pearse # 2018-10-29 #
# Combining 2018-2019 height data
# Elizabeth Simpson # 2019-07-31
setwd("~/Documents/projects/")

# Load data
data.18 <- read.csv("./raw_data/rhf_2018_traits_height.csv", as.is=TRUE)

# Magic matrix --> numeric conversion trick
subset <- t(as.matrix(data.18[,-1:-7]))
heights <- as.numeric(subset)

# Expand the original data to match (dropping unneeded columns)
expanded.data.18 <- data.18[rep(1:nrow(data.18), each=25),]
expanded.data.18 <- expanded.data.18[,1:7]

# Match up the two
expanded.data.18$height <- heights
expanded.data.18 <- na.omit(expanded.data.18)
expanded.data.18 <- expanded.data.18[,-7]
rownames(expanded.data.18) <- seq_along(1:nrow(expanded.data.18))

# load in 2019 data, match columns to 2018 data
data.19 <- read.csv("./functional_traits_rhf/raw_data/rhf_2019_traits_height.csv", as.is=TRUE)
data.19 <- data.19[,1:7]
colnames(data.19)[7] <- "height"

# put all height together 
all.height <- rbind(expanded.data.18, data.19)

#some height metrics
sp.ht.mean <- as.data.frame(with(all.height, tapply(height, Species, mean)))
sp.ht.max <- as.data.frame(with(all.height, tapply(height, Species, max)))
sp.ht.sd <- as.data.frame(with(all.height, tapply(height, Species, sd)))
