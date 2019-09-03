# Analyzing variance and scale of temperature and FT's at RHF
# Elizabeth Simpson # 2019-07-18 
setwd("~/Documents/projects/")

### data loading and setup
# load temp in deg C and subset to 2017-09-29 @ 0:00 to 2018-09-12 @ 0:00
temp <- read.csv("./functional_traits_rhf/clean_data/subset_formatted_temp_data.csv", as.is=TRUE)

# load cosine-transformed aspect & subset to match temp data
m.asp <- read.csv("./functional_traits_rhf/clean_data/m_aspect_cos.csv", as.is=TRUE)
m.asp <- subset(m.asp, Plot_id %in% temp$Plot_id)

### calculate temp metrics
mean.temp <- as.data.frame(with(temp, tapply(Temperature, Plot_id, function(x) mean(x,na.rm=TRUE))))
min.temp <-  as.data.frame(with(temp, tapply(Temperature, Plot_id, function(x) min(x,na.rm=TRUE))))
max.temp <- as.data.frame(with(temp, tapply(Temperature, Plot_id, function(x) max(x,na.rm=TRUE))))
var.temp <- as.data.frame(with(temp, tapply(Temperature, Plot_id, function(x) var(x,na.rm=TRUE))))

### Looking at how temperature varies across terrain
terr.temp <- cbind(m.asp, min.temp, mean.temp, max.temp, var.temp)
terr.temp <- terr.temp[,-1]
colnames(terr.temp) <- c("Plot_id", "m_cos_asp", "min_temp", "mean_temp", "max_temp", "var_temp")

#temperature min mean max
with(terr.temp, plot(mean_temp~m_cos_asp, pch=19, cex=2, ylim=c(-10,70), xlab="", ylab="", cex.axis=2))
with(terr.temp, points(min_temp~m_cos_asp, pch=19, cex=2, col="deepskyblue3"))
with(terr.temp, points(max_temp~m_cos_asp, pch=19, cex=2, col="chocolate2"))
mean.temp <- with(terr.temp, lm(mean_temp~m_cos_asp))
min.temp <- with (terr.temp, lm(min_temp~m_cos_asp))
max.temp <- with(terr.temp, lm(max_temp~m_cos_asp))
abline(mean.temp, lwd=3)
abline(min.temp, col="deepskyblue3", lwd=3)
abline(max.temp, col="chocolate2", lwd=3)

#variance in temperature across aspect
with(terr.temp, plot(var_temp~m_cos_asp, pch=19, cex=2, xlab="", ylab="", cex.axis=2))
var.temp <- with(terr.temp, lm(var_temp~m_cos_asp))
abline(var.temp, lwd=3)

# Making RHF height data long form
# Will Pearse # 2018-10-29 #
# Combining 2018-2019 height data
# Elizabeth Simpson # 2019-07-31

# Load data
data.18 <- read.csv("./functional_traits_rhf/raw_data/rhf_2018_traits_height.csv", as.is=TRUE)

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
data.19 <- read.csv("./functional_traits_rhf/raw_data/rhf_2019_traits_height_nearest_temp_plot.csv", as.is=TRUE)
data.19 <- data.19[,1:7]
colnames(data.19)[7] <- "height"

# put all height together 
all.ht <- expanded.data.18 #can included data.19 with rbind()
all.ht <- subset(all.ht, Plot_id %in% temp$Plot_id)

#############################################################################################################
#Calculate height metrics
mean.ht <- as.data.frame(with(all.ht, tapply(height, Plot_id, function(x) mean(x,na.rm=TRUE))))
max.ht <- as.data.frame(with(all.ht, tapply(height, Plot_id, function(x) max(x,na.rm=TRUE))))
min.ht <- as.data.frame(with(all.ht, tapply(height, Plot_id, function(x) min(x,na.rm=TRUE))))
var.ht <- as.data.frame(with(all.ht, tapply(height, Plot_id, function(x) var(x,na.rm=TRUE))))

#Figure out how to get one value of height for each species, use that in pez to calculate CWM an dispersion for the species you have....
terr.temp <- subset(terr.temp, Plot_id %in% as.integer(rownames(mean.ht)))
all.data <- cbind(terr.temp, min.ht, mean.ht, max.ht, var.ht)
colnames(all.data) <- c("Plot_id", "m_cos_asp", "min_temp", "mean_temp", "max_temp", "var_temp", "min_ht", "mean_ht", "max_ht", "var_ht")


with(all.data, plot(mean_ht~m_cos_asp, pch=19, cex=2, ylim=c(0, 165), xlab="", ylab="", cex.axis=2))
with(all.data, points(min_ht~m_cos_asp, pch=19, cex=2, col="deepskyblue3"))
with(all.data, points(max_ht~m_cos_asp, pch=19, cex=2, col="chocolate2"))
mean.ht <- with(all.data, lm(mean_ht~m_cos_asp))
min.ht <- with (all.data, lm(min_ht~m_cos_asp))
max.ht <- with(all.data, lm(max_ht~m_cos_asp))
abline(mean.ht, lwd=3)
abline(min.ht, col="deepskyblue3", lwd=3)
abline(max.ht, col="chocolate2", lwd=3)
#variance in temperature across aspect
with(all.data, plot(var_height~m_cos_asp, pch=19, xlab="Aspect", ylab="Temperature variance"))

#this is the ultimately important model for variance
with(all.data, plot(var_height~var_temp, pch=19, xlab="Variance in Temperature", ylab="Variance in Height"))
model <- with(all.data, lm(var_height~var_temp))
abline(model)

#versus the one for mean values....
with(all.data, plot(mean_height~mean_temp, pch=19, xlab="Mean Temperature", ylab="Mean Height"))
model <- with(all.data, lm(mean_height~mean_temp))
abline(model)

with(all.data, plot(mean_height~mean_temp, pch=19))
with(all.data, plot(max_height~mean_temp, pch=19))
with(all.data, plot(var_height~mean_temp, pch=19))
with(all.data, plot(var_temp~mean_temp, pch=19))
with(all.data, plot(var_height~mean_height, pch=19))

