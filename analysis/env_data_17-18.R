# Analyze RHF environmental data: soil temp (2017-18) & texture & terrain
# Elizabeth Simpson # 2019-07-29 | 2021-1-25

setwd("~/Documents/projects/functional_traits_rhf")

################# NOTES #######################################
### Temperature data - convert from .hobo file to .csv in HoboWare
### Use cleaning script to format the data for all years and plots
### Be aware that tempearture and texture data are not always collected at the same plots

################ LOAD PACKAGES ##################
library(tidyverse)
library(lubridate)
library(dplyr)
library(ggplot2)
library(soiltexture)

################# FUNCTIONS ####################

# Temperature correction function for hydrometer readings (+0.4 g/L for each degree above 20 C, - 0.4 g/L for each degree below 20 C)
temp.correct <- function(temp, hydro_rdg){
  return(hydro_rdg + ((temp-20)*0.4))
}

# Soil texture calculations
texture <- function(plot_ID, soil_wt, temp_calib, hydro_calib, temp_40, hydro_40, temp_120, hydro_120){
  clean_soil <- matrix(nrow=nrow(soil), ncol=4)
  colnames(clean_soil) <- c("plot_ID", "SAND", "CLAY", "SILT")
  clean_soil[,1] <- as.numeric(plot_ID)
  for(i in 1:nrow(soil)){
    B <- with(soil, temp.correct(temp_calib[i], hydro_calib[i]))
    A_40 <- with(soil, temp.correct(temp_40[i], hydro_40[i]))
    A_120 <- with(soil, temp.correct(temp_120[i], hydro_120[i]))
    silt_clay <- ((A_40-B)/soil_wt[i])*100
    clean_soil[i,2] <- as.numeric(100-silt_clay) #sand
    clean_soil[i,3] <- as.numeric(((A_120-B)/soil_wt[i])*100) #clay
    clean_soil[i,4] <- as.numeric(100 - (100-silt_clay) - (((A_120-B)/soil_wt[i])*100)) #silt
  }
  return(clean_soil)
}


################# LOAD DATA ####################
# Soil temperature
temp <- read.csv("./clean_data/17-18_temp_data.csv", as.is=TRUE)

# Soil texture
soil <- read.csv("./raw_data/rhf_2018_soil_texture.csv", as.is=TRUE)
soil <- na.omit(soil)

# Terrain - elevation in meters, aspect in radians, slope in degrees
terrain <- read.csv("clean_data/rhf_2018_terrain_measured_cleaned.csv", as.is=TRUE)

########### Subset temperature data ############
temp$datetime <- with(temp, paste0(year,"-",month,"-",day,"T",time))
temp$datetime <- as.POSIXct(temp$datetime,
                                    format = "%Y-%m-%dT%H:%M",
                                    tz = "MST")

# summarize date ranges to get time period completely covered by sensors - note this only works for one 'season' currently - example 2017-2018
# visually assess first and last date that covers the same range
date.range <- as.data.frame(temp %>%
                              group_by(plot) %>%
                              summarize(first.m=first(month),first.d=first(day),last.m=last(month),last.d=last(day), last.t=last(time)))

# subset data for 2017 and 2018 matched start and end. 
# NOTE: sensor at plot 3133 power reset / stopped on 2018-7-29 @ 13:30
temp.s <- subset(temp,
                   datetime >= as.POSIXct('2017-09-28 00:00',
                                          tz = "MST") &
                     datetime <= as.POSIXct('2018-09-12 00:00',
                                            tz = "MST"))

############## A big overview plot ############
qplot(x=datetime, y=temp, 
      data=temp.s,
      main="2017-2018 Daily Air Temperature at RHF sites",
      color=as.factor(plot)) +
      geom_point(shape=1)

###############################################
############# Bring data together #############

######## SOIL TEMPERATURE - 26 plots ##########
# Group temp data by plot
plot.temp <- group_by(temp.s, plot)
# Calculate mean temp from all raw points
temp.stats <- as.data.frame(dplyr::summarize(plot.temp, mean(temp)))
# Calculate SD as a measure of dispersion from the mean
temp.stats <- cbind(temp.stats, as.data.frame(dplyr::summarize(plot.temp, sd(temp)))[,2])
# Calculate min and max to calculate range for each plot
temp.stats <- cbind(temp.stats, as.data.frame(dplyr::summarize(plot.temp, max(temp)))[,2])
temp.stats <- cbind(temp.stats, as.data.frame(dplyr::summarize(plot.temp, min(temp)))[,2])
temp.stats <- cbind(temp.stats, as.data.frame(temp.stats[,4] - temp.stats[,5]))
colnames(temp.stats) <- c("plot", "mean", "sd", "max", "min", "range")

######### SOIL TEXTURE - 76 plots ##############
soil.text <- as.data.frame(with(soil, texture(plot_ID, soil_wt, temp_calib, hydro_calib, temp_40, hydro_40, temp_120, hydro_120)))
soil.text <- cbind(soil.text, TT.points.in.classes(tri.data = soil.text, class.sys = "USDA.TT", PiC.type="t"))
colnames(soil.text) <- c("plot", "SAND", "CLAY", "SILT", "class")

# plot onto USDA texture triangle 
TT.plot(class.sys = "USDA.TT", tri.data = soil.text, pch=20, cex=0.8, cex.axis=0.8, cex.lab=0.8, col="red", main="")

######## Add in TERRAIN - 78 plots #############
text.terr <- merge(terrain, soil.text, by.x = "Plot_id", by.y = "plot") # stop here to look at all 76 plots here for texture
env <- merge(text.terr, temp.stats, by.x = "Plot_id", by.y = "plot") # results in 25 plots with temp and texture for 2017-2018
soil.env <- env[,-2:-4]

write.csv(soil.env, "2018-soil-env.csv")

########## Exploratory plotting ################

# How soil classes distribute across Mean, SD, and Range of soil temperature at 25 core plots at RHF
par(mar=c(4,4.5,1,1))
par(mfrow=c(1,3))
with(env, boxplot(mean~class, ylab = "USDA Soil Classification", xlab = expression(paste('Mean Temperature (',~degree,'C)',sep='')), horizontal=TRUE))
with(env, boxplot(sd~class, ylab = "USDA Soil Classification", xlab = "SD of Mean Temperature", horizontal=TRUE))
with(env, boxplot(range~class, ylab = "USDA Soil Classification", xlab = expression(paste('Temperature Range (',~degree,'C)',sep='')), horizontal=TRUE))

# How soil classes distribute across terrain at 25 core plots at RHF
par(mar=c(4,4.5,1,1))
par(mfrow=c(1,3))
with(env, boxplot(elev_m~class, ylab = "USDA Soil Classification", xlab = "Elevation (meters)", horizontal=TRUE))
with(env, boxplot(aspect_rad~class, ylab = "USDA Soil Classification", xlab = "Aspect (S = -1, N = 1)", horizontal =TRUE))
with(env, boxplot(slope_deg~class, ylab = "USDA Soil Classification", xlab = expression(paste('Slope (',~degree,')',sep='')), horizontal=TRUE))

# How temperature stats vary across percent particle size
par(mar=c(4.5,4.5,1,1))
with(env, plot(mean~SAND, pch=19, col="goldenrod", xlim=c(0,80), xlab = "Percent particle size", ylab=expression(paste('Mean Temperature (',~degree,'C)',sep=''))))
with(env, points(mean~CLAY, pch=19, col="brown", xlim=c(0,80)))
with(env, points(mean~SILT, pch=19, col="darkslategray", xlim=c(0,80)))

with(env, plot(sd~SAND, pch=19, col="goldenrod", xlim=c(0,80), xlab = "Percent particle size", ylab="SD of Mean Temperature"))
with(env, points(sd~CLAY, pch=19, col="brown", xlim=c(0,80)))
with(env, points(sd~SILT, pch=19, col="darkslategray", xlim=c(0,80)))

with(env, plot(range~SAND, pch=19, col="goldenrod", xlim=c(0,80), xlab = "Percent particle size", ylab=expression(paste('Temperature Range (',~degree,'C)',sep=''))))
with(env, points(range~CLAY, pch=19, col="brown", xlim=c(0,80)))
with(env, points(range~SILT, pch=19, col="darkslategray", xlim=c(0,80)))
legend(60,70, legend=c("Sand", "Clay", "Silt"), pch=19,col=c("goldenrod", "brown", "darkslategray"), cex=0.9, y.intersp=1.5)

# How temp stats vary across terrain
par(mfrow=c(2,3))
with(env, plot(mean~elev_m, pch=19, ylim=c(-8,65), xlab = "Elevation (meters)", ylab=expression(paste('Mean Temperature (',~degree,'C)',sep=''))))
with(env, points(max~elev_m, pch=19, ylim=c(-8,65), col="red"))
with(env, points(min~elev_m, pch=19, ylim=c(-8, 65), col="blue"))

with(env, plot(mean~aspect_rad, pch=19, ylim=c(-8,65), xlab = "Aspect (S = -1, N = 1)", ylab=""))
with(env, points(max~aspect_rad, pch=19, ylim=c(-8,65), col="red"))
with(env, points(min~aspect_rad, pch=19, ylim=c(-8, 65), col="blue"))
abline(with(env, lm(max~aspect_rad)), col="red")
abline(with(env, lm(mean~aspect_rad)))

with(env, plot(mean~slope_deg, pch=19, ylim=c(-8,65), xlab = expression(paste('Slope (',~degree,')',sep='')), ylab="" ))
with(env, points(max~slope_deg, pch=19, ylim=c(-8,65), col="red"))
with(env, points(min~slope_deg, pch=19, ylim=c(-8, 65), col="blue"))

with(env, plot(sd~elev_m, pch=19, xlab = "Elevation (meters)", ylab="SD of Mean Temperature"))
with(env, plot(sd~aspect_rad, pch=19, xlab = "Aspect (S = -1, N = 1)", ylab=""))
abline(with(env, lm(sd~aspect_rad)))
with(env, plot(sd~slope_deg, pch=19, xlab = expression(paste('Slope (',~degree,')',sep='')), ylab=""))


