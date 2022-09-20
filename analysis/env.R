# Microenvironment at core RHF plots: soil temp (2017-18 and 2018-19), texture, & terrain
# Elizabeth Simpson # 2019-07-29 | 2021-1-25
# Set working directory
# setwd("~/Documents/projects/functional_traits_rhf")

library(soiltexture)
library(lubridate)
library(readxl)
library(weathermetrics)
library(tidyverse)

#############
# FUNCTIONS #
#############

### To calculate soil texture ###

# Temperature correction function for hydrometer readings (+0.4 g/L for each degree above 20 C, - 0.4 g/L for each degree below 20 C)
temp.correct <- function(temp, hydro_rdg){
  return(hydro_rdg + ((temp-20)*0.4))
}

# Soil texture calculations
texture <- function(soil, plot_ID, soil_wt, temp_calib, hydro_calib, temp_40, hydro_40, temp_120, hydro_120){
  clean_soil <- matrix(nrow=nrow(soil), ncol=3)
  for(i in 1:nrow(soil)){
    B <- with(soil, temp.correct(temp_calib[i], hydro_calib[i]))
    A_40 <- with(soil, temp.correct(temp_40[i], hydro_40[i]))
    A_120 <- with(soil, temp.correct(temp_120[i], hydro_120[i]))
    silt_clay <- ((A_40-B)/soil_wt[i])*100
    clean_soil[i,1] <- 100-silt_clay #sand
    clean_soil[i,2] <- ((A_120-B)/soil_wt[i])*100 #clay
    clean_soil[i,3] <- 100 - (100-silt_clay) - (((A_120-B)/soil_wt[i])*100) #silt
  }
  colnames(clean_soil) <- c("SAND", "CLAY", "SILT")
  output <- as.data.frame(cbind(plot_ID, clean_soil))
  output$SAND <- as.numeric(output$SAND)
  output$CLAY <- as.numeric(output$CLAY)
  output$SILT <- as.numeric(output$SILT)
  return(output)
}

# To format times in the temperature dataset - based on code from Michael Stemkovski (MS)
format.times <- function(x){ 
  #x <- all_data[4,1]
  hour <- as.character(hour(x))
  minute <- as.character(minute(x))
  if(nchar(hour) == 1) hour <- paste(0,hour,sep="")
  if(nchar(minute) == 1) minute <- paste(0,minute,sep="")
  return(paste(hour,minute,sep=":"))
}

#############
# LOAD DATA #
#############

### Topography - clean data from fractal paper
topo <- read.csv("./clean_data/clean_rhf_terrain_17-18.csv", as.is=TRUE)

### Soil texture
s.text <- read.csv("./raw_data/rhf_2018_soil_texture.csv", as.is=TRUE)

### Temperature - modified from code from MS
# Note, converted from .hobo file to .csv in HoboWare
temp.csvs <- list.files("./raw_data/17-18_temp/")[which(grepl("csv",list.files("raw_data/17-18_temp")))]
temp.csvs19 <- list.files("./raw_data/18-19_temp/")[which(grepl("csv",list.files("raw_Data/18-19_temp")))]

##############
# CLEAN DATA #

### Topography - remove extra first col and TPI from the dataframe
topo <- topo[,-c(1,10)]
topo$plot_id <- as.character(topo$plot_id)

### Soil texture
s.text <- na.omit(s.text)
text <- with(s.text, texture(s.text, plot_ID, soil_wt, temp_calib, hydro_calib, temp_40, hydro_40, temp_120, hydro_120))
text <- cbind(text, TT.points.in.classes(tri.data = text, class.sys = "USDA.TT", PiC.type="t"))
colnames(text) <- c("plot_id", "SAND", "CLAY", "SILT", "class")

### 2017-2018 temperature cleaning - based on code from MS
temp <- data.frame(datetime = rep(NA,10000*26),
                   temp = rep(NA,10000*26),
                   plot = rep(NA,10000*26))
top <- 1
for(i in temp.csvs){
  file <- read.csv(paste("raw_data/17-18_temp/",i,sep=""), header=F, as.is=TRUE)[,c(2,3)]
  colnames(file) <- c("datetime", "temp")  
  file <- file[-c(1,2),]
  file[,1] <- file[,1]
  file[,2] <- as.numeric(file[,2])
  file[,3] <- as.integer(substr(i,1,4)) 
  bottom <- top + nrow(file) - 1
  temp[top:bottom,] <- file
  top <- bottom + 1
}

temp <- temp[-which(is.na(temp$datetime)),]
temp[,1] <- mdy_hms(temp[,1], tz = "America/Denver") #if this throws an error, re install lubridate & restart R

# this next line takes a couple minutes to run
f_temp <- data.frame(plot = temp[,3],
                     datetime = temp[,1],
                     year = year(temp[,1]),
                     month = month(temp[,1]),
                     day = day(temp[,1]),
                     time = sapply(temp[,1], format.times),
                     temp = fahrenheit.to.celsius(temp[,2]))

f_temp <- na.omit(f_temp)

### 2018 - 2019 temperature cleaning
temp19 <- data.frame(datetime = rep(NA,10000*24),
                     temp = rep(NA,10000*24),
                     plot = rep(NA,10000*24))
top <- 1
for(i in temp.csvs19){
  file <- read.csv(paste("raw_data/18-19_temp/",i,sep=""), header=F, as.is=TRUE)[,c(2,3)]
  colnames(file) <- c("datetime", "temp")  
  file <- file[-c(1,2),]
  file[,1] <- file[,1]
  file[,2] <- as.numeric(file[,2])
  file[,3] <- as.integer(substr(i,1,4)) 
  bottom <- top + nrow(file) - 1
  temp19[top:bottom,] <- file
  top <- bottom + 1
}

temp19 <- temp19[-which(is.na(temp19$datetime)),]
temp19[,1] <- mdy_hms(temp19[,1], tz = "America/Denver") #if this throws an error, re install lubridate & restart R

# this next line takes a couple minutes to run
f_temp19 <- data.frame(plot = temp19[,3],
                       datetime = temp19[,1],
                       year = year(temp19[,1]),
                       month = month(temp19[,1]),
                       day = day(temp19[,1]),
                       time = sapply(temp19[,1], format.times),
                       temp = fahrenheit.to.celsius(temp19[,2]))

f_temp19 <- na.omit(f_temp19)

### Combine temp datas for full temperature data range 2017-2019
f_temp <- rbind(f_temp, f_temp19)

# Summarize date ranges to get time period completely covered by sensors - all ES code from hereon
date.range <- as.data.frame(f_temp %>%
                              group_by(plot) %>%
                              summarize(first.m=first(month),first.d=first(day), first.t=first(time), first.y=first(year),last.m=last(month),last.d=last(day), last.t=last(time), last.y=last(year)))

# Subset a full year of temperature data from 2017 -> 2018 matched to look at how microenvironment varies across topography
f_temp_year <- subset(f_temp,
                datetime >= as.POSIXct('2017-09-28 00:00:00',
                                        tz = "MST") &
                   datetime <= as.POSIXct('2018-09-28 00:00:00',
                                          tz = "MST"))

# Calculate statistics by plot, could also summarize at the monthly level first and then at plot
year_17.18 <- as.data.frame(group_by(f_temp_year, plot) %>%
                           summarize(mean(temp), sd(temp), max(temp), min(temp)))
colnames(year_17.18) <- c("plot_id", "mean", "sd", "max", "min")

# Subset date range for analyzing temp and 2019 cover data as the year prior to when the 2019 cover data was collected
# Many of the temperatures sensors power reset and otherwise didn't work within this date range leaving only twelve good ones
good_temp_19 <- c("1111","1122","1233","1311","1322","1333","2211","2222","2311","3133","3211","3333")

f_temp_19 <- f_temp[f_temp$plot %in% good_temp_19,]

f_temp_cvr19 <- subset(f_temp_19,
                       datetime >= as.POSIXct('2018-06-07 00:00:00',
                                              tz = "MST") &
                         datetime <= as.POSIXct('2019-06-07 00:00:00',
                                                tz = "MST"))

# Calculate statistics by plot, could also summarize at the monthly level first and then at plot
temp_cvr19_year <- as.data.frame(group_by(f_temp_cvr19, plot) %>%
                              summarize(mean(temp), sd(temp), max(temp), min(temp)))
colnames(temp_cvr19_year) <- c("plot_id", "mean", "sd", "max", "min")

###########################
# Put all env data together & export clean environmental data sheets

# Terrain and texture to use in analysis - 76 plots
env <- merge(topo, text, by.x="plot_id", by.y="plot_id")
write.csv(env, "./clean_data/texture-terrain-18.csv")

# Terrain, texture, and temperature (2017-18) - 25 plots 
env.t <- merge(env, year_17.18, by.x="plot_id", by.y="plot_id")
write.csv(env.t, "./clean_data/temp-year17-18_text_terr.csv")

# Terrain, texture, and temperature (year before 2019 cover survey) - 12 plots
env.cvr19.year <- merge(env, temp_cvr19_year, by.x="plot_id", by.y="plot_id")
write.csv(env.cvr19.year, "./clean_data/temp-CVR-year18-19_text_terr.csv")