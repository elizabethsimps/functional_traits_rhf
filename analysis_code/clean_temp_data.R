### Cleaning Right Hand Fork soil temperature data
# Michael Stemkovski
# modified to analyze for RHF by Elizabeth Simpson
setwd("~/Documents/projects/functional_traits_rhf")

library(lubridate)
library(readxl)
library(weathermetrics)

csvs <- list.files("raw_data/17-18_temp/")[which(grepl("csv",list.files("raw_data/17-18_temp")))]

all_data <- data.frame(datetime = rep(NA,10000*27),
                       temp = rep(NA,10000*27),
                       plot = rep(NA,10000*27))
top <- 1
for(i in csvs){
  file <- read.csv(paste("raw_data/17-18_temp/",i,sep=""), header=F, as.is=TRUE)[,c(2,3)]
  colnames(file) <- c("datetime", "temp")
  file <- file[-c(1,2),]
  file[,1] <- file[,1]
  file[,2] <- as.numeric(file[,2])
  file[,3] <- as.integer(substr(i,1,4))
  bottom <- top + nrow(file) - 1
  all_data[top:bottom,] <- file
  top <- bottom + 1
}
all_data <- all_data[-which(is.na(all_data$datetime)),]
all_data[,1] <- mdy_hms(all_data[,1], tz = "America/Denver")

format.times <- function(x){
  #x <- all_data[4,1]
  hour <- as.character(hour(x))
  minute <- as.character(minute(x))
  if(nchar(hour) == 1) hour <- paste(0,hour,sep="")
  if(nchar(minute) == 1) minute <- paste(0,minute,sep="")
  return(paste(hour,minute,sep=":"))
}

formatted_data <- data.frame(plot = all_data[,3],
                             year = year(all_data[,1]),
                             month = month(all_data[,1]),
                             day = day(all_data[,1]),
                             time = sapply(all_data[,1], format.times),
                             temp = fahrenheit.to.celsius(all_data[,2]))

formatted_data <- na.omit(formatted_data) #removes 63 rows

write.csv(formatted_data, "clean_data/17-18_temp_data.csv",row.names = F)