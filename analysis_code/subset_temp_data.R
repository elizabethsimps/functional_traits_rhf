# Subset temperature data change plot ID
# Elizabeth Simpson # 2019-07-29
setwd("~/Documents/projects/functional_traits_rhf")
library(tidyverse)
library(lubridate)

# load data
### Double check that this is correclty formatted ###
temp <- read.csv("./clean_data/17-18_temp_data.csv", as.is=TRUE)

temp$datetime <- with(temp, paste0(year,"-",month,"-",day,"T",time))
temp$datetime <- as.POSIXct(temp$datetime,
                                    format = "%Y-%m-%dT%H:%M",
                                    tz = "MST")

# summarize date ranges to get time period completely covered by sensors
date.range <- as.data.frame(temp %>%
                              group_by(plot) %>%
                              summarize(first.m=first(month),first.d=first(day), first.t=first(time),last.m=last(month),last.d=last(day), last.t=last(time)))

# subset data
temp.sub <- subset(temp,
                   datetime >= as.POSIXct('2017-09-28 00:00',
                                          tz = "MST") &
                     datetime <= as.POSIXct('2018-09-12 00:00',
                                            tz = "MST"))
temp.sub <- temp.sub[,1:6]

######################################################
### START ANALYSIS ... what are you interested in? ###
