# Subset temperature data change plot ID
# Elizabeth Simpson # 2019-07-29
setwd("~/Documents/projects/")
library(tidyverse)
library(lubridate)

# load data
temp <- read.csv("./functional_traits_rhf/clean_data/formatted_temp_data.csv", as.is=TRUE)

# reformat plot ID to match
temp <- temp %>% extract(Plotcode, c("first", "second", "third", "fourth"), "([[:digit:]]+)([[:digit:]]+)([[:digit:]]+)([[:digit:]]+)", remove = FALSE)
temp$Plot_id <- with(temp, paste0(first, second, third, fourth))
temp$Plot_id <- as.integer(temp$Plot_id)
temp <- temp[,-c(1:5)]
temp <- temp[c(6,1:5)]

temp$datetime <- with(temp, paste0(Year,"-",Month,"-",Day,"T",Time))
temp$datetime <- as.POSIXct(temp$datetime,
                                    format = "%Y-%m-%dT%H:%M",
                                    tz = "MST")
# subset data - 2009-2011
temp.sub <- subset(temp,
                    datetime >= as.POSIXct('2017-09-28 00:00',
                                                tz = "MST") &
                    datetime <= as.POSIXct('2018-09-12 00:00',
                                                tz = "MST"))

# test date range
date.range <- as.data.frame(temp.sub %>%
                              group_by(Plot_id) %>%
                              summarize(first.m=first(Month),first.d=first(Day), first.t=first(Time),last.m=last(Month),last.d=last(Day), last.t=last(Time)))

# make dataframe and save
temp.sub <- temp.sub[,1:6]
