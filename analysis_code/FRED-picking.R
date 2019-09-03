setwd("~/Documents/projects/functional_traits_rhf")
data <- read.csv("./raw_data/FRED/FRED2_20180518.csv", as.is=TRUE)

# reads back this error...
# Error in type.convert.default(data[[i]], as.is = as.is[i], dec = dec,  : invalid multibyte string at '<93>x<94> de<6e>otes that FRED includes duplicates of data. This may be due to the presence of both unpublished and published data, with additional individual replicates corresponding with previously published means, or data provided in multiple sources.'

# also, data loads in in a really strange format...