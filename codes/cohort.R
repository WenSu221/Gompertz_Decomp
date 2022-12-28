source("codes/functions.R")

library(data.table)
library(RColorBrewer)
library(ggplot2)

Ex <- fread(paste("data/HMD/","cExposures_1x1.txt",sep=""))[,list(Year,Age,Male)]

year.range <- unique(Ex$Year)

Dx <- fread(paste("data/HMD/","cMx_1x1.txt",sep=""))[,list(Year,Age,Male)]

Ex <- matrix(data = Ex$Male,nrow = 111,ncol = length(unique(Ex$Year)),
             dimnames = list(unique(Ex$Age),unique(Ex$Year)))

Dx <- matrix(data = Dx$Male,nrow = 111,ncol = length(unique(Dx$Year)),
             dimnames = list(unique(Dx$Age),unique(Dx$Year)))

age.interval <- c(0,110)

age <- (age.interval[1]:age.interval[2])

year.interval<-c(1835,1910)

death.counts <- as.matrix(Dx[(age.interval[1]+1):(age.interval[2]+1), 
                   as.character(c(year.interval[1]:year.interval[2]))])

expo.counts <- as.matrix(Ex[(age.interval[1]+1):(age.interval[2]+1), 
                  as.character(c(year.interval[1]:year.interval[2]))])

death.counts <- t(apply(death.counts,1,function(x){as.numeric(x)}))

expo.counts <- t(apply(expo.counts,1,function(x){as.numeric(x)}))

death.counts <- t(apply(death.counts,1,function(x){fifelse(is.na(x),0,as.numeric(x))}))

expo.counts <- t(apply(expo.counts,1,function(x){fifelse(is.na(x),0,as.numeric(x))}))

gomp <- Decomp.Siler(death.counts,expo.counts,
                     age.start = 0,age.end = 110)
