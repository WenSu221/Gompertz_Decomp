
###
### Decomposition of the Life expectancy changes
### into two components, the shift and the compression
###

##  This is primarily a replication of the work done by Bergeron-
##  Boucher et al., 2015. 

source("codes/functions.R")

library(data.table)
library(RColorBrewer)
library(ggplot2)

Dx <- fread("data/USA.Deaths_1x1.txt")[,list(Year,Age,Male)]

Dx <- matrix(data = Dx$Male,nrow = 111,ncol = length(unique(Dx$Year)),
             dimnames = list(unique(Dx$Age),unique(Dx$Year)))

Ex <- fread("data/USA.Exposures_1x1.txt")[,list(Year,Age,Male)]

Ex <- matrix(data = Ex$Male,nrow = 111,ncol = length(unique(Ex$Year)),
             dimnames = list(unique(Ex$Age),unique(Ex$Year)))

age.interval <- c(30,110)
age <- (age.interval[1]:age.interval[2])
year.interval<-c(as.numeric(colnames(Ex)[1]),
                 as.numeric(colnames(Ex)[ncol(Ex)]))

death.counts <- Dx[(age.interval[1]+1):(age.interval[2]+1), 
                   as.character(c(year.interval[1]:year.interval[2]))]
expo.counts <- Ex[(age.interval[1]+1):(age.interval[2]+1), 
                  as.character(c(year.interval[1]:year.interval[2]))]


# gomp <- Decomp.Gompertz(death.counts,expo.counts,
#                         age.start = 30,age.end = 110)

gomp <- Decomp.Makeham(death.counts,expo.counts,
                       age.start = 30,age.end = 110)

v <- gomp$deltaB
n <- 5

categories <- c(rep(1,2),rep(2:18,each=5))

c5.ex.make<- tapply(gomp$deltac  , categories , sum )
B5.ex.make<- tapply(gomp$deltaB , categories , sum )
M5.ex.make<- tapply(gomp$deltaM , categories , sum )
e5.ex.make<- tapply(gomp$deltae , categories , sum )

data <- data.table(label = c(rep("Background",18),
                             rep("Compression",18),
                             rep("Shift",18)),
                   value = c(c5.ex.make,
                             B5.ex.make,
                             M5.ex.make),
                   year = rep(c("1934-1935",
                                paste(seq(1936,2016,5),"-",
                                      seq(1940,2020,5)))
                              )
                   )

data$label <- factor(data$label,
                     levels = c("Background",
                                "Compression",
                                "Shift"))

ggplot(data, aes(x=year,y=value,fill=label))+
  geom_col()+
  scale_fill_manual(values = c("blue","navy","orange"))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1.2,
                                   hjust = 1.2),
        legend.position = "bottom")+
  labs(x="Year",y="Contribution",
       title = "Changes in Life expectancy for the US with Gompertz-Makeham Model,
       Male 1933-2020")+
  guides(fill=guide_legend("Components"))

ggsave("report/USA_Male,three_comp.pdf",
       width = 8,height = 6)
