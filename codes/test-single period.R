
###
### 
###

source("codes/functions.R")

library(data.table)
library(RColorBrewer)
library(ggplot2)

Dx <- fread("data/USA.Deaths_1x1.txt")[,list(Year,Age,Female)]

Dx <- matrix(data = Dx$Female,nrow = 111,ncol = length(unique(Dx$Year)),
             dimnames = list(unique(Dx$Age),unique(Dx$Year)))

Ex <- fread("data/USA.Exposures_1x1.txt")[,list(Year,Age,Female)]

Ex <- matrix(data = Ex$Female,nrow = 111,ncol = length(unique(Ex$Year)),
             dimnames = list(unique(Ex$Age),unique(Ex$Year)))

age.interval <- c(30,110)
age <- (age.interval[1]:age.interval[2])
year.interval<-c(as.numeric(colnames(Ex)[1]),
                 as.numeric(colnames(Ex)[ncol(Ex)]))

death.counts <- Dx[(age.interval[1]+1):(age.interval[2]+1), 
                   as.character(c(year.interval[1]:year.interval[2]))]
expo.counts <- Ex[(age.interval[1]+1):(age.interval[2]+1), 
                  as.character(c(year.interval[1]:year.interval[2]))]



gomp <- Decomp.Gompertz(death.counts,expo.counts,age.start = 30,age.end = 110)

v <- gomp$deltaB
n <- 5

categories <- c(rep(1,2),rep(2:18,each=5))

B5.ex.gomp<- tapply(gomp$deltaB , categories , sum )
M5.ex.gomp<- tapply(gomp$deltaM , categories , sum )
e5.ex.gomp<- tapply(gomp$deltae , categories , sum )

data <- data.table(label = c(rep("compression",18),rep("shift",18)),
                   value = c(B5.ex.gomp,M5.ex.gomp),
                   year = rep(c("1934-1935",
                                paste(seq(1936,2016,5),"-",
                                      seq(1940,2020,5)))
                              )
                   )

data$label <- factor(data$label,levels = c("shift","compression"))

ggplot(data, aes(x=year,y=value,fill=label))+
  geom_col()+
  scale_fill_manual(values = c("orange","navy"))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 45,vjust = 1.2,hjust = 1.2))+
  labs(x="Year",y="Contribution")+
  guides(fill=guide_legend("Components"))

ggsave("report/Replication-USA,two_comp.pdf",
       width = 6,height = 6)
