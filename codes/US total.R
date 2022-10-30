
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

#### Gompertz makeham #####

age.interval <- c(30,110)

age <- (age.interval[1]:age.interval[2])
year.interval<-c(as.numeric(colnames(Ex)[1]),
                 as.numeric(colnames(Ex)[ncol(Ex)]))

death.counts <- Dx[(age.interval[1]+1):(age.interval[2]+1), 
                   as.character(c(year.interval[1]:year.interval[2]))]
expo.counts <- Ex[(age.interval[1]+1):(age.interval[2]+1), 
                  as.character(c(year.interval[1]:year.interval[2]))]

gomp <- Decomp.Makeham(death.counts,expo.counts,
                       age.start = 30,age.end = 110)

#### Time trends ####

categories <- c(rep(1:29,each=3))

c5.ex.make<- tapply(gomp$deltac , categories , sum )
B5.ex.make<- tapply(gomp$deltaB , categories , sum )
M5.ex.make<- tapply(gomp$deltaM , categories , sum )
e5.ex.make<- tapply(gomp$deltae , categories , sum )

data <- data.table(label = c(rep("Background",29),
                             rep("Compression",29),
                             rep("Shift",29)),
                   value = c(c5.ex.make,
                             B5.ex.make,
                             M5.ex.make),
                   year = c(paste(seq(1933,2017,3),"-",
                                      seq(1936,2020,3),
                                  sep=""))
                   )

data$label <- factor(data$label,
                     levels = c("Background",
                                "Compression",
                                "Shift"))

ggplot(data, aes(x=year,y=value))+
  geom_col(aes(fill=label))+
  stat_summary(fun = sum,geom = "point")+
  scale_fill_manual(values = c("lightblue","navy","orange"))+
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

ggsave("report/USA_Male,three_comp.png",
       width = 8,height = 6)


#### Last ten years ####

data2 <- data.table(label = c(rep("Background",10),
                              rep("Compression",10),
                              rep("Shift",10)),
                    value = c(gomp$deltac[78:87],
                              gomp$deltaB[78:87],
                              gomp$deltaM[78:87]),
                    year = c(paste(seq(2010,2019,1),"-",
                                   seq(2011,2020,1),sep=""))
)

data2 <- data2[,year2 := fcase(year %in% c("2011-2012",
                                         "2012-2013",
                                         "2013-2014",
                                         "2014-2015"),
                             year2 = "2011-2015",
                             year %in% c("2015-2016",
                                         "2016-2017",
                                         "2017-2018",
                                         "2018-2019"),
                             year2 = "2015-2019",
                             year == "2019-2020",
                             year2 = "2019-2020")]

data2 <- data2[,list(value=sum(value)),
             by = list(year2,label)]

data2 <- data2[!is.na(year2),]

data2$label <- factor(data2$label,
                     levels = c("Background",
                                "Compression",
                                "Shift"))

ggplot(data2, aes(x=year2,y=value))+
  geom_col(aes(fill=label))+
  stat_summary(fun = sum,geom = "point")+
  scale_fill_manual(values = c("lightblue","navy","orange"))+
  scale_y_continuous(limits = c(-5,2))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1.2,
                                   hjust = 1.2),
        legend.position = "right")+
  labs(x="Year",y="Contribution",
       title = "Changes in Life expectancy for the US with Gompertz-Makeham Model,
       Male 2011-2020")+
  guides(fill=guide_legend("Components"))

ggsave("report/USA_Male,three_comp_2011-2020.png",
       width = 8,height = 6)


fwrite(data2[,value:=round(value,3)],"report/table total_HMD _Male 2011-2020.txt")

#### Life expectancy ####

ex <- data.table(
  e30 = gomp$ex,
  year = seq(1933,2020)
)

ggplot(ex,aes(year,e30))+
  geom_line()+
  scale_x_continuous(n.breaks = 10)+
  scale_y_continuous(n.breaks = 10)+
  theme_classic()

ggsave("report/time-trends_e30_Male.png",
       width = 6,height = 6)
