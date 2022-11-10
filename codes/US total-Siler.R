
####
#### Decomp with Siler
####

##  This is primarily a replication of the work done by Bergeron-
##  Boucher et al., 2015. 

source("codes/functions.R")

library(data.table)
library(RColorBrewer)
library(ggplot2)

cnty <- "USA"

Dx <- fread(paste("data/",cnty,".Deaths_1x1.txt",sep=""))[,list(Year,Age,Male)]

year.range <- unique(Dx$Year)

Dx <- matrix(data = Dx$Male,nrow = 111,ncol = length(unique(Dx$Year)),
             dimnames = list(unique(Dx$Age),unique(Dx$Year)))

Ex <- fread(paste("data/",cnty,".Exposures_1x1.txt",sep=""))[,list(Year,Age,Male)]

Ex <- matrix(data = Ex$Male,nrow = 111,ncol = length(unique(Ex$Year)),
             dimnames = list(unique(Ex$Age),unique(Ex$Year)))

age.interval <- c(0,110)

age <- (age.interval[1]:age.interval[2])
year.interval<-c(as.numeric(colnames(Ex)[1]),
                 as.numeric(colnames(Ex)[ncol(Ex)]))

death.counts <- Dx[(age.interval[1]+1):(age.interval[2]+1), 
                   as.character(c(year.interval[1]:year.interval[2]))]
expo.counts <- Ex[(age.interval[1]+1):(age.interval[2]+1), 
                  as.character(c(year.interval[1]:year.interval[2]))]

gomp <- Decomp.Siler(death.counts,expo.counts,
                     age.start = 0,age.end = 110)

#### Time trends ####

year.no <- length(year.range)%/%5
year.odd <- length(year.range)-year.no*5

categories <- c(rep(1,2),rep(2:18,each=5))

a5.ex.siler<- tapply(gomp$deltaa , categories , sum )
b5.ex.siler<- tapply(gomp$deltab , categories , sum )
c5.ex.siler<- tapply(gomp$deltac , categories , sum )
B5.ex.siler<- tapply(gomp$deltaB , categories , sum )
M5.ex.siler<- tapply(gomp$deltaM , categories , sum )
e5.ex.siler<- tapply(gomp$deltae , categories , sum )

data <- data.table(label = c(rep("infant mortlaiy initial",year.no+1),
                             rep("infant mortality decrease",year.no+1),
                             rep("Background",year.no+1),
                             rep("Compression",year.no+1),
                             rep("Shift",year.no+1)),
                   value = c(a5.ex.siler,
                             b5.ex.siler,
                             c5.ex.siler,
                             B5.ex.siler,
                             M5.ex.siler),
                   year = paste(c(year.range[year.odd],
                            seq(year.range[year.odd+5],
                                year.range[length(year.range)],5)))
)

data$label <- factor(data$label,
                     levels = c("infant mortlaiy initial",
                                "infant mortality decrease",
                                "Background",
                                "Compression",
                                "Shift"))

ggplot(data, aes(x=year,y=value))+
  geom_col(aes(fill=label))+
  stat_summary(fun = sum,geom = "point")+
  scale_fill_manual(values = c("pink","lightyellow","lightblue","navy","orange"))+
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

path_out <- paste("report/", cnty, "_Male,five_comp.pdf",sep="")

ggsave(path_out,
       width = 8,height = 6)


#### Last ten years ####

data2 <- data.table(label = c(rep("infant mortlaiy initial",9),
                              rep("infant mortality decrease",9),
                              rep("Background",9),
                              rep("Compression",9),
                              rep("Shift",9)),
                    value = c(gomp$deltaa[79:87],
                              gomp$deltab[79:87],
                              gomp$deltac[79:87],
                              gomp$deltaB[79:87],
                              gomp$deltaM[79:87]),
                    year = c(paste(seq(2010,2019,1),"-",
                                   seq(2011,2020,1)))
)

data2$label <- factor(data2$label,
                      levels = c("infant mortality decrease",
                                 "infant mortlaiy initial",
                                 "Background",
                                 "Compression",
                                 "Shift"))

ggplot(data2, aes(x=year,y=value))+
  geom_col(aes(fill=label))+
  stat_summary(fun = sum,geom = "point")+
  scale_fill_manual(values = c("pink","lightyellow","lightblue","navy","orange"))+
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

ggsave("report/USA_Male,five_comp_2011-2020.pdf",
       width = 8,height = 6)
