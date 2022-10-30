
###
### US by race & Ethnicity
### 

source("codes/functions.R")

library(data.table)
library(RColorBrewer)
library(ggplot2)

label <- c("blacks","native","nonhis white","total")

death <- fread(paste("data/deaths_", label[4], 
                     "_cleaned.txt",sep=""))
pop_smooth <- fread(paste("data/pop_", label[4], 
                          "_cleaned.txt",sep=""))

Dx <- death[Gender=="M",][,3:12]
Dx <- Dx[-c(1:30),]
Ex <- pop_smooth[gender=="M",][,3:12]
Ex <- Ex[-c(1:30)]

# Dx <- matrix(Data1$Deaths,nrow = 55, byrow = T,
#               dimnames = list(names(30:84),c(2001:2020)))
# Ex <- matrix(Data1$Population,nrow = 55, byrow = T,
#               dimnames = list(names(30:84),c(2001:2020)))

age.interval <- c(30,100)
age <- (age.interval[1]:age.interval[2])
year.interval<-c(as.numeric(colnames(Ex)[1]),
                 as.numeric(colnames(Ex)[ncol(Ex)]))

death.counts <- as.matrix(Dx)
expo.counts <- as.matrix(Ex)

gomp <- Decomp.Makeham(death.counts,expo.counts,age.start = 30,age.end = 100)

categories <- rep(1:9,each=1)

c5.ex.make<- tapply(gomp$deltac  , categories , sum )
B5.ex.make<- tapply(gomp$deltaB , categories , sum )
M5.ex.make<- tapply(gomp$deltaM , categories , sum )
e5.ex.make<- tapply(gomp$deltae , categories , sum )

data <- data.table(
  label = c(rep("Background",length(e5.ex.make)),
            rep("Compression",length(e5.ex.make)),
            rep("Shift",length(e5.ex.make))),
  value = c(c5.ex.make,
            B5.ex.make,
            M5.ex.make),
  year = paste(
    seq(2011,2019,1),"-",seq(2012,2020,1),sep = ""
  )
)

data <- data[,year2 := fcase(year %in% c("2011-2012",
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

data <- data[,list(value=sum(value)),
               by = list(year2,label)]

data <- data[!is.na(year2),]

data$label <- factor(data$label,levels = c("Background","Compression","Shift"))

ggplot(data, aes(x=year2,y=value))+
  geom_col(aes(fill=label))+
  stat_summary(fun = sum,geom = "point")+
  scale_fill_manual(values = c("lightblue","navy","orange"))+
  scale_y_continuous(limits = c(-5,2))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 45,vjust = 1.2,hjust = 1.2))+
  guides(fill=guide_legend("Components"))+
  labs(x="Year",y="Contribution",
       title = paste("Changes in Life expectancy for the US with Gompertz-Makeham Model,",
               label[4], "\n Male 2011-2020", sep=""))

ggsave(paste("report/USA_Male_", label[4], ",three_comp_2011-2020.png",sep=""),
       width = 8,height = 6)

fwrite(data[,value:=round(value,3)], paste("report/table", label[4],"_Male 2011-2020.txt"))
