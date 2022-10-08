
###
### US by race & Ethnicity
### 

source("codes/functions.R")

library(data.table)
library(RColorBrewer)
library(ggplot2)

Data <- fread("data/Blacks.txt",fill=T)
Data$Population <- as.numeric(Data$Population)

Data1 <- Data[Gender=="Male"&
                `Single-Year Ages Code`<85&
                `Single-Year Ages Code`>=30,]

Data1 <- Data1[,list(`Single-Year Ages Code`,
                   `Year Code`,Deaths,Population)]

Dx <- matrix(Data1$Deaths,nrow = 55, byrow = T,
              dimnames = list(names(30:84),c(1999:2020)))
Ex <- matrix(Data1$Population,nrow = 55, byrow = T,
              dimnames = list(names(30:84),c(1999:2020)))

age.interval <- c(30,84)
age <- (age.interval[1]:age.interval[2])
year.interval<-c(as.numeric(colnames(Ex)[1]),
                 as.numeric(colnames(Ex)[ncol(Ex)]))

death.counts <- Dx
expo.counts <- Ex

gomp <- Decomp.Gompertz(death.counts,expo.counts,age.start = 30,age.end = 110)

categories <- c(rep(1,1),rep(2:5,each=5))

B5.ex.gomp<- tapply(gomp$deltaB , categories , sum )
M5.ex.gomp<- tapply(gomp$deltaM , categories , sum )
e5.ex.gomp<- tapply(gomp$deltae , categories , sum )

data <- data.table(label = c(rep("compression",5),rep("shift",5)),
                   value = c(B5.ex.gomp,M5.ex.gomp),
                   year = rep(c("1999-2000",
                                paste(seq(2001,2016,5),"-",
                                      seq(2005,2020,5)))
                   )
)

data$label <- factor(data$label,levels = c("shift","compression"))

ggplot(data, aes(x=year,y=value,fill=label))+
  geom_col()+
  scale_fill_manual(values = c("orange","navy"))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 45,vjust = 1.2,hjust = 1.2))+
  guides(fill=guide_legend("Components"))
