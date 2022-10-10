
###
### US by race & Ethnicity
### 

source("codes/functions.R")

library(data.table)
library(RColorBrewer)
library(ggplot2)

death <- fread("data/deaths_nonhis white_cleaned.txt")
pop_smooth <- fread("data/pop_nonhis white_cleaned.txt")

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
  label = c(rep("background",length(e5.ex.make)),
            rep("compression",length(e5.ex.make)),
            rep("shift",length(e5.ex.make))),
  value = c(c5.ex.make,
            B5.ex.make,
            M5.ex.make),
  year = paste(
    seq(2011,2019,1),"-",seq(2012,2020,1),sep = ""
  )
)

data$label <- factor(data$label,levels = c("shift","compression","background"))

ggplot(data, aes(x=year,y=value))+
  geom_col(aes(fill=label))+
  stat_summary(fun = sum,geom = "point")+
  scale_fill_manual(values = c("orange","navy","lightblue"))+
  scale_y_continuous(limits = c(-5,2))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 45,vjust = 1.2,hjust = 1.2))+
  guides(fill=guide_legend("Components"))+
  labs(x="Year",y="Contribution",
       title = "Changes in Life expectancy for the US with Gompertz-Makeham Model,
        Native Male 2011-2020")

ggsave("report/USA_Male_native,three_comp_2011-2020.pdf",
       width = 8,height = 6)
