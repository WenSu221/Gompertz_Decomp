
###
### Gompertz contributions to life expectancy changes
### with Cohort
###

source("codes/functions.R")

library(data.table)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

Ex <- fread(paste("data/HMD/","DNK.cExposures_1x1.txt",sep=""))[,list(Year,Age,Male)]

Mx <- fread(paste("data/HMD/","DNK.cMx_1x1.txt",sep=""))[,list(Year,Age,Male)]

Ex$Male <- as.numeric(Ex$Male)
Ex$Male <- fifelse(is.na(Ex$Male),0,Ex$Male)

Mx$Male <- as.numeric(Mx$Male)
Mx$Male <- fifelse(is.na(Mx$Male),0,Mx$Male)

Ex <- matrix(data = Ex$Male,nrow = 111,ncol = length(unique(Ex$Year)),
             dimnames = list(unique(Ex$Age),unique(Ex$Year)))

Mx <- matrix(data = Mx$Male,nrow = 111,ncol = length(unique(Mx$Year)),
             dimnames = list(unique(Mx$Age),unique(Mx$Year)))

age.interval <- c(0,110)

age <- (age.interval[1]:age.interval[2])

year.interval<- range(fread(paste("data/HMD/","DNK.E0coh.txt",sep=""))[,list(Year)])

death.counts <- as.matrix(Mx[(age.interval[1]+1):(age.interval[2]+1), 
                   as.character(c(year.interval[1]:year.interval[2]))])

expo.counts <- as.matrix(Ex[(age.interval[1]+1):(age.interval[2]+1), 
                  as.character(c(year.interval[1]:year.interval[2]))])

death.counts <- death.counts*expo.counts

colnames(death.counts) <- c(year.interval[1]:year.interval[2])
colnames(expo.counts) <- c(year.interval[1]:year.interval[2])

fit <- Decomp.Siler(death.counts,expo.counts,
                     age.start = 0,age.end = 110)

n <- length(fit$deltaa)

table <- data.table(label = c(rep("infant mortlaiy initial",n),
                             rep("infant mortality decrease",n),
                             rep("Background",n),
                             rep("Compression",n),
                             rep("Shift",n)),
                    value = c(fit$deltaa,
                             fit$deltab,
                             fit$deltac,
                             fit$deltaB,
                             fit$deltaM),
                    year = rep((year.interval[1]+1):year.interval[2],5)
                    )

table$label <- factor(table$label,
                     levels = c("infant mortlaiy initial",
                                "infant mortality decrease",
                                "Background",
                                "Compression",
                                "Shift"))

table <- table %>% group_by(year = cut_width(year,5))

setDT(table)
table <- table[,list(value=sum(value)),by=.(label,year)]

ggplot(table, aes(x=year,y=value))+
  geom_col(aes(fill=label),width=0.8)+
  stat_summary(fun = sum,geom = "point")+
  scale_fill_manual(values = c("pink","lightblue","orange",
                               "navy","purple"))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1.2,
                                   hjust = 1.2),
        legend.position = "bottom")+
  labs(x="Year",y="Contribution",
       title = paste0("Changes in Cohort Life expectancy",
                      "\n with Siler Model"))+
  guides(fill=guide_legend("Components"))
