
###
### European Countries total
###

source("codes/functions.R")

library(data.table)
library(ggplot2)

cnty1 <- c("CHE","CHL","DNK","FIN","FRA",
           "DEU","IRL","JPN","NZL","NOR",
           "PRT","ESP","SWE","UK","USA")
cnty2 <- c("Switzerland","Chile","Denmark",
           "Finland","France","Germany",
           "Ireland","Japan","New Zealand",
           "Norway","Portugal","Spain",
           "Sweden","United Kingdom","United States")

Data <- c()

for(i in 1:length(cnty1)){
  cnty <- cnty1[i]
  
  col_name <- c("Year","Age","Male")
  
  Dx <- fread(paste("data/HMD/",cnty,".Deaths_1x1.txt",sep=""))[,..col_name]
  
  Dx <- Dx[Year>=2009,]
  
  names(Dx)[3] <- "V1"
  
  Dx <- matrix(data = Dx$V1,nrow = 111,ncol = length(unique(Dx$Year)),
               dimnames = list(unique(Dx$Age),unique(Dx$Year)))
  
  Ex <- fread(paste("data/HMD/",cnty,".Exposures_1x1.txt",sep=""))[,..col_name]
  
  Ex <- Ex[Year>=2009,]
  
  names(Ex)[3] <- "V1"
  
  Ex <- matrix(data = Ex$V1,nrow = 111,ncol = length(unique(Ex$Year)),
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
  
  n <- length(gomp$deltaa)
  
  data <- data.table(label = c(rep("infant mortlaiy initial",n),
                               rep("infant mortality decrease",n),
                               rep("Background",n),
                               rep("Compression",n),
                               rep("Shift",n)),
                     value = c(gomp$deltaa,
                               gomp$deltab,
                               gomp$deltac,
                               gomp$deltaB,
                               gomp$deltaM),
                     year = paste(seq(year.interval[1],year.interval[2]-1,1),
                                  "-",
                                  seq(year.interval[1]+1,year.interval[2],1))
  )
  
  data <- data[,Pop:=cnty2[i]]
  
  Data <- rbind(Data,data)
  
  data$label <- factor(data$label,
                       levels = c("infant mortlaiy initial",
                                  "infant mortality decrease",
                                  "Background",
                                  "Compression",
                                  "Shift"))
  
  ggplot(data, aes(x=year,y=value))+
    geom_col(aes(fill=label),width=0.8)+
    stat_summary(fun = sum,geom = "point")+
    scale_fill_manual(values = c("pink","lightyellow","lightblue",
                                 "navy","orange"))+
    theme_classic()+
    theme(panel.grid.major.x = element_line(colour = "grey"),
          axis.text.x = element_text(angle = 45,
                                     vjust = 1.2,
                                     hjust = 1.2),
          legend.position = "bottom")+
    labs(x="Year",y="Contribution",
         title = paste0("Changes in Life expectancy for the ",
                        cnty2[i],
                        "\n with Gompertz-Makeham Model,Male 2009-2021"))+
    guides(fill=guide_legend("Components"))
  
  
  path_out <- paste("report/", cnty, "_Male,five_comp.png",sep="")
  
  ggsave(path_out,
         width = 8,height = 6)
}

fwrite(Data,"data/EU-table_Male_five_comp.csv")
