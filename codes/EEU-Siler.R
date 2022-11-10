
###
### East European Countries
###

source("codes/functions.R")

library(data.table)
library(ggplot2)

cnty1 <- c("BLR","BGR","CZE","EST",
           "HUN","LVA","LTU","POL",
           "SVK","RUS")
cnty2 <- c("Belarus","Bulgaria","Czechia",
           "Estonia","Hungary","Latvia",
           "Lithuania","Poland","Slovakia",
           "Russia")

Data <- c()
ex <- c()

for(i in 1:(length(cnty1))){
  cnty <- cnty1[i]
  
  col_name <- c("Year","Age","Female")
  
  Dx <- fread(paste("data/HMD/",cnty,".Deaths_1x1.txt",sep=""))[,..col_name]
  
  Dx <- Dx[Year>=1959&Year<=2014,]
  
  names(Dx)[3] <- "V1"
  
  Dx <- matrix(data = Dx$V1,nrow = 111,ncol = length(unique(Dx$Year)),
               dimnames = list(unique(Dx$Age),unique(Dx$Year)))
  
  Ex <- fread(paste("data/HMD/",cnty,".Exposures_1x1.txt",sep=""))[,..col_name]
  
  Ex <- Ex[Year>=1959&Year<=2014,]
  
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
  
  categories <- c(rep(1:11,each=5))
  
  a5.ex.siler<- tapply(gomp$deltaa , categories , sum )
  b5.ex.siler<- tapply(gomp$deltab , categories , sum )
  c5.ex.siler<- tapply(gomp$deltac , categories , sum )
  B5.ex.siler<- tapply(gomp$deltaB , categories , sum )
  M5.ex.siler<- tapply(gomp$deltaM , categories , sum )
  e5.ex.siler<- tapply(gomp$deltae , categories , sum )
  
  data <- data.table(label = c(rep("infant mortlaiy initial",11),
                               rep("infant mortality decrease",11),
                               rep("Background",11),
                               rep("Compression",11),
                               rep("Shift",11)),
                     value = c(a5.ex.siler,
                               b5.ex.siler,
                               c5.ex.siler,
                               B5.ex.siler,
                               M5.ex.siler),
                     year = paste(seq(1959,2009,5),"-",
                                  seq(1964,2014,5))
  )
  
  data <- data[,Pop:=cnty2[i]]
  
  Data <- rbind(Data,data)
  
  ex <- cbind(ex, gomp$ex)
  colnames(ex)[i] <- cnty
  
  data$label <- factor(data$label,
                       levels = c("infant mortlaiy initial",
                                  "infant mortality decrease",
                                  "Background",
                                  "Compression",
                                  "Shift"))

  ggplot(data, aes(x=year,y=value))+
    geom_col(aes(fill=label))+
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
                        "\n with Gompertz-Makeham Model,Female 1959-2019"))+
    guides(fill=guide_legend("Components"))


  path_out <- paste("report/", cnty, "_Female,five_comp.png",sep="")

  ggsave(path_out,
         width = 8,height = 6)
}

fwrite(Data,"data/EEU-table_Female_five_comp.csv")
fwrite(as.data.table(ex),"data/EEU-table_Female_ex.csv")