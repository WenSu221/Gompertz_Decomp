
library(data.table)

DataDir <- "data/States"
File_List <- list.files(DataDir)

Data <- c()
selected_year <- c(2011:2020)

for(x in 1:51){
  name <- paste("data/States/",
                File_List[x],"/",
                File_List[x],"_mltper_1x1.txt",
                sep="")
  data <- fread(name,header = T,skip=1)
  data <- data[Year %in% selected_year,
               list(rep(File_List[x],
                        111*length(selected_year)),
                    Year,Age,mx)]
  Data <- rbind(Data,data)
}

names(Data)[1] <- "states"

Data2 <- dcast(Data,Year+Age~states,value.var = "mx",
               fun.aggregate = sum)
Data2$Age <- as.numeric(Data2$Age)
Data2$Age <- fifelse(is.na(Data2$Age),110,Data2$Age)

setorder(Data2,Year,Age)

fwrite(Data2,"data/mx_male.txt")              
