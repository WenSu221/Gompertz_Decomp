
####
####  US table by ethnicity
####


library(data.table)
library(ggplot2)

pop <- c("blacks","native","nonhis white","total")


table <- c()

for (i in 1:4){
  table_sub <- fread(paste("report/table", pop[i],"_Male 2011-2020.txt"))
  table_sub <- table_sub[,pop := pop[i]]
  table <- rbind(table,table_sub)
}

table2 <- fread("report/table total_HMD _Male 2011-2020.txt")

comp <- table[pop=="total",list(sum(value)),by=list(year2,label)]

comp <- comp[,`:=`(CDC=V1,V1=NULL)]

comp2 <- table2[,list(sum(value)),by=list(year2,label)]
comp2 <- comp2[,`:=`(HMD=V1,V1=NULL)]

comp <- comp[comp2,on=.(year2=year2,label=label)]

# comp <- comp[,year2 := fcase(year %in% c("2011-2012",
#                                          "2012-2013",
#                                          "2013-2014",
#                                          "2014-2015"),
#                              year2 = "2011-2015",
#                              year %in% c("2015-2016",
#                                          "2016-2017",
#                                          "2017-2018",
#                                          "2018-2019"),
#                              year2 = "2015-2019",
#                              year == "2019-2020",
#                              year2 = "2019-2020")]
# comp <- comp[,list(CDC = sum(CDC),HMD = sum(HMD)),
#              by = list(year2,label)]

comp[,diff := CDC-HMD]

comp[,list(CDC=sum(CDC),HMD=sum(HMD),
           diff=sum(diff)),by=.(year2)]

# table_pop <- table[,year2 := fcase(year %in% c("2011-2012",
#                                                "2012-2013",
#                                                "2013-2014",
#                                                "2014-2015"),
#                                    year2 = "2011-2015",
#                                    year %in% c("2015-2016",
#                                                "2016-2017",
#                                                "2017-2018",
#                                                "2018-2019"),
#                                    year2 = "2015-2019",
#                                    year == "2019-2020",
#                                    year2 = "2019-2020")]
# 
# table_pop <- table_pop[,list(Value=sum(value)),
#                        by=list(label,year2,pop)]


table_pop <- dcast(table,label+year2~pop,value.var = "value")

# table_pop2 <- table_pop[,`:=`(blacks=blacks/total*-1,native=native/total*-1,
#                 `nonhis white`=`nonhis white`/total*-1)]
# 
# table_pop2 <- table_pop2[table2,on=.(year2=year2,label=label)]
# 
# table_pop2[,`:=`(blacks=blacks*value*-1,native=native*value*-1,
#                 `nonhis white`=`nonhis white`*value*-1)]
