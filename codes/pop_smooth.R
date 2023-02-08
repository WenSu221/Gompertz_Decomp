
###
### Data Cleaning 
###

library(data.table)
library(ungroup)

table <- read.table("data/label_table.txt",sep=",",
                    header = T)

#### CDC population data ####

popA <- fread("data/US pop by states 2016-2020 Census.txt",
             fill=T)
popA <- popA[-c((nrow(popA)-35):nrow(popA))]
popA <- popA[table,on=list(States=fln)]

popB <- fread("data/US pop by states 2011-2015 Census.txt",
             fill=T)
popB <- popB[-c((nrow(popB)-35):nrow(popB))]
popB <- popB[table,on=list(States=fln)]

pop <- rbind(popA,popB)

#### Dis-aggregate ####

states <- unique(pop$abv)
year <- unique(pop$`Yearly July 1st Estimates`)
gender <- unique(pop$Gender)

age <- 0:85

weights <- c()

for (x in 1:length(states)){
  for (y in 1:length(year)){
    for (z in 1:length(gender)){
      states_sub <- states[x]
      year_sub <- year[y]
      gender_sub <- gender[z]
      pop_sub <- pop[abv %in% states_sub&
                       `Yearly July 1st Estimates` %in% year_sub&
                       Gender %in% gender_sub]$Population
      pop_sub <- pclm(age,pop_sub,nlast = 26)$fitted
      table_sub <- data.table(
        states = rep(states_sub,111),
        year = rep(year_sub,111),
        gender = rep(gender_sub,111),
        age = 0:110,
        count = pop_sub
      )
      weights <- rbind(weights,table_sub)
    }
  }
}

weights2 <- dcast(weights,year+gender+age~states,value.var = "count")

setorder(weights2,year)

fwrite(weights2,"data/pop_smooth_2011-2020.csv",sep=",")
