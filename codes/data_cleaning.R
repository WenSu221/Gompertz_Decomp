
###
### Data Cleaning
###

library(data.table)
library(ungroup)

#### Deaths ####

death <- fread("data/deaths_native 2011-2020.txt",
              fill=T)
death <- death[Notes=="",][,list(`Single-Year Ages Code`,
                             `Gender Code`,`Year Code`,
                             Deaths)]
death <- death[`Year Code`>=2011,]
names(death) <- c("Age","Gender","Year","Deaths")
death <- dcast(death,Gender+Age~Year,value.var = "Deaths")

fwrite(death,"data/deaths_native_cleaned.txt")

#### Population #####

pop <- fread("data/pop_native 2011-2020.txt",
             fill = T)
pop <- pop[Notes=="",][,list(`Single-Year Ages Code`,
                             `Gender Code`,
                             `Yearly July 1st Estimates Code`,
                             Population)]
names(pop) <- c("Age","Gender","Year","Population")
pop$Age <- as.numeric(pop$Age)
pop$Age <- fifelse(is.na(pop$Age),85,pop$Age)

age <- 0:85
year <- unique(pop$Year)
gender <- unique(pop$Gender)

weights <- c()

for (x in 1:length(year)){
  for (y in 1:length(gender)){
    
    year_sub <- year[x]
    gender_sub <- gender[y]
    
    pop_sub <- pop[Year==year_sub&Gender==gender_sub]$Population
    pop_sub <- pclm(age,pop_sub,nlast = 16)$fitted
    
    table_sub <- data.table(
      year = rep(year_sub,101),
      gender = rep(gender_sub,101),
      age = 0:100,
      count = pop_sub
    )
    weights <- rbind(weights,table_sub)
  }
}

pop_smooth <- dcast(weights,gender+age~year,value.var = "count")

fwrite(pop_smooth,"data/pop_native_cleaned.txt")
