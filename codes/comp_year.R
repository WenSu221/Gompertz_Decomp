
###
### Comparison by Year
###

library(ggplot2)
library(data.table)

#### Western European countries ####

covid <- c("2018 - 2019", "2019 - 2020", "2020 - 2021")

### Female ###

Data <- fread("data/EU-table_Female_five_comp.csv")

ggplot(Data,aes(x=Pop,y=value,fill=label))+
  geom_col()+
  facet_wrap(~year)+
  scale_fill_manual(values = c("lightblue","navy","pink",
                               "lightyellow","orange"))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1.2,
                                   hjust = 1.2),
        legend.position = "bottom")+
  labs(x="Populations",y="Contributions")
ggsave("report/comparison_among_EU_Female.png",width = 12,height = 6)

### Male ###

Data <- fread("data/EU-table_Male_five_comp.csv")

ggplot(Data[year%in%covid],aes(x=Pop,y=value,fill=label))+
  geom_col()+
  facet_wrap(~year)+
  scale_fill_manual(values = c("lightblue","navy","pink",
                               "lightyellow","orange"))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1.15,
                                   hjust = 1.2),
        text = element_text(face = "bold",size=12),
        legend.position = "bottom")+
  labs(x="Populations",y="Contributions")
ggsave("report/comparison_among_EU_Male_covid.png",width = 12,height = 6)

#### Eastern-european countries ####

### Female ###

Data <- fread("data/EEU-table_Female_five_comp.csv")

ggplot(Data,aes(x=Pop,y=value,fill=label))+
  geom_col()+
  facet_wrap(~year)+
  scale_fill_manual(values = c("lightblue","navy","pink",
                               "lightyellow","orange"))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1.15,
                                   hjust = 1.2),
        text = element_text(face = "bold",size=12),
        legend.position = "bottom")+
  labs(x="Populations",y="Contributions")
ggsave("report/comparison_among_EEU_Female.png",width = 12,height = 8)

### Male ###

Data <- fread("data/EEU-table_Male_five_comp.csv")

ggplot(Data,aes(x=Pop,y=value,fill=label))+
  geom_col()+
  facet_wrap(~year)+
  scale_fill_manual(values = c("lightblue","navy","pink",
                               "lightyellow","orange"))+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1.2,
                                   hjust = 1.2),
        legend.position = "bottom")+
  labs(x="Populations",y="Contributions")
ggsave("report/comparison_among_EEU_Male.png",width = 12,height = 8)


#### Life expectancy ####

Data <- fread("data/EEU-table_Female_ex.csv")
Data <- Data[,year:=1959:2014]
Data <- melt.data.table(Data,id.vars = "year")

ggplot(Data,aes(x=year,y=value))+
  geom_line()+
  facet_wrap(~variable)+
  theme_classic()+
  theme(panel.grid.major.x = element_line(colour = "grey"),
        legend.position = "bottom")+
  labs(x="Year",y="Life expectancy")
ggsave("report/comparison_among_EEU_Female_ex.png",width = 12,height = 8)
