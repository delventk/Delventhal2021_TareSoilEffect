#Figure S3 Differentially abundant rhizosphere OTUs identified to the genus level or higher associated with the two soil types for Fungi and Bacteria
#to get the csv files, first run code DiffAbudnace_bac.R and Diff_Abundance_fun.R

library(ggplot2)
library(magrittr)
library(patchwork)

#clear work space
rm(list = ls())

#read in objects
df <- read.csv("./Output/Tables/aldex_SOIL_effectsize_fun.csv")
db <- read.csv("./Output/Tables/aldex_SOIL_effectsize_bac.csv")

##Sort OTUs by effect size
#fungi
x <- tapply(df$effect, df$Genus, function(x) max(x))
x = sort(x, TRUE)
df$Genus = factor(as.character(df$Genus), levels=names(x))
#bacteria
x <- tapply(db$effect, db$Genus, function(x) max(x))
x = sort(x, TRUE)
db$Genus = factor(as.character(db$Genus), levels=names(x))


#############
#create plots

Fun2<- ggplot(df, aes(x=effect, y=Genus)) + geom_point(size=2) + 
  labs(x= "Effect Size", y = "Taxon") +
  theme_bw()

Bac2<-  
  ggplot(db, aes(x=effect, y=Genus)) + geom_point(size=2) + 
  labs(x= "Effect Size") +
  theme_bw() +
  theme(axis.title.y = element_blank())


boom <- Fun2 + Bac2
boom + plot_annotation(tag_levels = 'A')

ggsave("FigS3_SoilOrigin.jpeg", path = "./Output/Figures", width = 7, height = 5.5, units = "in")
