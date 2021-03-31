#Heat map of the top 25 most abundant genera across sample types, creates Figure 3

#libraries
library(microbiome)
library(phyloseq)
library(tidyverse)
library(magrittr)
library(vegan)
library(reshape2)
library(dplyr)
library(ggplot2)

rm(list = ls())

phy <- readRDS("./phy_object_fun/phy_fungi_prefilter.rds") #load phyloseq
phy %<>% transform_sample_counts(function(x)x/sum(x)) #convert to relative abundance  
phy %<>% subset_samples(., Sample.Type != "Bulk Soil") %>% prune_taxa(taxa_sums(.) > 0,.) #remove samples of bulk soil mixture

metareplace <- read.csv("./phy_object_fun/MetaHeatMap_M3.csv", row.names = 1) #replace meta file- important to porperly sort of the data for this analysis
sample_data(phy) <- sample_data(metareplace)

meta<- meta(phy)
meta$Sample.Group %>% tapply(.,.,length) #get amount of samples within each grroup, which is used later

#########################################

phy25 <-aggregate_top_taxa(phy, top = 28, 'Genus') #set to 27 because "unknown" and "other" are removed later
meta25<- meta(phy25) #make meta file
phy25 %<>% merge_samples(., meta25$Sample.Group) #merge samples by sample grouping variable
otucheck <- otu_table(phy25) %>% data.frame

#remove unidentified 
otucheck$Unknown <- NULL
otucheck$Other <- NULL
otucheck$unidentified <- NULL
otu_table(phy25)<- otu_table(otucheck, taxa_are_rows = FALSE)

#create a melted dataframe that can be added onto and used for building the heat map
dfm <- neat(abundances(phy25), method = "NMDS", distance = "bray")
dfm %<>% melt()
colnames(dfm)<- c("Taxa", "Sample.Group", "value")


#make new column that has the total number of samples for each group so that you can get an average later. use numbers from meta$Sample.Group %>% tapply(.,.,length) used earlier
dfm %<>% mutate(av.total = case_when(
  Sample.Group == "TS" ~ 10,
  Sample.Group == "SS AG" ~ 76,
  Sample.Group == "SS VG" ~ 74, 
  Sample.Group == "TS AG" ~ 74,
  Sample.Group == "TS VG" ~ 78,
  Sample.Group == "VG" ~ 1,
  Sample.Group == "AG" ~ 1,))

#new column that is value/ av.total
dfm %<>% mutate(Avg.Abundance = value/av.total)

#create new column of sample type
dfm %<>% mutate(Sample.Type = case_when(
  Sample.Group == "TS" ~ "Tare Soils",
  Sample.Group == "SS AG" ~ "Rhizosphere Soils",
  Sample.Group == "SS VG" ~ "Rhizosphere Soils", 
  Sample.Group == "TS AG" ~ "Rhizosphere Soils",
  Sample.Group == "TS VG" ~ "Rhizosphere Soils",
  Sample.Group == "VG" ~ "Inoculant Soils",
  Sample.Group == "AG" ~ "Inoculant Soils",))

dfm$Sqrt.Abundance <- sqrt(dfm$Avg.Abundance) #get squareroot average abundance

#order the levels
dfm$Sample.Group <- as.character(dfm$Sample.Group)
dfm$Sample.Group <- factor(dfm$Sample.Group, levels = c("AG", "VG", "TS AG", "SS AG", "TS VG", "SS VG", "TS"))


#plot
ggplot(data = dfm, mapping = aes(x = Sample.Group,
                                 y =  Taxa,
                                 fill = Sqrt.Abundance)) +
  geom_tile() +
  scale_fill_gradient(name = "Sqrt. Abundance",
                      low = "#FFFFFF",
                      high = "#012345") +
  facet_grid(~ Sample.Type, switch = "x", scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.placement = "outside",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"),
        legend.position="top")


#save
#ggsave("Fig3_HeatMap.jpeg", path = "./Output/Figures", width = 7, height = 6, units = "in" )

