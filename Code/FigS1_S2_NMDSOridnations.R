
library(phyloseq)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(patchwork)
library(vegan)
library(microbiome)

##########################
#NMDS all sample types
##########################


##########################
#Generate and save rds objects of the NMDS ordinations

#FUNGI
phy <- readRDS("./phy_object_fun/phy_fungi_allsamples.rds")
#distance and ordination object
phy.dist <- phy %>% phyloseq::distance("bray")
phy.ord <- metaMDS(phy.dist, 
                   k=2, 
                   trymax = 250, 
                   maxit = 999,
                   engine = "monoMDS")
saveRDS(phy.ord, "./output/RDS/fungi/phy.ord.all.rds")

#BACTERIA
rm(list = ls())
phy <- readRDS("./phy_object_bac/phy_bac_allsamples.rds")
phy.dist <- phy %>% phyloseq::distance("bray")
phy.ord <- metaMDS(phy.dist, k=2, 
                   trymax = 250, 
                   maxit = 999,
                   engine = "monoMDS")
saveRDS(phy.ord, "./output/RDS/bacteria/phy.ord.rds")


##########################
#Create NMDS ordinations/Figures


rm(list = ls())
theme_set(theme_bw())

phy_all_fun <- readRDS("./phy_object_fun/phy_fungi_allsamples.rds")
phy_all_bac <- readRDS("./phy_object_bac/phy_bac_allsamples.rds")

metafun <- meta(phy_all_fun)
metafun %<>% rename(Soil.Origin = Soil.Type)
metabac <- meta(phy_all_bac)
metabac %<>% rename(Soil.Origin = Soil.Type)

sample_data(phy_all_bac) <- sample_data(metabac)
sample_data(phy_all_fun)<- sample_data(metafun)

ord_all_fun <- readRDS("./Output/RDS/fungi/phy.ord.all.rds")
ord_all_bac <- readRDS("./Output/RDS/bacteria/phy.ord.rds")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pF <- plot_ordination(phy_all_fun, ord_all_fun, color= "Sample.Type", shape = "Soil.Origin") + 
  theme(legend.position = "none") + 
  geom_point(size=1.75) +
  scale_shape_manual(values=c(16, 15, 17 )) +
  scale_colour_manual(values= cbbPalette )

pB <-plot_ordination(phy_all_bac, ord_all_bac, color= "Sample.Type", shape = "Soil.Origin") + 
  geom_point(size=1.75) +
  scale_shape_manual(values=c(16, 15, 17 )) +
  scale_colour_manual(values=cbbPalette ) 

SampleTypes <-pF + pB + plot_annotation(tag_levels = 'A')
SampleTypes

ggsave("FigS1_NMDS_SampleTypes.jpeg", path = "./Output/Figures", width = 7, height = 3, units = "in")


########################
#Rhizosphere Samples only
########################

##########################
#Generate and save rds objects of the NMDS ordinations

rm(list = ls())

#FUNGI
phy <- readRDS("./phy_object_fun/phy_fungi_rhizo.rds")
phy.dist <- phy %>% phyloseq::distance("bray")
phy.ord <- metaMDS(phy.dist, 
                   k=2, 
                   trymax = 250, 
                   maxit = 999,
                   engine = "monoMDS")
stressplot(phy.ord)
saveRDS(phy.ord, "./output/RDS/fungi/phy.ord.rhizo.rds")

#BACTIER
rm(list = ls())

phy<- readRDS("./phy_object_bac/phy_bac_rhizo.rds")
phy.dist <- phy %>% phyloseq::distance("bray")
phy.ord <- metaMDS(phy.dist,
                   k=2, 
                   trymax = 250, 
                   maxit= 999, 
                   engine = "monoMDS")
saveRDS(phy.ord, "./output/RDS/bacteria/phy.ord.rhizo.rds")


##########################
#Create NMDS ordinations/Figures


#clear global environment
rm(list = ls())
(phyR_bac<- readRDS("./phy_object_bac/phy_bac_rhizo.rds"))
(phyR_fun <- readRDS("./phy_object_fun/phy_fungi_rhizo.rds"))

#rename the variable soil type to soil origin, for clarification purposes
metafun <- meta(phyR_fun)
metafun %<>% rename(Soil.Origin = Soil.Type)
metabac <- meta(phyR_bac)
metabac %<>% rename(Soil.Origin = Soil.Type)

sample_data(phyR_bac) <- sample_data(metabac)
sample_data(phyR_fun)<- sample_data(metafun)

ordR_bac <- readRDS("./Output/RDS/bacteria/phy.ord.rhizo.rds")
ordR_fun <- readRDS("./Output/RDS/fungi/phy.ord.rhizo.rds")

PRB <- plot_ordination(phyR_bac,ordR_bac,color= "Soil.Origin", shape = "Treatment") +
  #xlab("NMDS1(83.8%)") + ylab("NMDS2(37.4%)") + 
  geom_point(size=2)+
  scale_shape_manual(values=c(1, 16)) +
  scale_colour_manual(values= c("#E69F00","#0072B2"))

PRF <- plot_ordination(phyR_fun,ordR_fun,color= "Soil.Origin", shape = "Treatment") +
  #xlab("NMDS1(81%)") + ylab("NMDS2(30.2%)") + 
  geom_point(size=2) + 
  theme(legend.position = "none") +
  scale_shape_manual(values=c(1, 16)) +
  scale_colour_manual(values= c("#E69F00","#0072B2"))

plots <- PRF + PRB
plots + plot_annotation(tag_levels = 'A')

ggsave("FigS2_NMDS_Rhizosphere.jpeg", path = "./Output/Figures", width = 7, height = 3.5, units = "in")

