#Distance based redundancy analysis of rhizosphere samples constrained by sterilization treatment and bulk soil origin.

library(tidyverse)
library(magrittr)
library(phyloseq)
library(vegan)
library(ggthemes)
library(RColorBrewer)
library(dplyr)

#Set theme for figures and create color blind palette
theme_set(theme_bw())
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#################################
## FUNGI

# Load in data
fungi <- readRDS("./phy_object_fun/phy_fungi_rhizo.rds")

#fungi %<>% subset_samples(Soil.Type == "AG") %>% prune_taxa(taxa_sums(.) > 0,.)

meta_fun <- sample_data(fungi) %>% data.frame() #create meta file
meta_fun$Treatment %<>% factor()
meta_fun$Soil.Type %<>% factor()
meta_fun %<>% rename(
  Soil.Origin = Soil.Type)

# Get distance matrix
fungi.dist <- fungi %>% phyloseq::distance("bray") 

# db-RDA constrained by Site
dbrdafun <- dbrda(fungi.dist ~ Treatment + Soil.Origin, data = meta_fun)

scoresdbfun <- scores(dbrdafun) %>% data.frame #look at scores
db.sumfun <- summary(dbrdafun) #EXTRACT SUMMARY OF DB-RDA TO EXTRACT AXES DATA FROM. run summary(dbrda for full details)

dbscoresfun <- as.data.frame(scores(dbrdafun, display = "sites")) %>% #EXTRACT DB-RDA SCORES TO GRAPH
  cbind(meta_fun)

#################################
## BACTERIA

# Load in data
bac <- readRDS("./phy_object_bac/phy_bac_rhizo.rds")

meta_bac <- sample_data(bac) %>% data.frame()
meta_bac$Seed.Lot %<>% factor()
meta_bac$Treatment %<>% factor()
meta_bac$Soil.Type %<>% factor()
meta_bac %<>% rename(
  Soil.Origin = Soil.Type)

# Get distance matrix
bac.dist <- bac %>% phyloseq::distance("bray") 

# db-RDA constrained by Site
dbrdabac <- dbrda(bac.dist ~ Treatment + Soil.Origin, data = meta_bac)

scoresdbbac <- scores(dbrdabac) %>% data.frame

dbscoresbac <- as.data.frame(scores(dbrdabac, display = "sites")) %>% #EXTRACT DB-RDA SCORES TO GRAPH
  cbind(meta_bac)
db.sumbac <- summary(dbrdabac) #EXTRACT SUMMARY OF DB-RDA TO EXTRACT AXES DATA FROM

#################################
## CREATE PLOTS

FUN <- 
ggplot(dbscoresfun) +
  geom_point(mapping = aes(x = dbRDA1, y = dbRDA2, color = Treatment, shape = Soil.Origin), size = 2) +
  scale_colour_manual(values= c("#E69F00","#0072B2")) + 
  #scale_shape_manual(values=c(1, 16)) +
  xlab(paste("db-RDA1 (", round(100*(db.sumfun$cont$importance[2, "dbRDA1"]), digits = 2), "%)" )) + #
  ylab(paste("db-RDA2 (", round(100*(db.sumfun$cont$importance[2, "MDS2"]), digits = 2), "%)" )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  #theme(legend.position = "none") +
  xlim(-1.5, 1.5) +
  ylim(-2.5, 2)
 
  
BAC<- 
  ggplot(dbscoresbac, aes(x = dbRDA1, y = dbRDA2, color = Treatment, shape = Soil.Origin)) +
  geom_point(size= 2) +
  scale_colour_manual(values= c("#E69F00","#0072B2")) +
  xlab(paste("db-RDA1 (", round(100*(db.sumbac$cont$importance[2, "dbRDA1"]), digits = 2), "%)" )) +
  ylab(paste("db-RDA2 (", round(100*(db.sumbac$cont$importance[2, "dbRDA2"]), digits = 2), "%)" )) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  xlim(-2, 2) + 
  ylim(-5.5, 5)

library(patchwork)
plot <- FUN + BAC 

#to make it on top of the other with the legend centered
plot2 <- plot + plot_annotation(tag_levels = 'A')
plot3<- plot2 & theme(legend.position = "right")
plot3 + plot_layout(guides = "collect") 


#ggsave("dbRDA.jpeg", path = "./Output/Figures", width = 7, height = 3.5, units = "in")
