#Bacterial data: Code to filter phyloseq object and perform transformations.

#load in libraries
library(vegan)
library(phyloseq)
library(tidyverse)
library(magrittr)
library(goeveg)

#clear global environment
rm(list = ls())

##################################
##Create phyloseq objects for all samples
##################################

#read in phyloseq object
(phy <- readRDS("./phy_object_bac/phy_bac.rds"))

#Read in sample meta data
meta <- read.csv(file="./phy_object_bac/Meta.csv")
#make changes to have it match phy object
meta$Sample.ID <- meta$Sample.ID%>% gsub(" ",".",.)
#add new meta data to phylsoeq
meta %<>% column_to_rownames("Sample.ID")
meta$Seed.Lot <- as.factor(meta$Seed.Lot)
sample_data(phy) <- sample_data(meta)
#Look at meta data, sanity check that old matches
temp.meta<- sample_data(phy) %>% data.frame()

#remove samples and filter
phy %<>% subset_samples(., Replicates != "1") #these are samples that erroneously had the same sample code and need to be removed
phy %<>% subset_samples(., Sample.Type != "Control") #remove positive and negative controls
phy %<>% subset_samples(., sample_names(phy) != "38.TS.VG.3") #remove outlier
phy

#check sequencing depth 
minDepth <- 5000
data.frame(SeqDepth=sort(sample_sums(phy)), Treatment=sample_data(phy)$Treatment) %>%
  mutate(cutoff=SeqDepth>minDepth) %>%
  ggplot(aes(x=Treatment, y=SeqDepth)) +
  geom_violin() +
  geom_point(aes(color=cutoff),position=position_jitter(width=0.1)) +
  theme_classic()
check <- data.frame(SeqDepth=sort(sample_sums(phy)), Group = sample_data(phy)$Group, ID = sample_data(phy)$Plant.ID)

#' ### Remove samples below sequencing depth cutoff
#(phy %<>% prune_samples(sample_sums(.)>minDepth,.))


saveRDS(phy, "./phy_object_bac/phy_bac_prefilter.rds")

###############################################
## Transform and filter taxa
###############################################

phy<- readRDS("./phy_object_bac/phy_bac_prefilter.rds")

#relativize by rows, converting to relative abundances (to account for differences in sequencing depth)
phy %<>% transform_sample_counts(function(x)x/sum(x))
#log transformation: generalized log eq. = log(x +xmin)- log(xmin)
min(apply(otu_table(phy), 2 , FUN = function(x) {min(x[x > 0])})) #find min value
phy <- transform_sample_counts(phy, function(x){log10(x + 2.720792e-05) - log10(2.720792e-05)})

saveRDS(phy,"./phy_object_bac/phy_bac_allsamples.rds")

##################################
##Create phy object for rhizosphere samples only
##################################
phy<- readRDS("./phy_object_bac/phy_bac_allsamples.rds")

phy %<>% subset_samples(Sample.Type == "Rhizosphere Soil")
phy %<>% filter_taxa(., function(x) {sum(x>0) > 0.03*nsamples(.)}, TRUE) %>% prune_taxa(taxa_sums(.) > 0,.) #filter to retain OTUs present in at least 3% of all samples

saveRDS(phy,"./phy_object_bac/phy_bac_rhizo.rds")

