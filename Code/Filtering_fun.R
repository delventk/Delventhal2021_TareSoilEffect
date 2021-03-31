#Fungal data: Code to filter phyloseq object and perform transformations.

#load libraries
library(vegan)
library(phyloseq)
library(tidyverse)
library(magrittr)
library(microbiome)
library(goeveg)

#clear global environment
rm(list = ls())


##################################
##Create phyloseq objects for all samples
##################################

### Read in phyloseq object
(phy <- readRDS("./phy_object_fun/phy_fun.rds"))

##Read in meta data to replace old
meta <- read.csv(file="./phy_object_fun/Meta.csv")
meta$Sample.ID <- meta$Sample.ID%>% gsub(" ",".",.)
meta %<>% column_to_rownames("Sample.ID")
#add new meta data to phylsoeq
meta$Seed.Lot <- as.factor(meta$Seed.Lot)
sample_data(phy) <- sample_data(meta)


##subset data
phy %<>% subset_samples(., Replicates != "1") #these are samples that erroneously had the same sample code and need to be removed
phy %<>% subset_samples(., Sample.Type != "Control") #remove positive and negative PCR controls
phy %<>% subset_samples (., sample_names(phy) != "38.TS.VG.3") #remove outlier

#check sequencing depth
minDepth <- 5000
data.frame(SeqDepth=sort(sample_sums(phy)), Treatment=sample_data(phy)$Treatment) %>%
  mutate(cutoff=SeqDepth>minDepth) %>%
  ggplot(aes(x=Treatment, y=SeqDepth)) +
  geom_violin() +
  geom_point(aes(color=cutoff),position=position_jitter(width=0.1)) +
  theme_classic()
check <- data.frame(SeqDepth=sort(sample_sums(phy)), Group = sample_data(phy)$Group, ID = sample_data(phy)$Plant.ID)

#save the phyloseq object
saveRDS(phy,"./phy_object_fun/phy_fungi_prefilter.rds") #save object

##################################
##Filter taxa and transform
##################################
rm(list = ls())
phy <- readRDS("./phy_object_fun/phy_fungi_prefilter.rds")

###relativize and transform data
phy %<>% transform_sample_counts(function(x)x/sum(x)) #covert to relative abundance
##log transformation: generalized log eq = log(x +xmin)- log(xmin)
min(apply(otu_table(phy), 2 , FUN = function(x) {min(x[x > 0])})) #find min value
phy %<>% transform_sample_counts(function(x){log10(x + 2.171788e-05) - log10(2.171788e-05)})
otu<- otu_table(phy) %>% data.frame

##Check abundance curve
# GenLog <- racurve(otu, main= "Gen Log", ylog = FALSE, frequency = FALSE)
# relabun.genlog <- GenLog$rel.abund %>% data.frame
# tax <- tax_table(phy) %>% data.frame

saveRDS(phy,"./phy_object_fun/phy_fungi_allsamples.rds") #save object

##################################
##Create phy for rhizosphere samples only
##################################
rm(list = ls())

phy <- readRDS("./phy_object_fun/phy_fungi_allsamples.rds")
##filter data
phy %<>% subset_samples(Sample.Type == "Rhizosphere Soil")
phy %<>% filter_taxa(., function(x) {sum(x>0) > 0.03*nsamples(.)}, TRUE) %>% prune_taxa(taxa_sums(.) > 0,.) #filter to retain OTUs present in at least 3% of all samples

saveRDS(phy,"./phy_object_fun/phy_fungi_rhizo.rds")  #save object