#Compile MiSeq library and sample data into single phyloseq object 
#Also add taxonomy and remove non-target sequences
 
library(tidyverse)
library(magrittr)
library(foreach)
library(dada2)
library(Biostrings)
library(phyloseq)


  #read in denoised sequence table
  seqTab <- readRDS("output/dada/seqTab.rds")
  #extract sequences from OTU table 
  seqs <- getSequences(seqTab) %>% DNAStringSet
  #rename sequences
  OTU.names <- paste0("OTU.",1:ncol(seqTab))
  colnames(seqTab) <- OTU.names
  names(seqs) <- OTU.names
  
  #convert to phyloseq object
  otuTab <- seqTab %>% as.data.frame() %>% otu_table(taxa_are_rows = F)
  phy <- phyloseq(otuTab, refseq(seqs))
  
  #collapse sequences with >= 99% similarity
  cluster <- function(phy.in,method="single",dissimilarity=0.01){
    require(DECIPHER)
    clusters <- DistanceMatrix(refseq(phy.in), includeTerminalGaps = T, processors=NULL) %>%
      IdClusters(method=method, cutoff=dissimilarity, processors=NULL) %>%
      rownames_to_column("OTUid")
    for(i in unique(clusters$cluster)){
      foo <- clusters$OTUid[which(clusters$cluster==i)] 
      if(length(foo)>1){phy.in %<>% merge_taxa(foo)}
    }
    return(phy.in)
  }  
  phy %<>% cluster  

  
  #Assign taxonomy with bacteria database
  taxa <- assignTaxonomy(refseq(phy),"data/taxonomy_db/silva_nr_v132_train_set.fa.gz", multithread = T, minBoot=50)
  rownames(taxa) %<>% names

  #prune to have only bacteria
 
  Bact <- taxa[,1]=="Bacteria" 
  Bact[is.na(Bact)] <- F
  
  taxaBac <- taxa[Bact,]
  
  Chloro <- taxaBac[,4]!="Chloroplast"
  Chloro[is.na(Chloro)] <- T
  
  taxaBacChlor <- taxaBac[Chloro,]
  
  Mito <- taxaBacChlor[,5]!= "Mitochondria"
  Mito[is.na(Mito)] <- T
 
  taxaBacChlorMito <- taxaBacChlor[Mito,]

  taxa <- taxaBacChlorMito

  ###############   
### OUTPUTs ###
###############

  #set final OTU names
  otu.names <- paste0("OTU.",1:ntaxa(phy))
  names(otu.names) <- taxa_names(phy)[order(taxa_sums(phy),decreasing = T)]
  taxa_names(phy) <- otu.names[taxa_names(phy)]
  
###
  #Read in sample meta data
  meta <- read.csv(file="./data/Meta.csv") %>% mutate(ID=paste(index_i7_revComp,index_i5,sep="-")) %>% column_to_rownames("ID")
  #add meta data to phylsoeq
  sample_data(phy) <- sample_data(meta)
  
  #add taxonomy to phyloseq
  #was below script needed??
  #rownames(taxa) <- otu.names
 
   tax_table(phy) <- tax_table(taxa)
  
  #Update sample names
  sample_names(phy) <- sample_data(phy)$Sample.ID %>% gsub(" ",".",.)

#Save phyloseq object
saveRDS(phy,"output/phy.rds")

#Save components for possible manual inspection
otu_table(phy) %>% write.csv("output/OTU.table.csv")
tax_table(phy) %>% write.csv("output/taxonomy.table.csv")
