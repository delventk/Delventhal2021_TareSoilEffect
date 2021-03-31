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
  
  #Extract ITS region with ITSx
  #Create scratch directory
  dir.create("output/scratch")
  #write temp file
  writeXStringSet(seqs,file="output/scratch/tmp.fasta",width=1000)
  # ITSx parameters (these will depend on your data)
  ITSx.flags <- paste("-i output/scratch/tmp.fasta",
                        "-t 'all'",
                        "--preserve T",
                        "--complement T",
                        "--summary T",
                        "-o output/scratch/ITSx",
                        "--only_full T",
                        "-E 1e-2")
  system2("ITSx", args = ITSx.flags)
    
  #read in ITS2 sequences from ITSx output and filter out very short reads
  seqs.ITS2 <- readDNAStringSet("output/scratch/ITSx.ITS2.fasta") %>% .[.@ranges@width > 50]
  seqTab %<>% .[,names(seqs.ITS2)] 
  colnames(seqTab) <- seqs.ITS2 %>% as.character 
  seqTab %<>% collapseNoMismatch()
  
  #convert to phyloseq object
  seqs.ITS2.collapsed <- getSequences(seqTab) %>% DNAStringSet
  colnames(seqTab) %<>% openssl:::md5()
  names(seqs.ITS2.collapsed) <- seqs.ITS2.collapsed %>% as.character %>% openssl:::md5()
  otuTab <- seqTab %>% as.data.frame() %>% otu_table(taxa_are_rows = F)
  phy <- phyloseq(otuTab, refseq(seqs.ITS2.collapsed))
  
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
  
  #Identify and remove non-fungal reads with large UNITE database (all eukaryotes, including singletons)
  taxa.allEuk <- assignTaxonomy(refseq(phy),"data/taxonomy_db/sh_general_release_dynamic_s_all_02.02.2019.fasta", multithread = T, minBoot=50)
  rownames(taxa.allEuk) %<>% names
  #preliminary investigation showed that non-fungal Kingdom-level assignments without phylum-level assignments were likely to be fungal
  putativeFungi <- (taxa.allEuk[,1]=="k__Fungi" | is.na(taxa.allEuk[,2])) %>% .[.] %>% names
  phy %<>% prune_taxa(putativeFungi,.)
  
  #Re-assign taxonomy with fungi-only database
  taxa <- assignTaxonomy(refseq(phy),"data/taxonomy_db/sh_general_release_dynamic_02.02.2019.fasta", multithread = T, minBoot=50)
  rownames(taxa) %<>% names

###############
### OUTPUTs ###
###############

  #set final OTU names
  otu.names <- paste0("OTU.",1:ntaxa(phy))
  names(otu.names) <- taxa_names(phy)[order(taxa_sums(phy),decreasing = T)]
  taxa_names(phy) <- otu.names[taxa_names(phy)]

  #Read in sample meta data
  meta <- read.csv(file="./data/Meta.csv") %>% mutate(ID=paste(index_i7_revComp,index_i5,sep="-")) %>% column_to_rownames("ID")
  #add meta data to phylsoeq
  sample_data(phy) <- sample_data(meta)
  #add taxonomy to phyloseq
  rownames(taxa) <- otu.names
  tax_table(phy) <- tax_table(taxa)
  
  #Update sample names
  sample_names(phy) <- sample_data(phy)$Sample.ID %>% gsub(" ",".",.)

#Save phyloseq object
saveRDS(phy,"output/compiled/phy.rds")

#Save components for possible manual inspection
otu_table(phy) %>% write.csv("output/compiled/OTU.table.csv")
tax_table(phy) %>% write.csv("output/compiled/taxonomy.table.csv")
