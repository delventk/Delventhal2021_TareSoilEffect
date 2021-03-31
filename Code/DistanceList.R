#creates the .Rda objects (lists of distances between tare soil samples and rhizosphere samples, measured as a fucniton of community distance) to be used in permutation analysis (r code: Fig2_PermutationTest.R)

library(vegan)
library(phyloseq)
library(tidyverse)
library(magrittr)
library(spaa)
library(splitstackshape)

############################
#FUNGI
############################

##read in object

phy <- readRDS("./phy_object_fun/phy_fungi_allsamples.rds")
phy <- subset_samples(phy, Sample.Type != "Bulk Soil" & Sample.Type != "Bulk Soil Inoculant") #retain only rhizosphere and tare soil samples

##Make long-format pair-wise distance data frame

phy.dist <- phy %>% phyloseq::distance("bray")
list <-  dist2list(phy.dist)

##Filter long data frame (from step 1) to retain only rows with a tare soil sample in column one and only rhizosphere samples in column

#seperate sample ID into new rows to filter
newdf <- cSplit(list, "col", sep = ".", drop = FALSE)
newdf <- cSplit(newdf, "row", sep = ".", drop = FALSE)

#drop some columns and filter to retain only rhizospphere in col and tare soil in row 
df<- subset(newdf, select = -c(col_4, row_4))
df1 <- subset(df, col_3 != "NA")
df1$row_3 <- as.character(df1$row_3)
df1$row_3 %<>% replace_na("tare")
df2<- subset(df1, row_3 == "tare")
df2 %<>% subset(., select = -c(row_2, row_3))
#rename columns 
dist.list <- df2 %>% 
  rename(
    Rhizo = col,
    Tare = row,
    Distance = value,
    SL.R =col_1,
    Treatment.R = col_2,
    Soil.R= col_3,
    SL.T =row_1,
  )

# Add a column in the long data frame that is TRUE if the values (seed lots) in the columns add in steps 4 and 5 are the same and FALSE otherwise
dist.list$self <- dist.list$SL.R == dist.list$SL.T

#Save 
saveRDS(dist.list,file="./output/Rda/dist.list.fun.Rda")


############################
#BACTERIA
############################

#clear global environment
rm(list = ls())

phy <- readRDS("./phy_object_bac/phy_bac_allsamples.rds")
phy <- subset_samples(phy, Sample.Type != "Bulk Soil" & Sample.Type != "Bulk Soil Inoculant")

#Make long-format pair-wise distance data frame
phy.dist <- phy %>% phyloseq::distance("bray")
list <-  dist2list(phy.dist)

#Filter long data frame (from step 1) to retain only rows with a tare soil sample in column one and only rhizosphere samples in column

#seperate sample ID into new rows to filter
library(splitstackshape)
newdf <- cSplit(list, "col", sep = ".", drop = FALSE)
newdf <- cSplit(newdf, "row", sep = ".", drop = FALSE)

#drop some colums and filter to retain only rhizo in col and tare in row 
df<- subset(newdf, select = -c(col_4, row_4))
df1 <- subset(df, col_3 != "NA")
df1$row_3 <- as.character(df1$row_3)
df1$row_3 %<>% replace_na("tare")
df2<- subset(df1, row_3 == "tare")
df2 %<>% subset(., select = -c(row_2, row_3))
#rename columns 
dist.list <- df2 %>% 
  rename(
    Rhizo = col,
    Tare = row,
    Distance = value,
    SL.R =col_1,
    Treatment.R = col_2,
    Soil.R= col_3,
    SL.T =row_1,
  )

# Add a column in the long data frame that is TRUE if the values (seed lots) in the columns add in steps 4 and 5 are the same and FALSE otherwise
dist.list$self <- dist.list$SL.R == dist.list$SL.T

#Save 
saveRDS(dist.list,file="./output/Rda/dist.list.bac.Rda")