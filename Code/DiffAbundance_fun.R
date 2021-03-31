#FUNGI

#Differential abundance testing using the Aldex2 pipeline. Tables created here are used to create Table 3 and Supplementary Figure S3 of manuscript
#modular pipeline= aldex.clr --> aldex.ttest (for two groups)/ aldex.glm --> 
#aldex.effect --> aldex.plot


#libraries
library(phyloseq)
library(magrittr)
library(ALDEx2)
library(reshape2)

###########
#read in phy object and filter 
(phy <- readRDS("./phy_object_fun/phy_fungi_prefilter.rds"))
phy %<>% subset_samples(., Sample.Type == "Rhizosphere Soil")

phy %<>% filter_taxa(., function(x) {sum(x>0) > 0.03*nsamples(.)}, TRUE) %>% prune_taxa(taxa_sums(.) > 0,.)
sampledata <- sample_data(phy)%>% data.frame

#extract OTU table matrix from phyloseq object 
otu <- (as(otu_table(phy, taxa_are_rows = FALSE), 'matrix'))
otu %<>% t()

#Make vectors of explanatory variables to test #this is for aldex.glm()
treatment <- sample_data(phy)$Treatment %>% as.character()
soiltype <- sample_data(phy)$Soil.Type %>% as.character()

otu.clr.treatment <- aldex.clr(otu, conds = treatment, mc.samples = 128)
otu.clr.soiltype <- aldex.clr(otu, conds = soiltype, mc.samples = 128)

###########
###TREATMENT
###########

ttest.Treat <- aldex.ttest(otu.clr.treatment)
head(ttest.Treat)
#gives you p-values and such in a table. will need to filter. 
#effect size estimates
effect.Treat <-aldex.effect(otu.clr.treatment)
head(effect.Treat)
#combine ttest and effect sizes
ttest.effect.Treat <- data.frame(ttest.Treat, effect.Treat)

#make a table from the results
alpha = 0.05
sigtab = ttest.effect.Treat[which(ttest.effect.Treat$we.eBH < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)#dimensions of the table. First number is how many significant OTUs

#saveRDS(sigtab, "./output/Rda/aldexTreatment.Rda")
#write.csv(sigtab, "./output/Tables/AldexSigTreatment_fungi.csv")

###########
###SOILTYPE
###########

ttest.soil <- aldex.ttest(otu.clr.soiltype)
head(ttest.soil)
#gives you p-values and such in a table. will need to filter. 
#effect size estimates
effect.soil <-aldex.effect(otu.clr.soiltype)
head(effect.soil)
#combine ttest and effect sizes
ttest.effect.soil <- data.frame(ttest.soil, effect.soil)

##creating tables
#NOTE that for obtaining values of absolute effect size greater than one, this first step should be done for both postive and negative directions. 
sigsoil = ttest.effect.soil[which(ttest.effect.soil$effect < -1), ]
#sigsoil = ttest.effect.soil[which(ttest.effect.soil$effect > 1), ]
#sigsoil = ttest.effect.soil[which(ttest.effect.soil$we.eBH < alpha), ] #this instead gets p-values
sigsoil = cbind(as(sigsoil, "data.frame"), as(tax_table(phy)[rownames(sigsoil), ], "matrix"))
head(sigsoil)
dim(sigsoil)

#two CSvs were exported, one for positive and one for negative effect sizes greater than one. those data frames were combined in excell and saved as "aldex_SOIL_effectsize_fun.csv" which is used to create Fig S3 in code "FigS3_SoilOrigin.R"
