## Description of R Code found withing the "Code" directory for data analysis

Filtering_bac.R and Filtering_fun.R take the raw phyloseq object output by the bioinformatics/processing steps (found under the CodeProcessing directory) and applies filtering steps, including converting reads to relative abundances, filtering rare OTUs and applying transformations. 
Phyloseq objects produced from this code are used througout the subsequent analyses.

Fig1_dbRDA.R - code to perform a distance based redundancy analysis and create ordination

DistanceList.R and Fig2_PermutationTest.R are both used for the permutation test, which was used to determine if the mean Bray-Curtis dissimilarity between tare soil and rhizosphere samples paired from the same seed lot was lower than paired tare soil and rhizosphere samples from different seed lots. DistanceList.R creates the dataframe input to conduct the test and create a figure (Fig2_PermutationTest.R).

Fig3_HeatMap creates the heat map depicted in Figure 3

DiffAbundance_bac.R and Diff_abundace_fun.R are code for running a test of differential abundance using ALDEx2

the R code for the analysis and creation of supplementary figures is also included: FigS1_S2_NMDSOrdinations.R and FigS3_SoilOrigin.R

