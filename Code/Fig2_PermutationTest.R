#permutation test examines whether microbial communities in rhizosphere samples and tare soil samples from the same seed lot have a mean Bray-Curtis community dissimilarity that is lower than rhizosphere samples paired with tare soils from different seed lots
#Creates Fig2 - visual representation of the permutation test for fungi
#dist.list Rdas are generated in R code "DistanceList.R"


#read in libraries
library(tidyverse)
library(magrittr)
library(foreach)
library(ggplot2)
library(egg)

############################
#FUNGI
############################

#read in the distance matrix. (a distance matrix of all rhizospehre samples and tare soil samples was created and then filtered to retain only distances between "Rhizo" and "Tare" -- see R code "DistanceList.R". SL.T is the seed lot of the Tare soil and SL.R of Rhizoshere)
dist.list<- readRDS("./Output/Rda/dist.list.fun.Rda")

set.seed(27)
# get a test statistic for a difference in means
observed.stat <- dist.list %>% group_by(self) %>% summarize(mean=mean(Distance)) %$% mean %>% diff
#summarize() creates a new data frame, with a row for each combination of grouping variables. $ calls on the mean. then diff calls on difference between the means 

#Creates a little table identifying tare soil and its seed lot, in order to shuffle later
lookup<-dist.list %>% select(Tare,SL.T) %>% distinct %>% rename(SL.T.perm =SL.T)

# Run permutations
number.of.permutations <- 1000
perm.t <- foreach(i=1:number.of.permutations, .combine=c) %do% {
  #Shuffle the data under the null (no difference in means). reassign Tare Soil Seed lot (SL.T). reassign "Self" column
  lookup$SL.T.perm %<>% sample()
  dist.list.new <- dist.list %>% left_join(lookup)%>% mutate(self=(SL.R==SL.T.perm))
  #Get test statistic for permuted data
  #calculate difference in average distance between self and not self. 
  dist.list.new %>% group_by(self) %>% summarize(mean=mean(Distance)) %$% mean %>% diff %>% as.numeric()
}

# visualize distribution of the permutations
hist(perm.t, 50)
#plot the actual statistic over the distribution of permutation generated statistics
abline(v=observed.stat, col='red')

#Calculate the permutation p-value
#This is the number of times the permutation test generated test stats more extreme that the observed value
sum(abs(perm.t)>abs(observed.stat))/(1+number.of.permutations)

############################
#create figure

dist.list<- readRDS("./Output/Rda/dist.list.fun.Rda")
dist.list$SL.T <- as.factor(dist.list$SL.T)

ggplot(dist.list, aes(x= Distance, group = self, colour = self)) +
  geom_density() + 
  labs(colour="SL Match", shape = "Treatment") +
  xlab("Bray-Curtis Dissimilarity")+
  scale_x_continuous(breaks=seq(0.7,0.95,0.1)) +
  ylab("Density of Samples") +
  facet_grid(Treatment.R ~ SL.T)+
  theme_article() +
  theme(legend.position="bottom")

ggsave("Fig2_PermTest_fun.jpeg", path = "./Output/Figures", width = 9, height = 4.5, units = "in")


############################
#BACTERIA
############################

rm(list = ls()) #cLEAR global environment

#read in the distance matrix. (a distance matrix of all rhizospehre samples and tare soil samples was created and then filtered to retain only distances between "Rhizo" and "Tare" -- see table. SL.T is the seed lot of the Tare soil and SL.R of Rhizoshere)
dist.list<- readRDS("./output/Rda/dist.list.bac.Rda")

set.seed(27)
# get a test statistic for a difference in means
observed.stat <- dist.list %>% group_by(self) %>% summarize(mean=mean(Distance)) %$% mean %>% diff
#summarize() creates a new data frame, with a row for each combination of grouping variables. $ calls on the mean. then diff calls on difference between the means 

#Creates a little table identifying tare soil and its seed lot, in order to shuffle later
lookup<-dist.list %>% select(Tare,SL.T) %>% distinct %>% rename(SL.T.perm =SL.T)

# Run permutations
number.of.permutations <- 1000
perm.t <- foreach(i=1:number.of.permutations, .combine=c) %do% {
  #Shuffle the data under the null (no difference in means). reassign Tare Soil Seed lot (SL.T). reassign "Self" column
  lookup$SL.T.perm %<>% sample()
  dist.list.new <- dist.list %>% left_join(lookup)%>% mutate(self=(SL.R==SL.T.perm))
  #Get test statistic for permuted data calculate difference in average distance between self and not self. 
  dist.list.new %>% group_by(self) %>% summarize(mean=mean(Distance)) %$% mean %>% diff %>% as.numeric()
}

# visualize distribution of the permutations
hist(perm.t, 50)
#plot the actual statistic over the distribution of permutation generated statistics
abline(v=observed.stat, col='red')

#Calculate the permutation p-value
#This is the number of times the permutation test generated test stats more extreme that the observed value
sum(abs(perm.t)>abs(observed.stat))/(1+number.of.permutations)
