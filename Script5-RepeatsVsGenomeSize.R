#Script5: Impact of percentage of repeats on genome size

#Parts
#5.1: Impact of percentage of repeats on genome size (w/ A-mex)
#5.2: Impact of percentage of repeats on genome size (w/out A-mex)


#-------------------------PREPARE-ENVIRONMENT---------------------------------------#

#Open libraries
library(ggplot2)
library(ggsci)
library(tidyverse)  #for preserve order script
library("scales")

#Set working directory
setwd("/Users/tkosch/OneDrive - The University of Melbourne/Documents/UoM/Projects/Amphib-Comparative-Genomics/R_analyses/")


#-------------------------INPUT&FORMATTING------------------------------------------#

#Input data
file2 <-read.delim("RepeatMasker-percents_v5.csv", header = TRUE, stringsAsFactors = FALSE, quote = "", sep = ",")

#Change genome size to Gb
file2$length_Gb <- file2$Length.bp / 1000000000

#-------------------------PART5.1---------------------------------------------------#

#5.1: Impact of percentage of repeats on genome size (w/ A-mex)

#ScatterPlot
ggplot(data=file2, aes(x=Per.Bases.masked, y=length_Gb)) + 
  geom_point(color="blue") +
  theme_classic() +
  ylab ("Genome size (Gb)") +
  xlab ("Percent repeats") +
  theme(legend.position="none") +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  stat_regline_equation(label.y.npc = 0.8, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.87, aes(label = ..rr.label..)) 

#Linear regression
genomes.lm <- lm(length_Gb ~ Per.Bases.masked, data = file2)
summary(genomes.lm)
# Residual standard error: 4.408 on 30 degrees of freedom
# Multiple R-squared:  0.1636,	Adjusted R-squared:  0.1357 
# F-statistic: 5.867 on 1 and 30 DF,  p-value: 0.02167

#-------------------------PART5.2---------------------------------------------------#

#5.2: Impact of percentage of repeats on genome size (w/out A-mex)

#Exclude A-mex
file2 <- file2[!grepl("Ambystoma mexicanum",file2$Species),]

#Scatterplot
ggplot(data=file2, aes(x=Per.Bases.masked, y=length_Gb)) + 
  geom_point(color="blue") +
  theme_classic() +
  ylab ("Genome size (Gb)") +
  xlab ("Percent repeats") +
  theme(legend.position="none") +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  stat_regline_equation(label.y.npc = 0.8, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.87, aes(label = ..rr.label..)) 

#Linear regression
genomes.lm <- lm(length_Gb ~ Per.Bases.masked, data = file2)
summary(genomes.lm)
# Residual standard error: 1.102 on 29 degrees of freedom
# Multiple R-squared:  0.627,	Adjusted R-squared:  0.6141 
# F-statistic: 48.74 on 1 and 29 DF,  p-value: 0.0000001124