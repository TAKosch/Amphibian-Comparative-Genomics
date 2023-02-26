#Script4: Impact of genome quality on percentage of genome as chromosomes

#Parts
#4.1: Impact of contigN50 on percentage of genome as chromosomes
#4.2: Impact of scaffold count on percentage of genome as chromosomes
#4.3: Impact of %BUSCO on percentage of genome as chromosomes

#-------------------------PREPARE-ENVIRONMENT---------------------------------------#

#Open libraries
library("tidyverse")
library("ggplot2")
library("ggpmisc")
library("ggpubr")

#Set working directory
setwd("/Users/tkosch/OneDrive - The University of Melbourne/Documents/UoM/Projects/Amphib-Comparative-Genomics/R_analyses/")


#-------------------------INPUT&FORMATTING------------------------------------------#

#Input data
genomes <- read.csv("bbmap-busco6.csv", row.names = 1, , header = TRUE)

#Adjust to not use scientific notation in plots
options(scipen=999) 

#Pull relevant columns
genomes_sub <- genomes[ , c("n_scaffolds", "scaf_bp", "ctg_N50" , "readlength", "chr.assembly", 
                            "per_busco_complete", "unclassified_repeats", "classified_repeats", 
                            "per.gen.chrs")]  

#Exclude NA rows
genomes_sub_noNA <- genomes_sub[complete.cases(genomes_sub), ] 

#Create new cols with n_scaffolds converted to Kb and Mb
genomes_sub_noNA$n_scaffolds_kb <- genomes_sub_noNA$n_scaffolds / 1000
genomes_sub_noNA$n_scaffolds_Mb <- genomes_sub_noNA$n_scaffolds / 1000000
genomes_sub_noNA$scaf_bp_Gb <- genomes_sub_noNA$scaf_bp / 1000000000
genomes_sub_noNA$ctg_N50_Kb <- genomes_sub_noNA$ctg_N50 / 1000

#Log transform
genomes_sub_noNA$ctg_N50_Kb_log10 <- log10(genomes_sub_noNA$ctg_N50_Kb)
genomes_sub_noNA$n_scaffolds_kb_log10 <- log10(genomes_sub_noNA$n_scaffolds_kb)

#Exclude genomes not assembled to chromosome level
genomes_sub_noNA_noZero <- genomes_sub_noNA[!grepl("no",genomes_sub_noNA$chr.assembly),] 

#-------------------------PART4.1---------------------------------------------------#

#4.1: Impact of contigN50 on percentage of genome as chromosomes

#Scatterplot
ggplot(data=genomes_sub_noNA_noZero, aes(x=ctg_N50_Kb, y=per.gen.chrs)) + 
  geom_point(color="blue") +
  scale_x_log10() +
  stat_regline_equation(label.y.npc = 0.2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.27, aes(label = ..rr.label..)) +
  theme_classic() +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  ylab ("% assembled to chromosomes") +
  xlab (expression(Log[10]~contig~N50)) +
  theme(legend.position="none")

#Linear regression
genomes.lm <- lm(per.gen.chrs ~ ctg_N50_Kb_log10, data = genomes_sub_noNA_noZero)
summary(genomes.lm)
# Residual standard error: 8.193 on 13 degrees of freedom
# Multiple R-squared:  0.4492,	Adjusted R-squared:  0.4068 
# F-statistic:  10.6 on 1 and 13 DF,  p-value: 0.006253

#-------------------------PART4.2---------------------------------------------------#

#4.2: Impact of scaffold count on percentage of genome as chromosomes

#Scatterplot
ggplot(data=genomes_sub_noNA_noZero, aes(x=n_scaffolds, y=per.gen.chrs)) + 
  geom_point(color="blue") +
  scale_x_log10() +
  stat_regline_equation(label.y.npc = 0.2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.27, aes(label = ..rr.label..)) +
  theme_classic() +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  ylab ("% assembled to chromosomes") +
  xlab (expression(Log[10]~scaffold~count)) +
  theme(legend.position="none")

#Linear regression
genomes.lm <- lm(per.gen.chrs ~ n_scaffolds_kb_log10, data = genomes_sub_noNA_noZero)
summary(genomes.lm)
# Residual standard error: 6.791 on 13 degrees of freedom
# Multiple R-squared:  0.6215,	Adjusted R-squared:  0.5924 
# F-statistic: 21.35 on 1 and 13 DF,  p-value: 0.0004795


#-------------------------PART4.3---------------------------------------------------#

#4.3: Impact of %BUSCO on percentage of genome as chromosomes

#Scatterplot
ggplot(data=genomes_sub_noNA_noZero, aes(x=per_busco_complete, y=per.gen.chrs)) + 
  geom_point() +
  scale_x_log10() +
  stat_regline_equation(label.y.npc = 0.2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.27, aes(label = ..rr.label..)) +
  theme_classic() +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  ylab ("% assembled to chromosomes") +
  xlab ("BUSCO percent complete") +
  theme(legend.position="none")

#Linear regression
genomes.lm <- lm(per.gen.chrs ~ per_busco_complete, data = genomes_sub_noNA_noZero)
summary(genomes.lm)
# Residual standard error: 9.574 on 13 degrees of freedom
# Multiple R-squared:  0.2477,	Adjusted R-squared:  0.1899 
# F-statistic: 4.281 on 1 and 13 DF,  p-value: 0.05903