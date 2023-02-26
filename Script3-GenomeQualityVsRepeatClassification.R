#Script3: Impact of genome quality on repeat classification

#Parts
#3.1: Impact of contigN50 on %repeats unclassfied and classified
#3.2: Impact of scaffold count on %repeats unclassfied and classified
#3.3: Impact of %BUSCO on %repeats unclassfied and classified


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

#-------------------------PART3.1---------------------------------------------------#

#3.1: Impact of contigN50 on %repeats unclassfied and classified

#Scatterplot: unclassified
ggplot(data=genomes_sub_noNA, aes(x=ctg_N50_Kb, y=unclassified_repeats)) + 
  geom_point(color="blue") +
  scale_x_log10() +
  stat_regline_equation(label.y.npc = 0.8, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.87, aes(label = ..rr.label..)) +
  theme_classic() +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  ylab ("% unclassfied repeats") +
  xlab (expression(Log[10]~contig~N50)) +
  theme(legend.position="none")

#Scatterplot: classified
ggplot(data=genomes_sub_noNA, aes(x=ctg_N50_Kb, y=classified_repeats)) + 
  geom_point(color="blue") +
  scale_x_log10() +
  stat_regline_equation(label.y.npc = 0.8, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.87, aes(label = ..rr.label..)) +
  theme_classic() +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  ylab ("% classfied repeats") +
  xlab (expression(Log[10]~contig~N50)) +
  theme(legend.position="none")

#Linear regression: unclassified
genomes.lm <- lm(unclassified_repeats ~ ctg_N50_Kb_log10, data = genomes_sub_noNA)
summary(genomes.lm)
# Residual standard error: 7.463 on 29 degrees of freedom
# Multiple R-squared:  0.1824,	Adjusted R-squared:  0.1542 
# F-statistic: 6.471 on 1 and 29 DF,  p-value: 0.01655

#Linear regression: classified
genomes.lm <- lm(classified_repeats ~ ctg_N50_Kb_log10, data = genomes_sub_noNA)
summary(genomes.lm)
# Residual standard error: 10.13 on 29 degrees of freedom
# Multiple R-squared:  0.3063,	Adjusted R-squared:  0.2824 
# F-statistic:  12.8 on 1 and 29 DF,  p-value: 0.00124

#-------------------------PART3.2---------------------------------------------------#

#3.2: Impact of scaffold count on %repeats unclassfied and classified

#Scatterplot: unclassified
ggplot(data=genomes_sub_noNA, aes(x=n_scaffolds_kb_log10, y=unclassified_repeats)) + 
  geom_point(color="blue") +
  scale_x_log10() +
  stat_regline_equation(label.y.npc = 0.8, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.87, aes(label = ..rr.label..)) +
  theme_classic() +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  ylab ("% unclassfied repeats") +
  xlab (expression(Log[10]~scaffold~count)) +
  theme(legend.position="none")

#Scatterplot: classified
ggplot(data=genomes_sub_noNA, aes(x=n_scaffolds_kb_log10, y=classified_repeats)) + 
  geom_point(color="blue") +
  scale_x_log10() +
  stat_regline_equation(label.y.npc = 0.8, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.87, aes(label = ..rr.label..)) +
  theme_classic() +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  ylab ("% classfied repeats") +
  xlab (expression(Log[10]~scaffold~count)) +
  theme(legend.position="none")

#Linear regression: unclassified
genomes.lm <- lm(unclassified_repeats ~ n_scaffolds_kb_log10, data = genomes_sub_noNA)
summary(genomes.lm)
# Residual standard error: 7.985 on 29 degrees of freedom
# Multiple R-squared:  0.06407,	Adjusted R-squared:  0.03179 
# F-statistic: 1.985 on 1 and 29 DF,  p-value: 0.1695

#Linear regression: classified
genomes.lm <- lm(classified_repeats ~ n_scaffolds_kb_log10, data = genomes_sub_noNA)
summary(genomes.lm)
# Residual standard error: 10.79 on 29 degrees of freedom
# Multiple R-squared:  0.212,	Adjusted R-squared:  0.1848 
# F-statistic: 7.801 on 1 and 29 DF,  p-value: 0.009152

#-------------------------PART3.3---------------------------------------------------#

#3.3: Impact of %BUSCO on %repeats unclassfied and classified

#Scatterplot: unclassified
ggplot(data=genomes_sub_noNA, aes(x=per_busco_complete, y=unclassified_repeats)) + 
  geom_point(color="blue") +
  #  scale_x_log10() +
  stat_regline_equation(label.y.npc = 0.8, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.87, aes(label = ..rr.label..)) +
  theme_classic() +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  ylab ("% unclassfied repeats") +
  xlab ("% BUSCO complete") +
  theme(legend.position="none")

#Scatterplot: classified
ggplot(data=genomes_sub_noNA, aes(x=per_busco_complete, y=classified_repeats)) + 
  geom_point(color="blue") +
  #  scale_x_log10() +
  stat_regline_equation(label.y.npc = 0.8, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.87, aes(label = ..rr.label..)) +
  theme_classic() +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  ylab ("% classfied repeats") +
  xlab ("% BUSCO complete") +
  theme(legend.position="none")

#Linear regression: unclassified
genomes.lm <- lm(unclassified_repeats ~ per_busco_complete, data = genomes_sub_noNA)
summary(genomes.lm)
# Residual standard error: 7.372 on 29 degrees of freedom
# Multiple R-squared:  0.2022,	Adjusted R-squared:  0.1747 
# F-statistic: 7.351 on 1 and 29 DF,  p-value: 0.01115

#Linear regression: classified
genomes.lm <- lm(classified_repeats ~ per_busco_complete, data = genomes_sub_noNA)
summary(genomes.lm)
# Residual standard error: 11.3 on 29 degrees of freedom
# Multiple R-squared:  0.1368,	Adjusted R-squared:  0.107 
# F-statistic: 4.596 on 1 and 29 DF,  p-value: 0.04055

