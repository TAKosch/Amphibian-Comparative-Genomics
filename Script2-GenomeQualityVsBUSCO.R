#Script2: Impact of genome quality on %BUSCO complete

#Parts
#2.1: Impact of contigN50 on %BUSCO complete
#2.2: Impact of scaffold count on %BUSCO complete
#2.3: Impact of contigN50 on %BUSCO complete (colored by read length)
#2.4: Impact of scaffold count on %BUSCO complete (colored by read length)

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

#-------------------------PART2.1---------------------------------------------------#

#2.1: Impact of contigN50 on %BUSCO complete

#Scatterplot
ggplot(data=genomes_sub_noNA, aes(x=ctg_N50_Kb, y=per_busco_complete)) + 
  geom_point() + 
  scale_x_log10(breaks = c(0.1, 1.0, 10.0, 100.0, 1000.0), label = c("0.1", "1", "10", "100", "1000")) +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  stat_regline_equation(label.y.npc = 0.8, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.87, aes(label = ..rr.label..)) +
  ylab ("BUSCO percent complete") +
  xlab(expression(Log[10]~contig~N50~(Kb))) +
  theme_classic() +
  theme(legend.position="none") 

#Linear regression
genomes.lm <- lm(per_busco_complete ~ ctg_N50_Kb_log10, data = genomes_sub_noNA)
summary(genomes.lm)
# Residual standard error: 17.66 on 29 degrees of freedom
# Multiple R-squared:  0.6551,	Adjusted R-squared:  0.6432 
# F-statistic: 55.09 on 1 and 29 DF,  p-value: 0.00000003531

#-------------------------PART2.2---------------------------------------------------#

#2.2: Impact of scaffold count on %BUSCO complete

#Scatterplot
ggplot(data=genomes_sub_noNA, aes(x=n_scaffolds_kb, y=per_busco_complete)) + 
  geom_point() + 
  scale_x_log10() +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  stat_regline_equation(label.y.npc = 0.2, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y.npc = 0.27, aes(label = ..rr.label..)) +
  ylab ("BUSCO percent complete") +
  xlab ("Log10 Scaffold count (Kb)") +
  theme_classic() +
  theme(legend.position="none")

#Linear regression
genomes.lm <- lm(per_busco_complete ~ n_scaffolds_kb_log10, data = genomes_sub_noNA)
summary(genomes.lm)
# Residual standard error: 19.08 on 29 degrees of freedom
# Multiple R-squared:  0.5971,	Adjusted R-squared:  0.5832 
# F-statistic: 42.98 on 1 and 29 DF,  p-value: 0.0000003509

#-------------------------PART2.3---------------------------------------------------#

#2.3: Impact of contigN50 on %BUSCO complete (colored by read length)
ggplot(data=genomes_sub_noNA, aes(x=ctg_N50_Kb, y=per_busco_complete, color=readlength)) + 
  geom_point() + 
  scale_x_log10(breaks = c(0.1, 1.0, 10.0, 100.0, 1000.0), label = c("0.1", "1", "10", "100", "1000")) +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  ylab ("BUSCO percent complete") +
  xlab(expression(Log[10]~contig~N50~(Kb))) +
  theme_classic() +
  theme(legend.position="none") 

#-------------------------PART2.4---------------------------------------------------#

#2.4: Impact of scaffold count on %BUSCO complete (colored by read length)

#Scatterplot
ggplot(data=genomes_sub_noNA, aes(x=n_scaffolds_kb, y=per_busco_complete, color=readlength)) + 
  geom_point() + 
  scale_x_log10() +
  geom_smooth(method=lm, level=0.95, colour="black", linewidth=0.2, alpha=0.2) +
  ylab ("BUSCO percent complete") +
  xlab(expression(Log[10]~scaffold~count~(Kb))) +
  theme_classic() +
  theme(legend.position="none")
