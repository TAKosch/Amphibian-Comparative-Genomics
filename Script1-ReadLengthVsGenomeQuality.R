#Script1: Impact of read length on genome quality

#Parts
#1.1: Impact of read length on %BUSCO complete
#1.2: Impact of read length on contigN50
#1.3: Impact of read length on scaffold count
#1.4: Impact of read length on repeat classification
#1.5: Impact of read length on percentage of genome as chromosomes

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


#-------------------------PART1.1---------------------------------------------------#

#1.1: Impact of read length on %BUSCO complete

#Boxplot
ggplot(genomes_sub_noNA, aes(x=readlength, y=per_busco_complete)) + 
  geom_boxplot() + 
  theme_classic() +   
  ylab ("BUSCO complete") +
  xlab("Read type") 
#Is variance equal? NO

#T-test
t.test(per_busco_complete ~ readlength, var.equal = FALSE, data = genomes_sub_noNA)
#t = 2.8544, df = 13.448, p-value = 0.01316

#-------------------------PART1.2---------------------------------------------------#

#1.2: Impact of read length on contigN50

#Boxplot
ggplot(genomes_sub_noNA, aes(x=readlength, y=ctg_N50_Kb_log10)) + 
  geom_boxplot() + 
  theme_classic() +   
  ylab ("Congtig N50") +
  xlab("Read type")
#Is variance equal? YES

#T-test
t.test(ctg_N50_Kb_log10 ~ readlength, var.equal = TRUE, data = genomes_sub_noNA)
#t = 6.0471, df = 29, p-value = 0.0000014

#-------------------------PART1.3---------------------------------------------------#

#1.3: Impact of read length on scaffold count

#Boxplot
ggplot(genomes_sub_noNA, aes(x=readlength, y=n_scaffolds_kb_log10)) + 
  geom_boxplot() + 
  theme_classic() +   
  ylab (expression(Log[10]~scaffold~count)) +
  xlab("Read type")
#Is variance equal? YES

#T-test
t.test(n_scaffolds_kb_log10 ~ readlength, var.equal = TRUE, data = genomes_sub_noNA)
#t = -4.8559, df = 29, p-value = 0.00003786

#-------------------------PART1.4---------------------------------------------------#

#1.4: Impact of read length on repeat classification

#Box plot
ggplot(genomes_sub_noNA, aes(x=readlength, y=classified_repeats)) + 
  geom_boxplot() + 
  theme_classic() +   
  ylab ("% classfied repeats") +
  xlab("Read type")
#Is variance equal? NO

t.test(classified_repeats ~ readlength, var.equal = FALSE, data = genomes_sub_noNA)
#t = 3.5735, df = 26.241, p-value = 0.001393

#-------------------------PART1.5---------------------------------------------------#

#1.5: Impact of read length on percentage of genome as chromosomes

#Exclude genomes not assembled to chromosome level
genomes_sub_noNA_noZero <- genomes_sub_noNA[!grepl("no",genomes_sub_noNA$chr.assembly),] 

#Box plot
ggplot(genomes_sub_noNA, aes(x=readlength, y=per.gen.chrs)) + 
  geom_boxplot() + 
  theme_classic() +   
  ylab ("% assembled to chromosomes") +
  xlab("Read type")
#Is variance equal? NO

t.test(per.gen.chrs ~ readlength, var.equal = FALSE, data = genomes_sub_noNA)
#t = 3.6493, df = 28.484, p-value = 0.001047

