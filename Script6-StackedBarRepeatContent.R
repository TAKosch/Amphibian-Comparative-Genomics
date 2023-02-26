#Script6: Stacked bar plot of repeat content by species

#Parts
#6.1: Stacked bar plot of repeat content by species

#-------------------------PREPARE-ENVIRONMENT---------------------------------------#

#Open libraries
library(ggplot2)
library(ggsci)
library(tidyverse)
library("scales")

#Set working directory
setwd("/Users/tkosch/OneDrive - The University of Melbourne/Documents/UoM/Projects/Amphib-Comparative-Genomics/R_analyses/")


#-------------------------INPUT&FORMATTING------------------------------------------#

#Input data (sorted by %repeats, then repeat type)
file2 <-read.delim("RepeatMasker-percents_trunc_v9.csv", header = TRUE, stringsAsFactors = FALSE, quote = "", sep = ",") 


#-------------------------PART6.1---------------------------------------------------#

#Stacked bar plot (sorted by %repeats, then repeat type)
ggplot(file2, aes(fill=Repeat.type, y=Percent.of.genome, x=fct_inorder(Species))) + #"fct_inorder" preserves order
  geom_bar(position="stack", stat="identity") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  ylab("Percent of the genome") +
  xlab("Species") +
  labs(fill='Repeat types') +
  scale_fill_manual(values = c("bisque3","#E64B35FF", "#4DBBD5FF", "#00A087FF", 
                               "#3C5488FF", "#F39B7FFF"),
                    labels = c("Unclassified", "DNA transposons", "LINEs", "LTRs", 
                               "Simple", "SINEs"))

