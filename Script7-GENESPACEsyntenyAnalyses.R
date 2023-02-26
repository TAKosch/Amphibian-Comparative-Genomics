#Script7: GENESPACE synteny analyses

#Parts
#7.1: Run on all chr-level genomes
#7.2: Run on all chr-level frog genomes
#7.3: Run on all chr-level caecilian genomes

#-------------------------PREPARE-ENVIRONMENT---------------------------------------#

## Install 
#https://github.com/jtlovell/GENESPACE
## Need to git clone and compile MCScanX on system
## https://github.com/wyp1125/MCScanX#installation
## cd ~/Software
## git clone https://github.com/wyp1125/MCScanX.git
## cd MCScanX
## make

## If running on mac, we had some problems compiling matrix. 
# xcode-select --install
# also ran from an R project called Genome_alignment

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("Biostrings", "rtracklayer"), force = TRUE)
# 
# if (!requireNamespace("devtools", quietly = TRUE))
#   install.packages("devtools")
# devtools::install_github("jtlovell/GENESPACE", upgrade = F)

## matrix always defaulted to the previous version but seems to be running okay

##Orthofinder installation with anaconda
#conda install -c bioconda orthofinder
#Make executable by adding to path (add to bash profile)
#export PATH=/Users/tkosch/opt/anaconda3/pkgs/orthofinder-2.5.4-hdfd78af_0//bin:$PATH 

#Load genespace
library(GENESPACE)


#-------------------------INPUT&FORMATTING------------------------------------------#

#First need to format input files and directories by running bash scripts:
# 1. prepareGenSpGFF.slurm (prepares input directories)
# 2. prepareGenSpFASTA.slurm (prepares fasta file)

#Input files must be in defined database (ID=gid; e.g., Xtrop)
#/Genespace/rawGenomes/${ID}/${ID}/annotation
#These are created by "prepareGenSpGFF.slurm"

#-------------------------PART7.1---------------------------------------------------#

#7.1: Run on all chr-level genomes

#Set wd
runwd <- file.path("/Users/tkosch/MyDesktopFiles/R-Files/Genespace/Run-chr-files/Genespace")

#Set IDs
gids <- c("Ecoq","Bgarg","Bbufo","Epust","Rtemp","Pads","Lail","Lleish","Hboe","Xtrop","Xlae","Amex","Muni","Gser","Rbiv")
#Ordered by phylogeny

#Initialize run (this creates a bunch of sub-directories)
gpar <- init_genespace(
  genomeIDs = gids, 
  speciesIDs = gids, 
  versionIDs = gids, 
  ploidy = rep(1,1),  
  wd = runwd, 
  gffString = "gff", 
  pepString = "protein",  
  path2orthofinder = "orthofinder", 
  path2mcscanx = "/Users/tkosch/Applications/MCScanX-master/",  
  rawGenomeDir = file.path(runwd, "rawGenomes"),overwrite = TRUE)

#Format annotations
parse_annotations(
  gsParam = gpar, 
  gffEntryType = "gene", 
  gffIdColumn = "locus",
  gffStripText = "locus=", 
  headerEntryIndex = 1,
  headerSep = " ", 
  headerStripText = "locus=",
  troubleshoot = TRUE)


# #Run orthofinder from terminal (takes ~1 h to run)
# cd /Users/tkosch/MyDesktopFiles/R-Files/Genespace/Run-chr-files
# 
# #Run
# orthofinder -f /Users/tkosch/MyDesktopFiles/R-Files/Genespace/Run-chr-files/Genespace/peptide -t 4 -a 1 -X -o /Users/tkosch/MyDesktopFiles/R-Files/Genespace/Run-chr-files/Genespace/orthofinder
#Note: this took >30 minutes to run
#Had to increase computer file soft limit number to run
#ulimit -n 325


#Run synteny search (takes about 30 min to run)
gpar <- synteny(gsParam = gpar,overwrite = TRUE,overwriteGff = TRUE,overwriteHits = TRUE,minGenes4of = 10)

#Make synteny plot
ripdat <- plot_riparianHits(gpar,minGenes2plot = 10, refGenome="Rbiv", #reference genome=Rbiv
                            highlightRef="lightred", 
                            annotatePlot=TRUE, 
                            nGenomeLabChar=0, #plot w/out genome names
                            reorderChrs=TRUE,
                            labelTheseGenomes = NULL)

#-------------------------PART7.2---------------------------------------------------#

#7.2: Run on all chr-level frog genomes

#Set wd
runwd <- file.path("/Users/tkosch/MyDesktopFiles/R-Files/Genespace/All-frogs/Genespace")

#Set IDs
gids <- c("Bbufo","Bgarg","Ecoq","Epust","Rtemp","Pads","Lleish","Lail","Hboe","Xtrop","Xlae")

#Initialize run (this creates a bunch of sub-directories)
gpar <- init_genespace(
  genomeIDs = gids, 
  speciesIDs = gids, 
  versionIDs = gids, 
  ploidy = rep(1,1,1,1,1,1,1,1,1,1,2),  
  wd = runwd, 
  gffString = "gff", 
  pepString = "protein",  
  path2orthofinder = "orthofinder", 
  path2mcscanx = "/Users/tkosch/Applications/MCScanX-master/",   
  rawGenomeDir = file.path(runwd, "rawGenomes"),overwrite = TRUE)

#Format annotations
parse_annotations(
  gsParam = gpar, 
  gffEntryType = "gene", 
  gffIdColumn = "locus",
  gffStripText = "locus=", 
  headerEntryIndex = 1,
  headerSep = " ", 
  headerStripText = "locus=",
  troubleshoot = TRUE)

# #Run orthofinder from terminal
# cd /Users/tkosch/MyDesktopFiles/R-Files/Genespace/Run-chr-files
# 
# #Run
# orthofinder -f /Users/tkosch/MyDesktopFiles/R-Files/Genespace/All-frogs/Genespace/peptide -t 4 -a 1 -X -o /Users/tkosch/MyDesktopFiles/R-Files/Genespace/All-frogs/Genespace/orthofinder

#Run synteny search
gpar <- synteny(gsParam = gpar,overwrite = TRUE,overwriteGff = TRUE,overwriteHits = TRUE,minGenes4of = 10) 

#Make synteny plot
ripdat <- plot_riparianHits(gpar,minGenes2plot = 10, refGenome="Xlae", 
                            highlightRef="lightred", 
                            annotatePlot=TRUE, 
                            nGenomeLabChar=5, 
                            reorderChrs=TRUE, 
                            chrRectBuffer=2,
                            useOrder=FALSE,                 #change to TRUE to scale chrs by buscos
                            invertTheseChrs = data.frame(   #use to invert chrs to improve aesthetics
                              genome = c("Lail", "Lail","Lail","Lail","Lail","Lail","Lail",
                                         "Lleish", "Lleish", "Lleish", "Lleish", "Lleish", "Lleish", 
                                         "Pads", "Pads", "Pads", "Pads",
                                         "Rtemp", "Rtemp", "Rtemp", "Rtemp", "Rtemp", "Rtemp", "Rtemp",
                                         "Epust", "Epust", "Epust", "Epust", 
                                         "Ecoq", "Ecoq", "Ecoq", 
                                         "Bgarg", "Bgarg", "Bgarg", 
                                         "Bbufo", "Bbufo", "Bbufo", "Bbufo", "Bbufo", "Bbufo"),
                              chr = c(12, 1, 4, 2, 8, 3, 7, 
                                      1, 3, 10, 4, 5, 7, 
                                      3, 1, 4, 5, 
                                      1, 2, 3, 11, 7, 4, 13,
                                      2, 4, 3, 5,
                                      2, 11, 8,
                                      2, 4, 8,
                                      2, 3, 6, 8, 9, 4)))


#-------------------------PART7.3---------------------------------------------------#

#7.3: Run on all chr-level caecilian genomes

#Set wd
runwd <- file.path("/Users/tkosch/MyDesktopFiles/R-Files/Genespace/3-caecilians-run/Genespace")

#Set IDs
gids <- c("Gser","Muni","Rbiv") #This orders how the files are plotted: first genome is on the bottom, last in on the top

#Initialize run (this creates a bunch of sub-directories)
gpar <- init_genespace(
  genomeIDs = gids, 
  speciesIDs = gids, 
  versionIDs = gids, 
  ploidy = rep(1,1),  
  wd = runwd, 
  gffString = "gff", 
  pepString = "protein",  
  path2orthofinder = "orthofinder", 
  path2mcscanx = "/Users/tkosch/Applications/MCScanX-master/",   
  rawGenomeDir = file.path(runwd, "rawGenomes"),overwrite = TRUE)

#Format annotations
parse_annotations(
  gsParam = gpar, 
  gffEntryType = "gene", 
  gffIdColumn = "locus",
  gffStripText = "locus=", 
  headerEntryIndex = 1,
  headerSep = " ", 
  headerStripText = "locus=",
  troubleshoot = TRUE)

# #Run orthofinder from terminal
# cd /Users/tkosch/MyDesktopFiles/R-Files/Genespace/3-caecilians-run
# 
# #Run
# orthofinder -f /Users/tkosch/MyDesktopFiles/R-Files/Genespace/3-caecilians-run/Genespace/peptide -t 4 -a 1 -X -o /Users/tkosch/MyDesktopFiles/R-Files/Genespace/3-caecilians-run/Genespace/orthofinder

#Run synteny search
gpar <- synteny(gsParam = gpar,overwrite = TRUE,overwriteGff = TRUE,overwriteHits = TRUE,minGenes4of = 10)

#Make synteny plot
ripdat <- plot_riparianHits(gpar,minGenes2plot = 10, refGenome="Rbiv", 
                            highlightRef="lightred", 
                            annotatePlot=TRUE, 
                            nGenomeLabChar=0, #plot w/out genome names
                            reorderChrs=TRUE, 
                            #useOrder=FALSE, #scales by chr size
                            labelTheseGenomes = NULL,
                            invertTheseChrs = data.frame(   #use to invert chrs to improve aesthetics
                              genome = c("Muni", "Gser", "Gser", "Muni", "Gser", 
                                         "Muni", "Gser", "Muni", "Muni", "Gser", "Gser", "Gser"),
                              chr = c(2,2,3,9,7,7,5,5,6,15,8,18)))