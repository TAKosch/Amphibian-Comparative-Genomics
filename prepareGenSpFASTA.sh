#!/bin/bash

#                                     prepareGenSpFASTA

#SBATCH --job-name=prepGenSpFASTA
#SBATCH --mem=5M
#SBATCH --time=00:01:00  
#SBATCH --array=1-15 

#------------------------INFO---------------------------------#

#Create GFF input files to run Genespace

#Created by Tiffany Kosch (with assistance from Neil Young) on 17-Nov-2022

#Run from terminal as "sbatch prepareGenSpFASTA.slurm <array-name>"

#------------------------INPUT--------------------------------#

INPUT=(`head -n ${SLURM_ARRAY_TASK_ID} ${1} | tail -n 1`)

ID=$(echo ${INPUT[@]}  | cut -f 2 -d " ")  #Bgarg
ID2=$(echo ${INPUT[@]} | cut -f 3 -d " ")  #B-garg
DIR=/data/gpfs/projects/punim1525/Projects/Comparative-Genomics/TKosch_genome-edits/Busco_results/${ID2}_busco/run_tetrapoda_odb10/busco_sequences/single_copy_busco_sequences
DIR2=/data/scratch/projects/punim1525/GeneSpace/Practice/Genespace/rawGenomes/${ID}/${ID}/annotation

#------------------------RUN----------------------------------#

cd ${DIR}

awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-4); next} 1' *faa > ${ID}_protein.fa
mv ${ID}_protein.fa ${DIR2}