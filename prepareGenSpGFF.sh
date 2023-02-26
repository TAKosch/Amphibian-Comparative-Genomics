#!/bin/bash

#                                     prepareGenSpGFF

#SBATCH --job-name=prepGenSpGFF
#SBATCH --mem=5M
#SBATCH --time=00:01:00  
#SBATCH --array=1-15 

#------------------------INFO---------------------------------#

#Create GFF input files to run Genespace

#Created by Tiffany Kosch (with assistance from Neil Young) on 17-Nov-2022

#Run from terminal as "sbatch prepareGenSpGFF.slurm <array-name>"

#------------------------INPUT--------------------------------#

INPUT=(`head -n ${SLURM_ARRAY_TASK_ID} ${1} | tail -n 1`)

TSV=$(echo ${INPUT[@]} | cut -f 3 -d " ") #e.g. /full-path-to/X-trop_full_table.tsv
ID=$(echo ${INPUT[@]} | cut -f 2 -d " ")

#------------------------RUN----------------------------------#

#Make output directories
mkdir -p ./Genespace/rawGenomes/${ID}/${ID}/annotation

#Step1 (keep only complete, remove excess data, add extra cols)
awk '{if($2=="Complete") print $3,"busco","gene",$4,$5,$6,"locus="$1}' ${TSV} | sort -k1,1 > Genespace/rawGenomes/${ID}/${ID}/annotation/${ID}_gene.gff

#Step2
awk '{print $1}' Genespace/rawGenomes/${ID}/${ID}/annotation/${ID}_gene.gff | sort | uniq -c | awk '{if($1>10)print $2}' | sort | join - Genespace/rawGenomes/${ID}/${ID}/annotation/${ID}_gene.gff | awk 'OFS="\t" {print $1, $2, $3, $4, $5, $6, "", "", $7}' > ${ID}.tmp; mv ${ID}.tmp Genespace/rawGenomes/${ID}/${ID}/annotation/${ID}_gene.gff