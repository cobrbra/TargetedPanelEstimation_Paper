#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N analysis              
#$ -cwd                  
#$ -l h_rt=00:02:00 
#$ -l h_vmem=100G

# Initialise the environment modules
. /etc/profile.d/modules.sh

## RUN THIS FIRST
#export R_LIBS_USER=/exports/eddie/scratch/s1505825/R/library; . /etc/profile.d/modules.sh; module load igmm/apps/R/3.6.3; module load igmm/apps/BEDTools/2.27.1

R CMD BATCH Code.R
