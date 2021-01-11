#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N test              
#$ -cwd                  
#$ -l h_rt=00:01:00 

# Initialise the environment modules
. /etc/profile.d/modules.sh

export R_LIBS_USER=/exports/eddie/scratch/s1505825/R/library; . /etc/profile.d/modules.sh; module load igmm/apps/R/3.6.3; module load igmm/apps/BEDTools/2.27.1

Rscript BATCH test.R
