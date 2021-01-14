#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N analysis              
#$ -cwd                  
#$ -l h_rt=02:00:00 
#$ -l h_vmem=64G

# Initialise the environment modules
. /etc/profile.d/modules.sh

export R_LIBS_USER=/exports/eddie/scratch/s1505825/R/library; . /etc/profile.d/modules.sh; module load igmm/apps/R/3.6.3; module load igmm/apps/BEDTools/2.27.1

Rscript Code.R

# Note to push from Eddie requires 'git push origin HEAD:master'