#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N analysis              
#$ -cwd                  
#$ -l h_rt=00:02:00 
#$ -l h_vmem=100G

# Initialise the environment modules
. /etc/profile.d/modules.sh

rsetup

R CMD BATCH Code.R
