#!/bin/bash
#$ -l h_vmem=5G
#$ -pe smp 10
#$ -j yes
#$ -cwd
#$ -e /users/k1223459/brc_scratch/AD_Classifcation/AD_vs_control
#$ -V


# module

module load bioinformatics/R/3.5.0
module load compilers/gcc/8.1.0

for x in {102..1001}
do
qsub -N XGBoost_AD_vs_control$x -V -cwd -l h_rt=05:00:00 -l h_vmem=10G -pe smp 15 -b y Rscript --verbose /users/k1223459/brc_scratch/AD_Classifcation/AD_vs_control/AD_classifier_v5.0_LOGLOSS_gbtree_RFE.R $x
done
