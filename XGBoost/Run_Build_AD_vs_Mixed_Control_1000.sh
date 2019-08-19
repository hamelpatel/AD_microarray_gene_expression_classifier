#!/bin/bash
#$ -l h_vmem=5G
#$ -pe smp 10
#$ -j yes
#$ -cwd
#$ -e /users/k1223459/brc_scratch/AD_Classifcation 
#$ -V


# module

module load bioinformatics/R/3.5.0
module load compilers/gcc/8.1.0

#for x in {1}
#do
#qsub -N XGBoost_$x -V -cwd -l h_rt=05:00:00 -l h_vmem=10G -pe smp 10 -b y Rscript --verbose /users/k1223459/brc_scratch/AD_Classifcation/AD_vs_mixedcontrol/AD_specific_classifier_v5.0_LOGLOSS_gbtree_RFE.R $x
#done


#for x in "/users/k1223459/brc_scratch/AD_Classifcation/AD_vs_mixedcontrol/missing_data.txt"
#do
#qsub -N XGBoost_$x -V -cwd -l h_rt=05:00:00 -l h_vmem=10G -pe smp 10 -b y Rscript --verbose /users/k1223459/brc_scratch/AD_Classifcation/AD_vs_mixedcontrol/AD_specific_classifier_v5.0_LOGLOSS_gbtree_RFE.R $x
#done

while read x
do
qsub -N XGBoost_$x -V -cwd -l h_rt=20:00:00 -l h_vmem=20G -pe smp 10 -b y Rscript --verbose /users/k1223459/brc_scratch/AD_Classifcation/AD_vs_mixedcontrol/AD_specific_classifier_v5.0_LOGLOSS_gbtree_RFE.R $x
done < "/users/k1223459/brc_scratch/AD_Classifcation/AD_vs_mixedcontrol/missing_data.txt"
