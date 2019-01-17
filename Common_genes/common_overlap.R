##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                        COMMON PROBES                                                                   #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################


#
# NOTES - 
# 
# common reliably detected genes in datasets with both male + female samples

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### LOAD DATASETS- (reliably detected probes) #####

# AD
AD_E_GEOD_63060<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/1.AD/E-GEOD-63060/Clean_Data/E-GEOD-63060_QCd.RDS")
AD_E_GEOD_63061<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/1.AD/E-GEOD-63061/Clean_Data/E-GEOD-63061_QCd.RDS")
AD_E_GEOD_6613<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/1.AD/E-GEOD-6613/Clean_Data/E-GEOD-6613_QCd.RDS")

# PD
PD_E_GEOD_6613<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/2.PD/E-GEOD-6613/Clean_Data/E-GEOD-6613_QCd.RDS")
PD_E_GEOD_72267<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/2.PD/E-GEOD-72267/Clean_Data/E-GEOD-72267_QCd.RDS")

# MS
MS_E_GEOD_24427<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/3.MS/E-GEOD-24427/Clean_Data/E-GEOD-24427_QCd.RDS")
MS_E_GEOD_16214<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/3.MS/E-GEOD-16214/Clean_Data/E-GEOD-16214_QCd.RDS")
MS_E_GEOD_41890<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/3.MS/E-GEOD-41890/Clean_Data/E-GEOD-41890_QCd.RDS")

# Schizophrenia
SCZ_E_GEOD_38484<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/4.Schizophrenia/E-GEOD-38484/Clean_Data/E-GEOD-38484_QCd.RDS")
SCZ_E_GEOD_27383<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/4.Schizophrenia/E-GEOD-27383/Clean_Data/E-GEOD-27383_QCd.RDS")
SCZ_E_GEOD_38481<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/4.Schizophrenia/E-GEOD-38481/Clean_Data/E-GEOD-38481_QCd.RDS")

#Bipolar Disorder
BP_E_GEOD_46449<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/5.BD/E-GEOD-46449/Clean_Data/E-GEOD-46449_QCd.RDS")
BP_E_GEOD_23848<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/5.BD/E-GEOD-23848/Clean_Data/E-GEOD-23848_QCd.RDS")

# Cardiovascular Disease
CD_E_GEOD_46097<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/6.Cardiovascular_Disease/E-GEOD-46097/Clean_Data/E-GEOD-46097_QCd.RDS")
CD_E_GEOD_59867<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/6.Cardiovascular_Disease/E-GEOD-59867/Clean_Data/E-GEOD-59867_QCd.RDS")
CD_E_GEOD_12288<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/6.Cardiovascular_Disease/E-GEOD-12288/Clean_Data/E-GEOD-12288_QCd.RDS")

#Rheumatoid Arthritis
RA_E_GEOD_74143<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/7.Rheumatoid_Arthritis/E-GEOD-74143/Clean_Data/E-GEOD-74143_QCd.RDS")
RA_E_GEOD_54629<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/7.Rheumatoid_Arthritis/E-GEOD-54629/Clean_Data/E-GEOD-54629_QCd.RDS")
RA_E_GEOD_42296<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/7.Rheumatoid_Arthritis/E-GEOD-42296/Clean_Data/E-GEOD-42296_QCd.RDS")

#Chonic Obstructive Pulmonary Disease
COPD_E_GEOD_54837<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/8.Chronic_Obstructive_Pulmonary_Disease/E-GEOD-54837/Clean_Data/E-GEOD-54837_QCd.RDS")
COPD_E_GEOD_42057<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/8.Chronic_Obstructive_Pulmonary_Disease/E-GEOD-42057/Clean_Data/E-GEOD-42057_QCd.RDS")

#ALS
ALS_E_TABM_940<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/9.ALS/E_TABM-940/Clean_Data/E-TABM-940_QCd.RDS")

##### LOAD DATASETS- (full) #####

# AD
AD_E_GEOD_63060_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/1.AD/E-GEOD-63060/Clean_Data/E-GEOD-63060_QCd_FULL.RDS")
AD_E_GEOD_63061_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/1.AD/E-GEOD-63061/Clean_Data/E-GEOD-63061_QCd_FULL.RDS")
AD_E_GEOD_6613_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/1.AD/E-GEOD-6613/Clean_Data/E-GEOD-6613_QCd_FULL.RDS")

# PD
PD_E_GEOD_6613_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/2.PD/E-GEOD-6613/Clean_Data/E-GEOD-6613_QCd_FULL.RDS")
PD_E_GEOD_72267_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/2.PD/E-GEOD-72267/Clean_Data/E-GEOD-72267_QCd_FULL.RDS")

# MS
MS_E_GEOD_24427_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/3.MS/E-GEOD-24427/Clean_Data/E-GEOD-24427_QCd_FULL.RDS")
MS_E_GEOD_16214_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/3.MS/E-GEOD-16214/Clean_Data/E-GEOD-16214_QCd_FULL.RDS")
MS_E_GEOD_41890_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/3.MS/E-GEOD-41890/Clean_Data/E-GEOD-41890_QCd_FULL.RDS")

# Schizophrenia
SCZ_E_GEOD_38484_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/4.Schizophrenia/E-GEOD-38484/Clean_Data/E-GEOD-38484_QCd_FULL.RDS")
SCZ_E_GEOD_27383_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/4.Schizophrenia/E-GEOD-27383/Clean_Data/E-GEOD-27383_QCd_FULL.RDS")
SCZ_E_GEOD_48702_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/4.Schizophrenia/E-GEOD-48072/Clean_Data/E-GEOD-48072_QCd_FULL.RDS")


#Bipolar Disorder
BP_E_GEOD_46449_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/5.BD/E-GEOD-46449/Clean_Data/E-GEOD-46449_QCd_FULL.RDS")
BP_E_GEOD_23848_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/5.BD/E-GEOD-23848/Clean_Data/E-GEOD-23848_QCd_FULL.RDS")

# Cardiovascular Disease
CD_E_GEOD_46097_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/6.Cardiovascular_Disease/E-GEOD-46097/Clean_Data/E-GEOD-46097_QCd_FULL.RDS")
CD_E_GEOD_59867_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/6.Cardiovascular_Disease/E-GEOD-59867/Clean_Data/E-GEOD-59867_QCd_FULL.RDS")
CD_E_GEOD_12288_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/6.Cardiovascular_Disease/E-GEOD-12288/Clean_Data/E-GEOD-12288_QCd_FULL.RDS")

#Rheumatoid Arthritis
RA_E_GEOD_74143_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/7.Rheumatoid_Arthritis/E-GEOD-74143/Clean_Data/E-GEOD-74143_QCd_FULL.RDS")
RA_E_GEOD_54629_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/7.Rheumatoid_Arthritis/E-GEOD-54629/Clean_Data/E-GEOD-54629_QCd_FULL.RDS")
RA_E_GEOD_42296_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/7.Rheumatoid_Arthritis/E-GEOD-42296/Clean_Data/E-GEOD-42296_QCd_FULL.RDS")

#Chonic Obstructive Pulmonary Disease
COPD_E_GEOD_54837_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/8.Chronic_Obstructive_Pulmonary_Disease/E-GEOD-54837/Clean_Data/E-GEOD-54837_QCd_FULL.RDS")
COPD_E_GEOD_42057_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/8.Chronic_Obstructive_Pulmonary_Disease/E-GEOD-42057/Clean_Data/E-GEOD-42057_QCd_FULL.RDS")

#ALS
ALS_E_TABM_940_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/9.ALS/E_TABM-940/Clean_Data/E-TABM-940_QCd_FULL.RDS")

##### MERGE #####

common_probes<-as.data.frame(table(c(
  colnames(AD_E_GEOD_63060),
  colnames(AD_E_GEOD_63061),
  colnames(AD_E_GEOD_6613),
  colnames(PD_E_GEOD_6613),
  colnames(PD_E_GEOD_72267),
  colnames(MS_E_GEOD_24427),
  colnames(MS_E_GEOD_16214),
  colnames(MS_E_GEOD_41890),
  colnames(SCZ_E_GEOD_38484),
#  colnames(SCZ_E_GEOD_27383),
  colnames(SCZ_E_GEOD_38481),
#  colnames(BP_E_GEOD_46449),
  colnames(BP_E_GEOD_23848),
  colnames(CD_E_GEOD_46097),
  colnames(CD_E_GEOD_59867),
  colnames(CD_E_GEOD_12288),
  colnames(RA_E_GEOD_74143),
  colnames(RA_E_GEOD_54629),
  colnames(RA_E_GEOD_42296),
  colnames(COPD_E_GEOD_54837),
  colnames(COPD_E_GEOD_42057),
  colnames(ALS_E_TABM_940))))

#length

dim(common_probes)
head(common_probes)

table(common_probes$Freq)

common_probes<-common_probes[common_probes$Freq=="20",]
nrow(common_probes)-6

##### SAVE GENE LIST #####

setwd("/media/hamel/Workspace/Dropbox/Projects/AD-classification/4.Overlap_of_genes/")

saveRDS(common_probes$Var1, file="reliably_detetced_genes.RDS")

