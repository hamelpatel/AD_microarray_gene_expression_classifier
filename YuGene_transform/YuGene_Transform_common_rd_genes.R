README##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                          EXTRACT COMMON PROBES AND YUGENE TRANSFORM                                                    #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################


#
# NOTES - 
# extract reliably detected genes from all datasets
# 

##### SET PARAMETERS #####

rm(list=ls())

options=(stringAsFactors=FALSE)

##### LIBRARY #####
library(YuGene)
library(WGCNA) # for pca plot
library(factoextra) # for pca plot
library(lumi) # for density plot
library(org.Hs.eg.db) # for gene symbol to entrez mapping
library(reshape) # for probe detection plot 
library(gridExtra) # for probe detection plot 

##### SET DIRECTORY #####

work_dir="/media/hamel/Workspace/Dropbox/Projects/AD-classification/5.YuGene_Transform/common_rd_genes"

setwd(work_dir)

# create directory for full dataset
dir.create(paste(work_dir,"Yugene_transformed_data", sep="/"))
Yugene_transformed_data_dir=paste(work_dir,"Yugene_transformed_data", sep="/")

# create directory for PCA plots
dir.create(paste(work_dir,"PCA", sep="/"))
PCA_dir=paste(work_dir,"PCA", sep="/")

# create directory for Boxplots
dir.create(paste(work_dir,"Boxplots", sep="/"))
boxplots_dir=paste(work_dir,"Boxplots", sep="/")

# housekeeping genes
dir.create(paste(work_dir,"hk_gene_plots", sep="/"))
hk_gene_plots_dir=paste(work_dir,"hk_gene_plots", sep="/")

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
SCZ_E_GEOD_38481_full<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/1.Data/4.Schizophrenia/E-GEOD-38481/Clean_Data/E-GEOD-38481_QCd_FULL.RDS")

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

##### READ RD PROBES #####

RD_genes<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/4.Overlap_of_genes/reliably_detetced_genes.RDS")

##### YuGene Transform (RD) #####

yugene_transform<-function(dataset, probes_to_extract){
  #extract common probes
  dataset_subset<-dataset[colnames(dataset) %in% probes_to_extract]
  #yugene
  dataset_YT<-YuGene(t(dataset_subset[7:ncol(dataset_subset)]))
  #dataframe
  dataset_YT<-as.data.frame(t(dataset_YT[,]))
  #add diagnosis
  dataset_YT<-merge(dataset[1], dataset_YT, by="row.names")
  rownames(dataset_YT)<-dataset_YT$Row.names
  dataset_YT$Row.names<-NULL
  #print number of rows
  cat(c("number of genes: "), ncol(dataset_YT)-1)
  #return
  return(dataset_YT)
}

# run

AD_E_GEOD_63060_YT<-yugene_transform(AD_E_GEOD_63060_full, RD_genes)
AD_E_GEOD_63061_YT<-yugene_transform(AD_E_GEOD_63061_full, RD_genes)
AD_E_GEOD_6613_YT<-yugene_transform(AD_E_GEOD_6613_full, RD_genes)


PD_E_GEOD_6613_YT<-yugene_transform(PD_E_GEOD_6613_full, RD_genes)
PD_E_GEOD_72267_YT<-yugene_transform(PD_E_GEOD_72267_full, RD_genes)

MS_E_GEOD_24427_YT<-yugene_transform(MS_E_GEOD_24427_full, RD_genes)
MS_E_GEOD_16214_YT<-yugene_transform(MS_E_GEOD_16214_full, RD_genes)
MS_E_GEOD_41890_YT<-yugene_transform(MS_E_GEOD_41890_full, RD_genes)

SCZ_E_GEOD_38484_YT<-yugene_transform(SCZ_E_GEOD_38484_full, RD_genes)
SCZ_E_GEOD_27383_YT<-yugene_transform(SCZ_E_GEOD_27383_full, RD_genes)
SCZ_E_GEOD_38481_YT<-yugene_transform(SCZ_E_GEOD_38481_full, RD_genes)

BP_E_GEOD_46449_YT<-yugene_transform(BP_E_GEOD_46449_full, RD_genes)
BP_E_GEOD_23848_YT<-yugene_transform(BP_E_GEOD_23848_full, RD_genes)

CD_E_GEOD_46097_YT<-yugene_transform(CD_E_GEOD_46097_full, RD_genes)
CD_E_GEOD_59867_YT<-yugene_transform(CD_E_GEOD_59867_full, RD_genes)
CD_E_GEOD_12288_YT<-yugene_transform(CD_E_GEOD_12288_full, RD_genes)

RA_E_GEOD_74143_YT<-yugene_transform(RA_E_GEOD_74143_full, RD_genes)
RA_E_GEOD_54629_YT<-yugene_transform(RA_E_GEOD_54629_full, RD_genes)
RA_E_GEOD_42296_YT<-yugene_transform(RA_E_GEOD_42296_full, RD_genes)

COPD_E_GEOD_54837_YT<-yugene_transform(COPD_E_GEOD_54837_full, RD_genes)
COPD_E_GEOD_42057_YT<-yugene_transform(COPD_E_GEOD_42057_full, RD_genes)

ALS_E_TABM_940_YT<-yugene_transform(ALS_E_TABM_940_full, RD_genes)
# save

setwd(Yugene_transformed_data_dir)

saveRDS(AD_E_GEOD_63060_YT, file="AD_E_GEOD_63060_YT.RDS")
saveRDS(AD_E_GEOD_63061_YT, file="AD_E_GEOD_63061_YT.RDS")
saveRDS(AD_E_GEOD_6613_YT, file="AD_E_GEOD_6613_YT.RDS")

saveRDS(PD_E_GEOD_6613_YT, file="PD_E_GEOD_6613_YT.RDS")
saveRDS(PD_E_GEOD_72267_YT, file="PD_E_GEOD_72267_YT.RDS")
saveRDS(MS_E_GEOD_24427_YT, file="MS_E_GEOD_24427_YT.RDS")
saveRDS(MS_E_GEOD_16214_YT, file="MS_E_GEOD_16214_YT.RDS")
saveRDS(MS_E_GEOD_41890_YT, file="MS_E_GEOD_41890_YT.RDS")
saveRDS(SCZ_E_GEOD_38484_YT, file="SCZ_E_GEOD_38484_YT.RDS")
saveRDS(SCZ_E_GEOD_27383_YT, file="SCZ_E_GEOD_27383_YT.RDS")
saveRDS(SCZ_E_GEOD_38481_YT, file="SCZ_E_GEOD_38481_YT.RDS")
saveRDS(BP_E_GEOD_46449_YT, file="BP_E_GEOD_46449_YT.RDS")
saveRDS(BP_E_GEOD_23848_YT, file="BP_E_GEOD_23848_YT.RDS")
saveRDS(CD_E_GEOD_46097_YT, file="CD_E_GEOD_46097_YT.RDS")
saveRDS(CD_E_GEOD_59867_YT, file="CD_E_GEOD_59867_YT.RDS")
saveRDS(CD_E_GEOD_12288_YT, file="CD_E_GEOD_12288_YT.RDS")
saveRDS(RA_E_GEOD_74143_YT, file="RA_E_GEOD_74143_YT.RDS")
saveRDS(RA_E_GEOD_54629_YT, file="RA_E_GEOD_54629_YT.RDS")
saveRDS(RA_E_GEOD_42296_YT, file="RA_E_GEOD_42296_YT.RDS")
saveRDS(COPD_E_GEOD_54837_YT, file="COPD_E_GEOD_54837_YT.RDS")
saveRDS(COPD_E_GEOD_42057_YT, file="COPD_E_GEOD_42057_YT.RDS")
saveRDS(ALS_E_TABM_940_YT, file="ALS_E_TABM_940_YT.RDS")

##### CREATE SAMPLE-DISEASE MAPPING FILE #####

# change diagnosis
AD_E_GEOD_63060_YT[AD_E_GEOD_63060_YT$Diagnosis=="case",1]<-"AD_case"
AD_E_GEOD_63061_YT[AD_E_GEOD_63061_YT$Diagnosis=="case",1]<-"AD_case"
AD_E_GEOD_6613_YT[AD_E_GEOD_6613_YT$Diagnosis=="case",1]<-"AD_case"
PD_E_GEOD_6613_YT[PD_E_GEOD_6613_YT$Diagnosis=="case",1]<-"PD_case"
PD_E_GEOD_72267_YT[PD_E_GEOD_72267_YT$Diagnosis=="case",1]<-"PD_case"
MS_E_GEOD_24427_YT[MS_E_GEOD_24427_YT$Diagnosis=="case",1]<-"MS_case"
MS_E_GEOD_16214_YT[MS_E_GEOD_16214_YT$Diagnosis=="case",1]<-"MS_case"
MS_E_GEOD_41890_YT[MS_E_GEOD_41890_YT$Diagnosis=="case",1]<-"MS_case"
SCZ_E_GEOD_38484_YT[SCZ_E_GEOD_38484_YT$Diagnosis=="case",1]<-"SCZ_case"
SCZ_E_GEOD_27383_YT[SCZ_E_GEOD_27383_YT$Diagnosis=="case",1]<-"SCZ_case"
SCZ_E_GEOD_38481_YT[SCZ_E_GEOD_38481_YT$Diagnosis=="case",1]<-"SCZ_case"
BP_E_GEOD_46449_YT[BP_E_GEOD_46449_YT$Diagnosis=="case",1]<-"BP_case"
BP_E_GEOD_23848_YT[BP_E_GEOD_23848_YT$Diagnosis=="case",1]<-"BP_case"
CD_E_GEOD_46097_YT[CD_E_GEOD_46097_YT$Diagnosis=="case",1]<-"CD_case"
CD_E_GEOD_59867_YT[CD_E_GEOD_59867_YT$Diagnosis=="case",1]<-"CD_case"
CD_E_GEOD_12288_YT[CD_E_GEOD_12288_YT$Diagnosis=="case",1]<-"CD_case"
RA_E_GEOD_74143_YT[RA_E_GEOD_74143_YT$Diagnosis=="case",1]<-"RA_case"
RA_E_GEOD_54629_YT[RA_E_GEOD_54629_YT$Diagnosis=="case",1]<-"RA_case"
RA_E_GEOD_42296_YT[RA_E_GEOD_42296_YT$Diagnosis=="case",1]<-"RA_case"
COPD_E_GEOD_54837_YT[COPD_E_GEOD_54837_YT$Diagnosis=="case",1]<-"COPD_case"
COPD_E_GEOD_42057_YT[COPD_E_GEOD_42057_YT$Diagnosis=="case",1]<-"COPD_case"
ALS_E_TABM_940_YT[ALS_E_TABM_940_YT$Diagnosis=="case",1]<-"ALS_case"


AD_E_GEOD_63060_YT[AD_E_GEOD_63060_YT$Diagnosis=="control",1]<-"AD_control"
AD_E_GEOD_63061_YT[AD_E_GEOD_63061_YT$Diagnosis=="control",1]<-"AD_control"
AD_E_GEOD_6613_YT[AD_E_GEOD_6613_YT$Diagnosis=="control",1]<-"AD_control"

PD_E_GEOD_6613_YT[PD_E_GEOD_6613_YT$Diagnosis=="control",1]<-"PD_control"
PD_E_GEOD_72267_YT[PD_E_GEOD_72267_YT$Diagnosis=="control",1]<-"PD_control"
MS_E_GEOD_24427_YT[MS_E_GEOD_24427_YT$Diagnosis=="control",1]<-"MS_control"
MS_E_GEOD_16214_YT[MS_E_GEOD_16214_YT$Diagnosis=="control",1]<-"MS_control"
MS_E_GEOD_41890_YT[MS_E_GEOD_41890_YT$Diagnosis=="control",1]<-"MS_control"
SCZ_E_GEOD_38484_YT[SCZ_E_GEOD_38484_YT$Diagnosis=="control",1]<-"SCZ_control"
SCZ_E_GEOD_27383_YT[SCZ_E_GEOD_27383_YT$Diagnosis=="control",1]<-"SCZ_control"
SCZ_E_GEOD_38481_YT[SCZ_E_GEOD_38481_YT$Diagnosis=="control",1]<-"SCZ_control"
BP_E_GEOD_46449_YT[BP_E_GEOD_46449_YT$Diagnosis=="control",1]<-"BP_control"
BP_E_GEOD_23848_YT[BP_E_GEOD_23848_YT$Diagnosis=="control",1]<-"BP_control"
CD_E_GEOD_46097_YT[CD_E_GEOD_46097_YT$Diagnosis=="control",1]<-"CD_control"
CD_E_GEOD_59867_YT[CD_E_GEOD_59867_YT$Diagnosis=="control",1]<-"CD_control"
CD_E_GEOD_12288_YT[CD_E_GEOD_12288_YT$Diagnosis=="control",1]<-"CD_control"
RA_E_GEOD_74143_YT[RA_E_GEOD_74143_YT$Diagnosis=="control",1]<-"RA_control"
RA_E_GEOD_54629_YT[RA_E_GEOD_54629_YT$Diagnosis=="control",1]<-"RA_control"
RA_E_GEOD_42296_YT[RA_E_GEOD_42296_YT$Diagnosis=="control",1]<-"RA_control"
COPD_E_GEOD_54837_YT[COPD_E_GEOD_54837_YT$Diagnosis=="control",1]<-"COPD_control"
COPD_E_GEOD_42057_YT[COPD_E_GEOD_42057_YT$Diagnosis=="control",1]<-"COPD_control"
ALS_E_TABM_940_YT[ALS_E_TABM_940_YT$Diagnosis=="control",1]<-"ALS_control"

# check

table(AD_E_GEOD_63060_YT$Diagnosis)
table(AD_E_GEOD_63061_YT$Diagnosis)
table(AD_E_GEOD_6613_YT$Diagnosis)
table(PD_E_GEOD_6613_YT$Diagnosis)
table(PD_E_GEOD_72267_YT$Diagnosis)
table(MS_E_GEOD_24427_YT$Diagnosis)
table(MS_E_GEOD_16214_YT$Diagnosis)
table(MS_E_GEOD_41890_YT$Diagnosis)
table(SCZ_E_GEOD_38484_YT$Diagnosis)
table(SCZ_E_GEOD_27383_YT$Diagnosis)
table(SCZ_E_GEOD_38481_YT$Diagnosis)
table(BP_E_GEOD_46449_YT$Diagnosis)
table(BP_E_GEOD_23848_YT$Diagnosis)
table(CD_E_GEOD_46097_YT$Diagnosis)
table(CD_E_GEOD_59867_YT$Diagnosis)
table(CD_E_GEOD_12288_YT$Diagnosis)
table(RA_E_GEOD_74143_YT$Diagnosis)
table(RA_E_GEOD_54629_YT$Diagnosis)
table(RA_E_GEOD_42296_YT$Diagnosis)
table(COPD_E_GEOD_54837_YT$Diagnosis)
table(COPD_E_GEOD_42057_YT$Diagnosis)
table(ALS_E_TABM_940_YT$Diagnosis)

# add Dataset

AD_E_GEOD_63060_YT$Dataset<-"AD_GSE63060"
AD_E_GEOD_63061_YT$Dataset<-"AD_GSE63061"
AD_E_GEOD_6613_YT$Dataset<-"AD_E_GEOD_6613"
PD_E_GEOD_6613_YT$Dataset<-"PD_E_GEOD_6613"
PD_E_GEOD_72267_YT$Dataset<-"PD_E_GEOD_72267"
MS_E_GEOD_24427_YT$Dataset<-"MS_GSE24427"
MS_E_GEOD_16214_YT$Dataset<-"MS_E_GEOD_16214"
MS_E_GEOD_41890_YT$Dataset<-"MS_E_GEOD_41890"
SCZ_E_GEOD_38484_YT$Dataset<-"SCZ_GSE38484"
SCZ_E_GEOD_27383_YT$Dataset<-"SCZ_E_GEOD_27383"
SCZ_E_GEOD_38481_YT$Dataset<-"SCZ_GSE38481"
BP_E_GEOD_46449_YT$Dataset<-"BP_E_GEOD_46449"
BP_E_GEOD_23848_YT$Dataset<-"BP_GSE23848"
CD_E_GEOD_46097_YT$Dataset<-"CD_E_GEOD_46097"
CD_E_GEOD_59867_YT$Dataset<-"CD_GSE59867"
CD_E_GEOD_12288_YT$Dataset<-"CD_E_GEOD_12288"
RA_E_GEOD_74143_YT$Dataset<-"RA_E_GEOD_74143"
RA_E_GEOD_54629_YT$Dataset<-"RA_E_GEOD_54629"
RA_E_GEOD_42296_YT$Dataset<-"RA_E_GEOD_42296"
COPD_E_GEOD_54837_YT$Dataset<-"COPD_E_GEOD_54837"
COPD_E_GEOD_42057_YT$Dataset<-"COPD_E_GEOD_42057"
ALS_E_TABM_940_YT$Dataset<-"ALS_E_TABM_940"

# reorder column names
AD_E_GEOD_63060_YT<-AD_E_GEOD_63060_YT[c(1,1683,2:1682)]
AD_E_GEOD_63061_YT<-AD_E_GEOD_63061_YT[c(1,1683,2:1682)]
AD_E_GEOD_6613_YT<-AD_E_GEOD_6613_YT[c(1,1683,2:1682)]
PD_E_GEOD_6613_YT<-PD_E_GEOD_6613_YT[c(1,1683,2:1682)]
PD_E_GEOD_72267_YT<-PD_E_GEOD_72267_YT[c(1,1683,2:1682)]
MS_E_GEOD_24427_YT<-MS_E_GEOD_24427_YT[c(1,1683,2:1682)]
MS_E_GEOD_16214_YT<-MS_E_GEOD_16214_YT[c(1,1683,2:1682)]
MS_E_GEOD_41890_YT<-MS_E_GEOD_41890_YT[c(1,1683,2:1682)]
SCZ_E_GEOD_38484_YT<-SCZ_E_GEOD_38484_YT[c(1,1683,2:1682)]
SCZ_E_GEOD_27383_YT<-SCZ_E_GEOD_27383_YT[c(1,1683,2:1682)]
SCZ_E_GEOD_38481_YT<-SCZ_E_GEOD_38481_YT[c(1,1683,2:1682)]
BP_E_GEOD_46449_YT<-BP_E_GEOD_46449_YT[c(1,1683,2:1682)]
BP_E_GEOD_23848_YT<-BP_E_GEOD_23848_YT[c(1,1683,2:1682)]
CD_E_GEOD_46097_YT<-CD_E_GEOD_46097_YT[c(1,1683,2:1682)]
CD_E_GEOD_59867_YT<-CD_E_GEOD_59867_YT[c(1,1683,2:1682)]
CD_E_GEOD_12288_YT<-CD_E_GEOD_12288_YT[c(1,1683,2:1682)]
RA_E_GEOD_74143_YT<-RA_E_GEOD_74143_YT[c(1,1683,2:1682)]
RA_E_GEOD_54629_YT<-RA_E_GEOD_54629_YT[c(1,1683,2:1682)]
RA_E_GEOD_42296_YT<-RA_E_GEOD_42296_YT[c(1,1683,2:1682)]
COPD_E_GEOD_54837_YT<-COPD_E_GEOD_54837_YT[c(1,1683,2:1682)]
COPD_E_GEOD_42057_YT<-COPD_E_GEOD_42057_YT[c(1,1683,2:1682)]
ALS_E_TABM_940_YT<-ALS_E_TABM_940_YT[c(1,1683,2:1682)]

#check

head(ALS_E_TABM_940_YT)[1:5]
head(RA_E_GEOD_54629_YT)[1:5]

##### MERGE DATA #####

# create merge function

MyMerge       <- function(x, y){
  df            <- rbind(x,y[,colnames(x)])
  return(df)
}


merged_data_YT <- Reduce(MyMerge, list(AD_E_GEOD_63060_YT,
                                       AD_E_GEOD_63061_YT,
                                       AD_E_GEOD_6613_YT,
                                       PD_E_GEOD_6613_YT,
                                       PD_E_GEOD_72267_YT,
                                       MS_E_GEOD_24427_YT,
                                       MS_E_GEOD_16214_YT,
                                       MS_E_GEOD_41890_YT,
                                       SCZ_E_GEOD_38484_YT,
                                       SCZ_E_GEOD_27383_YT,
                                       SCZ_E_GEOD_38481_YT,
                                       BP_E_GEOD_46449_YT,
                                       BP_E_GEOD_23848_YT,
                                       CD_E_GEOD_46097_YT,
                                       CD_E_GEOD_59867_YT,
                                       CD_E_GEOD_12288_YT,
                                       RA_E_GEOD_74143_YT,
                                       RA_E_GEOD_54629_YT,
                                       RA_E_GEOD_42296_YT,
                                       COPD_E_GEOD_54837_YT,
                                       COPD_E_GEOD_42057_YT,
                                       ALS_E_TABM_940_YT))


head(merged_data_YT)[1:5]

dim(table(merged_data_YT$Dataset))

# save

setwd(Yugene_transformed_data_dir)
saveRDS(merged_data_YT, file="YuGene_t_merged_data_common_rd.RDS")

##### BOXPLOT + DENSITY PLOT #####

#boxplot(t(merged_data_YT[3:ncol(merged_data_YT)]))

setwd(boxplots_dir)

#plot
png("Boxplot_of_merged_data_full.png")
boxplot(t(merged_data_YT[3:ncol(merged_data_YT)]), xaxt="n", ylab="Expression", xlab="Samples")
axis(1, at=c(seq(1,2740,400), 2740))
dev.off()

tiff("Boxplot_of_merged_data_full.tif", width = 8, height = 8, unit="in",  res=300)
boxplot(t(merged_data_YT[3:ncol(merged_data_YT)]), xaxt="n", ylab="Expression", xlab="Samples", main="(d) Extracting common reliably detected genes boxplot")
axis(1, at=c(seq(1,2740,400), 2740))
dev.off()

# density plot

#plotDensity(t(merged_data_YT[3:ncol(merged_data_YT)]), logMode=F, addLegend=F)

#plot
png("DensityPlot_of_merged_data_full.png")
plotDensity(t(merged_data_YT[3:ncol(merged_data_YT)]), logMode=F, addLegend=F)
dev.off()

tiff("DensityPlot_of_merged_data_full.tif",  width = 8, height = 8, unit="in",  res=300)
plotDensity(t(merged_data_YT[3:ncol(merged_data_YT)]), logMode=F, addLegend=F, main="(c) Extracting common reliably detected genes boxplot")
dev.off()

##### PCA PLOT #####

setwd(PCA_dir)

pca_plot<-function(data, legend_position, shift){
  #run PCA
  pca<-prcomp(data[3:ncol(data)])
  # order of samples in expression data
  pheno<-data[1:2]
  # plot variance explained
  plot(fviz_eig(pca))
  # match order
  Diagnosis_pca_colour<-labels2colors(as.character(pheno$Diagnosis))
  Dataset_pca_colour<-labels2colors(as.character(pheno$Dataset))
  # pca plot - Diagnosis
  par(mar=c(5, 6, 6, 12), xpd=TRUE)
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Diagnosis",col="black", pch=21,bg=Diagnosis_pca_colour, bty='L')
  legend(legend_position, inset=c(shift,0), legend=unique(pheno$Diagnosis), fill=unique(Diagnosis_pca_colour), title="Diagnosis")
  # pca plot - Dataset
  par(mar=c(5, 6, 6, 12), xpd=TRUE)
  plot(pca$rotation[,1:2], main=" PCA plot coloured by Dataset",col="black", pch=21,bg=Dataset_pca_colour)
  legend(legend_position, inset=c(shift,0), legend=unique(pheno$Dataset), fill=unique(Dataset_pca_colour), title="Dataset")
}

pca_plot(merged_data_YT[1:10], 'bottomright', -0.3)

pdf("merged_data_PCA_plot.pdf")
pca_plot(merged_data_YT, 'bottomright', -0.6)
dev.off()

##### GENDER SPECIFIC PLOTS #####

# create gender sample mapping file

Gender <- Reduce(MyMerge, list(AD_E_GEOD_63060_full[2:3],
                               AD_E_GEOD_63061_full[2:3],
                               AD_E_GEOD_6613_full[2:3],
                               PD_E_GEOD_6613_full[2:3],
                               PD_E_GEOD_72267_full[2:3],
                               MS_E_GEOD_24427_full[2:3],
                               MS_E_GEOD_16214_full[2:3],
                               MS_E_GEOD_41890_full[2:3],
                               SCZ_E_GEOD_38484_full[2:3],
                               SCZ_E_GEOD_27383_full[2:3],
                               SCZ_E_GEOD_38481_full[2:3],
                               BP_E_GEOD_46449_full[2:3],
                               BP_E_GEOD_23848_full[2:3],
                               CD_E_GEOD_46097_full[2:3],
                               CD_E_GEOD_59867_full[2:3],
                               CD_E_GEOD_12288_full[2:3],
                               RA_E_GEOD_74143_full[2:3],
                               RA_E_GEOD_54629_full[2:3],
                               RA_E_GEOD_42296_full[2:3],
                               COPD_E_GEOD_54837_full[2:3],
                               COPD_E_GEOD_42057_full[2:3],
                               ALS_E_TABM_940_full[2:3]))


table(Gender[1])
table(Gender[2])

head(Gender)

# merge with data

merged_data_YT_with_gender<-merge(Gender, merged_data_YT, by="row.names")
rownames(merged_data_YT_with_gender)<-merged_data_YT_with_gender$Row.names
merged_data_YT_with_gender$Row.names<-NULL

head(merged_data_YT_with_gender)[1:5]

# get gene symbol list for chip

Gene_symbols_probes <- mappedkeys(org.Hs.egSYMBOL)

# Convert to a list
Gene_symbols <- as.data.frame(org.Hs.egSYMBOL[Gene_symbols_probes])

head(Gene_symbols)
dim(Gene_symbols)

# housekeeping genes

#Expressed in females only
XIST_probe_ID<-subset(Gene_symbols, symbol=="XIST")
XIST_probe_ID

#Expressed in Males only
PRKY_probe_ID<-subset(Gene_symbols, symbol=="PRKY")
PRKY_probe_ID

RPS4Y1_probe_ID<-subset(Gene_symbols, symbol=="RPS4Y1")
RPS4Y1_probe_ID

KDM5D_probe_ID<-subset(Gene_symbols, symbol=="KDM5D")
KDM5D_probe_ID

# HK genes expressed in all cells + males + females
MKRN1_probe_ID<-subset(Gene_symbols, symbol=="MKRN1")
MKRN1_probe_ID

ADIPOR1_probe_ID<-subset(Gene_symbols, symbol=="ADIPOR1")
ADIPOR1_probe_ID

BNIP3L_probe_ID<-subset(Gene_symbols, symbol=="BNIP3L")
BNIP3L_probe_ID

# merge all genes 
gene_list<-rbind(XIST_probe_ID,
                 PRKY_probe_ID,
                 RPS4Y1_probe_ID,
                 KDM5D_probe_ID,
                 MKRN1_probe_ID,
                 ADIPOR1_probe_ID,
                 BNIP3L_probe_ID)

gene_list

# create table of genes and state if KH or gender specific
gene_table<-read.table(text =
                         "Gene Expressed_in
                       ADIPOR1 All
                       MKRN1 All
                       RPS4Y1 Males", header=T)
gene_table


#create function to plot
plot_gender_specific_genes<-function(Expression_table, Gender, genes_to_extract, threshold, boxplot_title){
  #extract gene of interest
  Expression_table_gene_check<-as.data.frame(t(Expression_table[rownames(Expression_table)%in% genes_to_extract$gene_id,]))
  # change colnames TO GENE SYMBOL using genes to extract file
  for (x in 1:dim(Expression_table_gene_check)[2]){
    colnames(Expression_table_gene_check)[x]<-gene_list[genes_to_extract$gene_id==colnames(Expression_table_gene_check)[x],2]
  }
  # add in gender information
  Expression_table_gene_check_gender<-merge(Gender, Expression_table_gene_check, by="row.names")
  rownames(Expression_table_gene_check_gender)<-Expression_table_gene_check_gender$Row.names
  Expression_table_gene_check_gender$Row.names<-NULL
  #melt dataframe for plot
  Expression_table_gene_check_gender_melt<-melt(Expression_table_gene_check_gender, by=Gender)
  # change variable colun from factor to character
  Expression_table_gene_check_gender_melt$variable<-as.character(Expression_table_gene_check_gender_melt$variable)
  # order dataframe by variable
  Expression_table_gene_check_gender_melt<-Expression_table_gene_check_gender_melt[order(Expression_table_gene_check_gender_melt$variable),]
  # calculate user defined percentie threshold
  sample_quantiles<-apply(Expression_table, 2, quantile, probs=threshold)
  # mean of used defined threshold across samples
  mean_threshold=mean(sample_quantiles)
  #plot
  plot1<-qplot(variable, value, colour=get(colnames(Gender)), data = Expression_table_gene_check_gender_melt, geom = c("boxplot", "jitter")) + 
    geom_hline(yintercept = mean_threshold) +
    ggtitle(boxplot_title) +
    theme(text = element_text(size=20), axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x="Gene",y="Expression", colour = colnames(Gender)) 
  # 2nd legend
  plot2<-tableGrob(gene_table, rows=NULL)
  # plot
  grid.arrange(plot1, plot2,
               nrow=2,
               heights=c(3,1))
}

plot_gender_specific_genes(t(merged_data_YT_with_gender[5:ncol(merged_data_YT_with_gender)]), merged_data_YT_with_gender[1], gene_list, 0, "test")

# plot

setwd(hk_gene_plots_dir)

png("hk_gene_plot.png")
plot_gender_specific_genes(t(merged_data_YT_with_gender[5:ncol(merged_data_YT_with_gender)]), merged_data_YT_with_gender[2], gene_list, 0, "test")
dev.off()

##### SAVE #####

setwd(work_dir)
save.image("YuGene_transform_rd_genes.Rdata")
