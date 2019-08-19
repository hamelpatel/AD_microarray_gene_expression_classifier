#!/usr/bin/env Rscript

##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                 AD- SPECIFIC CLASSIFCATION MODEL                                                       #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

# following https://www.analyticsvidhya.com/blog/2016/03/complete-guide-parameter-tuning-xgboost-with-codes-python/

# Note: 
# 
# FIXED WEIGHTS - UPSAMPLE INDIVIDUAL CONTROLS 
# repeating to account for randomness in upsampling
#

# V5.0
# # 1. tune nrounds with default parameters (eta=1, max_depth=6, gamma=0, min_child_weight=1, subsample=1, colsample_bytree=1) + RFE
# # 2. tune max_depth (2:20, increasing by 1)
# # 3. tune min_child_weight (1:10, increasing by 1)
# # 4. tune gamma (0:10, increasing by 0.5)
# # 5. tune subsample (0.5:1, increasing by 0.1)
# # 6. tune colsample_bytree (0.5:1, increasing by 0.1)
# # 7. tune regularization (alpha - L1) (0:1, increasing by 1) 0, 0.001, 0.01, 0.1, 1, 10, 100
# # 8. tune regularization (lambda - L2) 0, 0.001, 0.01, 0.1, 1, 10, 100
# # 9. tune eta (0.01:0.2, increasing by 0.01) - nrounds cranked up to 10,000
# 
#
# using logloss

# with RFE

#export PATH=/home/hamel/R/x86_64-pc-linux-gnu-library/3.4

print(sessionInfo())

.libPaths()

##### SET SEED #####
#seed=2

# allow external arguments
args = commandArgs(trailingOnly=TRUE)
seed=args[1]

##### LOAD LIBRARIES ####

library(xgboost)
library(ggplot2)
library(caret)

##### SET DIRECTORIES #####

#linux
#master_work_dir="/media/hamel/Workspace/Dropbox/Projects/AD-classification/6.XGBoost/test_upsampling_and_weights_repeated/AD_specific_classifier_logloss_FINAL"

#windows
#master_work_dir="D:\\Dropbox\\Projects\\AD-classification\\6.XGBoost\\test_upsampling_and_weights_repeated\\AD_specific_classifier_logloss_FINAL"

# Rosalind
master_work_dir="/users/k1223459/brc_scratch/AD_Classifcation/AD_vs_mixedcontrol/results"

#create folder with seed

#linux
dir.create(paste(master_work_dir, seed, sep="/"))
work_dir=paste(master_work_dir, seed, sep="/")

#setwd(work_dir)

      
#windows
#dir.create(paste(master_work_dir, seed, sep="\\"))
#work_dir=paste(master_work_dir, seed, sep="\\")

setwd(work_dir)

##### LOAD DATA #####

#linux
#data<-readRDS("/media/hamel/Workspace/Dropbox/Projects/AD-classification/5.YuGene_Transform/common_rd_genes/Yugene_transformed_data/YuGene_t_merged_data_common_rd.RDS")

#windows
#data<-readRDS("D:\\Dropbox\\Projects\\AD-classification\\5.YuGene_Transform\\common_rd_genes\\Yugene_transformed_data\\YuGene_t_merged_data_common_rd.RDS")

#Rosalind
data<-readRDS("/users/k1223459/brc_scratch/AD_Classifcation/Data/YuGene_t_merged_data_common_rd.RDS")

# pool controls
table(data$Diagnosis)
data[grep("control", data$Diagnosis),1]<-"control"
table(data$Diagnosis)

##### TRAINING SET #####

table(data$Dataset)
table(data$Diagnosis)

#number of datasets
length(unique(c(data$Dataset)))

dim(data)

# train
training_data<-data[data$Dataset %in% c("AD_GSE63060", 
                                        "AD_E_GEOD_6613",
                                        "PD_E_GEOD_6613",
                                        "MS_E_GEOD_16214",
                                        "MS_E_GEOD_41890",
                                        "SCZ_GSE38484",
                                        "BP_E_GEOD_46449",
                                        "CD_E_GEOD_46097",
                                        "CD_E_GEOD_12288",
                                        "RA_E_GEOD_74143",
                                        "COPD_E_GEOD_54837",
                                        "ALS_E_TABM_940"),]

table(training_data$Dataset)
table(training_data$Diagnosis)
head(training_data)[1:5]

# split into AD and control group

training_data_AD<-subset(training_data, training_data$Diagnosis=="AD_case")

training_data_controls<-subset(training_data, training_data$Diagnosis!="AD_case")

dim(training_data_AD)
table(training_data_AD$Diagnosis)

dim(training_data_controls)
table(training_data_controls$Diagnosis)

# upsample controls to biggest sample number - 702
set.seed(seed)
training_data_controls<-caret::upSample(training_data_controls, as.factor(training_data_controls$Diagnosis))
dim(training_data_controls)
table(training_data_controls$Diagnosis)

#remove unwanted columns
training_data_controls$Class<-NULL

#check column number the same
ncol(training_data_controls)==ncol(training_data)

#add AD samples back in

all(colnames(training_data_controls)==colnames(training_data_AD))

training_data<-rbind(training_data_controls,
                     training_data_AD)


#hotcode - for xgboost
training_data_pheno_xgb<-training_data$Diagnosis
table(training_data_pheno_xgb)
training_data_pheno_xgb[grep("AD_case", training_data_pheno_xgb)]<-0
training_data_pheno_xgb[training_data_pheno_xgb!="0"]<-1

table(training_data_pheno_xgb)

##### TESTING SET #####

# train
testing_data<-data[data$Dataset %in% c("AD_GSE63061",
                                       "PD_E_GEOD_72267",
                                       "MS_GSE24427",
                                       "SCZ_E_GEOD_27383",
                                       "SCZ_GSE38481",
                                       "BP_GSE23848",
                                       "CD_GSE59867",
                                       "RA_E_GEOD_54629",
                                       "RA_E_GEOD_42296",
                                       "COPD_E_GEOD_42057"),]


#pheno - for caret
testing_data_pheno_caret<-as.factor(testing_data$Diagnosis)

#hotcode - for xgboost
testing_data_pheno_xgb<-testing_data$Diagnosis
table(testing_data_pheno_xgb)
testing_data_pheno_xgb[grep("AD_case", testing_data_pheno_xgb)]<-0
testing_data_pheno_xgb[testing_data_pheno_xgb!="0"]<-1

table(testing_data_pheno_xgb)

##### CONVERT TO xgboost FORMAT #####

training_data_xgb<-xgb.DMatrix(as.matrix(training_data[3:ncol(training_data)]), label=training_data_pheno_xgb)
testing_data_xgb<-xgb.DMatrix(as.matrix(testing_data[3:ncol(testing_data)]), label=testing_data_pheno_xgb)

##### DEFAULT MODEL PARAMATERS #####

set.seed(seed)

# default paramaters

params <- list(
  booster = "gbtree",
  objective = "binary:logistic",
  eta=0.3,
  max_depth=6,
  gamma=0,
  min_child_weight=1,
  subsample=1,
  colsample_bytree=1,
  scale_pos_weight=table(training_data_pheno_xgb)[1]/table(training_data_pheno_xgb)[2] # Control the balance of positive and negative weights, useful for unbalanced classes. A typical value to consider: sum(negative cases) / sum(positive cases)
)

##### TUNE  1 - FIND OPTIMAL NROUNDS - DEFAULT PARAMETERS #####

# 1st find optimum number of rounds - this will be further tuned again at the end.
# nrounds is number of trees in the model

set.seed(seed)
tune_results_1 <- xgb.cv(params = params,
                data = training_data_xgb,
                nfold= 10,
                nrounds=5000, # increase if plot does not converge
                verbose=T,
                showsd=T, # show SD of error
                eval_metric="logloss",
                #eval_metric="logloss", #
                stratified=T, # Stratification is the process of rearranging the data as to ensure each fold is a good representative of the whole. For example in a binary classification problem where each class comprises 50% of the data, it is best to arrange the data such that in every fold, each class comprises around half the instances.stratification is generally a better scheme, both in terms of bias and variance, when compared to regular cross-validation.
                early_stopping_rounds = 20
)


# find lowest nround

optimal_nround<-tune_results_1$best_iteration
optimal_nround

##### MODEL 1 #####

set.seed(seed)
# create initial model - basic model 

xgb_model_1 <- xgb.train(
  params = params,
  data = training_data_xgb,
  nrounds = optimal_nround # set to what was found session
)

xgb_model_1

# save model
setwd(work_dir)
xgb.save(xgb_model_1, "xgb_model_1")

##### TEST MODEL 1 #####

#test on training set
Model_1_train_prediction<-as.numeric(predict(xgb_model_1, training_data_xgb) > 0.5)
Model_1_train_prediction_results<-caret::confusionMatrix(as.factor(Model_1_train_prediction), as.factor(training_data_pheno_xgb))
Model_1_train_prediction_results

#test on validation set
Model_1_test_prediction<-as.numeric(predict(xgb_model_1, testing_data_xgb) > 0.5)
Model_1_test_prediction_results<-caret::confusionMatrix(as.factor(Model_1_test_prediction), as.factor(testing_data_pheno_xgb))
Model_1_test_prediction_results

#missclassified
Model_1_test_prediction_missclassifcation<-testing_data[1:2]
Model_1_test_prediction_missclassifcation$Prediction<-Model_1_test_prediction
# missclassifcation - disease 
table(Model_1_test_prediction_missclassifcation$Prediction, Model_1_test_prediction_missclassifcation$Diagnosis)
# missclassifcation - dataset
table(Model_1_test_prediction_missclassifcation$Prediction, Model_1_test_prediction_missclassifcation$Dataset)


##### RFE #####

# calculate importance

model_1_importance <- xgb.importance(feature_names = names(training_data[3:ncol(training_data)]), model = xgb_model_1)

head(model_1_importance)
dim(model_1_importance)


# paramters for loop
min_predictor<-10 #(lowest number of features to go to)
max_predictor<-nrow(model_1_importance)
predictor_interval<-2
predictor_counter=1

# empty list for results
RFE_results<-list()

for (x in rev(seq(min_predictor, max_predictor, predictor_interval))) {
  # counter for interation
  y=predictor_counter
  # print round
  print(paste("Using: ", x, "Features"))
  # extract features
  data<-training_data[head(model_1_importance, x)$Feature]
  # convert to xgb
  data_xgb<-xgb.DMatrix(as.matrix(data), label=training_data_pheno_xgb)
  # run xgboost cv
  seeds = set.seed(seed)
  RFE_results[[y]] <- xgb.cv(params = params,
                             data = data_xgb,
                             nfold= 10,
                             nrounds=5000, 
                             verbose=T,
                             showsd=T,
                             eval_metric="logloss", 
                             stratified=T,
                             early_stopping_rounds = 20)
  # add name
  names(RFE_results)[y]<-x
  # change y
  predictor_counter=predictor_counter+1
  
}

head(RFE_results)
names(RFE_results)


#find lowest logloss
RFE_results_lowest_logloss<-as.data.frame(unlist(lapply(RFE_results, function(x) min(x$evaluation_log[,4]))))
RFE_results_lowest_logloss

#nrounds
RFE_results_lowest_logloss_nr<-as.data.frame(unlist(lapply(RFE_results, function(x) x$best_iteration)))
RFE_results_lowest_logloss_nr

# merge with nrounds
RFE_results_lowest_logloss<-merge(RFE_results_lowest_logloss, RFE_results_lowest_logloss_nr, by="row.names")

#name
names(RFE_results_lowest_logloss)<-c("Features", "test_mlogloss_mean", "nrounds")
RFE_results_lowest_logloss$Features<-as.numeric(RFE_results_lowest_logloss$Features)

#order
RFE_results_lowest_logloss<-RFE_results_lowest_logloss[order(RFE_results_lowest_logloss$test_mlogloss_mean, RFE_results_lowest_logloss$Features),]
RFE_results_lowest_logloss

# create new training + testing data with optimum feature

RFE_training_data<-training_data[head(model_1_importance, as.numeric(RFE_results_lowest_logloss$Features[1]))$Feature]
RFE_training_data_xgb<-xgb.DMatrix(as.matrix(RFE_training_data), label=training_data_pheno_xgb)

RFE_testing_data<-testing_data[head(model_1_importance, as.numeric(RFE_results_lowest_logloss$Features[1]))$Feature]
RFE_testing_data_xgb<-xgb.DMatrix(as.matrix(RFE_testing_data), label=testing_data_pheno_xgb)

#otimum nrounds
optimal_nround_rfe<-as.numeric(RFE_results_lowest_logloss$nrounds[1])

##### MODEL 1 - RFE #####

set.seed(seed)
# create initial model - basic model 

xgb_model_1_rfe <- xgboost(
  params = params,
  data = RFE_training_data_xgb,
  nrounds = optimal_nround_rfe # set to what was found session
)

xgb_model_1_rfe

# save model
setwd(work_dir)
xgb.save(xgb_model_1_rfe, "xgb_model_1_rfe")

##### TEST MODEL 1 - rfe #####

#test on training set
Model_1_train_prediction_rfe<-as.numeric(predict(xgb_model_1_rfe, RFE_training_data_xgb) > 0.5)
Model_1_train_prediction_results_rfe<-caret::confusionMatrix(as.factor(Model_1_train_prediction_rfe), as.factor(training_data_pheno_xgb))
Model_1_train_prediction_results_rfe

#test on validation set
Model_1_test_prediction_rfe<-as.numeric(predict(xgb_model_1_rfe, RFE_testing_data_xgb) > 0.5)
Model_1_test_prediction_results_rfe<-caret::confusionMatrix(as.factor(Model_1_test_prediction_rfe), as.factor(testing_data_pheno_xgb))
Model_1_test_prediction_results_rfe

#missclassified
Model_1_test_prediction_missclassifcation_rfe<-testing_data[1:2]
Model_1_test_prediction_missclassifcation_rfe$Prediction<-Model_1_test_prediction_rfe
# missclassifcation - disease 
table(Model_1_test_prediction_missclassifcation_rfe$Prediction, Model_1_test_prediction_missclassifcation_rfe$Diagnosis)
# missclassifcation - dataset
table(Model_1_test_prediction_missclassifcation_rfe$Prediction, Model_1_test_prediction_missclassifcation_rfe$Dataset)


##### TUNE  2 - max_depth  #####

# max_depth:
# max depth of trees. used to control over-fitting. higher value allows model to learn relations specific to a sample

# copy default param
param2<-params

# paramters for loop
minimum_max_depth<-2
maximum_max_depth<-20
interval_max_depth<-1
max_depth_counter=1

# empty list for results
tune_results_2<-list()

for (x in seq(minimum_max_depth, maximum_max_depth, interval_max_depth)) {
  # counter for interation
  y=max_depth_counter
  # print round
  print(paste("round: ", y))
  # change reg param to x
  param2$max_depth<-x
  # run xgboost cv
  seeds = set.seed(seed)
  tune_results_2[[y]] <- xgb.cv(params = param2,
                                data = RFE_training_data_xgb,
                                nfold= 10,
                                nrounds=5000, 
                                verbose=T,
                                showsd=T,
                                eval_metric="logloss", 
                                stratified=T,
                                early_stopping_rounds = 20)
  # add name
  names(tune_results_2)[y]<-x
  # change y
  max_depth_counter=max_depth_counter+1
  
}

head(tune_results_2)

#find lowest logloss
tune_results_2_lowest_logloss<-as.data.frame(unlist(lapply(tune_results_2, function(x) min(x$evaluation_log[,4]))))
tune_results_2_lowest_logloss

#nrounds
tune_results_2_lowest_logloss_nr<-as.data.frame(unlist(lapply(tune_results_2, function(x) x$best_iteration)))
tune_results_2_lowest_logloss_nr

# merge with nrounds
tune_results_2_lowest_logloss<-merge(tune_results_2_lowest_logloss, tune_results_2_lowest_logloss_nr, by="row.names")

#name
names(tune_results_2_lowest_logloss)<-c("max_depth", "test_mlogloss_mean", "nrounds")
#convert to numeric
tune_results_2_lowest_logloss$max_depth<-as.numeric(tune_results_2_lowest_logloss$max_depth)
#order
tune_results_2_lowest_logloss<-tune_results_2_lowest_logloss[order(tune_results_2_lowest_logloss$test_mlogloss_mean, tune_results_2_lowest_logloss$max_depth),]
tune_results_2_lowest_logloss

##### MODEL 2 #####

# change param2 again
param2$max_depth<-as.numeric(tune_results_2_lowest_logloss$max_depth[1])
param2

set.seed(seed)
# create initial model - basic model 

xgb_model_2 <- xgboost(
  params = param2,
  data = RFE_training_data_xgb,
  nrounds = tune_results_2_lowest_logloss$nrounds[1] # set to what was found session
)

xgb_model_2

# save
setwd(work_dir)
xgb.save(xgb_model_2, "xgb_model_2")

##### TEST MODEL 2 #####

#test on training set
Model_2_train_prediction<-as.numeric(predict(xgb_model_2, RFE_training_data_xgb) > 0.5)
Model_2_train_prediction_results<-caret::confusionMatrix(as.factor(Model_2_train_prediction), as.factor(training_data_pheno_xgb))
Model_2_train_prediction_results

#test on validation set
Model_2_test_prediction<-as.numeric(predict(xgb_model_2, RFE_testing_data_xgb) > 0.5)
Model_2_test_prediction_results<-caret::confusionMatrix(as.factor(Model_2_test_prediction), as.factor(testing_data_pheno_xgb))
Model_2_test_prediction_results

#missclassified
Model_2_test_prediction_missclassifcation<-testing_data[1:2]
Model_2_test_prediction_missclassifcation$Prediction<-Model_2_test_prediction
# missclassifcation - disease 
table(Model_2_test_prediction_missclassifcation$Prediction, Model_2_test_prediction_missclassifcation$Diagnosis)
# missclassifcation - dataset
table(Model_2_test_prediction_missclassifcation$Prediction, Model_2_test_prediction_missclassifcation$Dataset)

##### TUNE  3 - min_child_weight  #####

# max_child:
# Defines the maximum sum of weights of all observations required in a child. Used to control over-fitting. Higher values prevent a model from learning relations which might be highly specific to the particular sample selected for a tree. 
# recommened to start with 1/sqrt(event rate)

# copy default param
param3<-param2
param3

# paramters for loop
minimum_min_child_weight<-1
maximum_min_child_weight<-10
interval_min_child_weight<-1
min_child_weight_counter=1

# empty list for results
tune_results_3<-list()

for (x in seq(minimum_min_child_weight, maximum_min_child_weight, interval_min_child_weight)) {
  # counter for interation
  y=min_child_weight_counter
  # print round
  print(paste("round: ", y))
  # change reg param to x
  param3$min_child_weight<-x
  # run xgboost cv
  seeds = set.seed(seed)
  tune_results_3[[y]] <- xgb.cv(params = param3,
                                data = RFE_training_data_xgb,
                                nfold= 10,
                                nrounds=5000, 
                                verbose=T,
                                showsd=T,
                                eval_metric="logloss", 
                                stratified=T,
                                early_stopping_rounds = 20)
  # add name
  names(tune_results_3)[y]<-x
  # change y
  min_child_weight_counter=min_child_weight_counter+1
  
}

head(tune_results_3)

#find lowest logloss
tune_results_3_lowest_logloss<-as.data.frame(unlist(lapply(tune_results_3, function(x) min(x$evaluation_log[,4]))))
tune_results_3_lowest_logloss

#nrounds
tune_results_3_lowest_logloss_nr<-as.data.frame(unlist(lapply(tune_results_3, function(x) x$best_iteration)))
tune_results_3_lowest_logloss_nr

# merge with nrounds
tune_results_3_lowest_logloss<-merge(tune_results_3_lowest_logloss, tune_results_3_lowest_logloss_nr, by="row.names")

#name
names(tune_results_3_lowest_logloss)<-c("min_child_weight", "test_mlogloss_mean", "nrounds")
#convert to numeric
tune_results_3_lowest_logloss$min_child_weight<-as.numeric(tune_results_3_lowest_logloss$min_child_weight)

#order
tune_results_3_lowest_logloss<-tune_results_3_lowest_logloss[order(tune_results_3_lowest_logloss$test_mlogloss_mean, tune_results_3_lowest_logloss$min_child_weight),]
tune_results_3_lowest_logloss

##### MODEL 3 #####

# change param3 again
param3$min_child_weight<-as.numeric(tune_results_3_lowest_logloss$min_child_weight[1])
param3

set.seed(seed)
# create initial model - basic model 

xgb_model_3 <- xgboost(
  params = param3,
  data = RFE_training_data_xgb,
  nrounds = tune_results_3_lowest_logloss$nrounds[1] # set to what was found session
)

xgb_model_3

# save
setwd(work_dir)
xgb.save(xgb_model_3, "xgb_model_3")

##### TEST MODEL 3 #####

#test on training set
Model_3_train_prediction<-as.numeric(predict(xgb_model_3, RFE_training_data_xgb) > 0.5)
Model_3_train_prediction_results<-caret::confusionMatrix(as.factor(Model_3_train_prediction), as.factor(training_data_pheno_xgb))
Model_3_train_prediction_results

#test on validation set
Model_3_test_prediction<-as.numeric(predict(xgb_model_3, RFE_testing_data_xgb) > 0.5)
Model_3_test_prediction_results<-caret::confusionMatrix(as.factor(Model_3_test_prediction), as.factor(testing_data_pheno_xgb))
Model_3_test_prediction_results

#missclassified
Model_3_test_prediction_missclassifcation<-testing_data[1:2]
Model_3_test_prediction_missclassifcation$Prediction<-Model_3_test_prediction
# missclassifcation - disease 
table(Model_3_test_prediction_missclassifcation$Prediction, Model_3_test_prediction_missclassifcation$Diagnosis)
# missclassifcation - dataset
table(Model_3_test_prediction_missclassifcation$Prediction, Model_3_test_prediction_missclassifcation$Dataset)

##### TUNE  4 - gamma  #####

# Gamma :
# A node is split only when the resulting split gives a positive reduction in the loss function. Gamma specifies the maximum loss reduction required to make a split.
# Makes the algorithm conservative. The values can vary depending on the loss function and should be tuned.
# It controls regularization (or prevents overfitting).
#Tune trick: Start with 0 and check CV error rate. If you see train error >>> test error, bring gamma into action. 
# Higher the gamma, lower the difference in train and test CV. If you have no clue what value to use, use gamma=5 and see the performance. 
# Remember that gamma brings improvement when you want to use shallow (low max_depth) trees.

# copy default param
param4<-param3
param4

# paramters for loop
minimum_gamma<-0
maximum_gamma<-10
interval_gamma<-0.25
gamma_counter=1

# empty list for results
tune_results_4<-list()

for (x in seq(minimum_gamma, maximum_gamma, interval_gamma)) {
  # counter for interation
  y=gamma_counter
  # print round
  print(paste("round: ", y))
  # change reg param to x
  param4$gamma<-x
  # run xgboost cv
  seeds = set.seed(seed)
  tune_results_4[[y]] <- xgb.cv(params = param4,
                                data = RFE_training_data_xgb,
                                nfold= 10,
                                nrounds=5000, 
                                verbose=T,
                                showsd=T,
                                eval_metric="logloss", 
                                stratified=T,
                                early_stopping_rounds = 20)
  # add name
  names(tune_results_4)[y]<-x
  # change y
  gamma_counter=gamma_counter+1
  
}

head(tune_results_4)

#find lowest logloss
tune_results_4_lowest_logloss<-as.data.frame(unlist(lapply(tune_results_4, function(x) min(x$evaluation_log[,4]))))
tune_results_4_lowest_logloss

#nrounds
tune_results_4_lowest_logloss_nr<-as.data.frame(unlist(lapply(tune_results_4, function(x) x$best_iteration)))
tune_results_4_lowest_logloss_nr

# merge with nrounds
tune_results_4_lowest_logloss<-merge(tune_results_4_lowest_logloss, tune_results_4_lowest_logloss_nr, by="row.names")

#name
names(tune_results_4_lowest_logloss)<-c("gamma", "test_mlogloss_mean", "nrounds")
#order
tune_results_4_lowest_logloss<-tune_results_4_lowest_logloss[order(tune_results_4_lowest_logloss$test_mlogloss_mean),]
tune_results_4_lowest_logloss

##### MODEL 4 #####

# change param4 again
param4$gamma<-as.numeric(tune_results_4_lowest_logloss$gamma[1])
param4

set.seed(seed)
# create initial model - basic model 

xgb_model_4 <- xgboost(
  params = param4,
  data = RFE_training_data_xgb,
  nrounds = tune_results_4_lowest_logloss$nrounds[1] # set to what was found session
)

xgb_model_4

# save
setwd(work_dir)
xgb.save(xgb_model_4, "xgb_model_4")

##### TEST MODEL 4 #####

#test on training set
Model_4_train_prediction<-as.numeric(predict(xgb_model_4, RFE_training_data_xgb) > 0.5)
Model_4_train_prediction_results<-caret::confusionMatrix(as.factor(Model_4_train_prediction), as.factor(training_data_pheno_xgb))
Model_4_train_prediction_results

#test on validation set
Model_4_test_prediction<-as.numeric(predict(xgb_model_4, RFE_testing_data_xgb) > 0.5)
Model_4_test_prediction_results<-caret::confusionMatrix(as.factor(Model_4_test_prediction), as.factor(testing_data_pheno_xgb))
Model_4_test_prediction_results

#missclassified
Model_4_test_prediction_missclassifcation<-testing_data[1:2]
Model_4_test_prediction_missclassifcation$Prediction<-Model_4_test_prediction
# missclassifcation - disease 
table(Model_4_test_prediction_missclassifcation$Prediction, Model_4_test_prediction_missclassifcation$Diagnosis)
# missclassifcation - dataset
table(Model_4_test_prediction_missclassifcation$Prediction, Model_4_test_prediction_missclassifcation$Dataset)

##### TUNE  5 - subsample  #####

# subsample :
# Denotes the fraction of observations to be randomly samples for each tree. Typical values: 0.5-1

# copy default param
param5<-param4
param5

# paramters for loop
minimum_subsample<-0.5
maximum_subsample<-1
interval_subsample<-0.05
subsample_counter=1

# empty list for results
tune_results_5<-list()

for (x in seq(minimum_subsample, maximum_subsample, interval_subsample)) {
  # counter for interation
  y=subsample_counter
  # print round
  print(paste("round: ", y))
  # change reg param to x
  param5$subsample<-x
  # run xgboost cv
  seeds = set.seed(seed)
  tune_results_5[[y]] <- xgb.cv(params = param5,
                                data = RFE_training_data_xgb,
                                nfold= 10,
                                nrounds=5000, 
                                verbose=T,
                                showsd=T,
                                eval_metric="logloss", 
                                stratified=T,
                                early_stopping_rounds = 20)
  # add name
  names(tune_results_5)[y]<-x
  # change y
  subsample_counter=subsample_counter+1
  
}

head(tune_results_5)

#find lowest logloss
tune_results_5_lowest_logloss<-as.data.frame(unlist(lapply(tune_results_5, function(x) min(x$evaluation_log[,4]))))
tune_results_5_lowest_logloss

#nrounds
tune_results_5_lowest_logloss_nr<-as.data.frame(unlist(lapply(tune_results_5, function(x) x$best_iteration)))
tune_results_5_lowest_logloss_nr

# merge with nrounds
tune_results_5_lowest_logloss<-merge(tune_results_5_lowest_logloss, tune_results_5_lowest_logloss_nr, by="row.names")

#name
names(tune_results_5_lowest_logloss)<-c("subsample", "test_mlogloss_mean", "nrounds")
#convert to numeric
tune_results_5_lowest_logloss$subsample<-as.numeric(tune_results_5_lowest_logloss$subsample)
#order
tune_results_5_lowest_logloss<-tune_results_5_lowest_logloss[order(tune_results_5_lowest_logloss$test_mlogloss_mean, tune_results_5_lowest_logloss$subsample),]
tune_results_5_lowest_logloss

##### MODEL 5 #####

# change param5 again
param5$subsample<-as.numeric(tune_results_5_lowest_logloss$subsample[1])
param5

set.seed(seed)
# create initial model - basic model 

xgb_model_5 <- xgboost(
  params = param5,
  data = RFE_training_data_xgb,
  nrounds = tune_results_5_lowest_logloss$nrounds[1] # set to what was found session
)

xgb_model_5

# save
setwd(work_dir)
xgb.save(xgb_model_5, "xgb_model_5")

##### TEST MODEL 5 #####

#test on training set
Model_5_train_prediction<-as.numeric(predict(xgb_model_5, RFE_training_data_xgb) > 0.5)
Model_5_train_prediction_results<-caret::confusionMatrix(as.factor(Model_5_train_prediction), as.factor(training_data_pheno_xgb))
Model_5_train_prediction_results

#test on validation set
Model_5_test_prediction<-as.numeric(predict(xgb_model_5, RFE_testing_data_xgb) > 0.5)
Model_5_test_prediction_results<-caret::confusionMatrix(as.factor(Model_5_test_prediction), as.factor(testing_data_pheno_xgb))
Model_5_test_prediction_results

#missclassified
Model_5_test_prediction_missclassifcation<-testing_data[1:2]
Model_5_test_prediction_missclassifcation$Prediction<-Model_5_test_prediction
# missclassifcation - disease 
table(Model_5_test_prediction_missclassifcation$Prediction, Model_5_test_prediction_missclassifcation$Diagnosis)
# missclassifcation - dataset
table(Model_5_test_prediction_missclassifcation$Prediction, Model_5_test_prediction_missclassifcation$Dataset)

##### TUNE  6 - colsample_bytree  #####

# colsample_bytree:
# Denotes the fraction of columns to be randomly samples for each tree.Typical values: 0.5-1

# copy default param
param6<-param5
param6

# paramters for loop
minimum_colsample_bytree<-0.5
maximum_colsample_bytree<-1
interval_colsample_bytree<-0.05
colsample_bytree_counter=1

# empty list for results
tune_results_6<-list()

for (x in seq(minimum_colsample_bytree, maximum_colsample_bytree, interval_colsample_bytree)) {
  # counter for interation
  y=colsample_bytree_counter
  # print round
  print(paste("round: ", y))
  # change reg param to x
  param6$colsample_bytree<-x
  # run xgboost cv
  seeds = set.seed(seed)
  tune_results_6[[y]] <- xgb.cv(params = param6,
                                data = RFE_training_data_xgb,
                                nfold= 10,
                                nrounds=5000, 
                                verbose=T,
                                showsd=T,
                                eval_metric="logloss", 
                                stratified=T,
                                early_stopping_rounds = 20)
  # add name
  names(tune_results_6)[y]<-x
  # change y
  colsample_bytree_counter=colsample_bytree_counter+1
  
}

head(tune_results_6)

#find lowest logloss
tune_results_6_lowest_logloss<-as.data.frame(unlist(lapply(tune_results_6, function(x) min(x$evaluation_log[,4]))))
tune_results_6_lowest_logloss

#nrounds
tune_results_6_lowest_logloss_nr<-as.data.frame(unlist(lapply(tune_results_6, function(x) x$best_iteration)))
tune_results_6_lowest_logloss_nr

# merge with nrounds
tune_results_6_lowest_logloss<-merge(tune_results_6_lowest_logloss, tune_results_6_lowest_logloss_nr, by="row.names")

#name
names(tune_results_6_lowest_logloss)<-c("colsample_bytree", "test_mlogloss_mean", "nrounds")
#numeric
tune_results_6_lowest_logloss$colsample_bytree<-as.numeric(tune_results_6_lowest_logloss$colsample_bytree)
#order
tune_results_6_lowest_logloss<-tune_results_6_lowest_logloss[order(tune_results_6_lowest_logloss$test_mlogloss_mean, tune_results_6_lowest_logloss$colsample_bytree),]
tune_results_6_lowest_logloss

##### MODEL 6 #####

# change param6 again
param6$colsample_bytree<-as.numeric(tune_results_6_lowest_logloss$colsample_bytree[1])
param6

set.seed(seed)
# create initial model - basic model 

xgb_model_6 <- xgboost(
  params = param6,
  data = RFE_training_data_xgb,
  nrounds = tune_results_6_lowest_logloss$nrounds[1] # set to what was found session
)

xgb_model_6

# save
setwd(work_dir)
xgb.save(xgb_model_6, "xgb_model_6")

##### TEST MODEL 6 #####

#test on training set
Model_6_train_prediction<-as.numeric(predict(xgb_model_6, RFE_training_data_xgb) > 0.5)
Model_6_train_prediction_results<-caret::confusionMatrix(as.factor(Model_6_train_prediction), as.factor(training_data_pheno_xgb))
Model_6_train_prediction_results

#test on validation set
Model_6_test_prediction<-as.numeric(predict(xgb_model_6, RFE_testing_data_xgb) > 0.5)
Model_6_test_prediction_results<-caret::confusionMatrix(as.factor(Model_6_test_prediction), as.factor(testing_data_pheno_xgb))
Model_6_test_prediction_results

#missclassified
Model_6_test_prediction_missclassifcation<-testing_data[1:2]
Model_6_test_prediction_missclassifcation$Prediction<-Model_6_test_prediction
# missclassifcation - disease 
table(Model_6_test_prediction_missclassifcation$Prediction, Model_6_test_prediction_missclassifcation$Diagnosis)
# missclassifcation - dataset
table(Model_6_test_prediction_missclassifcation$Prediction, Model_6_test_prediction_missclassifcation$Dataset)

##### TUNE  7 - alpha  #####

# alpha:

# L1 regularization term on weight (analogous to Lasso regression)
# Can be used in case of very high dimensionality so that the algorithm runs faster when implemented


# copy default param
param7<-param6
param7

# paramters for loop
minimum_alpha<-0
maximum_alpha<-1
interval_alpha<-0.1
alpha_counter=1

# empty list for results
tune_results_7<-list()

for (x in seq(minimum_alpha, maximum_alpha, interval_alpha)) {
  # counter for interation
  y=alpha_counter
  # print round
  print(paste("round: ", y))
  # change reg param to x
  param7$alpha<-x
  # run xgboost cv
  seeds = set.seed(seed)
  tune_results_7[[y]] <- xgb.cv(params = param7,
                                data = RFE_training_data_xgb,
                                nfold= 10,
                                nrounds=5000, 
                                verbose=T,
                                showsd=T,
                                eval_metric="logloss", 
                                stratified=T,
                                early_stopping_rounds = 20)
  # add name
  names(tune_results_7)[y]<-x
  # change y
  alpha_counter=alpha_counter+1
  
}

head(tune_results_7)

#find lowest logloss
tune_results_7_lowest_logloss<-as.data.frame(unlist(lapply(tune_results_7, function(x) min(x$evaluation_log[,4]))))
tune_results_7_lowest_logloss

#nrounds
tune_results_7_lowest_logloss_nr<-as.data.frame(unlist(lapply(tune_results_7, function(x) x$best_iteration)))
tune_results_7_lowest_logloss_nr

# merge with nrounds
tune_results_7_lowest_logloss<-merge(tune_results_7_lowest_logloss, tune_results_7_lowest_logloss_nr, by="row.names")

#name
names(tune_results_7_lowest_logloss)<-c("alpha", "test_mlogloss_mean", "nrounds")
#numeric
tune_results_7_lowest_logloss$alpha<-as.numeric(tune_results_7_lowest_logloss$alpha)
#order
tune_results_7_lowest_logloss<-tune_results_7_lowest_logloss[order(tune_results_7_lowest_logloss$test_mlogloss_mean, tune_results_7_lowest_logloss$alpha),]
tune_results_7_lowest_logloss

##### MODEL 7 #####

# change param7 again
param7$alpha<-as.numeric(tune_results_7_lowest_logloss$alpha[1])
param7

set.seed(seed)
# create initial model - basic model 

xgb_model_7 <- xgboost(
  params = param7,
  data = RFE_training_data_xgb,
  nrounds = tune_results_7_lowest_logloss$nrounds[1] # set to what was found session
)

xgb_model_7

# save
setwd(work_dir)
xgb.save(xgb_model_7, "xgb_model_7")

##### TEST MODEL 7 #####

#test on training set
Model_7_train_prediction<-as.numeric(predict(xgb_model_7, RFE_training_data_xgb) > 0.5)
Model_7_train_prediction_results<-caret::confusionMatrix(as.factor(Model_7_train_prediction), as.factor(training_data_pheno_xgb))
Model_7_train_prediction_results

#test on validation set
Model_7_test_prediction<-as.numeric(predict(xgb_model_7, RFE_testing_data_xgb) > 0.5)
Model_7_test_prediction_results<-caret::confusionMatrix(as.factor(Model_7_test_prediction), as.factor(testing_data_pheno_xgb))
Model_7_test_prediction_results

#missclassified
Model_7_test_prediction_missclassifcation<-testing_data[1:2]
Model_7_test_prediction_missclassifcation$Prediction<-Model_7_test_prediction
# missclassifcation - disease 
table(Model_7_test_prediction_missclassifcation$Prediction, Model_7_test_prediction_missclassifcation$Diagnosis)
# missclassifcation - dataset
table(Model_7_test_prediction_missclassifcation$Prediction, Model_7_test_prediction_missclassifcation$Dataset)

##### TUNE  8 - lambda  #####

# lambda:

# L1 regularization term on weight (analogous to Lasso regression)
# Can be used in case of very high dimensionality so that the algorithm runs faster when implemented


# copy default param
param8<-param7
param8

# paramters for loop
minimum_lambda<-0
maximum_lambda<-1
interval_lambda<-0.1
lambda_counter=1

# empty list for results
tune_results_8<-list()

for (x in seq(minimum_lambda, maximum_lambda, interval_lambda)) {
  # counter for interation
  y=lambda_counter
  # print round
  print(paste("round: ", y))
  # change reg param to x
  param8$lambda<-x
  # run xgboost cv
  seeds = set.seed(seed)
  tune_results_8[[y]] <- xgb.cv(params = param8,
                                data = RFE_training_data_xgb,
                                nfold= 10,
                                nrounds=5000, 
                                verbose=T,
                                showsd=T,
                                eval_metric="logloss", 
                                stratified=T,
                                early_stopping_rounds = 20)
  # add name
  names(tune_results_8)[y]<-x
  # change y
  lambda_counter=lambda_counter+1
  
}

#find lowest logloss
tune_results_8_lowest_logloss<-as.data.frame(unlist(lapply(tune_results_8, function(x) min(x$evaluation_log[,4]))))
tune_results_8_lowest_logloss

#nrounds
tune_results_8_lowest_logloss_nr<-as.data.frame(unlist(lapply(tune_results_8, function(x) x$best_iteration)))
tune_results_8_lowest_logloss_nr

# merge with nrounds
tune_results_8_lowest_logloss<-merge(tune_results_8_lowest_logloss, tune_results_8_lowest_logloss_nr, by="row.names")

#name
names(tune_results_8_lowest_logloss)<-c("lambda", "test_mlogloss_mean", "nrounds")
#numeric
tune_results_8_lowest_logloss$lambda<-as.numeric(tune_results_8_lowest_logloss$lambda)
#order
tune_results_8_lowest_logloss<-tune_results_8_lowest_logloss[order(tune_results_8_lowest_logloss$test_mlogloss_mean, tune_results_8_lowest_logloss$lambda),]
tune_results_8_lowest_logloss

##### MODEL 8 #####

# change param8 again
param8$lambda<-as.numeric(tune_results_8_lowest_logloss$lambda[1])
param8

set.seed(seed)
# create initial model - basic model 

xgb_model_8 <- xgboost(
  params = param8,
  data = RFE_training_data_xgb,
  nrounds = tune_results_8_lowest_logloss$nrounds[1] # set to what was found session
)

xgb_model_8

# save
setwd(work_dir)
xgb.save(xgb_model_8, "xgb_model_8")

##### TEST MODEL 8 #####

#test on training set
Model_8_train_prediction<-as.numeric(predict(xgb_model_8, RFE_training_data_xgb) > 0.5)
Model_8_train_prediction_results<-caret::confusionMatrix(as.factor(Model_8_train_prediction), as.factor(training_data_pheno_xgb))
Model_8_train_prediction_results

#test on validation set
Model_8_test_prediction<-as.numeric(predict(xgb_model_8, RFE_testing_data_xgb) > 0.5)
Model_8_test_prediction_results<-caret::confusionMatrix(as.factor(Model_8_test_prediction), as.factor(testing_data_pheno_xgb))
Model_8_test_prediction_results

#missclassified
Model_8_test_prediction_missclassifcation<-testing_data[1:2]
Model_8_test_prediction_missclassifcation$Prediction<-Model_8_test_prediction
# missclassifcation - disease 
table(Model_8_test_prediction_missclassifcation$Prediction, Model_8_test_prediction_missclassifcation$Diagnosis)
# missclassifcation - dataset
table(Model_8_test_prediction_missclassifcation$Prediction, Model_8_test_prediction_missclassifcation$Dataset)


##### TUNE  9 - eta  #####

# eta:

# Makes the model more robust by shrinking the weights on each step
# Typical final values to be used: 0.01-0.2

# copy default param
param9<-param8
param9

# paramters for loop
minimum_eta<-0.01
maximum_eta<-0.2
interval_eta<-0.01
eta_counter=1

# empty list for results
tune_results_9<-list()

for (x in seq(minimum_eta, maximum_eta, interval_eta)) {
  # counter for interation
  y=eta_counter
  # print round
  print(paste("round: ", y))
  # change reg param to x
  param9$eta<-x
  # run xgboost cv
  seeds = set.seed(seed)
  tune_results_9[[y]] <- xgb.cv(params = param9,
                                data = RFE_training_data_xgb,
                                nfold= 10,
                                nrounds=10000, 
                                verbose=T,
                                showsd=T,
                                eval_metric="logloss", 
                                stratified=T,
                                early_stopping_rounds = 20)
  # add name
  names(tune_results_9)[y]<-x
  # change y
  eta_counter=eta_counter+1
  
}

#find lowest logloss
tune_results_9_lowest_logloss<-as.data.frame(unlist(lapply(tune_results_9, function(x) min(x$evaluation_log[,4]))))
tune_results_9_lowest_logloss

#nrounds
tune_results_9_lowest_logloss_nr<-as.data.frame(unlist(lapply(tune_results_9, function(x) x$best_iteration)))
tune_results_9_lowest_logloss_nr

# merge with nrounds
tune_results_9_lowest_logloss<-merge(tune_results_9_lowest_logloss, tune_results_9_lowest_logloss_nr, by="row.names")

#name
names(tune_results_9_lowest_logloss)<-c("eta", "test_mlogloss_mean", "nrounds")
#numeric
tune_results_9_lowest_logloss$eta<-as.numeric(tune_results_9_lowest_logloss$eta)
#order
tune_results_9_lowest_logloss<-tune_results_9_lowest_logloss[order(tune_results_9_lowest_logloss$test_mlogloss_mean, tune_results_9_lowest_logloss$eta),]
tune_results_9_lowest_logloss

##### MODEL 9 #####

# change param9 again
param9$eta<-as.numeric(tune_results_9_lowest_logloss$eta[1])
param9

set.seed(seed)
# create initial model - basic model 

xgb_model_9 <- xgboost(
  params = param9,
  data = RFE_training_data_xgb,
  nrounds = tune_results_9_lowest_logloss$nrounds[1] # set to what was found session
)

xgb_model_9

# save
setwd(work_dir)
xgb.save(xgb_model_9, "xgb_model_9")

##### TEST MODEL 9 #####

#test on training set
Model_9_train_prediction<-as.numeric(predict(xgb_model_9, RFE_training_data_xgb) > 0.5)
Model_9_train_prediction_results<-caret::confusionMatrix(as.factor(Model_9_train_prediction), as.factor(training_data_pheno_xgb))
Model_9_train_prediction_results

#test on validation set
Model_9_test_prediction<-as.numeric(predict(xgb_model_9, RFE_testing_data_xgb) > 0.5)
Model_9_test_prediction_results<-caret::confusionMatrix(as.factor(Model_9_test_prediction), as.factor(testing_data_pheno_xgb))
Model_9_test_prediction_results

#missclassified
Model_9_test_prediction_missclassifcation<-testing_data[1:2]
Model_9_test_prediction_missclassifcation$Prediction<-Model_9_test_prediction
# missclassifcation - disease 
table(Model_9_test_prediction_missclassifcation$Prediction, Model_9_test_prediction_missclassifcation$Diagnosis)
# missclassifcation - dataset
table(Model_9_test_prediction_missclassifcation$Prediction, Model_9_test_prediction_missclassifcation$Dataset)

# ##### CHECK RESULTS #####
# Model_1_train_prediction_results
# Model_1_train_prediction_results_rfe
# Model_2_train_prediction_results
# Model_3_train_prediction_results
# Model_4_train_prediction_results
# Model_5_train_prediction_results
# Model_6_train_prediction_results
# Model_7_train_prediction_results
# Model_8_train_prediction_results
# Model_9_train_prediction_results
# 
# 
# Model_1_test_prediction_results
# Model_1_test_prediction_results_rfe
# Model_2_test_prediction_results
# Model_3_test_prediction_results
# Model_4_test_prediction_results
# Model_5_test_prediction_results
# Model_6_test_prediction_results
# Model_7_test_prediction_results
# Model_8_test_prediction_results
# Model_9_test_prediction_results
# 
# 
# table(Model_9_test_prediction_missclassifcation$Prediction, Model_9_test_prediction_missclassifcation$Diagnosis)
# # missclassifcation - dataset
# table(Model_9_test_prediction_missclassifcation$Prediction, Model_9_test_prediction_missclassifcation$Dataset)
# 
# # for results table
# cbind(
#   
#   as.data.frame(Model_1_test_prediction_results$byClass),
#   as.data.frame(Model_1_test_prediction_results_rfe$byClass),
#   as.data.frame(Model_2_test_prediction_results$byClass),
#   as.data.frame(Model_3_test_prediction_results$byClass),
#   as.data.frame(Model_4_test_prediction_results$byClass),
#   as.data.frame(Model_5_test_prediction_results$byClass),
#   as.data.frame(Model_6_test_prediction_results$byClass),
#   as.data.frame(Model_7_test_prediction_results$byClass),
#   as.data.frame(Model_8_test_prediction_results$byClass),
#   as.data.frame(Model_9_test_prediction_results$byClass)
#   
# )
# 
# 
##### SELECT BEST MODEL - default last model #####

best_model<-xgb_model_9
best_model

#test on training set
best_model_train_prediction<-as.numeric(predict(best_model, RFE_training_data_xgb) > 0.5)
best_model_train_prediction_results<-caret::confusionMatrix(as.factor(best_model_train_prediction), as.factor(training_data_pheno_xgb))
best_model_train_prediction_results

#test on validation set
best_model_test_prediction<-as.numeric(predict(best_model, RFE_testing_data_xgb) > 0.5)
best_model_test_prediction_results<-caret::confusionMatrix(as.factor(best_model_test_prediction), as.factor(testing_data_pheno_xgb))
best_model_test_prediction_results

#missclassified
best_model_test_prediction_missclassifcation<-testing_data[1:2]
best_model_test_prediction_missclassifcation$Prediction<-best_model_test_prediction
# missclassifcation - disease 
table(best_model_test_prediction_missclassifcation$Prediction, best_model_test_prediction_missclassifcation$Diagnosis)
# missclassifcation - dataset
table(best_model_test_prediction_missclassifcation$Prediction, best_model_test_prediction_missclassifcation$Dataset)

##### VARIBALE IMPORTANCE #####

# Gain is the improvement in accuracy brought by a feature to the branches it is on
# Cover measures the relative quantity of observations concerned by a feature.
# Frequency is a simpler way to measure the Gain. It just counts the number of times a feature is used in all generated trees. You should not use it (unless you know why you want to use it).


#best model xgb_model_1

importance <- xgb.importance(feature_names = names(RFE_training_data), model = best_model)
head(importance)
dim(importance)
barplot(importance$Gain)

# 

importanceRaw <- xgb.importance(feature_names = names(RFE_training_data), model = best_model, data = as.matrix(RFE_training_data), label = as.numeric(training_data_pheno_xgb))
importanceClean <- importanceRaw[,`:=`(Cover=NULL, Frequency=NULL)]
head(importanceClean)


##### PREDICTION PROBABILITY - TESTING DATA #####
# raw prediction for each sample
FINAL_test_prediction<-predict(best_model, RFE_testing_data_xgb)

head(FINAL_test_prediction)

#merge diagnosis + dataset information
FINAL_test_prediction_missclassifcation<-testing_data[1:2]
FINAL_test_prediction_missclassifcation$Prediction<-FINAL_test_prediction
head(FINAL_test_prediction_missclassifcation)

# add column if AD or not
FINAL_test_prediction_missclassifcation$Clinical<- ifelse(FINAL_test_prediction_missclassifcation$Diagnosis=="AD_case", "AD", "Non-AD")

head(FINAL_test_prediction_missclassifcation)
table(FINAL_test_prediction_missclassifcation$Clinical)

head(FINAL_test_prediction_missclassifcation)
table(FINAL_test_prediction_missclassifcation$Dataset)

# add column if illumina or AFFY
FINAL_test_prediction_missclassifcation$Platform<- ifelse(FINAL_test_prediction_missclassifcation$Dataset %in% c("AD_GSE63060", "BP_GSE23848", "SCZ_GSE38481", "SCZ_GSE38484"), "Illumina", "Affymetrix")

head(FINAL_test_prediction_missclassifcation)
table(FINAL_test_prediction_missclassifcation$Platform)

# susbet to non-AD
FINAL_test_prediction_missclassifcation_non_AD<-subset(FINAL_test_prediction_missclassifcation, FINAL_test_prediction_missclassifcation$Clinical=="Non-AD")
dim(FINAL_test_prediction_missclassifcation_non_AD)

# change BP to BD
FINAL_test_prediction_missclassifcation[FINAL_test_prediction_missclassifcation$Diagnosis=="BP_control",1]<-"BD_control"
FINAL_test_prediction_missclassifcation[FINAL_test_prediction_missclassifcation$Diagnosis=="BP_case",1]<-"BD_case"

# seperate by diseased and healthy
FINAL_test_prediction_missclassifcation_copy<-FINAL_test_prediction_missclassifcation
table(FINAL_test_prediction_missclassifcation_copy$Diagnosis)

FINAL_test_prediction_missclassifcation_copy[grep("control", FINAL_test_prediction_missclassifcation_copy$Diagnosis),1]<-"Controls"

FINAL_test_prediction_missclassifcation_copy[FINAL_test_prediction_missclassifcation_copy$Diagnosis=="BD_case",1]<-"BD"
FINAL_test_prediction_missclassifcation_copy[FINAL_test_prediction_missclassifcation_copy$Diagnosis=="AD_case",1]<-"AD"
FINAL_test_prediction_missclassifcation_copy[FINAL_test_prediction_missclassifcation_copy$Diagnosis=="CD_case",1]<-"CD"
FINAL_test_prediction_missclassifcation_copy[FINAL_test_prediction_missclassifcation_copy$Diagnosis=="COPD_case",1]<-"COPD"
FINAL_test_prediction_missclassifcation_copy[FINAL_test_prediction_missclassifcation_copy$Diagnosis=="MS_case",1]<-"MS"
FINAL_test_prediction_missclassifcation_copy[FINAL_test_prediction_missclassifcation_copy$Diagnosis=="PD_case",1]<-"PD"
FINAL_test_prediction_missclassifcation_copy[FINAL_test_prediction_missclassifcation_copy$Diagnosis=="RA_case",1]<-"RA"
FINAL_test_prediction_missclassifcation_copy[FINAL_test_prediction_missclassifcation_copy$Diagnosis=="SCZ_case",1]<-"SCZ"

table(FINAL_test_prediction_missclassifcation_copy$Diagnosis)

##### SAVE #####

setwd(work_dir)
save.image("AD_specific_classifier_logloss.Rdata")

# xgb_model_1<-xgb.load("xgb_model_1")
# xgb_model_1_rfe<-xgb.load("xgb_model_1_rfe")
# xgb_model_2<-xgb.load("xgb_model_2")
# xgb_model_3<-xgb.load("xgb_model_3")
# xgb_model_4<-xgb.load("xgb_model_4")
# xgb_model_5<-xgb.load("xgb_model_5")
# xgb_model_6<-xgb.load("xgb_model_6")
# xgb_model_7<-xgb.load("xgb_model_7")
# xgb_model_8<-xgb.load("xgb_model_8")
# xgb_model_9<-xgb.load("xgb_model_9")

# started at 11:44

Sys.time()


##### OUPUT FILE #####

setwd(work_dir)
#variable importance
write.csv(importance, file="importance.csv")
#raw prediction
write.csv(FINAL_test_prediction_missclassifcation_copy, file="raw_validation_prediction.csv")
#results
validation_results<-as.data.frame(Model_9_test_prediction_results$byClass)
# add number of features
validation_results<-rbind(validation_results, list(length(unique(importance$Feature))))
rownames(validation_results)[12]<-"number of features"
# add seed
colnames(validation_results)<-seed

#add tuning paramters
final_tuning_param<-as.data.frame(t(as.data.frame(param9)))
colnames(final_tuning_param)<-seed
final_tuning_param<-final_tuning_param

validation_results<-rbind(validation_results, final_tuning_param)

#add nrounds
validation_results<-rbind(validation_results, tune_results_9_lowest_logloss$nrounds[1])
rownames(validation_results)[24]<-"nrounds"

write.csv(validation_results, file="validation_results.csv")

Sys.time()

# 1. lowest test-logloss
# 2. optimal hyperparamters
# 3. Variable importance
# 4. probability for testing set
# 5. performance metrics (sensitivity + specificity)
# 6. gender 
