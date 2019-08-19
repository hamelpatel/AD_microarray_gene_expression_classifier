##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                 INTERPRETING REPEATED RESULTS                                                         #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################

# NOTES
# AD vs control - weights and repeated X1000
# AD vs mixedcontrol - bootstrap mixed control + weights - repeated x1000
# data on Rosalind
#
#

######## BOXPLOT OF EVALUATION METRIC #########

##### LOAD DATA #####

#linux
setwd("/media/hamel/Workspace/Dropbox/Projects/AD-classification/6.XGBoost/Upsampling_and_weights_repeated_FINAL")

#windows
#setwd("D:\\Dropbox\\Projects\\AD-classification\\6.XGBoost\\Upsampling_and_weights_repeated_FINAL")

AD_vs_Control<-read.csv("AD_vs_control_combined_results.csv", check.names = F, row.names = 1, stringsAsFactors = F)
AD_vs_MixedControl<-read.csv("AD_vs_mixedcontrol_combined_results.csv", check.names = F, row.names = 1, stringsAsFactors = F)

head(AD_vs_Control[1:5])
head(AD_vs_MixedControl[1:5])

##### REMOVE ROWS 13 and 14 #####

AD_vs_Control<-AD_vs_Control[-c(13,14),]
AD_vs_MixedControl<-AD_vs_MixedControl[-c(13,14),]

AD_vs_Control[1]
AD_vs_MixedControl[1]

# make numeric

head(str(AD_vs_Control))
AD_vs_Control_new<-as.data.frame(apply(AD_vs_Control, 2, function(x) as.numeric(x)))
rownames(AD_vs_Control_new)<-rownames(AD_vs_Control)
head(str(AD_vs_Control_new))
AD_vs_Control_new[1]

head(str(AD_vs_MixedControl))
AD_vs_MixedControl_new<-as.data.frame(apply(AD_vs_MixedControl, 2, function(x) as.numeric(x)))
rownames(AD_vs_MixedControl_new)<-rownames(AD_vs_MixedControl)
AD_vs_MixedControl_new[1]


##### T-Test #####

t.test(AD_vs_Control_new[1,], AD_vs_MixedControl_new[1,])

##### CALCULATE AVERAGE + 95% CI #####

AD_vs_Control_ave<-as.data.frame(apply(AD_vs_Control_new, 1, function(x) mean(x)))
AD_vs_Control_ave[2]<-apply(AD_vs_Control_new, 1, function(x) quantile(x, 0.025))
AD_vs_Control_ave[3]<-apply(AD_vs_Control_new, 1, function(x) quantile(x, 0.975))
colnames(AD_vs_Control_ave)<-c("Average", "0.25_Quantile", "0.975_Quantile")

AD_vs_Control_ave

AD_vs_MixedControl_ave<-as.data.frame(apply(AD_vs_MixedControl_new, 1, function(x) mean(x)))
AD_vs_MixedControl_ave[2]<-apply(AD_vs_MixedControl_new, 1, function(x) quantile(x, 0.025))
AD_vs_MixedControl_ave[3]<-apply(AD_vs_MixedControl_new, 1, function(x) quantile(x, 0.975))
colnames(AD_vs_MixedControl_ave)<-c("Average", "0.25_Quantile", "0.975_Quantile")

AD_vs_MixedControl_ave

##### WRTIE FILE #####

write.table(AD_vs_Control_ave, file="AD_vs_Control_ave.txt", quote=F, sep=",")
write.table(AD_vs_MixedControl_ave, file="AD_vs_MixedControl_ave.txt", quote=F, sep=",")


##### PLOTS #####

library(reshape2)
library(ggplot2)

#format data
#columns wanted - sens, spec, balanced acc, PPV, NPV

AD_vs_Control_plot_data<-as.data.frame(t(AD_vs_Control_new[c(1, 2, 11, 3, 4),]))
AD_vs_Control_plot_data<-melt(AD_vs_Control_plot_data)
AD_vs_Control_plot_data$Model<-"AD vs Healthy Controls"
head(AD_vs_Control_plot_data)

AD_vs_MixedControl_plot_data<-as.data.frame(t(AD_vs_MixedControl_new[c(1, 2, 11, 3, 4),]))
AD_vs_MixedControl_plot_data<-melt(AD_vs_MixedControl_plot_data)
AD_vs_MixedControl_plot_data$Model<-"AD vs Mixed Controls"
head(AD_vs_MixedControl_plot_data)

#merge
plot_data<-rbind(AD_vs_Control_plot_data, AD_vs_MixedControl_plot_data)

# percentage
plot_data$value<-(plot_data$value*100)

#colnames
colnames(plot_data)<-c("Evaluation metric", "Percentage (%)", "Classifcation Model")

ggplot(plot_data, aes(x=`Evaluation metric`, y=`Percentage (%)`, fill=`Classifcation Model`)) + 
  geom_boxplot() +
  ggtitle("Comparison of classification models performances") +
  theme(plot.title = element_text(hjust = 0.5))


tiff("Comparison of classification models performances.tiff", width = 8, height = 8, unit="in",  res=300)
ggplot(plot_data, aes(x=`Evaluation metric`, y=`Percentage (%)`, fill=`Classifcation Model`)) + 
  geom_boxplot() +
  ggtitle("Comparison of classification models performances") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("Comparison of classification models performances.png", width = 8, height = 8, unit="in",  res=300)
ggplot(plot_data, aes(x=`Evaluation metric`, y=`Percentage (%)`, fill=`Classifcation Model`)) + 
  geom_boxplot() +
  ggtitle("Comparison of classification models performances") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


##### CUI #####
# CUI -  AD_vs_control

head(AD_vs_Control)[1:5]

AD_vs_Control_CUI_pos<-list()
AD_vs_Control_CUI_neg<-list()

for (x in 1:1000) {
  AD_vs_Control_CUI_pos<-c(AD_vs_Control_CUI_pos, as.numeric(AD_vs_Control[3,x])*as.numeric(AD_vs_Control[1,x]))
  AD_vs_Control_CUI_neg<-c(AD_vs_Control_CUI_neg, as.numeric(AD_vs_Control[4,x])*as.numeric(AD_vs_Control[2,x]))
}

AD_vs_Control_CUI_pos<-unlist(AD_vs_Control_CUI_pos)
AD_vs_Control_CUI_neg<-unlist(AD_vs_Control_CUI_neg)

summary(AD_vs_Control_CUI_pos)
summary(AD_vs_Control_CUI_neg)

quantile(AD_vs_Control_CUI_pos, 0.025)
quantile(AD_vs_Control_CUI_pos, 0.975)

quantile(AD_vs_Control_CUI_neg, 0.025)
quantile(AD_vs_Control_CUI_neg, 0.975)


# CUI -  AD_vs_MixedControl

head(AD_vs_MixedControl)[1:5]

AD_vs_MixedControl_CUI_pos<-list()
AD_vs_MixedControl_CUI_neg<-list()

for (x in 1:1000) {
  AD_vs_MixedControl_CUI_pos<-c(AD_vs_MixedControl_CUI_pos, as.numeric(AD_vs_MixedControl[3,x])*as.numeric(AD_vs_MixedControl[1,x]))
  AD_vs_MixedControl_CUI_neg<-c(AD_vs_MixedControl_CUI_neg, as.numeric(AD_vs_MixedControl[4,x])*as.numeric(AD_vs_MixedControl[2,x]))
}

AD_vs_MixedControl_CUI_pos<-unlist(AD_vs_MixedControl_CUI_pos)
AD_vs_MixedControl_CUI_neg<-unlist(AD_vs_MixedControl_CUI_neg)

summary(AD_vs_MixedControl_CUI_pos)
summary(AD_vs_MixedControl_CUI_neg)

quantile(AD_vs_MixedControl_CUI_pos, 0.025)
quantile(AD_vs_MixedControl_CUI_pos, 0.975)

quantile(AD_vs_MixedControl_CUI_neg, 0.025)
quantile(AD_vs_MixedControl_CUI_neg, 0.975)


########## BOXPLOT OF RAW PREDICTION #####

AD_vs_Control_pred<-read.csv("AD_vs_control_combined_prediction.csv", check.names = F, row.names = 1, stringsAsFactors = F)
AD_vs_MixedControl_pred<-read.csv("AD_vs_mixedcontrol_combined_prediction.csv", check.names = F, row.names = 1, stringsAsFactors = F)

#remove unwanted cols
AD_vs_Control_pred$Dataset<-NULL
AD_vs_Control_pred$Platform<-NULL

AD_vs_MixedControl_pred$Dataset<-NULL
AD_vs_MixedControl_pred$Platform<-NULL

# reorder
AD_vs_Control_pred<-AD_vs_Control_pred[c(1,3,2,4:1002)]
AD_vs_MixedControl_pred<-AD_vs_MixedControl_pred[c(1,3,2,4:1002)]

#check
head(AD_vs_Control_pred)[1:10]
head(AD_vs_MixedControl_pred)[1:10]

# reshape data

AD_vs_Control_pred_plot_data<-AD_vs_Control_pred
AD_vs_Control_pred_plot_data$SampleID<-rownames(AD_vs_Control_pred_plot_data)
AD_vs_Control_pred_plot_data<-melt(AD_vs_Control_pred_plot_data)

head(AD_vs_Control_pred_plot_data)
dim(AD_vs_Control_pred_plot_data)

AD_vs_MixedControl_pred_plot_data<-AD_vs_MixedControl_pred
AD_vs_MixedControl_pred_plot_data$SampleID<-rownames(AD_vs_MixedControl_pred_plot_data)
AD_vs_MixedControl_pred_plot_data<-melt(AD_vs_MixedControl_pred_plot_data)

head(AD_vs_MixedControl_pred_plot_data)
dim(AD_vs_MixedControl_pred_plot_data)


AD_vs_MixedControl_pred_plot_data<-AD_vs_MixedControl_pred
AD_vs_MixedControl_pred_plot_data$SampleID<-rownames(AD_vs_MixedControl_pred_plot_data)
AD_vs_MixedControl_pred_plot_data<-melt(AD_vs_MixedControl_pred_plot_data)

head(AD_vs_MixedControl_pred_plot_data)
dim(AD_vs_MixedControl_pred_plot_data)


#merge
AD_vs_Control_pred_plot_data$Model<-"AD vs Healthy Controls"
AD_vs_MixedControl_pred_plot_data$Model<-"AD vs Mixed Controls"

data_for_prediction_plot<-rbind(AD_vs_Control_pred_plot_data, AD_vs_MixedControl_pred_plot_data)
head(data_for_prediction_plot)

# plot

ggplot(AD_vs_MixedControl_pred_plot_data, aes(x=Diagnosis, y=value, color=Diagnosis)) + 
  geom_boxplot() +
  ggtitle("(b) AD-specific classification model prediction probability") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(AD_vs_Control_pred_plot_data, aes(x=Diagnosis, y=value, color=Diagnosis)) + 
  geom_boxplot() +
  ggtitle("(a) AD classification model prediction probability") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# plot to tiff

tiff("Classification Prediction Probability by Clinical Diagnosis1.tiff", width = 8, height = 8, unit="in",  res=300)
ggplot(AD_vs_Control_pred_plot_data, aes(x=Diagnosis, y=value, color=Diagnosis)) + 
  geom_boxplot() +
  ggtitle("(a) AD vs Healthy Control Classifcation Model Prediction Probability") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

tiff("Classification Prediction Probability by Clinical Diagnosis2.tiff", width = 8, height = 8, unit="in",  res=300)
ggplot(AD_vs_MixedControl_pred_plot_data, aes(x=Diagnosis, y=value, color=Diagnosis)) + 
  geom_boxplot() +
  ggtitle("(b) AD vs Mixed Control Classifcation Model Prediction Probability") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()#######

png("Classification Prediction Probability by Clinical Diagnosis1.png", width = 8, height = 8, unit="in",  res=300)
ggplot(AD_vs_Control_pred_plot_data, aes(x=Diagnosis, y=value, color=Diagnosis)) + 
  geom_boxplot() +
  ggtitle("(a) AD classification model prediction probability") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

png("Classification Prediction Probability by Clinical Diagnosis2.png", width = 8, height = 8, unit="in",  res=300)
ggplot(AD_vs_MixedControl_pred_plot_data, aes(x=Diagnosis, y=value, color=Diagnosis)) + 
  geom_boxplot() +
  ggtitle("(b) AD-specific classification model prediction probability") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()#######

###### AVERAGE AUC ######
library(ROCR)

#copy raw predcition - AD vs control

AD_vs_control_AUC_pred<-AD_vs_Control_pred

head(AD_vs_control_AUC_pred)[1:5]

AD_vs_control_AUC_pred$Diagnosis<-ifelse(AD_vs_control_AUC_pred$Clinical=="Non-AD", 1, 0)

head(AD_vs_control_AUC_pred)[1:5]
table(AD_vs_control_AUC_pred$Diagnosis)

AD_vs_control_AUC<-list()

for (x in 3:1002) {
  best_model_ROC_pred<-prediction(as.vector(t(AD_vs_control_AUC_pred[x])), as.vector(t(AD_vs_control_AUC_pred[1])))
  best_model_ROC_perf<-performance( best_model_ROC_pred, measure="auc")
  AD_vs_control_AUC<-c(AD_vs_control_AUC, best_model_ROC_perf@y.values)
}

AD_vs_control_AUC<-unlist(AD_vs_control_AUC)
summary(AD_vs_control_AUC)

quantile(AD_vs_control_AUC, 0.025)
quantile(AD_vs_control_AUC, 0.975)

#copy raw predcition - AD_vs_MixedControl_AUC

AD_vs_MixedControl_AUC_pred<-AD_vs_MixedControl_pred

head(AD_vs_MixedControl_AUC_pred)[1:5]

AD_vs_MixedControl_AUC_pred$Diagnosis<-ifelse(AD_vs_MixedControl_AUC_pred$Clinical=="Non-AD", 1, 0)

head(AD_vs_MixedControl_AUC_pred)[1:5]
table(AD_vs_MixedControl_AUC_pred$Diagnosis)

AD_vs_MixedControl_AUC<-list()

for (x in 3:1002) {
  best_model_ROC_pred<-prediction(as.vector(t(AD_vs_MixedControl_AUC_pred[x])), as.vector(t(AD_vs_MixedControl_AUC_pred[1])))
  best_model_ROC_perf<-performance( best_model_ROC_pred, measure="auc")
  AD_vs_MixedControl_AUC<-c(AD_vs_MixedControl_AUC, best_model_ROC_perf@y.values)
}

AD_vs_MixedControl_AUC<-unlist(AD_vs_MixedControl_AUC)
summary(AD_vs_MixedControl_AUC)

quantile(AD_vs_MixedControl_AUC, 0.025)
quantile(AD_vs_MixedControl_AUC, 0.975)

######### FEATURES ############################

AD_vs_Control_feat<-read.csv("AD_vs_control_combined_features.csv", check.names = F, row.names = 1, stringsAsFactors = F)
AD_vs_MixedControl_feat<-read.csv("AD_vs_mixedcontrol_combined_features.csv", check.names = F, row.names = 1, stringsAsFactors = F)

head(AD_vs_Control_feat)[1:5]
head(AD_vs_MixedControl_feat)[1:5]

# count
AD_vs_Control_feat_count<-as.data.frame(table(as.vector(t(AD_vs_Control_feat))))
AD_vs_MixedControl_feat_count<-as.data.frame(table(as.vector(t(AD_vs_MixedControl_feat))))

colnames(AD_vs_Control_feat_count)[1]<-"EntrezID"
colnames(AD_vs_MixedControl_feat_count)[1]<-"EntrezID"

#sort
AD_vs_Control_feat_count<-AD_vs_Control_feat_count[order(-AD_vs_Control_feat_count$Freq),]
AD_vs_MixedControl_feat_count<-AD_vs_MixedControl_feat_count[order(-AD_vs_MixedControl_feat_count$Freq),]

AD_vs_Control_feat_count
AD_vs_MixedControl_feat_count

#write
write.table(AD_vs_Control_feat_count, file="AD_vs_Control_feat_count.txt", quote=F, row.names = F)
write.table(AD_vs_MixedControl_feat_count, file="AD_vs_MixedControl_feat_count.txt", quote=F, row.names = F)


###### AVERAGE PRECITION PER SAMPLE OVER 1000 #####

# AD vs control

head(AD_vs_Control_pred)[1:5]

AD_vs_Control_pred_average<-AD_vs_Control_pred
AD_vs_Control_pred_average$Diagnosis<-NULL
AD_vs_Control_pred_average$Clinical<-NULL

AD_vs_Control_pred_average<-apply(AD_vs_Control_pred_average, 1, function(x) mean(x))
head(AD_vs_Control_pred_average)
length(AD_vs_Control_pred_average)


AD_vs_Control_pred_average_prediction<-as.numeric(AD_vs_Control_pred_average > 0.5)

AD_vs_Control_pred_average_prediction_results<-caret::confusionMatrix(as.factor(AD_vs_Control_pred_average_prediction), as.factor(as.vector(t(AD_vs_control_AUC_pred[1]))))

AD_vs_Control_pred_average_prediction_results


# AD vs Mixed Control

head(AD_vs_MixedControl_pred)[1:5]

AD_vs_MixedControl_pred_average<-AD_vs_MixedControl_pred
AD_vs_MixedControl_pred_average$Diagnosis<-NULL
AD_vs_MixedControl_pred_average$Clinical<-NULL

AD_vs_MixedControl_pred_average<-apply(AD_vs_MixedControl_pred_average, 1, function(x) mean(x))
head(AD_vs_MixedControl_pred_average)
length(AD_vs_MixedControl_pred_average)


AD_vs_MixedControl_pred_average_prediction<-as.numeric(AD_vs_MixedControl_pred_average > 0.5)

AD_vs_MixedControl_pred_average_prediction_results<-caret::confusionMatrix(as.factor(AD_vs_MixedControl_pred_average_prediction), as.factor(as.vector(t(AD_vs_MixedControl_AUC_pred[1]))))

AD_vs_MixedControl_pred_average_prediction_results
