setwd("C:/Users/test/OneDrive - webmail.hzau.edu.cn/DB_research/Mollah Sir/Paper Draft/Results File ML/")
#################################################################################

d1= read.csv("GSE23878_hub_protein.csv", header = TRUE)
d2= read.csv("GSE110224_hub_protein.csv", header = TRUE)
d3= read.csv("GSE9348_hub_protein.csv", header = TRUE)
d4= read.csv("GSE35279_hub_protein.csv", header = TRUE)
ddc = rbind(d1,d2,d3,d4)

library(BBmisc)
dataN = normalize(ddc[, -1])
X = ddc[,1]
dd = cbind(X, dataN)
set.seed(111)
sam1 = sort(sample(nrow(dd), nrow(dd)*0.8))
data_train <- dd[sam1,]
data_test <- dd[-sam1,]

tiff(file="heatmap_all.tiff")
heatmap(as.matrix(dataN))
dev.off()


dataI = read.csv("GSE23878_hub_protein.csv", header = TRUE)
dataIN = normalize(dataI[, -1])
X = dataI[,1]
data_ind = cbind(X,dataIN) 

#########TCGA#############
TCGA=read.csv("TCGA_colon.csv")

TCGAN = log(TCGA[, -1]+1)
X = TCGA[ ,1]
TCGA = cbind(X,TCGAN) 

###########################
#            RF           #
###########################

library(randomForest)
library(caret)
library(ROSE)

data_train$X = as.factor(data_train$X)
class(data_train$X)
set.seed(100)
rfFit <-  randomForest(data_train$X ~ CXCL8+MMP7+CA4+ADH1C+GUCA2A+GUCA2B+CEMIP+ZG16+CLCA4+MS4A12+CLDN1, 
                       data = data_train, ntree = 500, norm.votes=FALSE)


########## Train data set ###########
Pred_RF_train <- predict(rfFit, data_train[2:12], type = "class")

rf_tr=table(as.factor(data_train$X), Pred_RF_train)
TP = rf_tr[1,1]; TN = rf_tr[2,2]; FN = rf_tr[1,2]; FP = rf_tr[2,1];
TPR = (TP/(TP+FN)); TPR; #SENSITIVITY
TNR = (TN/(TN+FP)); TNR; #specificity
FDR = (FP/(FP+TP)); FDR; #false discovery rate
ACC = (TP+TN)/(TP+TN+FP+FN); ACC;
MCC = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)); MCC

rfOuttrain = cbind(TPR, TNR, FDR, MCC, ACC)

########## test data set ###########
Pred_RF_test <- predict(rfFit, data_test[2:12], type = "class")

rf_te=table(as.factor(data_test$X), Pred_RF_test)
TP = rf_te[1,1]; TN = rf_te[2,2]; FN = rf_te[1,2]; FP = rf_te[2,1];
TPR = (TP/(TP+FN)); TPR; #SENSITIVITY
TNR = (TN/(TN+FP)); TNR; #specificity
FDR = (FP/(FP+TP)); FDR; #false discovery rate
ACC = (TP+TN)/(TP+TN+FP+FN); ACC;
MCC = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)); MCC

rfOuttest = cbind(TPR, TNR, FDR, MCC, ACC)

########## ind data set ###########
Pred_RF_ind <- predict(rfFit, data_ind[2:12], type = "class")

rf_in=table(as.factor(data_ind$X), Pred_RF_ind)
TP = rf_in[1,1]; TN = rf_in[2,2]; FN = rf_in[1,2]; FP = rf_in[2,1];
TPR = (TP/(TP+FN)); TPR; #SENSITIVITY
TNR = (TN/(TN+FP)); TNR; #specificity
FDR = (FP/(FP+TP)); FDR; #false discovery rate
ACC = (TP+TN)/(TP+TN+FP+FN); ACC;
MCC = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)); MCC

rfOutind = cbind(TPR, TNR, FDR, MCC, ACC)

########## TCGA data set ###########

Pred_RF_TCGA <- predict(rfFit, TCGA[2:12], type = "class")

rf_TCGA=table(as.factor(TCGA$X), Pred_RF_TCGA)
TP = rf_TCGA[1,1]; TN = rf_TCGA[2,2]; FN = rf_TCGA[1,2]; FP = rf_TCGA[2,1];
TPR = (TP/(TP+FN)); TPR; #SENSITIVITY
TNR = (TN/(TN+FP)); TNR; #specificity
FDR = (FP/(FP+TP)); FDR; #false discovery rate
ACC = (TP+TN)/(TP+TN+FP+FN); ACC;
MCC = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)); MCC

rfOutiTCGA = cbind(TPR, TNR, FDR, MCC, ACC)

######### RF ROC Curve ###########

library(pROC)
tiff("ROC_RF4.tiff", width = 5, height = 4.5, res = 300, units = "in", compression = c("lzw"))

roc(as.numeric(as.factor(data_test$X)), as.numeric(as.factor(Pred_RF_test)),plot=TRUE,print.auc=TRUE,
    col="blue",lwd = 3,print.auc.y=0.55,legacy.axes=TRUE,main="ROC curve for RF Classifier")
roc(as.numeric(as.factor(data_train$X)), as.numeric(as.factor(Pred_RF_train)),plot=TRUE,print.auc=TRUE,
    col="green",lwd = 3,print.auc.y=0.6,legacy.axes=TRUE,add = TRUE)
roc(as.numeric(as.factor(data_ind$X)), as.numeric(as.factor(Pred_RF_ind)),plot=TRUE,print.auc=TRUE,
    col="red",lwd = 3,print.auc.y=0.5,legacy.axes=TRUE,add = TRUE)
roc(as.numeric(as.factor(TCGA$X)), as.numeric(as.factor(Pred_RF_TCGA)),plot=TRUE,print.auc=TRUE,
    col="black",lwd = 3,print.auc.y=0.45,legacy.axes=TRUE,add = TRUE)
legend("bottomright",legend=c("Train data","Test data", "Independent data","TCGA data"),col=c("green","blue","red","black"),lwd=3)
dev.off()


