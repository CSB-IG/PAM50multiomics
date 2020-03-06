#!/usr/bin/env Rscript

################ libraries #############################
library(data.table)
library(glmnet)
library(doParallel)
library(methods)
library(caret)

################ ARRANGE INPUT #########################
#registerDoParallel(cores=3)
args = commandArgs(trailingOnly=TRUE)
subtype=args[1]
gen=args[2]
subtype=fread(paste("../",subtype,sep=""),sep='\t')
#fread loads a data frame with feature names as additional column factor
nombres=subtype$V1
subtype$V1=NULL
#index of training set
i=round(ncol(subtype)*0.8)
subtype=t(as.matrix(subtype))
#correct feature names
colnames(subtype)=nombres
#separate training & testing data
training=subtype[1:i,]
testing=subtype[(i+1):nrow(subtype),]
#choose k folds
k=5
if(nrow(training)<100){k=3}
registerDoParallel(cores=k)
print("input done")

################ MODEL FIT #############################
cvfit=function(gen,k){
#function repeat cross-validation
#print(dim(training[,colnames(training)!=gen]))
#print(length(training[,colnames(training)==gen]))
cv.glmnet(
	y =training[,colnames(training)==gen],
	x =training[,colnames(training)!=gen],
	lambda=10^ seq (3, -3 , length=50),
	nfolds=k,
	parallel=T,
	intercept=F,
	alpha=0.5,
	standardize=T,
#weights control shrinkage per omic
	weights=rep(1,nrow(training)))}

#base fit
print("start main fit")
#seed(123)
model=cvfit(gen,k)

#choose the lambda with the smallest MSE
print("start repeated cv")
cvm=sapply(1:99,function(i) cvfit(gen,k)$cvm)
cvm=cbind(cvm,model$cvm)
bestLambda=model$lambda[which.min(rowMeans(cvm))]

################ BEST LAMBDA PREDICTORS ################
coefs=coef(model,s=bestLambda)
#keep only predicting predictors 
coefs1=matrix(coefs[which(coefs!=0)],ncol=1)
rownames(coefs1)=rownames(coefs)[which(coefs!=0)]
#get the error for testing dataset
ajuste=RMSE(predict(model,testing[,colnames(testing)!=gen]),testing[,colnames(testing)==gen])
colnames(coefs1)=paste("lambda",bestLambda,"RMSE",ajuste,sep=':')
write.table(coefs1,
	    file=paste(gen,args[1],"coefs",sep='.'),
	    quote=F,
	    sep='\t')
