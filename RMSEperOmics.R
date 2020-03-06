#!/usr/bin/env Rscript

################ LIBRARIES #############################
library(data.table)
library(caret)#RMSE function

################ ARRANGE INPUT #########################
#repeat per subtype!!!!
args = commandArgs(trailingOnly=TRUE)
subtype=args[1]
interacs=fread("slctdPrdctrs.tsv")#output of elasticNetToSif.R
interacs=interacs[interacs$subtype==subtype,]

subtype=fread(""paste("../",subtype,sep=""),sep='\t'"")
#fread loads a data frame with feature names as additional column factor
names=subtype$V1
subtype$V1=NULL
#index of training set
i=round(ncol(subtype)*0.8)
subtype=t(as.matrix(subtype))
#correct feature names
colnames(subtype)=names
#separate training & testing data
training=subtype[1:i,]
testing=subtype[(i+1):nrow(subtype),]

pam50=read.table("pam50.tsv",header=T)#table with pam50 ids from ensembl,hgnc & class

################ TO GET THE RMSE PER OMIC #########################
rmse_omic=function(gen,omic){
#keep only the predictors of the gen of interest
 model=interacs[interacs$pam50==gen,]
 subset=testing[,colnames(testing)%in%model$predictor]
#if only 1 predictor was selected 
 if(nrow(model)==1){
    return(RMSE(colSums(t(subset)*model$coef),testing[,colnames(testing)==gen]))}
 #if more than 1 predictors were selected
 subset=subset[,order(match(colnames(subset),model$predictor))]
 #index of predictors from the omic of interest
 i=which(substr(colnames(subset),1,1)==omic)
 if(omic=="all"){i=seq(1,ncol(subset))}
 #if no predictors of the omic of interest were selected
 if(length(i)==0){return(NA)}
 predis=colSums(t(subset[,i])*model$coef[i])
 obs=testing[,colnames(testing)==gen]
#actual RMSE calculation 
 RMSE(predis,obs)}

################ RMSE PER OMIC #########################
#repeat per subtype!!!!
omicsContri=t(pbsapply(pam50$ensembl_gene_id,function(y) 
              sapply(c("c","E","h","all"),function(x) rmse_omic(y,x))))
rownames(omicsContri)=pam50$ensembl_gene_id
write.table(omicsContri,
    file=paste(interacs$subtype[1],"omicsContri","tsv",sep='.'),
    sep='\t',
    quote=F,
    row.names=F)
