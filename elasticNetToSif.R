################ libraries #############################
library(ggplot2)
library(biomaRt)

bestModels=read.table("bestModels.tsv",sep='\t')

############# build sifs for non null models############  
bestModels=bestModels[bestModels$dev.ratio>0,]
#read files from best models
files=apply(bestModels,1,function(x) paste(x[2],x[1],"coefs",sep='.'))
coefs=sapply(files,readLines)
#write them nicely
names(coefs)=gsub(".coefs","",files)
coefs=sapply(coefs,function(x) t(do.call(cbind,strsplit(x,"\t"))))
coefs=do.call(rbind,lapply(1:length(files),function(x)
 cbind(names(coefs)[x],coefs[[x]])))
coefs=coefs[coefs[,2]!="x",]
coefs=coefs[coefs[,2]!="(Intercept)",]
coefs=coefs[order(coefs[,1]),]
#get hgnc symbols from ensembl_ids for pam50 genes
temp1=table(sapply(strsplit(coefs[,1],".",fixed=T),function(x) x[1]))
ids=unlist(sapply(1:length(temp1),function(x) 
	rep(as.character(pam50$hgnc_symbol)[pam50$ensembl_gene_id==names(temp1)[x]],
		temp1[x])))
temp=cbind(coefs,ids)
#get hgnc symbols from ensembl_ids for predictors
mart=useEnsembl("ensembl",
	dataset="hsapiens_gene_ensembl",
	host="http://apr2019.archive.ensembl.org")
myannot=getBM(attributes = c("ensembl_gene_id", "hgnc_id","hgnc_symbol"),
	filters = "ensembl_gene_id", 
	values=unique(as.character(temp[,3])),
	mart=mart)
temp=temp[order(temp[,3]),]
temp1=table(temp[grep("ENSG",temp[,3]),3])
ids=unlist(sapply(1:length(temp1),function(x) 
	rep(myannot$hgnc_symbol[myannot$ensembl_gene_id==names(temp1)[x]],temp1[x])))
temp=cbind(coefs,ids)
#fit all into a long table with model per subtype per gene
coefs1=lapply(unique(bestModels$subtype),function(x) temp[grep(x,temp[,2]),c(1,3,4)])
coefs1=do.call(rbind,lapply(1:5,function(x) cbind(names(coefs1)[x],coefs1[[x]])))
colnames(coefs1)=c("subtype","pam50","predictor","coef","pam50Symbol","predictorSymbol")
write.table(coefs1,"slctdPrdctrs.tsv",sep='\t',quote=F,row.names=F)
