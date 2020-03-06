##############################################################################
library(biomaRt)
library(rentrez)
library(data.table)
library(pbapply)

interacs=fread("slctdPrdctrs.tsv")#output of elasticNetToSif.R
pam50=read.table("pam50.tsv",header=T)#table with pam50 ids from ensembl,hgnc & class

########################### CpGs ###########################
methy=fread("HumanMethylation450_15017482_v1-2.csv",sep=',',header=T,fill=T)#use Illumina file
methy=methy[,c(1,22,24)]
methy=methy[methy[,2]!="",]
methy=do.call(rbind,apply(methy,1,function(x) 
	cbind(x[1],
		unlist(strsplit(x[2],";")),
		unlist(strsplit(x[3],";")))))
writr.table(methy,
	"MapMethy.tsv",
	sep='\t',
	quote=F,
	row.names=F)

withCpG=do.call(rbind,lapply(unique(interacs$pam50Symbol),function(x)
	interacs[interacs$pam50Symbol==x,c(1,3:5)][which(interacs$predictor[interacs$pam50Symbol==x]%in%methy[methy[,2]==x,1]),]))

########################### miRNAs ###########################

library(multiMiR)#https://www.bioconductor.org/packages/devel/bioc/vignettes/multiMiR/inst/doc/multiMiR.html
#Searching mirecords, mirtarbase, tarbase,diana_microt,elmmo, microcosm, miranda, mirdbpictar, pita, targetscan, pharmaco_mir ...
pam50=read.table("../ini/pam50.tsv",header=T)

#get all validated miRs interactions for PAM50 genes
miRtargets=get_multimir(target=pam50$ensembl_gene_id,
	summary=F,
	table="all",
	legacy.out=F)
miRtargetsV=select(miRtargets,
	keys="validated",
	columns=columns(miRtargets),
	keytype="type")
miRtargetsP=select(miRtargets,
	keys="predicted",
	columns=columns(miRtargets),
	keytype="type")
#miRtargets=unique(rbind(miRtargetsV1[,3:4],miRtargetsP[,3:4]))
mirIDs=fread(miR.ids.map.tsv)#ftp://mirbase.org/pub/mirbase/CURRENT/aliases.txt.gz
withMIR=do.call(rbind,lapply(unique(interacs$pam50Symbol),function(x) 
	interacs[interacs$pam50Symbol==x,][interacs$predictor[interacs$pam50Symbol==x]%in%mirIDs$precursor[mirIDs$mature%in%miRtargets$mature_mirna_id[miRtargets$target_symbol==x]],]))

########################### TFs ###########################

library(tftargets)#https://github.com/slowkow/tftargets

#transform TF target lists to tables 
load("TFtargets.RData")
tfs=list()
tfs$TRED=do.call(rbind,lapply(1:length(TRED),function(x) 
	cbind(names(TRED)[x],TRED[[x]])))
temp=lapply(1:length(Neph2012),function(x) lapply(1:length(Neph2012[[x]]),function(y) 
	cbind(names(Neph2012[[x]])[y],Neph2012[[x]][[y]])))
temp=lapply(temp,function(x) do.call(rbind,x[sapply(x,ncol)==2]))
tfs$Neph2012=do.call(rbind,sapply(1:length(temp),function(x) 
	cbind(names(Neph2012)[x],temp[[x]])))
#for each TF source

#get hgnc_id for data annotated with entrez id
mart=useEnsembl("ensembl",
	dataset="hsapiens_gene_ensembl",
	host="http://jan2019.archive.ensembl.org")
#mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl",version=95)
myannot=getBM(
 attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
 mart=mart)
myannot=myannot[!is.na(myannot$entrezgene),]
myannot=myannot[!duplicated(myannot$entrezgene),]
tfs$TRED=tfs$TRED[tfs$TRED[,2]%in%myannot$entrezgene,]
tfs$TRED=tfs$TRED[order(tfs$TRED[,2]),]
temp=table(tfs$TRED[,2])
ids=unlist(sapply(1:length(temp),function(x)
 rep(as.character(myannot$hgnc_symbol)[which(as.character(myannot$entrezgene)==names(temp)[x])],
 	temp[x])))
tfs$TRED[,2]=unlist(ids)
tfs$TRED=cbind("",tfs$TRED)
tfs=do.call(rbind,lapply(1:6,function(x) cbind(names(tfs)[x],tfs[[x]])))
tfs$V3=gsub("E2F-","E2F",tfs$V3)
tfs$V3=gsub("c-","C",tfs$V3)

#keep only the TF of intereset
TFtargets=tfs[tfs$V4%in%pam50$hgnc_symbol,]
interacs=interacs[interacs$predictorSymbol%in%TFtargets[,3],]
#true TF interactions found
withTF=apply(interacs,1,function(x) 
	TFtargets[TFtargets[,3]==x[6]&TFtargets[,4]==x[5],])
#format nicely
interacs=interacs[sapply(withTF,length)>0,]
withTF=withTF[sapply(withTF,length)>0]
withTF=sapply(withTF,function(x) matrix(x,ncol=4))
withTF=lapply(withTF,function(x) 
	cbind(x[1,1],paste(x[,2],collapse=", "),x[1,3],x[1,4]))
withTF=cbind(interacs[,1],do.call(rbind,withTF))
i=paste(withTF$V3,withTF$V4)
withTF=lapply(unique(i),function(x) withTF[i==x,])
withTF=t(sapply(withTF,function(x) 
	cbind(paste(x$subtype,collapse=", "),x$V1[1],x$V2[1],x$V3[1],x$V4[1])))


######################################################################
#######how many times terms are comentioned in literature
#####################################################################
#both on the same paper
set_entrez_key(k)#ncbi account to submit 10 queries per second

comention=pbsapply(seq(1,12112,10),function(i) 
	apply(interacs[i:(i+9),],1,function(x) {
		reque=entrez_search(db = "pubmed",
		 term = paste(x[2]," AND ",x[3],
		 	collapse=" "));
		Sys.sleep(0.1);
		return(reque)}))


cuentas=as.numeric(apply(comention,c(1,2),function(x) unlist(x[[1]][2])))
comen=cbind(interacs[,2:3],cuentas[1:12112])#further than interacs nrow, search was NA AND NA
comen=comen[comen[,2]!="0",]
colnames(comen)[3]="comention"
write.table(comen,"comention.tsv",sep='\t',quote=F,row.names=F)

#each mentioned
query=unique(unlist(comen[,3:4]))
mention=pbsapply(seq(1,10746,7),function(i)
	sapply(query[i:(i+6)],function(x){
	reque=entrez_search(db = "pubmed", term = x);
	Sys.sleep(0.1);
	return(reque)}))
#gives a matrix of output lines(5)*searched terms(7) rows per length(seq(1,10746,7)) columns
cuentas=unlist(apply(mention,2,function(x) x[seq(2,35,5)]))[1:length(query)]#2 is the index of id counts in the columns
#further than query len, search was NA
mention=cbind(query,cuentas)

continComent=function(interac){
	members=as.character(interac)
	a=interac$comention
	b=mention$cuentas[mention$query==members[1]]-a
	c=mention$cuentas[mention$query==members[2]]-a
	d=30151833-b-c-a
	mat=matrix(c(a,b,c,d),ncol=2,nrow=2)
	pval=fisher.test(mat,alternative="greater")$p.val
	return(c(mat,pval))}
fshr=sapply(1:nrow(comen),function(x) continComent(comen))
fshr=cbind(fshr,p.adjust(fshr[,5],"fdr")
colnames(fshr)=c("comention","pam50mention","predimention","neither","p.val","fdr")
comention=cbind(comention,fshr)
write.table(comention,"comention.tsv",sep='\t',quote=F,row.names=F)

