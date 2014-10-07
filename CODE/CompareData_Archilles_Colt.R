rm(list = ls())
# Data preparation

library(affy)
library(CePa)
############## Project Archilles with CCLE EXP and Cell Viability
A<-read.gct("~/shRNA_database//Project_Achilles/Data//Achilles_102lines_gene_solutions.gct")
id_exprLayer <- "syn1757082" 
layer_expr <- loadEntity(id_exprLayer)
expr.ccle <- exprs(layer_expr$objects$eSet_expr)
name1.archilles<-sapply(strsplit(rownames(A),"_"),function(x){x[[1]]})
rownames(A)<-name1.archilles
commonColNames.archilles <- intersect(colnames(A), colnames(expr.ccle))

dataSet.archilles <- list(feature = expr.ccle[,commonColNames.archilles],response = A[,commonColNames.archilles])
colnames(dataSet.archilles$feature)<-sapply(strsplit(colnames(dataSet.archilles$feature),"_"),function(x){x[[1]]})
colnames(dataSet.archilles$response)<-sapply(strsplit(colnames(dataSet.archilles$response),"_"),function(x){x[[1]]})
############## Colt with Sanger EXP and Cell Viability
E<-read.delim("~/shRNA_database/COLT/Data/GARP-score.txt")
id_exprLayer <- "1742878"
layer_expr <- loadEntity(id_exprLayer)
expr.sanger <- exprs(layer_expr$objects$eSet_expr)

name1.colt<-E$Gene.name
a<-which(name1.colt == "")
newE<-E[-a,-c(1:4)]
name1.colt.new<-name1.colt[-a]
b<-duplicated(name1.colt.new)
d<-which(b!=1)
newEE<-newE[d,]
rownames(newEE)<-name1.colt.new[d]

name.duplicate<-c()
K<-c()
for(k in 1:length(which(b==1))){  
  name.duplicate<-c(name.duplicate, name1.colt.new[which(b==1)[k]])
  K<-rbind(K,apply(newE[which(name1.colt.new == name1.colt.new[which(b==1)[k]]),],2,mean))
}
rownames(K)<-name.duplicate
finalE<-rbind(newEE,K)

name2.colt<-gsub("\\.","",toupper(make.names(colnames(finalE))))
colnames(finalE)<-name2.colt

commonColNames.colt <- intersect(name2.colt, colnames(expr.sanger))
dataSet.colt <- list(feature = expr.sanger[,commonColNames.colt],response = finalE[,commonColNames.colt])
###################

######### input variable data between archilles(ccle) and colt(sanger) should be compared  (same features)

common.feature<-intersect(rownames(dataSet.archilles$feature),rownames(dataSet.colt$feature))
common.responses<-intersect(rownames(dataSet.archilles$response),rownames(dataSet.colt$response))
common.samples<-intersect(colnames(dataSet.archilles$response),colnames(dataSet.colt$response))

tmp.archilles<-dataSet.archilles$response[common.responses,common.samples]
tmp.colt     <-as.matrix(dataSet.colt$response[common.responses,common.samples])

tmp1.archilles<-dataSet.archilles$feature[common.feature,common.samples]
tmp1.colt     <-as.matrix(dataSet.colt$feature[common.feature,common.samples])


cc<-c()
for(k in 1:nrow(tmp.archilles)){
  cc<-c(cc,cor(tmp.archilles[k,],tmp.colt[k,],method = "spearman"))
}

cc1<-c()
for(k in 1:nrow(tmp1.archilles)){
  cc1<-c(cc1,cor(tmp1.archilles[k,],tmp1.colt[k,],method = "spearman"))
}

hist(cc,30)
hist(cc1,30)


