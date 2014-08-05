rm(list = ls())
# Data preparation
library(CePa)
A<-read.gct("~/shRNA_database//Project_Achilles/Data//Achilles_102lines_gene_solutions.gct")
B<-read.gct("~/shRNA_database//Project_Achilles/Data//20110303_achilles2_PMAD_adjFC.rnai.gct")
C<-read.delim("~/shRNA_database/Project_Achilles//Data/Achilles_102lines_shRNA_table.txt")
D<-read.delim("~/shRNA_database/Project_Achilles//Data/Achilles_v2.0_SampleInfo_small.txt.original.txt")

# gene symbol in RNAi screen
name<-sapply(strsplit(rownames(A),"_"),function(x){x[[1]]})
name.A<-strsplit(rownames(A),"_")
name<-c()
for(k in 1:length(name.A)){
  name<-c(name,name.A[[k]][1])
}

name.duplicate = name[which(duplicated(name)==1)]

which(name == name.duplicate[2])

# cell line sample name
name1<-colnames(A)

# expression and copy level in CCLE
library(predictiveModeling)
library(synapseClient)
id_exprLayer <- "syn1757082" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- (layer_expr$objects$eSet_expr)

id_copyLayer <- "syn1757086"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- (layer_copy$objects$eSet_copy)

# # Gistic score in CCLE : Corrupted
# id_copyLayer1 <- "syn2367236"     
# layer_copy1 <- loadEntity(id_copyLayer1)
# eSet_cgh <- (layer_copy1$objects$eSet_copy)

id_hybridLayer <-  "syn1757084" 
layer_hybrid <- loadEntity(id_hybridLayer)
eSet_hybrid <- layer_hybrid$objects$eSet_hybrid

featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr,copy = eSet_copy, mut = eSet_hybrid))    

common.sample<-intersect(colnames(featureData),colnames(A))

##### final input and response data matching with sample names in CCLE
data.input<-featureData[,common.sample]
data.response<-A[,common.sample]

data.input.exp<-exprs(eSet_expr)[,common.sample]
data.input.copy<-exprs(eSet_copy)[,common.sample]
data.input.mut<-exprs(eSet_hybrid)[,common.sample]

tissue.sample<-sapply(strsplit(common.sample,"_"),function(x){x[[2]]})

table(tissue.sample)

tissue.interest<-c("LUNG","OVARY","PANCREAS")

tissue.interest[1]

a<-which(tissue.sample == tissue.interest[2])

data.input.exp.interest<-data.input.exp[,a]
data.input.copy.interest<-data.input.copy[,a]
data.input.mut.interest<-data.input.mut[,a]

M<-apply(data.input.mut.interest,1,sum)
name.M<-names(which(M==length(a)))

data.response.target<-data.response[,a]



tempFun<-function(k){
  kk<-match(name[k],rownames(data.input.exp.interest))
  if(length(kk)>0){
    CC<-cor(data.response.target[k,],data.input.exp.interest[kk,],method = "spearman")
  }else{
    CC<-NA
  }
  return(CC)
}

library(parallel)
cc<-mclapply(1:length(name),function(x)tempFun(x),mc.cores=10)
CC<-do.call("c",cc)
hist(CC,30)


Z <- (CC-mean(CC,na.rm = T))/(sd(CC,na.rm = T))
Pval<-2*pnorm(-abs(Z))
name[which(Pval<=0.01)]
#########

gene<-intersect(name[which(Pval<=0.01)],hits.RNAi)

b1<-match(gene[3],name)
b2<-match(gene[3],rownames(data.input.exp.interest))


plot(data.response.target[b1,],data.input.exp.interest[b2,])

plot(rank(data.response.target[b1,]),rank(data.input.exp.interest[b2,]))
cor(rank(data.response.target[b1,]),rank(data.input.exp.interest[b2,]))

BBB<-read.delim("~/OV_2014/OV_Druggable_20140725.tsv")
hits.RNAi<-data.matrix(BBB$Gene_symbol[which(BBB$Untreated.median<=75)])

intersect(name[which(Pval<=0.05)],hits.RNAi)

mutFunction<-function(synID){
  aa<-synGet(synID)
  a<-read.delim(aa@filePath,skip = 2)  
  
  MAT_mut<-matrix(NA,nrow = length(unique(a$Hugo_Symbol)),ncol = length(unique(a$Tumor_Sample_Barcode)))
  rownames(MAT_mut)<-unique(a$Hugo_Symbol)
  colnames(MAT_mut)<-unique(make.names(a$Tumor_Sample_Barcode))
  
  for(k in 1:length(a$Tumor_Sample_Barcode)){
    MAT_mut[a$Hugo_Symbol[k],a$Tumor_Sample_Barcode[k]]<-as.character(as.matrix(a$Variant_Classification[k]))
  }
  return(MAT_mut)
}


M09 <- mutFunction("syn2363322") # HNSCC IlluminaGA_DNASeq 

#synID = "syn2319946"
expFunction<-function(synID){
  aa<-synGet(synID)
  a<-read.table(aa@filePath, sep="\t",header=TRUE,comment="",row.names=1)
  A<-strsplit(rownames(a),"\\|")
  genes<-c()
  for(k in 1:length(A)){
    genes<-c(genes,A[[k]][1])
  }
  aaa<-as.matrix(a)
  MAT<-log2(aaa+1)
  rownames(MAT)<-genes
  return(MAT)
}

E09<-expFunction("syn2319946")

name1<-substr(colnames(M09),1,16)
Name1<-substr(colnames(E09),1,16)
num1<-match(Name1,name1)

b<-which(!is.na(num1))
match(Name1[b],name1)
Name2<-Name1[b]
match(Name2,name1)

exp.HNSCC<-E09[-grep("\\?",rownames(E09)),b]
mut.HNSCC<-M09[,match(Name2,name1)]

countNA<-function(x){
  length(which(!is.na(x)))
}
FEQ<-apply(mut.HNSCC,1,countNA)
N<-FEQ[which(FEQ>=50)]

covFunction<-function(x){
  var(x,na.rm = T)/mean(x,na.rm = T)
}

V<-apply(exp.HNSCC,1,var)
COV<-apply(exp.HNSCC,1,covFunction)

v<-(which(COV>=1 & V>=2))
v<-(which(V>=2))
e.HNSCC<-t(apply(exp.HNSCC[v,],1,rank))
# rank transformed

rownames(e.HNSCC)<-rownames(exp.HNSCC)[v]
library("infotheo")

miCompute1<-function(k){
  a1<-match(NN[k,1],rownames(e.HNSCC))
  a2<-match(NN[k,2],rownames(e.HNSCC))
  if(is.na(a1) | is.na(a2)){
    return(NA)
  }else{
    return(mutinformation(as.matrix(e.HNSCC[a1,]),as.matrix(e.HNSCC[a2,])))
  }
}

MI1<-mclapply(1:nrow(NN),function(x)miCompute1(x),mc.cores= 12)
mi1<-do.call("c",MI1)


cmiCompute1<-function(k){
  a1<-match(NN[k,1],rownames(e.HNSCC))
  a2<-match(NN[k,2],rownames(e.HNSCC))
  
  a3<-match(names(N)[1],rownames(mut.HNSCC))
  b.mut<-which(!is.na(mut.HNSCC[a3,]))
  b.wt<-which(is.na(mut.HNSCC[a3,]))
  
  if(is.na(a1) | is.na(a2)){
    return(NA)
  }else{
    return(cbind(mutualinfo_ap(as.matrix(e.HNSCC[a1,b.mut]),as.matrix(e.HNSCC[a2,b.mut])),mutualinfo_ap(as.matrix(e.HNSCC[a1,b.wt]),as.matrix(e.HNSCC[a2,b.wt]))))
  }
}


CMI1<-mclapply(1:nrow(NN),function(x)cmiCompute1(x),mc.cores= 12)
cmi1<-do.call("rbind",CMI1)

pcmi1<-(cmi1[,1]-cmi1[,2])
ncmi1<-(cmi1[,2]-cmi1[,1])
dcmi1<-abs(cmi1[,1]-cmi1[,2])

K<-which(dcmi1>=0.4)
a3<-match(names(N)[1],rownames(mut.HNSCC))
b.mut<-which(!is.na(mut.HNSCC[a3,]))
b.wt<-which(is.na(mut.HNSCC[a3,]))

par(ask=T)
for(k in 1:length(K)){
  a1<-match(NN[K[k],1],rownames(e.HNSCC))
  a2<-match(NN[K[k],2],rownames(e.HNSCC))
  plot(e.HNSCC[a1,b.mut],e.HNSCC[a2,b.mut],pch =19, col = "red")
  points(e.HNSCC[a1,b.wt],e.HNSCC[a2,b.wt],pch = 19,col = "blue")
  rbind(cor(e.HNSCC[a1,b.wt],e.HNSCC[a2,b.wt],method="spearman"),
        cor(e.HNSCC[a1,b.mut],e.HNSCC[a2,b.mut],method="spearman"),
        cor(e.HNSCC[a1,],e.HNSCC[a2,],method="spearman"))
}