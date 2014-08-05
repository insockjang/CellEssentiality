# Data preparation
library(CePa)
A<-read.gct("~/shRNA_database//Project_Achilles/Data//Achilles_102lines_gene_solutions.gct")
B<-read.gct("~/shRNA_database//Project_Achilles/Data//20110303_achilles2_PMAD_adjFC.rnai.gct")
C<-read.delim("~/shRNA_database/Project_Achilles//Data/Achilles_102lines_shRNA_table.txt")
D<-read.delim("~/shRNA_database/Project_Achilles//Data/Achilles_v2.0_SampleInfo_small.txt.original.txt")

# gene symbol in RNAi screen
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
id_exprLayer <- "syn1757082" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- exprs(layer_expr$objects$eSet_expr)

id_copyLayer <- "syn1757086"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- exprs(layer_copy$objects$eSet_copy)



a1<-match(name1,colnames(eSet_expr))
b1<-match(name1,colnames(eSet_copy))



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