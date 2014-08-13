rm(list = ls())
# Data preparation
library(CePa)
A<-read.gct("~/shRNA_database//Project_Achilles/Data//Achilles_102lines_gene_solutions.gct")
B<-read.gct("~/shRNA_database//Project_Achilles/Data//20110303_achilles2_PMAD_adjFC.rnai.gct")
C<-read.delim("~/shRNA_database/Project_Achilles//Data/Achilles_102lines_shRNA_table.txt")
D<-read.delim("~/shRNA_database/Project_Achilles//Data/Achilles_v2.0_SampleInfo_small.txt.original.txt")

# gene symbol in RNAi screen
name<-sapply(strsplit(rownames(A),"_"),function(x){x[[1]]})
name.duplicate = name[which(duplicated(name)==1)]


# cell line sample name
name1<-colnames(A)

# expression and copy level in CCLE
library(predictiveModeling)
library(synapseClient)

id_exprLayer <- "syn1757082" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- exprs(layer_expr$objects$eSet_expr)

## CGH gistic score in CCLE
ent <- synGet("syn2341702")
ccle <- read.delim(ent@filePath)
dim(ccle)
CGH.annot<-ccle[,c(1:3)]
CGH.gistic<-ccle[,-c(1:3)]

tmp1<-substring(colnames(CGH.gistic),12)
tmp2<-make.names(tmp1)
colnames(CGH.gistic)<-tmp2


id_hybridLayer <-  "syn1757084" 
layer_hybrid <- loadEntity(id_hybridLayer)
eSet_hybrid <- exprs(layer_hybrid$objects$eSet_hybrid)


common.sample<-intersect(intersect(colnames(eSet_expr),colnames(CGH.gistic)),colnames(A))

##### final input and response data matching with sample names in CCLE
data.response<-A[,common.sample]
data.input.exp<-eSet_expr[,common.sample]
data.input.copy<-CGH.gistic[,common.sample]

tissue.sample<-sapply(strsplit(common.sample,"_"),function(x){x[[2]]})

table(tissue.sample)

tissue.interest<-c("LUNG","OVARY","PANCREAS")

tissue.interest[2]

a<-which(tissue.sample == tissue.interest[2])

data.input.exp.interest<-data.input.exp[,a]
data.input.copy.interest<-data.input.copy[,a]
data.response.target<-data.response[,a]

# find condition which we might be interested
# First, I would like to choose CGH gain through targeted samples and later, I will choose CGH lost.

M<-apply(data.input.copy.interest,1,mean)
Z1 <- (M-mean(M,na.rm = T))/(sd(M,na.rm = T))
Pval1<-2*pnorm(-abs(Z1))

aa1<-which(Pval1<=0.05 & Z1>0)
aa2<-which(Pval1<=0.05 & Z1<0)

CGH.gained.genes<-CGH.annot[aa1,1]
CGH.loss.genes<-CGH.annot[aa2,1]

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


P.threshold = 0.05
Z <- (CC-mean(CC,na.rm = T))/(sd(CC,na.rm = T))
Pval<-2*pnorm(-abs(Z))
correlated.gene<-name[which(Pval<=P.threshold)]



#########
BBB<-read.delim("~/OV_2014/OV_Druggable_20140725.tsv")
hits.RNAi<-data.matrix(BBB$Gene_symbol[which(BBB$Untreated.median<=75)])

gene<-intersect(name[which(Pval<=P.threshold)],hits.RNAi)

b1<-match(gene[1],name)
b2<-match(gene[1],rownames(data.input.exp.interest))
plot(data.response.target[b1,],data.input.exp.interest[b2,])
plot(rank(data.response.target[b1,]),rank(data.input.exp.interest[b2,]))
cor(rank(data.response.target[b1,]),rank(data.input.exp.interest[b2,]))


# Now we have candidate pairs between 
save(correlated.gene,CGH.gained.genes,CGH.loss.genes,file ="~/Synthetic_Lethality_Prediction/OV/targetedGene_CGH.Rdata")