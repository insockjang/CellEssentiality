## prepare data for predictive model : CCLE and Project Archilles (training) and Sanger and Colt (testing)
## 1. matching sampels between CCLE exp and Project Archilles
## 2. matching samples between Sanger exp and Colt
## 3. matching features between CCLE exp and Sanger exp
## 4. train

## 1. matching samples between CCLE exp and Project Archilles
library(CePa)
A<-read.gct("/home/ijang/shRNA_database//Project_Achilles/Data//Achilles_102lines_gene_solutions.gct")

# gene symbol in RNAi screen
name<-sapply(strsplit(rownames(A),"_"),function(x){x[[1]]})
# cell line sample name
name1<-colnames(A)

# expression and copy level in CCLE
library(predictiveModeling)
library(synapseClient)
id_exprLayer <- "syn1757082" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- exprs(layer_expr$objects$eSet_expr)

common.sample.name<-intersect(name1,colnames(eSet_expr))

EXP.ccle<-eSet_expr[,common.sample.name]

## 2. matching samples between Sanger exp and Colt
#get COlT Data
coltData <- synGet("syn2582519")
coltDataGarp <- read.table(coltData@filePath, sep=" ", header=TRUE)
coltAnnot <- coltDataGarp[,1:4]
coltDataGarp <- coltDataGarp[,-c(1:4)]
coltCells <- make.names(toupper(colnames(coltDataGarp)))
coltCells <- gsub(".","",coltCells,fixed=TRUE)
colnames(coltDataGarp) <- coltCells

#expression data
id_exprLayer <- "1742878" 
layer_expr <- loadEntity(id_exprLayer)
eSet_Sanger <- exprs(layer_expr$objects$eSet_expr)


inBoth <- intersect(colnames(coltDataGarp), colnames(eSet_Sanger))

EXP.sanger<-eSet_Sanger[,inBoth]

# validating Cell Viability from Colt database
final.CV.Colt<-coltDataGarp[,inBoth]
rownames(final.CV.Colt)<-coltAnnot$Gene.name

## 3. matching features between CCLE exp and Sanger exp
common.feature.name<-intersect(rownames(EXP.ccle),rownames(EXP.sanger))

final.EXP.ccle<-EXP.ccle[common.feature.name,]
final.EXP.sanger<-EXP.sanger[common.feature.name,]


## 4. build predictive model with training input and response data (CCLE exp and Project Archilles)
# training response dataset from Project Archilles
final.CV.ProjA<-A[,common.sample.name]
rownames(final.CV.ProjA)<-name

# testing response dataset from Colt

source("~/PredictiveModel_pipeline/R5/crossValidatePredictiveModel1.R")
source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")

X<-final.EXP.ccle

VCV_Func<-function(x){  
  Y<-final.CV.ProjA[x,]
  bsModel <- myEnetModel1$new()        
  bsModel$customTrain(t(X), Y, alpha = 1, nfolds = 10)  
  return(bsModel)  
}

tmp.lasso<-mclapply(1:nrow(final.CV.ProjA),function(x)VCV_Func(x),mc.cores= 10)
Pred<-lapply(1:length(tmp),function(x){return(tmp[[x]]$customPredict(t(final.EXP.sanger)))})
Pred.Lasso<-do.call("cbind",Pred)


VCV_Func_enet<-function(x){  
  Y<-final.CV.ProjA[x,]
  bsModel <- myEnetModel1$new()        
  bsModel$customTrain(t(X), Y, alpha = 0.005, nfolds = 10)  
  return(bsModel)  
}

tmp.enet<-mclapply(1:nrow(final.CV.ProjA),function(x)VCV_Func_enet(x),mc.cores= 15)
Pred1<-lapply(1:length(tmp.enet),function(x){return(tmp.enet[[x]]$customPredict(t(final.EXP.sanger)))})
Pred.ENet<-do.call("cbind",Pred1)

VCV_Func_ridge<-function(x){  
  Y<-final.CV.ProjA[x,]
  bsModel <- myEnetModel1$new()        
  bsModel$customTrain(t(X), Y, alpha = 0, nfolds = 10)  
  return(bsModel)  
}

tmp.ridge<-mclapply(1:nrow(final.CV.ProjA),function(x)VCV_Func_ridge(x),mc.cores= 15)
Pred2<-lapply(1:length(tmp.ridge),function(x){return(tmp.ridge[[x]]$customPredict(t(final.EXP.sanger)))})
Pred.Ridge<-do.call("cbind",Pred2)

colnames(Pred.ENet)<-colnames(Pred.Lasso)<-rownames(final.CV.ProjA)
# save(tmp.lasso,tmp.ridge,tmp.enet, Pred.Lasso,Pred.Ridge, Pred.ENet, file = "~/shRNA_database/predictiveModel.Rdata")

save(tmp.lasso, Pred.Lasso,file = "~/shRNA_database/predictiveModel_Lasso.Rdata")
save(tmp.ridge, Pred.Ridge,file = "~/shRNA_database/predictiveModel_Ridge.Rdata")
save(tmp.enet, Pred.ENet,file = "~/shRNA_database/predictiveModel_ENet.Rdata")



common.CV.feature<-intersect(name,rownames(final.CV.ProjA))

COR<-lapply(1:length(common.CV.feature),function(x){
  if(var(Pred.ENet[,common.CV.feature[x]])==0){
    return(NA)
  }else{
    return(cor((Pred.ENet[,common.CV.feature[x]]),as.numeric(final.CV.Colt[which(name == common.CV.feature[x]),]),method="spearman"))
  }
}
)


COR<-c()
for(x in 1:length(common.CV.feature)){
  if(var(Pred.ENet[,common.CV.feature[x]])==0){
    COR<-c(COR,NA)
  }else{
    if(length(which(name == common.CV.feature[x]))==1){
      COR<-c(COR,cor((Pred.ENet[,common.CV.feature[x]]),as.numeric(final.CV.Colt[which(name == common.CV.feature[x]),]),method="spearman"))
    }else{
      COR<-c(COR,cor((Pred.ENet[,common.CV.feature[x]]),as.numeric(final.CV.Colt[which(name == common.CV.feature[x])[1],]),method="spearman"))
    }
    
  }
}
  
hist(COR,50)

######### correlation of correlation in order to check how two CVs are related
name.unique<-name[-which(duplicated(name)==1)]
final.CV.Colt.unique<-final.CV.Colt[-which(duplicated(name)==1),]
rownames(final.CV.Colt.unique)<-name.unique

CV1<-final.CV.ProjA[common.CV.feature,]
CV2<-final.CV.Colt.unique[common.CV.feature,]

corCV1<-cor(t(CV1),method = "spearman")
corCV2<-cor(t(CV2),method = "spearman")

CC<-lapply(1:nrow(corCV1),function(x){cor(corCV1[x,],corCV2[x,],method = "spearman")})
CorCor<-do.call("c",CC)
hist(CorCor,50)