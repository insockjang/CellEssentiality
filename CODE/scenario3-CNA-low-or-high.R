library(synapseClient)
library(parallel)
synapseLogin()
setwd("~/CellEssentiality/")
load("CCLE-Achilles.RData")

options(stringsAsFactors=F)
dataResponseMap <- data.frame(response=as.character(rownames(data.response)), gene=as.character(name))


cellEssFolder <- synGet("syn2582478")

##grab GISTIC data for CCLE
ent <- synGet("syn2341702")
ccle <- read.delim(ent@filePath)
dim(ccle)

CGH.annot<-ccle[,c(1:3)]
CGH.gistic<-ccle[,-c(1:3)]

rm(ccle)

tmp1<-substring(colnames(CGH.gistic),12)
tmp2<-make.names(tmp1)
colnames(CGH.gistic)<-tmp2

#filter out to those only in expression data
CGH.gistic<-CGH.gistic[,colnames(data.input.exp)]


#look at correlation in all up/wt candidates
tempUpFun<-function(k){
  geneSym <- dataResponseMap[dataResponseMap$response ==k,"gene"]
  kk<-match(geneSym,rownames(data.input.exp.interest))
  cvk <- match(geneSym, rownames(CGH.up))
  if(length(kk)>0 & length(cvk)>0){
    up <-names(which(CGH.up[cvk,]==1))
    print(up)
    wt <- names(which(CGH.up[cvk,]!=1))
    CCCGH<-cor(data.response.target[k,up],data.input.exp.interest[kk,up],
               method = "spearman")
    CCwt <- cor(data.response.target[k,wt], data.input.exp.interest[kk,wt], 
                method="spearman")
    
  }else{
    CCwt<-NA
    CCCGH <- NA
  }
  return(c(CCwt, CCCGH))
}

#look at correlation in all down candidates
tempDownFun<-function(k){
  geneSym <- dataResponseMap[dataResponseMap$response ==k,"gene"]
  kk<-match(geneSym,rownames(data.input.exp.interest))
  cvk <- match(geneSym, rownames(CGH.down))
  if(length(kk)>0 & length(cvk)>0){
    down <-names(which(CGH.down[cvk,]==1))
    wt <- names(which(CGH.down[cvk,]!=1))
    CCCGH<-cor(data.response.target[k,down],data.input.exp.interest[kk,down],
               method = "spearman")
    CCwt <- cor(data.response.target[k,wt], data.input.exp.interest[kk,wt], 
                method="spearman")
    
  }else{
    CCwt<-NA
    CCCGH <- NA
  }
  return(c(CCwt, CCCGH))
}


#tissue.interest<-c("LUNG","OVARY","PANCREAS","LARGE")
#for(tiss in tissue.interest){
  tiss <- "OVARY"
  a<-which(tissue.sample == tiss)
  
  data.input.exp.interest<-data.input.exp[,a]
  data.input.copy.interest<-data.input.copy[,a]
  data.input.mut.interest<-data.input.mut[,a]
  ##select tissue samples of interest
  data.response.target<-data.response[,a]
  
  
  data.index <- grep(tiss,colnames(CGH.gistic))
  
  CGH.tiss <- CGH.gistic[,data.index]
  rownames(CGH.tiss) <- CGH.annot$Gene.Symbol
  CGH.up <- apply(CGH.tiss,2,function(x){ifelse(x>0,1,0)})
  CGH.down <- apply(CGH.tiss,2,function(x){ifelse(x<0,1,0)})
  CGH.none <- apply(CGH.tiss,2,function(x){ifelse(x==0,1,0)})
  numUp <- apply(CGH.up,1,sum)
  numDown <- apply(CGH.down,1,sum)
  numNone <- apply(CGH.none,1,sum)
  CGHframe <- data.frame(CGH.annot,numUp, numDown, numNone)
  
  attach(CGHframe)

  noCGHdown <- as.character(CGHframe[numDown==0,"Gene.Symbol"])
  noCGHup <- as.character(CGHframe[numUp ==  0,"Gene.Symbol"])

  detach(CGHframe)
  
  noCGHdownDR <- as.character(dataResponseMap[dataResponseMap$gene %in% noCGHdown,"response"])
  noCGHupDR <- as.character(dataResponseMap[dataResponseMap$gene %in% noCGHup, "response"])
  
  ccNoDown <- mclapply(noCGHdownDR, tempUpFun, mc.cores=10)
  ccNoUp <- mclapply(noCGHupDR, tempDownFun, mc.cores=10)
  
  CCNoDown <- data.frame(do.call(rbind,ccNoDown))
  colnames(CCNoDown) <- c("WT","CNA")
  rownames(CCNoDown) <- noCGHdownDR
  CCNoUp <- data.frame(do.call(rbind,ccNoUp))
  colnames(CCNoUp) <- c("WT","CNA")
  rownames(CCNoUp) <- noCGHupDR

  CCNoDown <- data.frame(CCNoDown, cordiff = CCNoDown$CNA - CCNoDown$WT)
  CCNoUp <- data.frame(CCNoUp, cordiff = CCNoUp$CNA - CCNoUp$WT)
  
  CCNoUp[order(abs(CCNoUp$cordiff),decreasing=TRUE),]
  CCNoDown[order(abs(CCNoDown$cordiff),decreasing=TRUE),][1:50,]
#}




