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


tempUpS3pval<-function(k){
  geneSym <- dataResponseMap[dataResponseMap$response ==k,"gene"]
  kk<-match(geneSym,rownames(data.input.exp.interest))
  print(kk)
  cvk <- match(geneSym, rownames(CGH.up))
  print(cvk)
  if(length(kk)>0 & length(cvk)>0 & !is.na(kk) & !is.na(cvk)){
    up <-names(which(CGH.up[cvk,]==1))
    print(up)
  #  wt <- names(which(CGH.up[cvk,]!=1))
    CCCGH<-cor(data.response.target[k,up],data.input.exp.interest[kk,up],
               method = "spearman")
    CCfull <- cor(data.response.target[k,], data.input.exp.interest[kk,], 
                method="spearman")
    if(length(up)>2){
      pvalFull <- cor.test(data.response.target[k,], data.input.exp.interest[kk,], method=
                             "spearman", alternative="two.sided")$p.value
      pvalUp <- cor.test(data.response.target[k,up], data.input.exp.interest[kk,up],
                          method = "spearman", alternative="two.sided")$p.value
    }
    else{pvalUp <- NA
         pvalFull <- NA
    }
    
  }else{
    CCfull<-NA
    pvalUp <- NA
    CCCGH <- NA
    pvalFull <- NA
  }
  return(c(CCfull, pvalFull, CCCGH, pvalUp))
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




#look at correlation in all down candidates - compare to full correlation
#
tempDownS3pval<-function(k){
  geneSym <- dataResponseMap[dataResponseMap$response ==k,"gene"]
  kk<-match(geneSym,rownames(data.input.exp.interest))
  cvk <- match(geneSym, rownames(CGH.down))
  if(length(kk)>0 & length(cvk)>0){
    down <-names(which(CGH.down[cvk,]==1))
    wt <- names(which(CGH.down[cvk,]!=1))
    CCCGH<-cor(data.response.target[k,down],data.input.exp.interest[kk,down],
               method = "spearman")
    CCfull <- cor(data.response.target[k,], data.input.exp.interest[kk,], 
                method="spearman")

    if(length(down)>2){
      pvalFull <- cor.test(data.response.target[k,], data.input.exp.interest[kk,], method=
                             "spearman", alternative="two.sided")$p.value
      pvalDown <- cor.test(data.response.target[k,up], data.input.exp.interest[kk,up],
                         method = "spearman", alternative="two.sided")$p.value
    }
    else{pvalDown <- NA
         pvalFull <- NA
    }
    
    
  }else{
    CCfull<-NA
    pvalDown <- NA
    CCCGH <- NA
    pvalFull <- NA
  }
  return(c(CCfull,pvalFull, CCCGH, pvalDown))
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
  
  #ccNoDown <- mclapply(noCGHdownDR, tempUpFun, mc.cores=10)
  #ccNoUp <- mclapply(noCGHupDR, tempDownFun, mc.cores=10)
  
  ccNoDown <- mclapply(noCGHdownDR, tempUpS3pval, mc.cores=10)
  ccNoUp <- mclapply(noCGHupDR, tempDownS3pval, mc.cores=10)

  CCNoDown <- data.frame(do.call(rbind,ccNoDown))
  colnames(CCNoDown) <- c("Full", "pvalFull","CGHup", "pvalUp")
  rownames(CCNoDown) <- noCGHdownDR
  CCNoUp <- data.frame(do.call(rbind,ccNoUp))
  colnames(CCNoUp) <- c("Full", "pvalFull", "CGHdown", "pvalDown")
  rownames(CCNoUp) <- noCGHupDR

  CCNoDown <- data.frame(CCNoDown, cordiff = CCNoDown$CGHup - CCNoDown$Full)
  CCNoUp <- data.frame(CCNoUp, cordiff = CCNoUp$CGHdown - CCNoUp$Full)
  
#  CCNoUp <- CCNoUp[order(abs(CCNoUp$cordiff),decreasing=TRUE),]
#  CCNoDown <- CCNoDown[order(abs(CCNoDown$cordiff),decreasing=TRUE),]

  CCNoUp <- data.frame(gene=dataResponseMap[rownames(CCNoUp), "gene"],CCNoUp)
  CCNoDown <-  data.frame(gene=dataResponseMap[rownames(CCNoDown), "gene"],CCNoDown)

CCHalfDown <- na.omit(merge(CCNoUp, CGHframe, by.x="gene", by.y="Gene.Symbol"))
CCHalfUp <- na.omit(merge(CCNoDown, CGHframe, by.x="gene", by.y="Gene.Symbol"))

CCHalfDown <- CCHalfDown[CCHalfDown$pvalDown < 0.1,]
CCHalfUp <- CCHalfUp[CCHalfUp$pvalUp < 0.1,]

write.table(CCHalfDown, "Scenario3-half-down-half-wt-pvalue.txt", quote=F, sep="\t", row.name=T)
write.table(CCHalfUp, "Scenario3-half-up-half-wt.pvalue.txt", quote=F, sep="\t", row.name=T)

#}




