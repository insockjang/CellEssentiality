rm(list = ls())

load("CCLE-Achilles.RData")

#load tissues of interest
tissue.interest<-c("LUNG","OVARY","PANCREAS","LARGE")


#BBB<-read.delim("/home/ijang/OV_2014/OV_Druggable_20140725.tsv")
#hits.RNAi<-data.matrix(BBB$Gene_symbol[which(BBB$Untreated.median<=75)])

CClist <- list()
highConfidenceCorrelations <- list()
interestingHits <- list()
highFrequencyMuts <- list()
for(tiss in tissue.interest){

  a<-which(tissue.sample == tiss)

  data.input.exp.interest<-data.input.exp[,a]
  data.input.copy.interest<-data.input.copy[,a]
  data.input.mut.interest<-data.input.mut[,a]

  M<-apply(data.input.mut.interest,1,sum)
  name.M<-names(which(M==length(a)))

  ##select tissue samples of interest
  data.response.target<-data.response[,a]

  tempFun<-function(k){
    kk<-match(name[k],rownames(data.input.exp.interest))
    if(length(kk)>0){
      CC<-cor(data.response.target[k,],data.input.exp.interest[kk,],
              method = "spearman")
      }else{
            CC<-NA
      }
    return(CC)
  }

  library(parallel)
  cc<-mclapply(1:length(name),function(x){tempFun(x)},mc.cores=10)
  CC<-do.call("c",cc)
  CClist[[tiss]]<-CC

  pdf(paste("correlation-histogram-", tiss, ".pdf", sep=""))
  hist(CC,30, main = paste("Distribution of Correlations for", tiss))
  dev.off()

  Z <- (CC-mean(CC,na.rm = T))/(sd(CC,na.rm = T))
  Pval<-2*pnorm(-abs(Z))
  highConfidenceCorrelations[[tiss]] <- name[which(Pval<=0.01)]
  #########

  gene<-intersect(name[which(Pval<=0.01)],hits.RNAi)
  interestingHits[[tiss]] <- gene
  numMuts <- apply(data.input.mut.interest,1,sum)
  threshold <- floor(0.9 *ncol(data.input.mut.interest))
  highFrequencyMuts[[tiss]]<- numMuts[numMuts > threshold]
  
#  b1<-match(gene[3],name)
#  b2<-match(gene[3],rownames(data.input.exp.interest))

#  plot(data.response.target[b1,],data.input.exp.interest[b2,])

#  plot(rank(data.response.target[b1,]),rank(data.input.exp.interest[b2,]))
#  cor(rank(data.response.target[b1,]),rank(data.input.exp.interest[b2,]))


#  intersect(name[which(Pval<=0.05)],hits.RNAi)

}



save.image("test02-laderas.RData")

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

##grab CCLE gistic results
ent <- synGet("syn2341702")
ccle <- read.delim(ent@filePath)
dim(ccle)


permuteCorrCandidates <- function(pvalvec, data.input.exp,data.response,namevec,numPermutes=10000){
  pvalOut <- list()
  corrOut <- list()
  for(pv in pvalvec){
    pvInd <- match(pv,namevec)
    print(pvInd)
    corrval <- cor(data.input.exp[pv,], data.response[pvInd,], method="spearman")
    corrDistribution <- mclapply(1:numPermutes,function(x){
      sampIdx <- sample(rownames(data.response),1)
      cor(data.input.exp[pv,], data.response[sampIdx,], method="spearman")
    }, mc.cores=10)
    corrOut[[pv]]<- corrval
    pvalOut[[pv]] <- length(which(corrDistribution > corrval))/numPermutes
  }
  pvalOut <- unlist(pvalOut)
  corrOut <- unlist(corrOut)
  return(data.frame(pvalOut,corrOut))
}

a<-which(tissue.sample == "OVARY")

data.input.exp.interest<-data.input.exp[,a]
data.input.copy.interest<-data.input.copy[,a]
data.input.mut.interest<-data.input.mut[,a]

M<-apply(data.input.mut.interest,1,sum)
name.M<-names(which(M==length(a)))

##select tissue samples of interest
data.response.target<-data.response[,a]

#ovaryCandidates <- intersect(highConfidenceCorrelations[["OVARY"]], rownames(data.response.target))

pvalOut <- permuteCorrCandidates(highConfidenceCorrelations[["OVARY"]],data.input.exp.interest, data.response.target,name)

save(highConfidenceCorrelations, data.input.exp.interest, data.response.target, name,file="permutation-analysis.RData")


