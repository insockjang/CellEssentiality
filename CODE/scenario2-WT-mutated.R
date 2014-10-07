##look for those pairs that have the highest correlation difference between WT and mutated
load("CCLE-Achilles.RData")

#concentrate on OVARY tissue first.
tiss <- "OVARY"  
a<-which(tissue.sample == tiss)

data.input.exp.interest<-data.input.exp[,a]
data.input.copy.interest<-data.input.copy[,a]
data.input.mut.interest<-data.input.mut[,a]

M<-apply(data.input.mut.interest,1,sum)
name.M<-names(which(M==length(a)))

numMuts <- apply(data.input.mut.interest,1,sum)
threshold <- 3
highFreqMuts <-which(numMuts > threshold)



##select tissue samples of interest
data.response.target<-data.response[,a]

mutNames <- rownames(data.input.mut.interest)

mutMatch <- lapply(mutNames, function(x){match(x, name)})
mutMatch <- do.call("c",mutMatch)
mutMatch <- na.omit(mutMatch)

##modified version of tempFun from test02.R
#this function finds the indices of WT and mutated across samples
#and calculates correlation in expression and shRNA CV in these subsets
tempFun<-function(k){
  kk<-match(name[k],rownames(data.input.exp.interest))
  cvk <- match(name[k], rownames(data.input.mut.interest))
  if(length(kk)>0){
    muts <-which(data.input.mut.interest[cvk,]==1)
    wt <- which(data.input.mut.interest[cvk,]!=1)
    CCmut<-cor(data.response.target[k,muts],data.input.exp.interest[kk,muts],
               method = "spearman")
#    pvalMut <- cor.test(data.response.target[k,muts], data.input.exp.interest[kk,muts],
#                        method = "spearman")$p.value
    CCwt <- cor(data.response.target[k,wt], data.input.exp.interest[kk,wt], 
                method="spearman")
#    pvalWt <- cor.test(data.response.target[k,wt], data.input.exp.interest[kk,wt], method=
#                         "spearman")$p.value
    
  }else{
    CCwt<-NA
 #   pvalMut <- NA
    CCmut <- NA
  #  pvalWt <- NA
  }
#  return(c(CCwt,pvalWt, CCmut, pvalMut))
  return(c(CCwt,CCmut))
}


tempFunPval<-function(k){
  print(k)
  kk<-match(name[k],rownames(data.input.exp.interest))
  cvk <- match(name[k], rownames(data.input.mut.interest))
  if(length(kk)>0 & !is.na(kk)){
    muts <-which(data.input.mut.interest[cvk,]==1)
    wt <- which(data.input.mut.interest[cvk,]!=1)

    CCmut<-cor(data.response.target[k,muts],data.input.exp.interest[kk,muts],
               method = "spearman")
    
    CCwt <- cor(data.response.target[k,wt], data.input.exp.interest[kk,wt], 
                method="spearman")
    
    if(length(muts)>2 & length(wt)>2){
        pvalWt <- cor.test(data.response.target[k,wt], data.input.exp.interest[kk,wt], method=
                             "spearman", alternative="two.sided")$p.value
        pvalMut <- cor.test(data.response.target[k,muts], data.input.exp.interest[kk,muts],
                        method = "spearman", alternative="two.sided")$p.value
    }
    else{pvalMut <- NA
      pvalWt <- NA
    }
    
  }else{
    CCwt<-NA
       pvalMut <- NA
    CCmut <- NA
      pvalWt <- NA
  }
    return(c(CCwt,pvalWt, CCmut, pvalMut))
  #return(c(CCwt,CCmut))
}


tempFunPvalFull<-function(k){
  print(k)
  kk<-match(name[k],rownames(data.input.exp.interest))
  cvk <- match(name[k], rownames(data.input.mut.interest))
  if(length(kk)>0 & !is.na(kk)){
    muts <-which(data.input.mut.interest[cvk,]==1)
    #wt <- which(data.input.mut.interest[cvk,]!=1)
    
    CCmut<-cor(data.response.target[k,muts],data.input.exp.interest[kk,muts],
               method = "spearman")
    
    CCfull <- cor(data.response.target[k,], data.input.exp.interest[kk,], 
                method="spearman")
    
    if(length(muts)>2){
      pvalFull <- cor.test(data.response.target[k,], data.input.exp.interest[kk,], method=
                           "spearman", alternative="two.sided")$p.value
      pvalMut <- cor.test(data.response.target[k,muts], data.input.exp.interest[kk,muts],
                          method = "spearman", alternative="two.sided")$p.value
    }
    else{pvalMut <- NA
         pvalFull <- NA
    }
    
  }else{
    CCfull<-NA
    pvalMut <- NA
    CCmut <- NA
    pvalFull <- NA
  }
  return(c(CCfull,pvalFull, CCmut, pvalMut))
  #return(c(CCwt,CCmut))
}




library(parallel)


cc<-lapply(mutMatch[1:50],function(x){tempFunPval(x)})


cc<-mclapply(mutMatch,function(x){tempFunPval(x)},mc.cores=10)



CC<-do.call(rbind,cc)
CC <- data.frame(CC)
rownames(CC) <- name[mutMatch]
colnames(CC)<- c("WT","pvalMut","Mutated", "pvalWt")
ccFullOvary <- CClist[["OVARY"]]
names(ccFullOvary) <- name



#CC <- na.omit(CC)

##build data frame with correlations and size of each subset
CC <- data.frame(CC,full=ccFullOvary[rownames(CC)],numMuts=numMuts[rownames(CC)],numWT=ncol(data.input.mut.interest)-numMuts[rownames(CC)])

corrange <- CC$full -CC$WT

CC <- data.frame(CC, corrange)

CCfiltered <- na.omit(CC)
CCfiltered[CCfiltered$pvalMut < 0.1,]
CCfiltered[CCfiltered$pvalWt < 0.05,]


cordiff <- CCfiltered$corrange
Z <- (cordiff-mean(cordiff,na.rm = T))/(sd(cordiff,na.rm = T))
Pval<-2*pnorm(-abs(Z))
qval <- p.adjust(Pval, method="BH")
CCfiltered <- data.frame(CCfiltered, pvalue=Pval, qvalue=qval)
write.table(CCfiltered, "Scenario2-Mut-WT-corrdiff.txt", quote=F, sep="\t", row.name=T)

#visualize correlations on histogram
p1 <- hist(CCfiltered$corrange)
p2 <- hist(CCfiltered$Mutated)
pdf("distribution-of-correlations-full-minus-wt.pdf")
plot(p1, col=rgb(0,0,1,1/4), xlim=c(-1,1), main="Distribution of Correlations of Full-WT (Blue)")
dev.off()

pdf("distribution-of-correlations-mutated.pdf")
plot(p2, col=rgb(1,0,0,1/4), xlim=c(-1,1), main="Distribution of Correlations of Mutated (Pink)")
dev.off()

#look at difference between WT and mutated samples
pdf("range-correlations-btw-mut-and-wt.pdf")
hist(corrange, main="Ovary Correlation Differences btw WT and Mutated samples")
dev.off()
highCorRange <- sort(corrange[abs(corrange)>.75])
intersect(highConfidenceCorrelations[["OVARY"]],names(highCorRange))
write.table(CC[names(highCorRange),], "highest-cor-wt-mut-candidates.txt", quote=F, sep="\t")


cc <- mclapply(mutMatch, function(x){tempFunPvalFull(x)}, mc.cores=2)
CC<-do.call(rbind,cc)
CC <- data.frame(CC)
rownames(CC) <- name[mutMatch]
colnames(CC)<- c("Full","pvalFull","Mutated", "pvalMut")
ccFullOvary <- CClist[["OVARY"]]
names(ccFullOvary) <- name

CC <- data.frame(CC,full=ccFullOvary[rownames(CC)],numMuts=numMuts[rownames(CC)],numWT=ncol(data.input.mut.interest)-numMuts[rownames(CC)])

corrange <- CC$Full -CC$Mutated

CC <- data.frame(CC, corrange)

CCfiltered <- CC[CC$numMuts > 3,]
CCfiltered <- na.omit(CCfiltered)
write.table(CCfiltered[CCfiltered$pvalMut < 0.1,], "scenario2-significantly-cor-half-mutated-genes.txt", quote=F, sep="\t")

