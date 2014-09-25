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
    CCwt <- cor(data.response.target[k,wt], data.input.exp.interest[kk,wt], 
                method="spearman")
    
  }else{
    CCwt<-NA
    CCmut <- NA
  }
  return(c(CCwt, CCmut))
}

library(parallel)

cc<-mclapply(mutMatch,function(x){tempFun(x)},mc.cores=10)
CC<-do.call(rbind,cc)
CC <- data.frame(CC)
rownames(CC) <- name[mutMatch]
colnames(CC)<- c("WT","Mutated")
ccFullOvary <- CClist[["OVARY"]]
names(ccFullOvary) <- name



#CC <- na.omit(CC)

##build data frame with correlations and size of each subset
CC <- data.frame(CC,full=ccFullOvary[rownames(CC)],numMuts=numMuts[rownames(CC)],numWT=ncol(data.input.mut.interest)-numMuts[rownames(CC)])

corrange <- CC$full -CC$WT

CC <- data.frame(CC, corrange)

CCfiltered <- na.omit(CC)

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