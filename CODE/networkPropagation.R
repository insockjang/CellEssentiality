# Network propagation
source("~/NetworkAnalysis_RandomWalk/networkPropagate_vector.R")
load("~/NetworkAnalysis_RandomWalk/STRING/AdjMat_90percent.Rdata")

parallelNetworkPropagate<-function(k){
  NP<-networkPropgate_vector(f.0[,k],0.7,MAT)
  return(NP)
}
library(parallel)


# provide list of genes to be tested
load("~/Synthetic_Lethality_Prediction/OV/targetedGene_CGH.Rdata")
geneList1<-intersect(colnames(MAT),correlated.gene)
geneList2<-intersect(colnames(MAT),CGH.gained.genes)
geneList3<-intersect(colnames(MAT),CGH.loss.genes)
load("~/Synthetic_Lethality_Prediction/OV/targetedGene_MUT.Rdata")
geneList4<-intersect(colnames(MAT),correlated.gene)
geneList5<-intersect(colnames(MAT),MUT.genes)

f.0<-matrix(0,nrow = nrow(MAT),ncol = 5)
rownames(f.0)<-rownames(MAT)
f.0[geneList1,1]<-1 
f.0[geneList2,2]<-1 
f.0[geneList3,3]<-1 
f.0[geneList4,4]<-1 
f.0[geneList5,5]<-1 

test<-mclapply(1:ncol(f.0),function(x)parallelNetworkPropagate(x),mc.cores=10)
Test<-do.call("cbind",test)
colnames(Test)<-c("NP.CGH.apair.gene","NP.CGH.gain.gene","NP.CGH.loss.gene","NP.MUT.apair.gene","NP.MUT.gene")
save(Test,f.0,file = "~/Synthetic_Lethality_Prediction/OV/NetProp_CGH_MUT.Rdata")
### 
