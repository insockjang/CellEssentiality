# pairs between no CGH vs. gained CGH
##look for those pairs that have the highest correlation difference between WT and mutated
aa<-synGet("syn2744664")
load(aa@filePath)

# load("~/shRNA_database//Project_Achilles/CCLE-Achilles.RData")
#concentrate on OVARY tissue first.
tiss <- "PANCREAS" # "OVARY", "LARGE", "PANCREAS" (23,15,11) 
a<-which(tissue.sample == tiss)

rownames(data.response)<-name

data.input.exp.interest<-data.input.exp[,a]
data.input.copy.interest<-data.input.copy[,a]
data.input.mut.interest<-data.input.mut[,a]
data.response.target<-data.response[,a]


# correlation between exrpession and cell viability
common.name.exp<-intersect(name,rownames(data.input.exp.interest))
K<-lapply(1:length(common.name.exp),function(x){kk<-cor.test(data.input.exp.interest[common.name.exp[x],],data.response.target[common.name.exp[x],])
                                                return(c(kk$estimate,kk$p.value))})
KK<-do.call("rbind",K)
rownames(KK)<-common.name.exp
colnames(KK)<-c("cor.estimate","cor.p.value")

tmpFunGain<-function(k){#} in 1:nrow(data.input.mut.interest)){
  D<-data.input.copy.interest[k,]
  d0<-which(D<=1 & D>=-1)
  d1<-which(D>1)
  if(length(d0)>2 & length(d1)>2){
    tmp0<-lapply(1:length(common.name.exp),function(x){
      kk<-cor.test(data.input.exp.interest[common.name.exp[x],d0],data.response.target[common.name.exp[x],d0])
      return(c(kk$estimate,kk$p.value))
    }
    )
    tmp1<-lapply(1:length(common.name.exp),function(x){
      kk<-cor.test(data.input.exp.interest[common.name.exp[x],d1],data.response.target[common.name.exp[x],d1])
      return(c(kk$estimate,kk$p.value))
    }
    )
    Tmp0<-do.call("rbind",tmp0)
    Tmp1<-do.call("rbind",tmp1)
    rownames(Tmp0)<-rownames(Tmp1)<-common.name.exp
    colnames(Tmp0)<-c("WT.cor.estimate","WT.cor.p.value")    
    colnames(Tmp1)<-c("Gain.cor.estimate","Gain.cor.p.value")    
    
    return(cbind(Tmp0,Tmp1))
  }else{
    #     TMP[[k]]<-NA    
    return(NA)
  }  
}

countGain<-function(k){#} in 1:nrow(data.input.mut.interest)){
  D<-data.input.copy.interest[k,]
  d0<-which(D<=1 & D>=-1)
  d1<-length(which(D>1))
  return(d1)
}

countLoss<-function(k){#} in 1:nrow(data.input.mut.interest)){
  D<-data.input.copy.interest[k,]
  d0<-which(D<=1 & D>=-1)
  d1<-length(which(D< -1))
  return(d1)
}


lossCount <-lapply(1:nrow(data.input.copy.interest),function(x){countLoss(x)})
names(lossCount)<- rownames(data.input.copy.interest)
lossCount <- unlist(lossCount)
lossCount <- data.frame(gene=names(lossCount),lossCount)

gainCount <- lapply(1:nrow(data.input.copy.interest), function(x){countGain(x)})
names(gainCount)<- rownames(data.input.copy.interest)
gainCount <- unlist(gainCount)
gainCount <- data.frame(gene=names(gainCount), gainCount)

tmpFunLoss<-function(k){#} in 1:nrow(data.input.mut.interest)){
  D<-data.input.copy.interest[k,]
  d0<-which(D<=1 & D>=-1)
  d1<-which(D< -1)
  if(length(d0)>2 & length(d1)>2){
    tmp0<-lapply(1:length(common.name.exp),function(x){
      kk<-cor.test(data.input.exp.interest[common.name.exp[x],d0],data.response.target[common.name.exp[x],d0])
      return(c(kk$estimate,kk$p.value))
    }
    )
    tmp1<-lapply(1:length(common.name.exp),function(x){
      kk<-cor.test(data.input.exp.interest[common.name.exp[x],d1],data.response.target[common.name.exp[x],d1])
      return(c(kk$estimate,kk$p.value))
    }
    )
    Tmp0<-do.call("rbind",tmp0)
    Tmp1<-do.call("rbind",tmp1)
    rownames(Tmp0)<-rownames(Tmp1)<-common.name.exp
    colnames(Tmp0)<-c("WT.cor.estimate","WT.cor.p.value")    
    colnames(Tmp1)<-c("Loss.cor.estimate","Loss.cor.p.value")    
    return(cbind(Tmp0,Tmp1))
  }else{
    return(NA)
  }  
}

cc.gain<-mclapply(1:nrow(data.input.copy.interest),function(x){tmpFunGain(x)},mc.cores=5)
cc.loss<-mclapply(1:nrow(data.input.copy.interest),function(x){tmpFunLoss(x)},mc.cores=5)

S.gain<-c()
S.loss<-c()
for(k in 1:nrow(data.input.copy.interest)){
  D<-data.input.copy.interest[k,]
  d0<-which(D<=1 & D>=-1)
  d1<-which(D> 1)
  if(length(d0)>3 & length(d1)>3){
    S.gain<-c(S.gain,k)
  }
  
  
  d0<-which(D<=1 & D>=-1)
  d2<-which(D< -1)
  
  if(length(d0)>3 & length(d2)>3){
    S.loss<-c(S.loss,k)
  }
  
}

name.apair.gain<-rownames(data.input.copy.interest)[S.gain]
name.apair.loss<-rownames(data.input.copy.interest)[S.loss]

for(k in 1:nrow(data.input.copy.interest)){
  if(!is.na(cc.loss[[k]])){
    #     C<-cbind(cc[[k]],KK)
    C<-cc.gain[[k]]
    
    k3<-which(C[,2]<=0.1 & C[,4]<=0.1)    
    E<-C[k3,]    
    e<-(C[k3,1]-C[k3,3])
    FF<-E[which(e>1 | e< -1),]
    
    N<-rownames(FF)
    for(kk in 1:nrow(FF)){
      MM<-N[kk]
      par(ask=TRUE)
      r1<-rank(data.input.exp.interest[MM,])
      r2<-rank(data.response.target[MM,])
      plot(r1,r2,col = "blue",pch =19)
      points(r1[data.input.copy.interest[k,]> 1],r2[data.input.copy.interest[k,]> 1],col = "red",pch = 19)
      plot(data.input.exp.interest[MM,],data.response.target[MM,],pch = 19,col = "blue",main = MM)
      points(data.input.exp.interest[MM,data.input.copy.interest[k,]> 1],data.response.target[MM,data.input.copy.interest[k,]> 1],col = "red",pch = 19)      
    }   
  }
}



HowMany.gain<-matrix(0,nrow = length(S.gain),ncol = 1)
HowMany.loss<-matrix(0,nrow = length(S.loss),ncol = 1)
for(k in 1:length(S.gain)){
  
  C<-cc.gain[[S.gain[k]]]  
  k3<-which(C[,2]<=0.1 & C[,4]<=0.1)    
  E<-C[k3,]    
  e<-(C[k3,1]-C[k3,3])
  
  HowMany.gain[k,1]<-length(which(e>1 | e< -1))
  print(k)
}
rownames(HowMany.gain)<-name.apair.gain
for(k in 1:length(S.loss)){
  
  C<-cc.loss[[S.loss[k]]]  
  k3<-which(C[,2]<=0.1 & C[,4]<=0.1)    
  E<-C[k3,]    
  e<-(C[k3,1]-C[k3,3])
  
  HowMany.loss[k,1]<-length(which(e>1 | e< -1))
  print(k)
}
rownames(HowMany.loss)<-name.apair.loss

save(KK,cc.gain,cc.loss,HowMany.gain,HowMany.loss,name.apair.gain,name.apair.loss,file = "~/shRNA_database//Project_Achilles/PANCREAS_scenario3.Rdata")
