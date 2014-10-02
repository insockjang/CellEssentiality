##look for those pairs that have the highest correlation difference between WT and mutated
aa<-synGet("syn2744664")
load(aa@filePath)

# load("~/shRNA_database//Project_Achilles/CCLE-Achilles.RData")
#concentrate on OVARY tissue first.
tiss <- "OVARY" # "OVARY", "LARGE", "PANCREAS" (23,15,11) 
a<-which(tissue.sample == tiss)

rownames(data.response)<-name

data.input.exp.interest<-data.input.exp[,a]
data.input.copy.interest<-data.input.copy[,a]
data.input.mut.interest<-data.input.mut[,a]
data.response.target<-data.response[,a]



common.name.exp<-intersect(name,rownames(data.input.exp.interest))
K<-lapply(1:length(common.name.exp),function(x){kk<-cor.test(data.input.exp.interest[common.name.exp[x],],data.response.target[common.name.exp[x],])
                                                return(c(kk$estimate,kk$p.value))})


KK<-do.call("rbind",K)
rownames(KK)<-common.name.exp
colnames(KK)<-c("cor.estimate","cor.p.value")

tmpFun<-function(k){#} in 1:nrow(data.input.mut.interest)){
  D<-data.input.mut.interest[k,]
  d0<-which(D==0)
  d1<-which(D==1)
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
    colnames(Tmp1)<-c("MT.cor.estimate","MT.cor.p.value")    
#     TMP[[k]]<-Tmp
    return(cbind(Tmp0,Tmp1))
  }else{
#     TMP[[k]]<-NA    
    return(NA)
  }  
}

cc<-mclapply(1:nrow(data.input.mut.interest),function(x){tmpFun(x)},mc.cores=10)

S<-c()
for(k in 1:nrow(data.input.mut.interest)){
  if(!is.na(cc[[k]])){
    S<-c(S,k)
  }
}
    
name.apair<-rownames(data.input.mut.interest)[S]

for(k in 1:nrow(data.input.mut.interest)){
  if(!is.na(cc[[k]])){
#     C<-cbind(cc[[k]],KK)
    C<-cc[[k]]
    
    k3<-which(C[,2]<=0.1 & C[,4]<=0.1)    
    E<-C[k3,]    
    e<-(C[k3,1]-C[k3,3])
    FF<-E[which(e>1 | e< -1),]
    
    N<-rownames(FF)
    for(kk in 1:nrow(FF)){
      MM<-N[kk]
      par(ask=TRUE)
      plot(data.input.exp.interest[MM,],data.response.target[MM,],pch = 19,col = "blue",main = MM)
      points(data.input.exp.interest[MM,data.input.mut.interest[k,]==1],data.response.target[MM,data.input.mut.interest[k,]==1],col = "red",pch = 19)      
    }   
  }
}



HowMany<-matrix(0,nrow = length(S),ncol = 1)
for(k in 1:length(S)){
  
#   C<-cbind(cc[[S[k]]],KK)
  C<-cc[[S[k]]]
  
    k3<-which(C[,2]<=0.1 & C[,4]<=0.1)    
    E<-C[k3,]    
    e<-(C[k3,1]-C[k3,3])
      
    
    HowMany[k,1]<-length(which(e>1 | e< -1))
  print(k)
}
rownames(HowMany)<-name.apair

save(KK,cc,HowMany,name.apair,file = "~/shRNA_database//Project_Achilles/scenario2.Rdata")
