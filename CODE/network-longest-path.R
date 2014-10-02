

##Use reshape to find pairs in pathMat
library(reshape2)
##subset reachability matrix
mutPairs <- pathMat[color.red,c(color.orange,color.gold)]
##melt reachability matrix
mutPairsMelt<- melt(mutPairs)
##remove those pairs that are unreachable (infinite)
mutPairsMelt<- mutPairsMelt[!is.infinite(mutPairsMelt$value),]
dim(mutPairsMelt)

mutPairsMelt <- mutPairsMelt[duplicated(mutPairsMelt),]
mutPairsMelt$Var1 <- as.character(mutPairsMelt$Var1)
dim(mutPairsMelt)
mutPairsList <- by(mutPairsMelt,mutPairsMelt$Var1,function(x){return(x)})

maxMutPairs <- NULL
for(mut in mutPairsList){
  maxMutPairs <- rbind(maxMutPairs, mut[which.max(mut$value),])
}

save(hid1,g,gg,pathMat, maxMutPairs,file="network-analysis-ovarian.Rdata")
