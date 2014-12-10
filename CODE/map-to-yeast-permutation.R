#load("/gluster/home/ijang/shRNA_database/Project_Achilles/PANCREAS_scenario3.Rdata")

homeDir <- "/gluster/home/tladeras/CellEssentiality/"

#load in orthology files
yeastOrtho <- read.delim(paste(homeDir,"yeast-human-orthologs.txt", sep=""))
colnames(yeastOrtho)
hist(yeastOrtho$X..Identity.with.respect.to.query.gene)
hist(yeastOrtho$X..Identity.with.respect.to.Yeast.gene)

#grab only those inviable synthetic lethal interactions
YeastInt <- read.delim(paste(homeDir,"BIOGRID-ORGANISM-Saccharomyces_cerevisiae-3.2.117.tab2.txt", sep=""))
YeastSynthLeth <- YeastInt[YeastInt$Experimental.System == "Synthetic Lethality",]
YeastSynthLeth <- YeastSynthLeth[grep("inviable", YeastSynthLeth$Phenotypes),]
rm(YeastInt)

numPermutes <- 10000
yeastOrthoF <- yeastOrtho[yeastOrtho$X..Identity.with.respect.to.query.gene > 40,]
permuteCandidates <- function(ind, numPairs=22){
    print(ind)
    YeastX <- yeastOrthoF[sample(1:nrow(yeastOrthoF), numPairs),"Yeast.Ensembl.Gene.ID"]
    YeastY <- yeastOrtho[sample(1:nrow(yeastOrthoF),numPairs),"Yeast.Ensembl.Gene.ID"]
    yeastCand <- data.frame(YeastX, YeastY)
    
    attach(yeastCand)
    map1 <- yeastCand[YeastX %in% YeastSynthLeth$Systematic.Name.Interactor.A & YeastX %in% YeastSynthLeth$Systematic.Name.Interactor.B,]
    map2 <- yeastCand[YeastY %in% YeastSynthLeth$Systematic.Name.Interactor.A & YeastY %in% YeastSynthLeth$Systematic.Name.Interactor.B,]
    detach(yeastCand)
    outMap <- rbind(map1,map2)
    #print(outMap)
    
    return(nrow(outMap[!duplicated(outMap),]))
}

library(parallel)

#run Permutation 1
#pull random 22 genes from 

obs <- 22
permuteResults <- mclapply(1:numPermutes, function(x,...){permuteCandidates(x)}, 
	mc.cores=5, numPairs=obs)
distYeast <- unlist(permuteResults)
hist(distYeast)

#return pvalue
pval <- length(which(distYeast >= obs))/numPermutes

png("ovary-permutations-yeast-orthologs.png")
hist(distYeast, xlab="Number of Synthetic Pairs Picked", main ="Permutations of 22 random yeast orthologs\n mapped to synthetic lethal interactions")
dev.off()


#Permutation 2 for scenario 2 Ovarian
#Fix mutated genes as GeneX
HM <- HowMany[,1]

mapToYeast <- function(pairFrame){
  firstPair <- yeastOrtho[yeastOrtho$Associated.Gene.Name %in% pairFrame$geneX & yeastOrtho$X..Identity.with.respect.to.query.gene > 40,]
  
  firstPairMerge <- merge(pairFrame, firstPair, by.x=1, by.y=3)
  colnames(firstPairMerge)[5] <- "YeastX"
  
  secondPair <- yeastOrtho[yeastOrtho$Associated.Gene.Name %in% pairFrame$geneY & yeastOrtho$X..Identity.with.respect.to.query.gene > 40,]
  secondPairMerge <- merge(firstPairMerge, secondPair, by.x =2, by.y  =3)
  colnames(secondPairMerge)[12] <- "YeastY"  
  return(secondPairMerge[,c("geneX","geneY", "YeastX", "YeastY")])
}

#permuteFromHumanGenes takes an index as function, 
#HM (HowMany as vector, and cc (correlation results)
permuteFromHumanGenes <- function(ind,HM,cc){
   ccMut <- cc[names(HM)]
   genes <- rownames(ccMut[[1]])

   #build new set of pairs by sampling from genes randomly with
   #same degree for each gene X
   geneY <- unlist(lapply(HM,function(x){sample(genes,size=x,replace=FALSE)}))
   geneX <- lapply(list(HM), function(x){rep(names(x),times=x)})

   #build pairFrame (input for mapToYeast()
   pairFrame <- data.frame(geneX=geneX[[1]], geneY)
   yeastPairs <- mapToYeast(pairFrame)
   yeastPairs <- yeastPairs[!duplicated(yeastPairs),]
   
   #look for synthetic lethal pairs
   attach(yeastPairs)
   map1 <- yeastPairs[YeastX %in% YeastSynthLeth$Systematic.Name.Interactor.A & YeastX %in% YeastSynthLeth$Systematic.Name.Interactor.B,]
   map2 <- yeastPairs[YeastY %in% YeastSynthLeth$Systematic.Name.Interactor.A & YeastY %in% YeastSynthLeth$Systematic.Name.Interactor.B,]
   detach(yeastPairs)
   outMap <- rbind(map1,map2)
   #print(outMap)
   
   #return number of unique pairs found
   return(nrow(outMap[!duplicated(outMap),]))
   
}

#run Permutation 2 results on Ovary Scenario 2
load("/gluster/home/ijang/shRNA_database/Project_Achilles/OVARY_scenario2.Rdata")
numPermutes <- 10000
permuteResultsHuman <- mclapply(1:numPermutes, function(x){permuteFromHumanGenes(x,HM,cc)}, mc.cores=5)
distResHuman <- unlist(permuteResultsHuman)

obsPairs <- 22
#return p-value
length(which(distResHuman >= obsPairs))/10000


png("Mutated-fixed-ovarian-permutations-scenario2.png")
hist(distResHuman, xlab="Number of Synthetic Pairs Picked", main ="Permutations with Mutated gene fixed, other gene random")
dev.off()


#run results on OVARY scenario 3
load("/gluster/home/ijang/shRNA_database/Project_Achilles/OVARY_scenario3.Rdata")

HMgain <- HowMany.gain[,1]

#take out entries in cc.gain that have no matrix
ccCandGain <- unlist(lapply(cc.gain, function(x){is.matrix(x)}))
ccFilteredGain <- cc.gain[ccCandGain]
#name gene X
names(ccFilteredGain )<- name.apair.gain

#Permutation 2 results for Scenario 3 ovary
permuteResultsHumanGain <- mclapply(1:numPermutes, function(x){permuteFromHumanGenes(x,HMgain,ccFilteredGain)}, mc.cores=5)
obs <- 80
distResHumanGain <- unlist(permuteResultsHumanGain)
#return p-value
length(which(distResHumanGain >= obs))/numPermutes

#Permutation 1 results for Scenario 1 ovaary
permuteResultsGain <- mclapply(1:numPermutes, function(x,...){permuteCandidates(x)}, mc.cores=5, numPairs=obs)
distYeastGain <- unlist(permuteResultsGain)
#return p-value
length(which(distYeastGain > obs))/numPermutes

#run Permutation 1 and 2 for scenario 2 LARGE intestine
load("/gluster/home/ijang/shRNA_database/Project_Achilles/LARGE_scenario2.Rdata")
HM <- HowMany[,1]
ccCand <- unlist(lapply(cc, function(x){is.matrix(x)}))
ccFiltered <- cc[ccCand]
names(ccFiltered )<- name.apair

#Permutation 2 results for scenario 2 LARGE intestine
numPermutes <- 10000
permuteResultsHumanGain <- mclapply(1:numPermutes, function(x){permuteFromHumanGenes(x,HM,ccFiltered)}, mc.cores=5)
obs <- 39
distResHumanGain <- unlist(permuteResultsHumanGain)
#return p-value
length(which(distResHumanGain >= obs))/numPermutes

#permutation 1 results for scenario 2 LARGE intestine
obs <- 39  
permuteResults <- mclapply(1:numPermutes, function(x,...){permuteCandidates(x)}, mc.cores=5, numPairs=obs)
distYeast <- unlist(permuteResults)
hist(distYeast)

#return p-value
length(which(distYeast >= obs))/numPermutes



