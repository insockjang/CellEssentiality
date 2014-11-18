load("/gluster/home/ijang/shRNA_database/Project_Achilles/PANCREAS_scenario3.Rdata")

outFileGain <- "mapped-scenario3-PAN-gain-candidates.txt"
outFileLoss <- "mapped-scenario3-PAN-loss-candidates.txt"
options(stringsAsFactors=FALSE)

#find entries in cc list that have candidates
ccCandGain <- unlist(lapply(cc.gain, function(x){is.matrix(x)}))
ccFilteredGain <- cc.gain[ccCandGain]
names(ccFilteredGain )<- name.apair.gain

ccCandLoss <- unlist(lapply(cc.loss, function(x){is.matrix(x)}))
ccFilteredLoss <- cc.loss[ccCandLoss]
names(ccFilteredLoss)<- name.apair.loss



#filter genes for interactions
namfilteredListGain <- lapply(names(ccFilteredGain), function(x){
  test <- data.frame(ccFilteredGain[[x]])
  attach(test)
  cands <- test[Gain.cor.p.value < 0.1 & WT.cor.p.value < 0.1 & abs(Gain.cor.estimate-WT.cor.estimate)> 1,]
  detach(test)
  geneX <- rep(x, nrow(cands))
  geneY <- rownames(cands)
  out <- data.frame(geneX, geneY,cands)
  out
})

filteredListLoss <- lapply(names(ccFilteredLoss), function(x){
  test <- data.frame(ccFilteredLoss[[x]])
  attach(test)
  cands <- test[Loss.cor.p.value < 0.1 & WT.cor.p.value < 0.1 & abs(Loss.cor.estimate-WT.cor.estimate)> 1,]
  detach(test)
  geneX <- rep(x, nrow(cands))
  geneY <- rownames(cands)
  out <- data.frame(geneX, geneY,cands)
  out
})



filteredFrameGain <- do.call(rbind, filteredListGain)
#write.table(filteredFrameGain,"scenario3-OV-gain-filtered-interactions.txt", sep="\t", quote=F)

filteredFrameLoss <- do.call(rbind, filteredListLoss)
#write.table(filteredFrameGain,"scenario3-OV-loss-filtered-interactions.txt", sep="\t", quote=F)




yeastOrtho <- read.delim("yeast-human-orthologs.txt")
colnames(yeastOrtho)
hist(yeastOrtho$X..Identity.with.respect.to.query.gene)
hist(yeastOrtho$X..Identity.with.respect.to.Yeast.gene)

#try mapping geneX (name.apair) first
dim(yeastOrtho[yeastOrtho$Associated.Gene.Name %in% name.apair.gain & 
                 yeastOrtho$X..Identity.with.respect.to.query.gene < 30,])

#grab only those inviable synthetic lethal interactions
YeastInt <- read.delim("BIOGRID-ORGANISM-Saccharomyces_cerevisiae-3.2.117.tab2.txt")
YeastSynthLeth <- YeastInt[YeastInt$Experimental.System == "Synthetic Lethality",]
YeastSynthLeth <- YeastSynthLeth[grep("inviable", YeastSynthLeth$Phenotypes),]
rm(YeastInt)

#save.image("OV_scenario3_yeastMap.RData")

firstPairGain <- yeastOrtho[yeastOrtho$Associated.Gene.Name %in% filteredFrameGain$geneX & yeastOrtho$X..Identity.with.respect.to.query.gene > 40,]

firstPairMergeGain <- merge(filteredFrameGain, firstPairGain, by.x=1, by.y=3)
colnames(firstPairMergeGain)[9] <- "YeastX"

secondPairGain <- yeastOrtho[yeastOrtho$Associated.Gene.Name %in% filteredFrameGain$geneY & yeastOrtho$X..Identity.with.respect.to.query.gene > 40,]
secondPairMergeGain <- merge(firstPairMergeGain, secondPairGain, by.x =2, by.y  =3)
colnames(secondPairMergeGain)[16] <- "YeastY"

attach(secondPairMergeGain)
map1 <- secondPairMergeGain[YeastX %in% YeastSynthLeth$Systematic.Name.Interactor.A & YeastX %in% YeastSynthLeth$Systematic.Name.Interactor.B,]
map2 <- secondPairMergeGain[YeastY %in% YeastSynthLeth$Systematic.Name.Interactor.A & YeastY %in% YeastSynthLeth$Systematic.Name.Interactor.B,]
detach(secondPairMergeGain)
outMapGain <- rbind(map1, map2)


firstPairLoss <- yeastOrtho[yeastOrtho$Associated.Gene.Name %in% filteredFrameLoss$geneX & yeastOrtho$X..Identity.with.respect.to.query.gene > 40,]

firstPairMergeLoss <- merge(filteredFrameLoss, firstPairLoss, by.x=1, by.y=3)
colnames(firstPairMergeLoss)[9] <- "YeastX"

secondPairLoss <- yeastOrtho[yeastOrtho$Associated.Gene.Name %in% filteredFrameLoss$geneY & yeastOrtho$X..Identity.with.respect.to.query.gene > 40,]
secondPairMergeLoss <- merge(firstPairMergeLoss, secondPairLoss, by.x =2, by.y  =3)
colnames(secondPairMergeLoss)[16] <- "YeastY"

attach(secondPairMergeLoss)
map1 <- secondPairMergeLoss[YeastX %in% YeastSynthLeth$Systematic.Name.Interactor.A & YeastX %in% YeastSynthLeth$Systematic.Name.Interactor.B,]
map2 <- secondPairMergeLoss[YeastY %in% YeastSynthLeth$Systematic.Name.Interactor.A & YeastY %in% YeastSynthLeth$Systematic.Name.Interactor.B,]
detach(secondPairMergeLoss)
outMapLoss <- rbind(map1, map2)



#write.table(outMap,"mapped-scenario3-yeast-interaction-candidates.txt", quote=F, sep="\t")

#map only unique interactions
attach(outMapGain)
simpleMapGain <- data.frame(geneX, geneY, YeastX, YeastY)
detach(outMapGain)
gainScenario3OV <- outMapGain[!duplicated(simpleMapGain),]
write.table(gainScenario3OV, outFileGain, quote=F, sep="\t")

attach(outMapLoss)
simpleMapLoss <- data.frame(geneX, geneY, YeastX, YeastY)
detach(outMapLoss)
lossScenario3OV <- outMapLoss[!duplicated(simpleMapLoss),]
write.table(lossScenario3OV, outFileLoss, quote=F, sep="\t")


