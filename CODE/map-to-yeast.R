load("/gluster/home/ijang/shRNA_database/Project_Achilles//OVARY_scenario2.Rdata")
options(stringsAsFactors=FALSE)

#build list of pairs
ccCand <- unlist(lapply(cc, function(x){is.matrix(x)}))
ccFiltered <- cc[ccCand]
names(ccFiltered )<- name.apair

#filter genes for interactions
filteredList <- lapply(names(ccFiltered), function(x){
  test <- data.frame(ccFiltered[[x]])
  attach(test)
  cands <- test[MT.cor.p.value < 0.1 & WT.cor.p.value < 0.1 & abs(MT.cor.estimate-WT.cor.estimate)> 1,]
  detach(test)
  geneX <- rep(x, nrow(cands))
  geneY <- rownames(cands)
  out <- data.frame(geneX, geneY,cands)
  out
})

filteredFrame <- do.call(rbind, filteredList)
write.table(filteredFrame,"scenario2-IJ-filtered-interactions.txt", sep="\t", quote=F)

yeastOrtho <- read.delim("yeast-human-orthologs.txt")
colnames(yeastOrtho)
hist(yeastOrtho$X..Identity.with.respect.to.query.gene)
hist(yeastOrtho$X..Identity.with.respect.to.Yeast.gene)

#try mapping geneX (name.apair) first
dim(yeastOrtho[yeastOrtho$Associated.Gene.Name %in% name.apair & 
                 yeastOrtho$X..Identity.with.respect.to.query.gene < 30,])

#grab only those inviable synthetic lethal interactions
YeastInt <- read.delim("BIOGRID-ORGANISM-Saccharomyces_cerevisiae-3.2.117.tab2.txt")
YeastSynthLeth <- YeastInt[YeastInt$Experimental.System == "Synthetic Lethality",]
YeastSynthLeth <- YeastSynthLeth[grep("inviable", YeastSynthLeth$Phenotypes),]
rm(YeastInt)

save.image("OV_scenario2_yeastMap.RData")

firstPair <- yeastOrtho[yeastOrtho$Associated.Gene.Name %in% filteredFrame$geneX & yeastOrtho$X..Identity.with.respect.to.query.gene > 40,]

firstPairMerge <- merge(filteredFrame, firstPair, by.x=1, by.y=3)
colnames(firstPairMerge)[9] <- "YeastX"

secondPair <- yeastOrtho[yeastOrtho$Associated.Gene.Name %in% filteredFrame$geneY & yeastOrtho$X..Identity.with.respect.to.query.gene > 40,]
secondPairMerge <- merge(firstPairMerge, secondPair, by.x =2, by.y  =3)
colnames(secondPairMerge)[16] <- "YeastY"

attach(secondPairMerge)
map1 <- secondPairMerge[YeastX %in% YeastSynthLeth$Systematic.Name.Interactor.A & YeastX %in% YeastSynthLeth$Systematic.Name.Interactor.B,]
map2 <- secondPairMerge[YeastY %in% YeastSynthLeth$Systematic.Name.Interactor.A & YeastY %in% YeastSynthLeth$Systematic.Name.Interactor.B,]

outMap <- rbind(map1, map2)

write.table(outMap,"mapped-scenario1-yeast-interaction-candidates.txt", quote=F, sep="\t")

