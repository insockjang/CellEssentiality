##prepare data for analysis in different scenarios
##code is from In Sock's test02.R
##load Achilles data

library(CePa)
A<-read.gct("/home/ijang/shRNA_database//Project_Achilles/Data//Achilles_102lines_gene_solutions.gct")
B<-read.gct("/home/ijang/shRNA_database//Project_Achilles/Data//20110303_achilles2_PMAD_adjFC.rnai.gct")
C<-read.delim("/home/ijang/shRNA_database/Project_Achilles//Data/Achilles_102lines_shRNA_table.txt")
D<-read.delim("/home/ijang/shRNA_database/Project_Achilles//Data/Achilles_v2.0_SampleInfo_small.txt.original.txt")

# gene symbol in RNAi screen
name<-sapply(strsplit(rownames(A),"_"),function(x){x[[1]]})
name.A<-strsplit(rownames(A),"_")
name<-c()
for(k in 1:length(name.A)){
  name<-c(name,name.A[[k]][1])
}

name.duplicate = name[which(duplicated(name)==1)]

which(name == name.duplicate[2])

# cell line sample name
name1<-colnames(A)

# expression and copy level in CCLE
library(predictiveModeling)
library(synapseClient)
synapseLogin()
id_exprLayer <- "syn1757082" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- (layer_expr$objects$eSet_expr)

id_copyLayer <- "syn1757086"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- (layer_copy$objects$eSet_copy)

# # Gistic score in CCLE : Corrupted
# id_copyLayer1 <- "syn2367236"     
# layer_copy1 <- loadEntity(id_copyLayer1)
# eSet_cgh <- (layer_copy1$objects$eSet_copy)

id_hybridLayer <-  "syn1757084" 
layer_hybrid <- loadEntity(id_hybridLayer)
eSet_hybrid <- layer_hybrid$objects$eSet_hybrid

featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr,copy = eSet_copy, mut = eSet_hybrid))    

common.sample<-intersect(colnames(featureData),colnames(A))

##### final input and response data matching with sample names in CCLE
data.input<-featureData[,common.sample]
data.response<-A[,common.sample]

data.input.exp<-exprs(eSet_expr)[,common.sample]
data.input.copy<-exprs(eSet_copy)[,common.sample]
data.input.mut<-exprs(eSet_hybrid)[,common.sample]

tissue.sample<-sapply(strsplit(common.sample,"_"),function(x){x[[2]]})

table(tissue.sample)

save.image("CCLE-Achilles.RData")


#get COlT Data
coltData <- synGet("syn2582519")
coltDataGarp <- read.table(coltData@filePath, sep=" ", header=TRUE)
coltAnnot <- coltDataGarp[,1:4]
coltDataGarp <- coltDataGarp[,-c(1:4)]
coltCells <- make.names(toupper(colnames(coltDataGarp)))
coltCells <- gsub(".","",coltCells,fixed=TRUE)
colnames(coltDataGarp) <- coltCells


##get Sanger copy Data
id_copyLayer <- "1742880"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <-  layer_copy@objOwn$objects$eSet_copy

#mutation data
id_oncomapLayer <- "syn1742882"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$eSet_oncomap

id_hybridLayer <- "syn2368559"  
layer_hybrid <- synGet(id_hybridLayer,load = T)
load(layer_hybrid@filePath)

#expression data
id_exprLayer <- "1742878" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$eSet_expr

featureSanger <- createAggregateFeatureDataSet(list(expr = eSet_expr,copy = eSet_copy, mut=eSet_oncomap))

coltCellLines <- read.delim("coltCellLines.txt")
coltCellLines$cellLine <- gsub("-", "", coltCellLines$cellLine, fixed=TRUE)
coltCellLines$cellLine <- gsub(" ", "", coltCellLines$cellLine, fixed=TRUE)
coltCellLines$cellLine <- gsub(".","",coltCellLines$cellLine, fixed=TRUE)
coltCellLines$cellLine <- toupper(coltCellLines$cellLine)

inBoth <- intersect(colnames(coltDataGarp), colnames(featureSanger))

grep("OV",colnames(coltDataGarp))