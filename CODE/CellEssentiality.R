#################################################################################################################################
##### Cell Essentiality #####
# Clear Screen
cat("\014")  

# Clear workspace
rm(list=ls())

# Set working and library paths
setwd('/gluster/home/tperumal/Work/Cell Essentiality/')
.libPaths('/gluster/home/tperumal/mylibs/')

# Load external libraries
library('CePa')
library('MRCE')
library('Biobase')
library('pracma')
library('synapseClient')
synapseLogin()

# CellEssentiality_MRCE <- function(
  CANCER_CENSUS = TRUE;
  METHOD = 'single'; #options are single; cv
  
  # options if method is 'single'
  L1 = 1; 
  L2 = 1; 
  
  # options if method is 'cv'
  LVEC1 = c(1,0.8,0.5); 
  LVEC2 = c(1,0.8,0.5);
  
  RESULTS_DIR = '.';
  FILE_NAME = 'MODEL.RData';
#   ){   
  
  #################################################################################################################################
  ##### Load Data #####
  # Load Colt Data #
  COLT_SCORE <- read.table('Colt/processed_GARP_Score.txt',sep=' ',header=T)
  COLT_PVAL <- read.table('Colt/processed_GARP_Pval.txt',sep=' ',header=T)
  
  # Load Achilles Data #
  ACHILLES_SCORE <- read.gct("Project_Achilles/Achilles_102lines_gene_solutions.gct")
  
  # Load CCLE Data #
  id_exprLayer <- "syn1757082" 
  layer_expr <- loadEntity(id_exprLayer)
  CCLE_EXPR <- exprs(layer_expr$objects$eSet_expr)
  colnames(CCLE_EXPR) <- sapply(colnames(CCLE_EXPR),function(x){return(strsplit(x,'_')[[1]][1])})
  
  # Load Sanger Data #
  id_exprLayer <- "1742878" 
  layer_expr <- loadEntity(id_exprLayer)
  SANGER_EXPR <- t(as.data.frame(layer_expr$objects$eSet_expr))
  
  #  Retrive and Save Cancer Census GeneList #
  Cancer_census <- read.csv('Census_allFri Sep  5 05-50-14 2014.csv')
  CENSUS_GENES <- unique(Cancer_census$Gene.Symbol)
  #################################################################################################################################
  
  
  #################################################################################################################################
  ##### Sample pre-processing #####
  
  #### Processing Sample Names ####
  # COLT
  colnames(COLT_SCORE) <- gsub('\\.|\\_','',toupper(colnames(COLT_SCORE)))
  colnames(COLT_PVAL) <- gsub('\\.|\\_','',toupper(colnames(COLT_PVAL)))
  
  # Achilles
  colnames(ACHILLES_SCORE) <- unlist(sapply(colnames(ACHILLES_SCORE),function(x){strsplit(x,'_')[[1]][1]}))
  
  # CCLE
  colnames(CCLE_EXPR) <- unlist(sapply(colnames(CCLE_EXPR),function(x){strsplit(x,'_')[[1]][1]}))
  
  # Sanger
  colnames(SANGER_EXPR) <- gsub('\\.|\\_','',toupper(colnames(SANGER_EXPR)))
  
  #### Processing Gene Names ####
  # COLT
  rownames(COLT_PVAL) <- COLT_PVAL$TARGETGENEREF
  rownames(COLT_SCORE) <- COLT_SCORE$REFSEQ
  
  COLT_SCORE <- COLT_SCORE[!duplicated(COLT_SCORE$GENENAME),]
  COLT_PVAL <- COLT_PVAL[rownames(COLT_SCORE),]
  
  COLT_SCORE <- COLT_SCORE[,-which(colnames(COLT_SCORE) %in% c('REFSEQ','GENEID','DESCRIPTION'))]
  COLT_PVAL <- COLT_PVAL[,-which(colnames(COLT_PVAL) %in% c('TARGETGENEREF','NCBIGENEID',
                                                                           'GENEDESCRIPTION','SPECIES','TARGETREGION'))]
  
  rownames(COLT_PVAL) <- COLT_PVAL$GENENAME 
  rownames(COLT_SCORE) <- COLT_SCORE$GENENAME
  
  COLT_PVAL <- COLT_PVAL[,-which(colnames(COLT_PVAL) %in% c('GENENAME'))]
  COLT_SCORE <- COLT_SCORE[,-which(colnames(COLT_SCORE) %in% c('GENENAME'))]
  
  COLT_PVAL <- data.matrix(COLT_PVAL)
  COLT_SCORE <- data.matrix(COLT_SCORE)
  
  # Replace pval > 0.05 with 0 and remove rows with zero sum
  COLT_SCORE[which(COLT_PVAL > 0.05,arr.ind=T)] <- 0
  COLT_PVAL <- COLT_PVAL[rowSums(COLT_SCORE)!=0,]
  COLT_SCORE <- COLT_SCORE[rowSums(COLT_SCORE)!=0,]
  
  # Achilles
  rownames(ACHILLES_SCORE) <- sapply(rownames(ACHILLES_SCORE),function(x){return(strsplit(x,'_')[[1]][1])})
  #################################################################################################################################
  
  
  #################################################################################################################################
  #### Match Colt and Sanger samples ####
  SANGER_EXPR <- SANGER_EXPR[,intersect(colnames(COLT_SCORE),colnames(SANGER_EXPR))]
  COLT_SCORE <- COLT_SCORE[,intersect(colnames(COLT_SCORE),colnames(SANGER_EXPR))]
  
  #### Match Achilles and CCLE samples ####
  CCLE_EXPR <- CCLE_EXPR[,intersect(colnames(ACHILLES_SCORE),colnames(CCLE_EXPR))]
  ACHILLES_SCORE <- ACHILLES_SCORE[,intersect(colnames(ACHILLES_SCORE),colnames(CCLE_EXPR))]
  #################################################################################################################################
  
  
  
  #################################################################################################################################
  #### Remove zero perturbation samples and genes ####
  CCLE_EXPR <- CCLE_EXPR[rowSums(CCLE_EXPR) != 0,]
  SANGER_EXPR <- SANGER_EXPR[rowSums(SANGER_EXPR) != 0,]
  COLT_SCORE <- COLT_SCORE[rowSums(COLT_SCORE) != 0,]
  ACHILLES_SCORE <- ACHILLES_SCORE[rowSums(ACHILLES_SCORE) != 0,]
  #################################################################################################################################
  
  
  
  #################################################################################################################################
  #### Retain Common Perturbation between Colt and Achilles ####
  Perturbations <- intersect(rownames(COLT_SCORE),rownames(ACHILLES_SCORE))
  COLT_SCORE <- COLT_SCORE[Perturbations,]
  ACHILLES_SCORE <- ACHILLES_SCORE[Perturbations,]
  
  #### Retain Common Observations between CCLE and Sanger ####
  Observations <- intersect(rownames(CCLE_EXPR),rownames(SANGER_EXPR))
  CCLE_EXPR <- CCLE_EXPR[Observations,]
  SANGER_EXPR <- SANGER_EXPR[Observations,]
  #################################################################################################################################
  
  
  
  #################################################################################################################################
  #### Arrange and save Colt, Achilles, CCLE and Sanger Data ####
  COLT_SCORE <- COLT_SCORE[c(intersect(Perturbations,Observations),
                                       setdiff(Perturbations,Observations)),]
  ACHILLES_SCORE <- ACHILLES_SCORE[c(intersect(Perturbations,Observations),
                                     setdiff(Perturbations,Observations)),]
  CCLE_EXPR <- CCLE_EXPR[c(intersect(Perturbations,Observations),
                           setdiff(Observations,Perturbations)),]
  SANGER_EXPR <- SANGER_EXPR[c(intersect(Perturbations,Observations),
                               setdiff(Observations,Perturbations)),]
  #################################################################################################################################
  
  
  
  #################################################################################################################################
  #### Pick Cancer Census Gene Sets for Prediction ####
  if (CANCER_CENSUS == 1){
    ACHILLES_SCORE <- ACHILLES_SCORE[intersect(CENSUS_GENES,rownames(ACHILLES_SCORE)),]
    COLT_SCORE <- COLT_SCORE[intersect(CENSUS_GENES,rownames(COLT_SCORE)),]    
  }  
  #################################################################################################################################
  
  
  #################################################################################################################################
  #### Save Data ####
  write.table(COLT_SCORE,file="COLT_SCORE",sep="\t",quote=F)
  write.table(ACHILLES_SCORE,file="Achilles_Score",sep="\t",quote=F)
  write.table(CCLE_EXPR,file="CCLE_EXPR",sep="\t",quote=F)
  write.table(SANGER_EXPR,file="SANGER_EXPR",sep="\t",quote=F)
  #################################################################################################################################
  