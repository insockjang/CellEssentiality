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

source('./CellEssentiality_MRCE.R')
##########################################################################################################
#### Parameter Scanning for MRCE ####
LAMBDA = c(1000,100,10,1) #seq(10,0,-2)
ERRORS <- array(0,dim=c(length(LAMBDA),length(LAMBDA),4))
for (L1 in LAMBDA)
  for (L2 in LAMBDA){
    MODEL <- CellEssentiality_MRCE(METHOD = 'single', L1 = L1, 
                                   L2 = L2, FILE_NAME = paste('MODEL','census','single',L1,L2,'.RData',sep='_'))
    ERRORS[which(L1==LAMBDA), which(L2==LAMBDA),1:4] <- c(norm(MODEL$PERROR_ACHILLES_CCLE,'F'),
                                                          norm(MODEL$VERROR_ACHILLES_CCLE,'F'),
                                                          norm(MODEL$PERROR_COLT_SANGER,'F'),
                                                          norm(MODEL$VERROR_COLT_SANGER,'F'))     
  }