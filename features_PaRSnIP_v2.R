library(data.table)
library(doParallel)
library(foreach)

args = commandArgs(trailingOnly=TRUE)

#setwd('/home/local/QCRI/rmall/Proteins/Protein_Crystallization/scripts/')
processes <- 4
registerDoParallel(processes)
source('PaRSnIP_v2.R')

#Get list of protein sequences
alns <- read.fasta(args[1])
N <- length(alns$id)

#Scratch path
SCRATCH.path <- './SCRATCH-1D_1.2/bin/run_SCRATCH-1D_predictors.sh'
DISORDER.path <- './DISOPRED/run_disopred.pl'

#Create feature matrix along with true labels
full_feature_matrix <- foreach (j = 1:N, .inorder = TRUE, .combine = 'rbind') %dopar%
{
  aln.ali <- alns$ali[j,]
  output_prefix <- tempfile( pattern = paste0("tmp_",j,"_"),tmpdir="/tmp",
                             fileext = "" )
  temp_features <- PaRSnIP.calc.features.RSAhydro.test(aln.ali,SCRATCH.path,DISORDER.path,
                                                       output_prefix,n.cores=1);
  temp_features
}

full_feature_matrix <- as.matrix(full_feature_matrix);

if (dim(full_feature_matrix)[2] == 1)
{
  full_feature_matrix <- t(full_feature_matrix)
}

write.csv(full_feature_matrix,"features.csv",row.names = T)
