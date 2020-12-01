
library("getopt")
options(stringAsfactors = F, useFancyQuotes = F)

args <- commandArgs(trailingOnly = TRUE)

option_specification = matrix(c('input1', 'i1', 2, 'character',
								'input2', 'i2', 2, 'character',
								'output1', 'o1', 2, 'character'
                               ),
                              byrow=TRUE, ncol=4);

options = getopt(option_specification);


options(warn=-1)

library(sva, verbose = F)
library(pamr, verbose = F)
library(limma, verbose = F)

expMat <- as.matrix(read.table(file = options$input1,
                               sep = "\t",
                               header = T,
                               quote = "",
                               row.names = 1))
sampleMat <- read.table(file = options$input2,
                                  sep = "\t",
                                  header = F)


var <- apply(expMat, 1, var)
curExpMat <- expMat[-which(var == 0), ]

sampleMat <- as.data.frame(sampleMat)
colnames(sampleMat) <- c("sampleID", "batch")

tt <- ComBat(dat = as.matrix(curExpMat), batch = sampleMat$batch, mod = NULL, 
             par.prior = TRUE, prior.plots = FALSE,
             mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam"))

write.table(tt, file = options$output1, sep = "\t", quote = F, row.names = T, col.names = T)








