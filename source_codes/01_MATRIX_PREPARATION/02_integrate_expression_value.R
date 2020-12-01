

############################# integrate stringtie results
library("getopt")
options(stringAsfactors = F, useFancyQuotes = F)

args <- commandArgs(trailingOnly = TRUE)
option_specification <- matrix(c('input1', 'i1', 2, 'character',
                                 'output1', 'o1', 2, 'character'),
                               byrow=TRUE, ncol=4)

options = getopt(option_specification);

######### For Gene #########
############################
resDic <- paste0(options$input1, "/")
resDir <- list.files(resDic)

#### check gene ID
geneID <- c()
for( i in 1:length(resDir)){
  
  curDir <- paste0(resDic, resDir[i])
  curRes <- as.matrix(read.table(file = curDir, sep = "\t", quote = "", header = T))
  
  geneID <- c(geneID, curRes[ , 1])
  
  cat(i, "\n")
}
geneID <- unique(geneID)

#### extract gene expression 
geneExp <- matrix(0, nrow = length(geneID), ncol = length(resDir))
rownames(geneExp) <- geneID
colnames(geneExp) <- resDir

for( i in 1:length(resDir)){
  
  curDir <- paste0(resDic, resDir[i])
  curRes <- as.matrix(read.table(file = curDir, sep = "\t", quote = "", header = T))
  rownames(curRes) <- curRes[ , 1]
  
  geneExp[intersect(rownames(geneExp), curRes[ , 1]), i] <- curRes[intersect(rownames(geneExp), curRes[ , 1]), "TPM"] 
  
  cat(i, "\n")
}

colnames(geneExp) <- sub("_gene.gtf", "", colnames(geneExp))
rownames(geneExp) <- sub("gene:", "", rownames(geneExp))

write.table(geneExp, file = options$output1, sep = "\t", quote = F, col.names = T, row.names = T)
