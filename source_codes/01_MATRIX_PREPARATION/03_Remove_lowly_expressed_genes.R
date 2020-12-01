


library("getopt")

options(stringAsfactors = F, useFancyQuotes = F)

args <- commandArgs(trailingOnly = TRUE)
option_specification <- matrix(c('input1', 'i1', 2, 'character',
                                 'input2', 'i2', 2, 'character',
                                 'input3', 'i3', 2, 'character',
                                 'output1', 'o1', 2, 'character'),
                               byrow=TRUE, ncol=4)

options = getopt(option_specification);


minExpValue <- options$input2
minExpSamNum <- options$input3

extractExpGene <- function(x){
  
  expSamNum <- length(which(as.numeric(x) >= minExpValue))
  return(expSamNum)
  
}

geneExp <- as.matrix(read.table(file = options$input1,
                                 sep = "\t",
                                 quote = "",
                                 comment.char = "",
                                 header = T,
                                 row.names = 1))

expGeneID <- rownames(geneExp[which(apply(geneExp, 1, extractExpGene) >= minExpSamNum), ])
expGeneExp <- geneExp[expGeneID, ]


write.table(expGeneExp,
            file = options$output1,
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = T)
