
library("getopt")

options(stringAsfactors = F, useFancyQuotes = F)

args <- commandArgs(trailingOnly = TRUE)
option_specification <- matrix(c('input1', 'i1', 2, 'character',
                                 'output1', 'o1', 2, 'character'),
                               byrow=TRUE, ncol=4)

options = getopt(option_specification);

options(warn=-1)
geneExp <- as.matrix(read.table(file = options$input1,
                                 sep = "\t",
                                 quote = "",
                                 comment.char = "",
                                 header = T,
                                 row.names = 1))

class(geneExp) <- "numeric"
geneExp <- log2(geneExp + 1)

write.table(geneExp,
            file = options$output1,
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = T)
