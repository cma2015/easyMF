######################spatial course analysis#########################################
# Dynamic Transcriptome Landscape of Maize Embryo and Endosperm Development
# Setup R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})
# We need to not crash galaxy with an UTF8 error on German LC settings.
#loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
# Import required libraries
library("getopt")
options(stringAsfactors = F, useFancyQuotes = F)
# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Get options using the spec as defined by the enclosed list
# Read the options from the default: commandArgs(TRUE)

# input1: input the prepared gene expression matrix
# input2: input the amplitude matrix
# input3: input the pattern matrix
# input4: input the developmental information of the samples
# input5: input the the gene description file
# input6: input the cpu numbers to parallelly compute
# output1: significant genes of each module
# pdf1: 
option_specification = matrix(c('input1', 'i1', 2, 'character',
                                'input2', 'i2', 2, 'character',
                                'input3', 'i3', 2, 'character',
                                'input4', 'i4', 2, 'character',
                                'output1', 'o1', 2, 'character',
                                'pdf1', 'o2', 2, 'character'
),
byrow=TRUE, ncol=4);
# Parse options
options(warn=-1)
options = getopt(option_specification);

patternMarkers <- function (Amatrix = NA, scaledPmatrix = FALSE, Pmatrix = NA, 
                            threshold = "All", lp = NA, full = FALSE){
  if (scaledPmatrix == FALSE) {
    if (!is.na(Pmatrix)) {
      pscale <- apply(Pmatrix, 1, max)
      Amatrix <- sweep(Amatrix, 2, pscale, FUN = "*")
    }
    else (warning("P values must be provided if not already scaled"))
  }
  Arowmax <- t(apply(Amatrix, 1, function(x) x/max(x)))
  if(length(which(is.na(Arowmax[,1]) == T))){
    Arowmax[which(is.na(Arowmax[,1]) == T), ] <- 0
  }
  pmax <- apply(Amatrix, 1, max)
  ssranks <- matrix(NA, nrow = nrow(Amatrix), ncol = ncol(Amatrix), 
                    dimnames = dimnames(Amatrix))
  ssgenes <- matrix(NA, nrow = nrow(Amatrix), ncol = ncol(Amatrix), 
                    dimnames = NULL)
  nP = dim(Amatrix)[2]
  if (!is.na(lp)) {
    if (length(lp) != dim(Amatrix)[2]) {
      warning("lp length must equal the number of columns of the Amatrix")
    }
    sstat <- apply(Arowmax, 1, function(x) sqrt(t(x - lp) %*% 
                                                  (x - lp)))
    ssranks[order(sstat), i] <- 1:length(sstat)
    ssgenes[, i] <- names(sort(sstat, decreasing = FALSE))
  }
  else {
    for (i in 1:nP) {
      lp <- rep(0, dim(Amatrix)[2])
      lp[i] <- 1
      sstat <- apply(Arowmax, 1, function(x) sqrt(t(x - 
                                                      lp) %*% (x - lp)))
      ssranks[order(sstat), i] <- 1:length(sstat)
      ssgenes[, i] <- names(sort(sstat, decreasing = FALSE))
    }
  }
  if (threshold == "cut") {
    geneThresh <- sapply(1:nP, function(x) min(which(ssranks[ssgenes[, 
                                                                     x], x] > apply(ssranks[ssgenes[, x], ], 1, min))))
    ssgenes.th <- sapply(1:nP, function(x) ssgenes[1:geneThresh[x], 
                                                   x])
  }
  if (threshold == "All") {
    pIndx <- apply(ssranks, 1, which.min)
    gBYp <- lapply(sort(unique(pIndx)), function(x) names(pIndx[pIndx == 
                                                                  x]))
    ssgenes.th <- lapply(1:nP, function(x) ssgenes[which(ssgenes[, 
                                                                 x] %in% gBYp[[x]]), x])
  }
  if (full) {
    return(list(PatternMarkers = ssgenes.th, PatternRanks = ssranks))
  }
  else {
    return(PatternMarkers = ssgenes.th)
  }
}


.pasteVector <- function(input){
  output <- c()
  for(i in 1:length(input)){
    output <- paste0(output, ";", input[i])
  }
  output <- substr(x = output, 2, nchar(output))
  output
}

.intersectGene <- function(A, B){
  intersect(A, B)
}

prePareGO <- function(dataset = "zmays_eg_gene"){
  # if(!require(biomaRt)){
  #   source("https://bioconductor.org/biocLite.R")
  #   biocLite("biomaRt")
  # }
  # 
  # if(!require(topGO)){
  #   source("https://bioconductor.org/biocLite.R")
  #   biocLite("topGO")
  # }
  require(biomaRt)
  require(topGO)
  mart <- useMart(biomart = "plants_mart", dataset = dataset, host = 'plants.ensembl.org')
  GTOGO <- getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
  GTOGO <- GTOGO[GTOGO$go_id != '', ]
  geneID2GO <- by(GTOGO$go_id, GTOGO$ensembl_gene_id, function(x) as.character(x))
  return(list(GTOGO = GTOGO, geneID2GO = geneID2GO))
}

runTopGO <- function(geneID, prePareGO,  statistic = "fisher", algorithm = "elim",
                     topNodes = 200, plot = TRUE){
  library(ggplot2)
  # if(!require(biomaRt)){
  #   source("https://bioconductor.org/biocLite.R")
  #   biocLite("biomaRt")
  # }
  # 
  # if(!require(topGO)){
  #   source("https://bioconductor.org/biocLite.R")
  #   biocLite("topGO")
  # }
  # mart <- useMart(biomart = "plants_mart", dataset = dataset, host = 'plants.ensembl.org')
  # GTOGO <- getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
  # GTOGO <- GTOGO[GTOGO$go_id != '', ]
  # geneID2GO <- by(GTOGO$go_id, GTOGO$ensembl_gene_id, function(x) as.character(x))
  require(biomaRt)
  require(topGO)
  GTOGO <- prePareGO$GTOGO
  geneID2GO <- prePareGO$geneID2GO
  all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
  int.genes <- geneID
  int.genes <- intersect(int.genes, names(geneID2GO))
  int.genes <- factor(as.integer(all.genes %in% int.genes))
  names(int.genes) = all.genes
  
  go.obj.BP <- new("topGOdata", ontology='BP',
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   gene2GO = geneID2GO)
  
  go.obj.MF <- new("topGOdata", ontology='MF',
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   gene2GO = geneID2GO)
  
  go.obj.CC <- new("topGOdata", ontology='CC',
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   gene2GO = geneID2GO)
  
  ##########retrieve the gene list related to a GO ID######################
  allGO.BP <- genesInTerm(object = go.obj.BP)
  allGO.MF <- genesInTerm(object = go.obj.MF)
  allGO.CC <- genesInTerm(object = go.obj.CC)
  
  #########retrive the significant GO terms
  results.BP <- runTest(go.obj.BP, algorithm = algorithm, statistic = statistic)
  results.tab.BP <- GenTable(object = go.obj.BP, elimFisher = results.BP,
                             topNodes = topNodes)
  gene.BP <- genesInTerm(object = go.obj.BP, whichGO = results.tab.BP$GO.ID)
  inter.gene.BP <- lapply(X = gene.BP, FUN = .intersectGene, B = geneID)
  inter.gene.BP <- unlist(lapply(X = inter.gene.BP, FUN = .pasteVector))
  results.tab.BP$significantGene <- inter.gene.BP
  
  if(length(which(results.tab.BP$elimFisher == "< 1e-30")) != 0){
    results.tab.BP[which(results.tab.BP$elimFisher == "< 1e-30"), ]$elimFisher <- 1e-30
  }
  
  
  results.MF <- runTest(go.obj.MF, algorithm = algorithm, statistic = statistic)
  results.tab.MF <- GenTable(object = go.obj.MF, elimFisher = results.MF, 
                             topNodes = topNodes)
  gene.MF <- genesInTerm(object = go.obj.MF, whichGO = results.tab.MF$GO.ID)
  inter.gene.MF <- lapply(X = gene.MF, FUN = .intersectGene, B = geneID)
  inter.gene.MF <- unlist(lapply(X = inter.gene.MF, FUN = .pasteVector))
  results.tab.MF$significantGene <- inter.gene.MF
  if(length(which(results.tab.MF$elimFisher == "< 1e-30")) != 0){
    results.tab.MF[which(results.tab.MF$elimFisher == "< 1e-30"), ]$elimFisher <- 1e-30
  }
  
  results.CC <- runTest(go.obj.CC, algorithm = algorithm, statistic = statistic)
  results.tab.CC <- GenTable(object = go.obj.CC, elimFisher = results.CC, 
                             topNodes = topNodes)
  gene.CC <- genesInTerm(object = go.obj.CC, whichGO = results.tab.CC$GO.ID)
  inter.gene.CC <- lapply(X = gene.CC, FUN = .intersectGene, B = geneID)
  inter.gene.CC <- unlist(lapply(X = inter.gene.CC, FUN = .pasteVector))
  results.tab.CC$significantGene <- inter.gene.CC
  if(length(which(results.tab.CC$elimFisher == "< 1e-30")) != 0){
    results.tab.CC[which(results.tab.CC$elimFisher == "< 1e-30"), ]$elimFisher <- 1e-30
  }
  
  
  if(plot){
    df <- data.frame(Category = c(rep("BP", topNodes), rep("CC", topNodes), rep("MF", topNodes)), 
                     x = c(results.tab.BP$Significant, results.tab.CC$Significant, 
                           results.tab.MF$Significant),
                     y = c(-log10(as.numeric(results.tab.BP$elimFisher)), 
                           -log10(as.numeric(results.tab.CC$elimFisher)), 
                           -log10(as.numeric(results.tab.MF$elimFisher))),
                     size = c(-log10(as.numeric(results.tab.BP$elimFisher)),
                              -log10(as.numeric(results.tab.CC$elimFisher)), 
                              -log10(as.numeric(results.tab.MF$elimFisher)))
    )
    
    kk <- ggplot(data = df, aes(x = x, y = y)) + 
      geom_point(aes(color = Category, size = size)) + 
      scale_size_continuous(range = c(2,10)) + 
      labs(x = "The number of significant genes", y = "The adjusted p-values for each GO term")
    print(kk)
  }
  
  results <- list(BP = results.tab.BP, CC = results.tab.CC, MF = results.tab.MF)
  results
}


SpatialCourseAnalysis <- function(){
  
  
  gene_expression <- as.matrix(read.table(options$input1, header = T, stringsAsFactors = F))
  Amplitude <- as.matrix(read.table(options$input2, header = T, sep = "\t", stringsAsFactors = F))
  Pattern <- as.matrix(read.table(options$input3, header = T, sep = "\t", stringsAsFactors = F))
  plants_ensembl <- as.matrix(read.table("/galaxy/tools/TAMF/test_data/plants_ensemble.txt", header = F, stringsAsFactors = F, sep = "\t"))          
  dataSet <- plants_ensembl[as.numeric(options$input4),1]
  
  
  library(CoGAPS)
  patternmarker <- patternMarkers(Amatrix = Amplitude, Pmatrix = t(Pattern), scaledPmatrix = T, threshold = "cut")
  tempsc <- t(apply(gene_expression[unlist(patternmarker),], 1, scale))
  colnames(tempsc) <- colnames(gene_expression)
  
  library(pheatmap)
  gene_expression_scale <- t(apply(gene_expression[unlist(patternmarker),], 1, scale))
  ph=pheatmap(gene_expression[unlist(patternmarker),], scale = "row", show_rownames = F, show_colnames = F, cluster_rows = T,cluster_cols = T, cutree_rows = 10, cutree_cols = 10, silent = T)
  row_cluster=cutree(ph$tree_row,k=10)
  col_cluster=cutree(ph$tree_col,k=10)
  a1 <- sapply(names(row_cluster), function(ii){return(grep("[.]",ii))})
  for(j in which(a1>0)){
    names(row_cluster)[j] <- unlist(strsplit(names(row_cluster)[j], split = "[.]"))[1]
  }
  clustermaptissues <- matrix(, 10, 10)
  for(cl in 1:10){
    clustermaptissues[cl,] <- sapply(1:10, function(kk){
      return(mean(na.exclude(gene_expression_scale[names(which(row_cluster == kk)), cl])))
    })
  }
  a1 <- apply(clustermaptissues, 1, which.max)
  cluster_result <- list()
  for(cl in 1:ncol(gene_expression)){
    cluster_result[[colnames(gene_expression)[cl]]] <- names(which(row_cluster == a1[cl]))
  }
  gene_expression_plot <- gene_expression
  colnames(gene_expression_plot) <- colnames(gene_expression)[apply(clustermaptissues, 2, which.max)]
  pdf(options$pdf1)
  ph=pheatmap(gene_expression[unlist(patternmarker),], scale = "row", show_rownames = F, show_colnames = T, cluster_rows = T,cluster_cols = T, cutree_rows = 10, cutree_cols = 10)
  dev.off()
  prePareGO <- prePareGO(dataset = dataSet)
  GSEA <- NULL
  for(cl in 1:length(cluster_result)){
    ap = runTopGO(cluster_result[[cl]],prePareGO)
    BP <- cbind(domain = "BP", ap$BP)
    CC <- cbind(domain = "CC", ap$CC)
    MF <- cbind(domain = "MF", ap$MF)
    BP <- BP[which((as.numeric(BP[,"Significant"]) >= 5) & (as.numeric(BP[,"elimFisher"]) < 0.05) ), ]
    CC <- CC[which((as.numeric(CC[,"Significant"]) >= 5) & (as.numeric(CC[,"elimFisher"]) < 0.05) ), ]
    MF <- MF[which((as.numeric(MF[,"Significant"]) >= 5) & (as.numeric(MF[,"elimFisher"]) < 0.05) ), ]
    result <- cbind(names(cluster_result)[cl], rbind(BP, CC, MF))
    GSEA <- rbind(GSEA, result)
  }
  colnames(GSEA)[1] <- c("tissue")
  GSEA <- as.matrix(GSEA)
  write.table(GSEA, options$output1, quote = F, sep = "\t", row.names = F)
  
}
SpatialCourseAnalysis()
#################################################################################################

