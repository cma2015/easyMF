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
                                'input5', 'i5', 2, 'character',
                                'input6', 'i6', 2, 'character',
                                'output1', 'o1', 2, 'character',
                                'output2', 'o2', 2, 'character',
                                'pdf1', 'o3', 2, 'character'
),
byrow=TRUE, ncol=4);
# Parse options
options(warn=-1)
options = getopt(option_specification);

#############################################################
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
  sstatMat <- matrix(NA, nrow = nrow(Amatrix), ncol = ncol(Amatrix), 
                     dimnames = NULL)
  rownames(sstatMat) <- rownames(Amatrix)
  
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
      sstatMat[,i] <- sstat
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
    return(list(PatternMarkers = ssgenes.th, PatternRanks = ssranks, PatternScore = sstatMat, PatternGenes = ssgenes))
  }
  else {
    return(PatternMarkers = ssgenes.th)
  }
}

expressionMarkers <- function(Amplitude = NULL, Pattern = NULL, geneExp = NULL, cpu = 1, PatternTC = NULL, PCC = 0.6, pvalue = 1e-03){
  expressionMarkers <- list()
  for(k in 1:ncol(Amplitude)){
    options(stringsAsFactors=F)
    options(scipen=999)
    suppressPackageStartupMessages(library(doParallel))
    suppressPackageStartupMessages(library(foreach))
    cl <-  makeCluster(cpu)
    registerDoParallel(cl)
    resMatList <- foreach(i = 1 : length(PatternTC$PatternMarkers[[k]])) %dopar%{#length(PatternTC[[k]])
      result <- NULL
      # for(j in 1:length(PatternTC)){
      a <- cor.test(as.numeric(Pattern[,k]), as.numeric(geneExp[PatternTC$PatternMarkers[[k]][i],]))
      result <- cbind(result,cbind(a$estimate, a$p.value))
      # }
      return(result)
    }
    stopCluster(cl)
    # corMattemp <- matrix(unlist(resMatList), ncol = length(PatternTC)*2, byrow = T)
    corMattemp <- do.call(rbind, resMatList)
    rownames(corMattemp) <- PatternTC$PatternMarkers[[k]]
    a1 <- which((corMattemp[,1] > PCC) & (corMattemp[,2] < pvalue))
    expressionMarkers[[paste0("Metagene",k)]] <- PatternTC$PatternMarkers[[k]][a1]
  }
  return(expressionMarkers)
}
#############################################################################
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


runTopGO <- function(geneID, prePareGO,  statistic = "fisher", algorithm = "elim",
                     topNodes = 200, plot = TRUE){
  library(ggplot2)
  
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
  }
  
  results <- list(BP = results.tab.BP, CC = results.tab.CC, MF = results.tab.MF)
  results
}


SpatialCourseAnalysis <- function(){
  
  
  gene_expression <- as.matrix(read.table(options$input1, header = T, stringsAsFactors = F))
  Amplitude <- as.matrix(read.table(options$input2, header = T, sep = "\t", stringsAsFactors = F))
  Pattern <- as.matrix(read.table(options$input3, header = T, sep = "\t", stringsAsFactors = F))
  Tissue <- as.matrix(read.table(options$input4, header = T, stringsAsFactors = F))
  require(biomaRt)
  mart <- useMart(biomart = "plants_mart", dataset =  options$input5, host = "plants.ensembl.org")
  GTOGO <- getBM(attributes = c("ensembl_gene_id", "go_id", "description","name_1006"), mart = mart)
  geneAnnota <- GTOGO[!duplicated(GTOGO[,1]), c(1,3)]
  rownames(geneAnnota) <- geneAnnota[,1]
  rownames(Tissue) <- Tissue[,1]
  cpu <- as.numeric(options$input6)

  
  
  pdf(options$pdf1)
  library(pheatmap)
  pheatmap(Pattern, labels_row = Tissue[rownames(Pattern),2])
  dev.off()
  ###############################################################
  TC <- patternMarkers(Amatrix = Amplitude, Pmatrix = t(Pattern), threshold = "cut", full = T)
  if(class(TC$PatternMarkers) == "matrix"){
    pp = TC$PatternMarkers
    TC$PatternMarkers <- list()
    for(i in 1:ncol(pp)){
      TC$PatternMarkers[[i]] <- pp[,i]
    }
  }
  geneExp_sig_gene <- expressionMarkers(Amplitude = Amplitude, Pattern = Pattern, geneExp = gene_expression
                                        , cpu = cpu, PatternTC = TC)
  
  
  
  geneID2GO <- by(GTOGO$go_id, GTOGO$ensembl_gene_id, function(x) as.character(x))
  prePareGO <- list(GTOGO = GTOGO, geneID2GO = geneID2GO)
  
  GSEA <- list()
  for(cl in 1:length(geneExp_sig_gene)){
    ap = runTopGO(geneExp_sig_gene[[cl]],prePareGO)
    BP <- cbind(domain = "BP", ap$BP)
    CC <- cbind(domain = "CC", ap$CC)
    MF <- cbind(domain = "MF", ap$MF)
    BP <- BP[which((as.numeric(BP[,"Significant"]) >= 5) & (as.numeric(BP[,"elimFisher"]) < 0.05) ), ]
    CC <- CC[which((as.numeric(CC[,"Significant"]) >= 5) & (as.numeric(CC[,"elimFisher"]) < 0.05) ), ]
    MF <- MF[which((as.numeric(MF[,"Significant"]) >= 5) & (as.numeric(MF[,"elimFisher"]) < 0.05) ), ]
    result <- rbind(BP, CC, MF)
    if(nrow(result) == 0){
      result[1,] <- NA
    }
    result <- cbind("Metagene"=paste0("Metagene",cl),result)
    GSEA[[names(geneExp_sig_gene)[cl]]] <- result
  }
  
  result <- do.call(rbind, GSEA)
  colnames(result) <- c("Metagene","Domain","GO","Term","Annotated","Significant","Expected","elimFisher","SignificantGene")
  
  sigStage_gene_annota <- NULL
  for(i in 1:length(geneExp_sig_gene)){
    sigStage_gene_annota <- rbind(sigStage_gene_annota,
                                  cbind("Metagene"=paste0("Metagene",i),
                                        geneID=geneAnnota[unlist(geneExp_sig_gene[[i]]),1],
                                        geneDescription=geneAnnota[unlist(geneExp_sig_gene[[i]]),2]))
  }
  
  write.table(sigStage_gene_annota, file = options$output1, quote = F, sep = "\t", row.names = F)
  write.table(result, file = options$output2, quote = F, sep = "\t", row.names = F)
  
}
SpatialCourseAnalysis()
#################################################################################################

