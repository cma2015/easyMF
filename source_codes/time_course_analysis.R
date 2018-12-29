######################time course analysis#########################################
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
# pdf1: visulization of the speficic module along the stages
option_specification = matrix(c('input1', 'i1', 2, 'character',
                                'input2', 'i2', 2, 'character',
                                'input3', 'i3', 2, 'character',
                                'input4', 'i4', 2, 'character',
                                'input5', 'i5', 2, 'character',
                                'input6', 'i6', 2, 'character',
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

TimeCourseAnalysis <- function(){
  
  # time_course <- as.matrix(read.table("/home/malab12/softwear/galaxy/tools/TAMF/test_data/default_maize_high_quality_gene_expression_timecourse.txt", header = T, stringsAsFactors = F))
  
  # Stage <- as.matrix(read_xlsx("/home/malab12/research/Galaxy/MF/design/polyploid_wheat.xlsx"))
  gene_expression <- as.matrix(read.table(options$input1, header = T, stringsAsFactors = F))
  Amplitude <- as.matrix(read.table(options$input2, header = T, stringsAsFactors = F))
  Pattern <- as.matrix(read.table(options$input3, header = T, stringsAsFactors = F))
  Stage <- as.matrix(read.table(options$input4, header = T, stringsAsFactors = F))
  require(biomaRt)
  plants_ensembl <- as.matrix(read.table("/galaxy/tools/TAMF/test_data/plants_ensemble.txt", header = F, stringsAsFactors = F, sep = "\t"))          
  dataSet <- plants_ensembl[as.numeric(options$input5),1]
  mart <- useMart(biomart = "plants_mart", dataset = dataSet, host = "plants.ensembl.org")
  GTOGO <- getBM(attributes = c("ensembl_gene_id", "go_id", "namespace_1003", "go_linkage_type", "description"), mart = mart)
  geneAnnota <- GTOGO[!duplicated(GTOGO[,1]), c(1,5)]
  rownames(geneAnnota) <- geneAnnota[,1]
  rownames(Stage) <- Stage[,1]
  Stage_specific <- Stage[rownames(Pattern),]
  cpu <- as.numeric(options$input6)
  
  
  if(ncol(Pattern) > 30){
    num <- 2
  }else{
    num <- ncol(Pattern)
  }
  
  library(CoGAPS)
  patternmarker <- patternMarkers(Amatrix = Amplitude[,1:num], Pmatrix = t(Pattern[,1:num]), threshold = "cut")
  aa <- cbind(Pattern[,1:num],as.numeric(t(Stage_specific[,2])))
  stage <- unique(aa[,ncol(aa)])
  bb <- matrix(,length(stage),ncol(aa) - 1)
  rownames(bb) <- stage
  colnames(bb) <- paste0("Pattern",1:(ncol(aa)-1) )
  for(j in 1:ncol(bb)){
    for(i in 1:length(stage)){
      bb[i,j] <- mean(as.numeric(aa[which((aa[,ncol(aa)] == stage[i]) ),j]))
    }
  }
  pdf(options$pdf1)
  plot(bb[,1], pch = 20, type = "o", xaxt="n", xlab = "Stage", ylab = "pattern value", ylim = range(bb))
  axis(1,c(1:nrow(bb)),labels=rownames(bb))
  for(k in 2:ncol(bb)){
    points(bb[,k], pch = 20, col = k, type="o")
  }
  dev.off()
  
  sigStage_gene_annota <- NULL
  ###############################
  # correlation
  ###NMF
  T1 <- Sys.time()
  time_course_sig_gene <- list()
  for(k in 1:length(patternmarker)){
    t1 <- Sys.time()
    options(stringsAsFactors=F)
    options(scipen=999)
    suppressPackageStartupMessages(library(doParallel))
    suppressPackageStartupMessages(library(foreach))
    cl <-  makeCluster(cpu)
    registerDoParallel(cl)
    resMatList <- foreach(i = 1 : length(patternmarker[[k]])) %dopar%{
      result <- NULL
      for(j in 1:length(patternmarker)){
        a <- cor.test(Pattern[,j], gene_expression[patternmarker[[k]][i],])
        result <- cbind(result,cbind(a$estimate, a$p.value))
      }
      return(result)
    }
    stopCluster(cl)
    t2 <- Sys.time()
    print(t2 - t1)
    corMattemp <- matrix(unlist(resMatList), ncol = length(patternmarker)*2, byrow = T)
    a1 <- which((corMattemp[,k*2-1] > 0.9) & (corMattemp[,k*2] < 1e-05))
    # for(l in setdiff(1:ncol(PCATC), k)){
    #   a1 <- intersect(a1, which((corMatPCAtemp[,l*2-1] < 0) & (corMatPCAtemp[,l*2] < 1e-05)))
    # }
    time_course_sig_gene[[paste0("Pattern",k)]] <- patternmarker[[k]][a1]
  }
  T2 <- Sys.time()
  print(T2 - T1)
  
  if(length(setdiff(unlist(time_course_sig_gene),geneAnnota[,1]))){
    a3 <- cbind(setdiff(unlist(time_course_sig_gene),geneAnnota[,1]), "")
    colnames(a3) <- colnames(geneAnnota)
    geneAnnota <- rbind(geneAnnota, a3)
    rownames(geneAnnota) <- geneAnnota[,1]
  }
  
  for(i in 1:length(time_course_sig_gene)){
    if(length(time_course_sig_gene[[paste0("Pattern",i)]])){
      sigStage_gene_annota <- rbind(sigStage_gene_annota, cbind(i, time_course_sig_gene[[paste0("Pattern",i)]],geneAnnota[time_course_sig_gene[[paste0("Pattern",i)]],2]))
    }else{
      sigStage_gene_annota <- rbind(sigStage_gene_annota, cbind(i,"",""))
      
    }
  }
  rownames(sigStage_gene_annota) <- 1:nrow(sigStage_gene_annota)
  colnames(sigStage_gene_annota) <- c("Pattern", "geneName", "geneDiscription")
  write.table(sigStage_gene_annota, file = options$output1, quote = F, sep = "\t", row.names = F)
}
TimeCourseAnalysis()
#################################################################################################


