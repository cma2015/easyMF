###################### single cell analysis #########################################
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

# input1: code of matrix source from the user upload or fetch from GEO or default
# input2: code of fetch source
# input3: description of the matrix source
# pdf1: Statistics analysis conducted on samples to filter low quality samples
# output1: 
option_specification = matrix(c('input1', 'i1', 2, 'integer',
                                'input2', 'i2', 2, 'character',
                                'input3', 'i3', 2, 'character',
                                'input4', 'i4', 2, 'character',
                                'input5', 'i5', 2, 'character',
                                'input6', 'i6', 2, 'character',
                                'input7', 'i7', 2, 'character',
                                'input8', 'i8', 2, 'character',
                                'input9', 'i9', 2, 'character',
                                'input10', 'i10', 2, 'character',
                                'pdf1', 'o1', 2, 'character',
                                'output1', 'o2', 2, 'character'
                       
),
byrow=TRUE, ncol=4);
# Parse options
options(warn=-1)
options = getopt(option_specification);
PrepDR <- function(
  object,
  genes.use = NULL,
  use.imputed = FALSE,
  assay.type="RNA"
) {
  if (length(object@var.genes) == 0 && is.null(x = genes.use)) {
    stop("Variable genes haven't been set. Run MeanVarPlot() or provide a vector
          of genes names in genes.use and retry.")
  }
  if (use.imputed) {
    data.use <- t(x = scale(x = t(x = object@imputed)))
  } else {
    data.use <- GetAssayData(object, assay.type = assay.type,slot = "scale.data")
  }
  genes.use <- SetIfNull(x = genes.use, default = object@var.genes)
  genes.use <- unique(x = genes.use[genes.use %in% rownames(x = data.use)])
  genes.var <- apply(X = data.use[genes.use, ], MARGIN = 1, FUN = var)
  genes.use <- genes.use[genes.var > 0]
  genes.use <- genes.use[!is.na(x = genes.use)]
  data.use <- data.use[genes.use, ]
  return(data.use)
}
SetIfNull <- function(x, default) {
  if(is.null(x = x)){
    return(default)
  } else {
    return(x)
  }
}
# MF decomposition
MFdecomposition <- function(geneExp, cpu = 2, methods = "PCA"){
  
  if(methods == "PCA"){
    
    sampleNum <- ncol(geneExp)#update the sample number
    sampleName <- colnames(geneExp)
    geneNum <- nrow(geneExp)
    #PCA analysis on the genes
    geneExp <- t(geneExp)  ##make sure samples in rows and genes in columes
    genePCA <- prcomp(geneExp)
    GSPNum <- length(genePCA$sdev)#principal component number
    PCEV <- sapply(1:GSPNum, function(PCIdx){genePCA$sdev[PCIdx]^2/sum(genePCA$sdev^2)})#expalined variance for each PC
    PCCEV <- sapply(1:GSPNum, function(PCIdx){sum(PCEV[1:PCIdx])})#cumulative explained variance
    if(ncol(geneExp) < 15){
      selectGSPNum <- GSPNum
      selectPCIdx <- 1:GSPNum
      PCEVPlot <- PCEV[selectPCIdx]
      PCCEVPlot <- sapply(1:selectGSPNum, function(selectGSPNumIdx){sum(PCEVPlot[1:selectGSPNumIdx])})
    }else{
      #choose PC for further analysis based on internal consistency alpha coefficient > 0.7
      ###alpha coefficient
      pcaScores <- genePCA$x
      alphaVec <- NULL
      for(clflag in 1:floor(GSPNum/100)){
        t1 <- Sys.time()
        options(stringsAsFactors=F)
        options(scipen=999)
        suppressPackageStartupMessages(library(doParallel))
        suppressPackageStartupMessages(library(foreach))
        cl <-  makeCluster(cpu)
        registerDoParallel(cl)
        resMatList <- foreach(PCIdx = ((clflag - 1)*100+1) : (100 * clflag)) %dopar%{
          # for(PCIdx in 1:GSPNum){ ##for each PC
          variance <- sapply(1:ncol(geneExp), function(geneIdx){var(geneExp[,geneIdx]) * genePCA$rotation[geneIdx,PCIdx]^2})#item variance
          sumVariance <- sum(variance)#total variance
          alphaVec <- (ncol(geneExp)/(ncol(geneExp) - 1)) * (1 - sumVariance/(var(pcaScores[,PCIdx])))
          return(alphaVec)
          # cat("Alpha of component", PCIdx,": ", alphaVec[PCIdx],"\n")
        }
        stopCluster(cl)
        alphaVec <- c(alphaVec,unlist(resMatList))
        print(range(alphaVec))
        if(range(alphaVec)[1] < 0.5){
          break
        }
        t2 <- Sys.time()
        print(t2 - t1)
      }
      selectGSPNum <- length(which(alphaVec > 0.7))
      selectPCIdx <- which(alphaVec > 0.7)
      PCEVPlot <- PCEV[selectPCIdx]
      PCCEVPlot <- sapply(1:selectGSPNum, function(selectGSPNumIdx){sum(PCEVPlot[1:selectGSPNumIdx])})
    }
    ########################PDF3
    ##done PCA analysis on the genes
    #save PC for further analysis
    # Pattern <- genePCA$rotation[,selectPCIdx]
    # Amplitude <- genePCA$x[,selectPCIdx]
    Amplitude <- genePCA$rotation[,selectPCIdx]
    Pattern <- genePCA$x[,selectPCIdx]
    colnames(Pattern) <- paste0("NC",1:selectGSPNum)
    colnames(Amplitude) <- paste0("NC",1:selectGSPNum)
    
    return(list(Amplitude = Amplitude, Pattern = Pattern, sdev = genePCA$sdev[selectPCIdx]))
  }
  
  
  if(methods == "ICA"){
    library("ica")
    if(ncol(geneExp) < 15){
      t1 <- Sys.time()
      geneICA <- icafast(geneExp, nc = 2)
      t2 <- Sys.time()
      print(t2 - t1)
      residule <- mean(as.matrix(abs(geneExp - geneICA$S %*% t(geneICA$M))))
      cl = cl + 1
    }else{
      max_cl <- min(10, ncol(geneExp))
      residule <- matrix(,max_cl-1, 1)
      dresidule <- matrix(,max_cl-2, 1)
      t1 <- Sys.time()
      geneICA <- icafast(geneExp, nc = 2)
      t2 <- Sys.time()
      print(t2 - t1)
      residule[1,1] <- mean(as.matrix(abs(geneExp - geneICA$S %*% t(geneICA$M))))
      t1 <- Sys.time()
      geneICA <- icafast(geneExp, nc = 3)
      t2 <- Sys.time()
      print(t2 - t1)
      residule[2,1] <- mean(as.matrix(abs(geneExp - geneICA$S %*% t(geneICA$M))))
      dresidule[1,1] <- abs(residule[2,1] - residule[1,1])
      Pattern <- geneICA$M
      Amplitude <- geneICA$S
      # print(residule)
      for(cl in 4:max_cl){
        t1 <- Sys.time()
        geneICA <- icafast(geneExp, nc = cl)
        t2 <- Sys.time()
        print(t2 - t1)
        residule[cl-1,1] <- mean(as.matrix(abs(geneExp - geneICA$S %*% t(geneICA$M))))
        dresidule[cl-2,1] <- abs(residule[cl-1,1] - residule[cl-2,1])
        print(dresidule)
        if(dresidule[cl-2,1] <= dresidule[cl-3,1]){
          break
        }else{
          Pattern <- geneICA$M
          Amplitude <- geneICA$S
        }
      }
    }
    
    if(ncol(Pattern) == 2){
      t1 <- Sys.time()
      geneICA <- icafast(geneExp, nc = 3)
      t2 <- Sys.time()
      print(t2 - t1)
      Pattern <- geneICA$M
      Amplitude <- geneICA$S
      cl = 4  
    }
    #save PC for further analysis
    # Pattern <- genePCA$rotation[,selectPCIdx]
    # Amplitude <- genePCA$x[,selectPCIdx]
    colnames(Pattern) <- paste0("NC",1:(cl-1))
    rownames(Pattern) <- colnames(geneExp)
    colnames(Amplitude) <- paste0("NC",1:(cl-1))
    return(list(Amplitude = Amplitude, Pattern = Pattern, dresidule = dresidule))
  }
  
  if(methods == "NMF"){
    library(bignmf)
    if(ncol(geneExp) < 15){
      cl = ncol(geneExp)
      t1 <- Sys.time()
      geneNMF <- bignmf(geneExp, r = cl, max.iteration = 200, stop.condition = 1e-05) # geneNMF@fit@
      t2 <- Sys.time()
      print(t2 - t1)
      residule <- sum(abs(sum(geneExp - geneNMF$W %*% geneNMF$H)))
      Pattern <- t(geneNMF$H)
      Amplitude <- geneNMF$W
      cl = cl + 1
    }else{
      max_cl <- min(10, ncol(geneExp))
      residule <- matrix(,max_cl-1, 1)
      t1 <- Sys.time()
      geneNMF <- bignmf(geneExp, r = 2, max.iteration = 200, stop.condition = 1e-05) # geneNMF@fit@
      t2 <- Sys.time()
      print(t2 - t1)
      residule[1,1] <- sum(abs(sum(geneExp - geneNMF$W %*% geneNMF$H)))
      Pattern <- t(geneNMF$H)
      Amplitude <- geneNMF$W
      print(residule)
      for(cl in 3:max_cl){
        t1 <- Sys.time()
        geneNMF <- bignmf(geneExp, r = cl, max.iteration = 200, stop.condition = 1e-05) # geneNMF@fit@
        t2 <- Sys.time()
        print(t2 - t1)
        residule[cl-1,1] <- sum(abs(sum(geneExp - geneNMF$W %*% geneNMF$H)))
        print(residule)
        if(residule[cl-1,1] >= residule[cl-2,1]){
          break
        }else{
          Pattern <- t(geneNMF$H)
          Amplitude <- geneNMF$W
        }
      }
    }
    if(ncol(Pattern) == 2){
      t1 <- Sys.time()
      geneICA <- icafast(geneExp, nc = 3)
      t2 <- Sys.time()
      print(t2 - t1)
      Pattern <- geneICA$M
      Amplitude <- geneICA$S
      cl = 4  
    }
    #save PC for further analysis
    # Pattern <- genePCA$rotation[,selectPCIdx]
    # Amplitude <- genePCA$x[,selectPCIdx]
    colnames(Pattern) <- paste0("NC",1:(cl-1))
    rownames(Pattern) <- colnames(geneExp)
    colnames(Amplitude) <- paste0("NC",1:(cl-1))
    rownames(Amplitude) <- rownames(geneExp)
    write.table(Amplitude, file = options$output1, quote = FALSE, sep = "\t")
    write.table(Pattern, file = options$output2, quote = FALSE, sep = "\t")
    return(list(Amplitude = Amplitude, Pattern = Pattern, residule = residule))
  }
}

SingleCellAnalysis <- function(){
  
  gene_expression <- as.matrix(read.table(options$input1, header = T, stringsAsFactors = F))
 
  SpecScore <- as.matrix(read.table(options$input4, header = T, stringsAsFactors = F))
  # gene_expression <- as.matrix(read.table("/home/malab12/Downloads/GSE116614_RAW/single_cell.txt"))
  # cell_na_num <- apply(gene_expression, 1, function(ll){
  #   return(length(which(is.na(ll) == T)))
  # })
  # gene_expression <- gene_expression[which(cell_na_num == 0),]
  # cell_zero_num <- apply(gene_expression, 1, function(ll){
  #   return(length(which(ll > 0)))
  # })
  # gene_expression <- gene_expression[which(cell_zero_num > 0 ),]
  # gene_expression_sd <- apply(gene_expression, 1, sd)
  # gene_expression <- gene_expression[which(gene_expression_sd > 0 ),]
  # 
  disturb_code <- as.numeric(options$input5)
  if(disturb_code == 1){
    plast <- read.table(options$input6, stringsAsFactors = F, header = T)[,1]
    gene_expression <- gene_expression[setdiff(rownames(gene_expression), plast),]
  }
  
  library(Seurat)
  gene_expression_seurat <- CreateSeuratObject(gene_expression, normalization.method = "LogNormalize", do.scale = TRUE, do.center = F, display.progress = F)
  gene_expression_seurat <- FindVariableGenes(object = gene_expression_seurat, do.plot = F, display.progress = F)
  # length(genes@var.genes)
  data.use <- PrepDR(object = gene_expression_seurat, genes.use = gene_expression_seurat@var.genes, use.imputed = F, assay.type = "RNA")
  
  cpu <- as.numeric(options$input7)
  
  if(as.numeric(options$input8) == 4){
    Amplitude <- as.matrix(read.table(options$input9, header = T, stringsAsFactors = F))
    Pattern <-  as.matrix(read.table(options$input10, header = T, stringsAsFactors = F))
  }else{
    methods <- c("PCA", "ICA", "NMF")
    PCA_sc <- MFdecomposition(geneExp = data.use, cpu = 2, methods = methods[as.numeric(options$input8)])
    # pca_obj <- RunPCA(object = genes, pc.genes = genes@var.genes, pcs.compute = 100, do.print = F)
    Amplitude <- PCA_sc$Amplitude
    Pattern <- PCA_sc$Pattern
  }
  
  # 
  library(methods)
  pca.obj <- new(Class = "dim.reduction", gene.loadings = Amplitude, cell.embeddings = Pattern, sdev = 1:ncol(Amplitude), key = "PC")
  
  eval(expr = parse(text = paste0("gene_expression_seurat@dr$", "pca","<- pca.obj")))
  
  # ica_obj = JackStraw(pca_obj, num.pc = length(PCA_sc$sdev))
  # # PCA_obj <- ProjectPCA(object = ica_obj)
  # ica_obj <- JackStrawPlot(object = ica_obj, PCs = 1:length(which(alphaVec > 0.65)), nCol = 6)
  # informative <- ica_obj@dr$pca@jackstraw@overall.p.values[,2]
  
  gene_expression_seurat <- FindClusters(gene_expression_seurat, dims.use = 1:ncol(Amplitude), resolution = 0.6, force.recalc = T)
  tsne_obj <- RunTSNE(gene_expression_seurat, dims.use = 1:ncol(Amplitude))
 
  
  gene_expression_cluster <- GetClusters(object = gene_expression_seurat)
  gene_expression_cluster[,2] <- as.numeric(gene_expression_cluster[,2])
  a1 <- GetClusters(tsne_obj)
  
  # 
  # save(list = ls(), file = "/home/malab12/research/Galaxy/MF/gene_expression2.RData")
  # save(list = ls(), file = "/home/malab12/research/Galaxy/MF/temp.RData")
  
  # allmarker <- FindAllMarkers(object = tsne_obj)
  
  # result <- list()
  # for(cl in 0:(max(gene_expression_cluster[,2]) - 1)){
  #   print(cl)
  #   a1 <- FindMarkers(object = tsne_obj, ident.1 = cl)
  #   print(rownames(a1)[1])
  #   result[[paste0(cl)]] <- list(a1 = a1)
  # }
  # 
  # aa <- NULL
  # for(cl in 0:(max(gene_expression_cluster[,2]) - 1)){
  #   aa <- rbind(aa, cbind(result[[paste0(cl)]]$a1, cluster = cl, gene = rownames(result[[paste0(cl)]]$a1))  )
  # }
  # 
  # gene1 <- read.table("/home/malab12/research/Galaxy/MF/gene1.txt", stringsAsFactors = F, header = F)[,1]
  
  # t1 <- Sys.time()
  # specscore <- getAllSpec(data = gene_expression_seurat@scale.data, header = colnames(gene_expression)
  #                         , cuts = 2, distshape = 2)
  # t2 <- Sys.time()
  # print(t2 - t1)
  # source("/home/malab12/research/Galaxy/MF/spec.R")
  source("/galaxy/tools/TAMF/test_data/identity.r")
  # SpecScore <- as.matrix(read.table("/home/malab12/research/Galaxy/MF/specscore.txt", header = T, stringsAsFactors = F))
  ci <- list()
  ci[[1]] <- SpecScore
  t1 <- Sys.time()
  markers <- getMarkerList(ci, 20, rownames(gene_expression_seurat@raw.data))
  ICIScore <- getIdentity(data = gene_expression_seurat@raw.data, ci, markers, returnSig=FALSE, universe = rownames(gene_expression_seurat@raw.data))
  t2 <- Sys.time()
  print(t2 - t1)
  # save(list = ls(), file = "/home/malab12/research/Galaxy/MF/temp3.RData")
  score_mat <- ICIScore[[2]]
  score_result <- t(sapply(1:nrow(score_mat), function(ll){
    return(cbind(rownames(score_mat)[ll], colnames(score_mat)[which.max(score_mat[ll,])], score_mat[ll,which.max(score_mat[ll,])]))
  }))
  score_result <- as.data.frame(score_result, stringsAsFactors = F)
  rownames(score_result) <- score_result[,1]
  score_result <- score_result[,-1]

  cluster_name <- NULL
  for(l in 1:max(gene_expression_cluster[,2])){
    temp <- table((score_result[gene_expression_cluster[which(gene_expression_cluster[,2] == l),1],1]))
    cluster_name <- c(cluster_name, names(which.max(temp)))
    print(table((score_result[gene_expression_cluster[which(gene_expression_cluster[,2] == l),1],1])))
  }
  # 
  current.cluster.ids <- c(0:(max(gene_expression_cluster[,2]) - 1))
  new.cluster.ids <-cluster_name
  # new.cluster.ids <- c("Columella", "Hair Cells", "Non-hair", "QC", "Endodermis", "Cortex", "Endodermis"
  #                      , "Cortex", "Non-hair", "QC", "Hair Cells","Stele","Non-hair")
  tsne2_obj <- tsne_obj
  tsne2_obj@ident <- plyr::mapvalues(x = tsne2_obj@ident, from = current.cluster.ids, to = new.cluster.ids)
  pdf(file = options$pdf1)
  par(mfcol =c(1,2))
  TSNEPlot(tsne_obj, do.label = TRUE, label.size = 6)
  TSNEPlot(object = tsne2_obj, do.label = TRUE, label.size = 6)
  dev.off()
  
  cellType <- gene_expression_cluster
  for(cl in 0:(max(gene_expression_cluster[,2]) - 1)){
    cellType[which(cellType[,2] == (cl + 1)), 2] <- cluster_name[cl + 1]
  }
  colnames(cellType) <- c("cell_name", "type")
  write.table(cellType, options$output1, quote = F, sep = "\t", row.names = F)
  
  
}
SingleCellAnalysis()
#################################################################################################

