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
option_specification = matrix(c('input1', 'i1', 2, 'character',
                                'input2', 'i2', 2, 'character',
                                'input3', 'i3', 2, 'character',
                                'input4', 'i4', 2, 'character',
                                'input5', 'i5', 2, 'character',
                                'input6', 'i6', 2, 'character',
                                'input7', 'i7', 2, 'character',
                                'input8', 'i8', 2, 'character',
                                'input9', 'i9', 2, 'character',
                                'input10', 'i10', 2, 'character',
                                'input11', 'i11', 2, 'character',
                                'pdf1', 'o1', 2, 'character',
                                'output1', 'o2', 2, 'character'
),
byrow=TRUE, ncol=4);
# Parse options
options(warn=-1)
options = getopt(option_specification);

# Number of permuations used to calculate significance
RAND_ITER=1000

# The minimal Spec score that can be used as a marker
# MIN_USEFUL_SPEC= 0.15
MIN_USEFUL_SPEC= 0.15
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
#####################################################################################
#
# Spec
#
# Calculate information scores
#
# Based on "Measuring cell identity in noisy biological systems", Birnbaum and Kussell, 
#
# Spec score is calculate as followed:
#
# I(level) = 1 + sum(P(type|level)logn(P(type|level)
# Spec(type) = sum(I(level)P(level|type))
#
#  Idan Efroni (ie10@nyu.edu)
#  Dec  14th, 2014


#####################################################################################
# specy
#
# Calculate spec score for a given gene
# 
# Arguments
#    gene - vector of gene expression vector for samples
#    header - vector of sample codes
#    binsize - size of the bin
# Returns
#    list of Spec scores for the gene, and the mean expression/binsize for each score


#######################################################################
#
# getAllSpec
#
# accepts a data table with the first row (Group) marks the group codes
#
# Arguments
#   data - expression matrix
#   header - vector of sample codes
#   medianfilter - filter gene whose median is above this
#    cuts - number of bins (0 to not use binning) (l)
#    distshape - max bin to be an acceptable background bin (u)
#
# Returns
#    The Spec data structure
#

getAllSpec <- function(data, header, cpu, medianfilter=0, cuts= FALSE, distshape=0) {
  
  # Settings
  float_prec <- 5
  if(distshape==0) {
    distshape = ceiling(cuts*0.33)
  }
  if(medianfilter==0) {
    medianfilter=max(data)
  }
  
  # Get Spec values row by row (gene by gene)
  specs <- matrix(nrow=nrow(data), ncol= length(unique(as.vector(header))), data=0)
  specMeanValues <- matrix(nrow=nrow(data), ncol= length(unique(as.vector(header))), data=0)
  
  t1 <- Sys.time()
  options(stringsAsFactors=F)
  options(scipen=999)
  suppressPackageStartupMessages(library(doParallel))
  suppressPackageStartupMessages(library(foreach))
  cl <-  makeCluster(cpu)
  registerDoParallel(cl)
  resMatList <- foreach(gene = 1 : nrow(data) ) %dopar%{#(length(genes[,1]) - 1) 
    # for (gene in 1:nrow(data)) {
    
    
    specy <- function(gene, header, binsize) {
      ntypes <- length(unique(as.vector(header)))
      
      bylist=list()
      bylist[[1]]=header
      # gene <- gene[sample(1:8)]
      # bylist[[1]]=names(gene)
      meanValues_t = aggregate(gene, bylist,median)
      meanValues= meanValues_t[,2]
      names(meanValues)=meanValues_t[,1]
      
      bins <- cut(gene, breaks=c(min(gene)-1, binsize), labels=F)
      bins[is.na(bins)] <- 2 # For values above hmax
      
      tab <- as.matrix(table(as.data.frame(cbind(header, bins))))
      
      # If not all bins are represented, fill a 0 matrix instead
      if (ncol(tab) < 2) {
        binned <- matrix(0,ntypes,2)
        for (i in colnames(tab)) {
          binned[,as.numeric(i)] <- tab[,i]
        }
      } else {
        binned <- tab
      }
      
      pleveltype <- binned/rowSums(binned)
      pleveltype[is.na(pleveltype)] <- 0
      ptypelevel <- t(pleveltype)/colSums(pleveltype)
      ptypelevel[is.na(ptypelevel)] <- 0
      logNptypelevel <- log(ptypelevel, ntypes)
      logNptypelevel[is.infinite(logNptypelevel)] = 0 # For log0 values
      #logNptypelevel[is.na(logNptypelevel)] = 0 # For log0 values
      ilevel <- 1 + rowSums(ptypelevel*logNptypelevel)
      ilevel[is.na(ilevel)] <- 0
      
      spec <- pleveltype%*%ilevel
      # Information is in the negative if the info is in the bottom bins
      infosign = as.integer(aggregate (gene, by=list(header), function(x) { mean(x)})[,2] > mean(gene)) -1 
      # we do not use "absent" information for the time being
      spec[infosign<0]=0
      
      return (list(spec, meanValues/binsize))
      
    }
    #######################################################################3
    #
    #  optimizebinsize
    #
    #  Identify the optimal binsize for a given gene
    #
    #  Arguments
    #    gene - vector of gene expression vector for samples
    #    header - vector of sample codes
    #    cuts - number of bins (0 to not use binning) (l)
    #    distshape - max bin to be an acceptable background bin (u)
    #  Returns:
    #	optimal binsize
    optimizebinsize <- function(gene, header, cuts=0, distshape=0) {
      
      background_cutoff=0
      
      if(cuts) {
        # gene should conform to a (very rough) power-law dist
        histo <- (table(cut(gene, breaks=seq(min(gene),max(gene),(max(gene)-min(gene))/cuts)))) >= ceiling(length(header)/cuts)
        if(length(which(histo))>0) {
          background_cutoff = max(which(histo))
        } else {
          background_cutoff=length(histo)
        }
      }
      
      # If gene does not conform, abort
      if(background_cutoff <= distshape ) {
        ntypes <- length(unique(as.vector(header)))
        minbin = min(gene) + (max(gene)-min(gene))*0.1
        maxbin = max(gene) - (max(gene)-min(gene))*0.1
        if(cuts) {# if using distribution shape filtering, just go for the cuts
          minbin = min(gene)+(max(gene)-min(gene))/cuts*(background_cutoff-1.5)
          maxbin = min(gene)+(max(gene)-min(gene))/cuts*(background_cutoff+0.5)
        }
        
        binsizes <- seq(minbin, maxbin, (maxbin-minbin)/50)
        infosum = matrix(nrow=length(binsizes), ncol=ntypes, data=0)
        rownames(infosum) = binsizes
        for(i in 1:length(binsizes)) {
          infosum[i,]= specy(gene, header, binsizes[i])[[1]][,1]
        }
        # return maximizie information for the highest sample
        binsizes[ match(max(infosum),apply(infosum,1,max))]
      } else {
        NA
      }
    }
    
    
    binsize <- optimizebinsize(data[gene,],header, cuts, distshape)
    if(!is.na(binsize) && median(data[gene,])<medianfilter) {
      specg <- specy(gene = data[gene,], header = header, binsize = binsize)
      # specs[gene,] = specg[[1]]
      return(specg[[1]])
    } else {
      # specs[gene,] = rep(0, ncol(specs))
      return(rep(0, ncol(specs)))
    }
  }
  stopCluster(cl)
  t2 <- Sys.time()
  print(t2 - t1)
  specs <- t(do.call(cbind, resMatList))
  
  rownames(specs) <- rownames(data)
  colnames(specs) <- unique(as.vector(header))
  
  
  list(round(specs, float_prec))
}
# MF decomposition
MFdecomposition <- function(geneExp, cpu = 1, methods = "PCA"){
  
  if(methods == "PCA"){
    
    #prepare the final used gene expression data
    sampleNum <- ncol(geneExp)#update the sample number
    sampleName <- colnames(geneExp)
    geneNum <- nrow(geneExp)
    #PCA analysis on the genes
    geneExp <- t(geneExp)  ##make sure samples in rows and genes in columes
    genePCA <- prcomp(geneExp)
    
    GSPNum <- length(genePCA$sdev)#principal component number
    PCEV <- sapply(1:GSPNum, function(PCIdx){genePCA$sdev[PCIdx]^2/sum(genePCA$sdev^2)})#expalined variance for each PC
    PCCEV <- sapply(1:GSPNum, function(PCIdx){sum(PCEV[1:PCIdx])})#cumulative explained variance
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
    ##done PCA analysis on the genes
    #save PC for further analysis
    # Pattern <- genePCA$rotation[,selectPCIdx]
    # Amplitude <- genePCA$x[,selectPCIdx]
    Amplitude <- genePCA$rotation[,selectPCIdx]
    Pattern <- genePCA$x[,selectPCIdx]
    colnames(Pattern) <- paste0("Metagene",1:selectGSPNum)
    rownames(Pattern) <- rownames(geneExp)
    colnames(Amplitude) <- paste0("Metagene",1:selectGSPNum)
    rownames(Amplitude) <- colnames(geneExp)
    
    return(list(Amplitude = Amplitude, Pattern = Pattern, sdev = genePCA$sdev[selectPCIdx]))
  }
  
  
  if(methods == "ICA"){
    library("ica")
    max_cl <- min(30, ncol(geneExp))
    residule <- matrix(,max_cl-1, 1)
    dresidule <- matrix(,max_cl-2, 1)
    set.seed(1)
    geneICA <- icafast(geneExp, nc = 2)
    residule[1,1] <- mean(as.matrix(abs(geneExp - geneICA$S %*% t(geneICA$M))))
    Pattern <- geneICA$M
    Amplitude <- geneICA$S
    set.seed(1)
    geneICA <- icafast(geneExp, nc = 3)
    residule[2,1] <- mean(as.matrix(abs(geneExp - geneICA$S %*% t(geneICA$M))))
    dresidule[1,1] <- residule[1,1] - residule[2,1]
    print(residule[1:2,])
    if(dresidule[1,1] < 0){
      pdf(file = options$pdf1)
      plot(residule[1:2,], type = "l", lwd = 2, xlab = "Optimal metagenes", ylab = "Residule", xaxt="n", main = "ICA statistics")
      axis(side = 1,at = 1:2
           , labels = c(2:3))       
      dev.off()
      rownames(Amplitude) <- rownames(geneExp)
      colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
      rownames(Pattern) <- colnames(geneExp)
      colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
      
    }else{
      # print(residule)
      Pattern <- geneICA$M
      Amplitude <- geneICA$S
      for(cl in 4:max_cl){
        t1 <- Sys.time()
        set.seed(seed)
        geneICA <- icafast(geneExp, nc = cl)
        t2 <- Sys.time()
        print(t2 - t1)
        residule[cl-1,1] <- mean(as.matrix(abs(geneExp - geneICA$S %*% t(geneICA$M))))
        dresidule[cl-2,1] <- residule[cl-2,1] - residule[cl-1,1]
        print(dresidule[1:(cl-1),])
        if(dresidule[cl-2,1] <= dresidule[cl-3,1]){
          break
        }else{
          Pattern <- geneICA$M
          Amplitude <- geneICA$S
        }
      }
    } 
    #save PC for further analysis
    # Pattern <- genePCA$rotation[,selectPCIdx]
    # Amplitude <- genePCA$x[,selectPCIdx]
    rownames(Amplitude) <- rownames(geneExp)
    colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
    rownames(Pattern) <- colnames(geneExp)
    colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
    return(list(Amplitude = Amplitude, Pattern = Pattern, dresidule = dresidule))
  }
  
  if(methods == "NMF"){
    library(bignmf)
    max_cl <- min(30, ncol(geneExp))
    residule <- matrix(,max_cl-1, 1)
    dresidule <- matrix(,max_cl-2, 1)
    geneNMF <- bignmf(geneExp, r = 2, max.iteration = 200, stop.condition = 1e-05)
    residule[1,1] <- sum(abs(sum(geneExp - geneNMF$W %*% geneNMF$H)))
    Pattern <- t(geneNMF$H)
    Amplitude <- geneNMF$W
    geneNMF <- bignmf(geneExp, r = 3, max.iteration = 200, stop.condition = 1e-05)
    residule[2,1] <- sum(abs(sum(geneExp - geneNMF$W %*% geneNMF$H)))
    dresidule[1,1] <- residule[1,1] - residule[2,1]
    print(residule[1:2,])
    if(dresidule[1,1] < 0){
      pdf(file = options$pdf1)
      plot(residule[1:2,], type = "l", lwd = 2, xlab = "Optimal metagenes", ylab = "Residule", xaxt="n", main = "NMF statistics")
      axis(side = 1,at = 1:2
           , labels = c(2:3))       
      dev.off()
      rownames(Amplitude) <- rownames(geneExp)
      colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
      rownames(Pattern) <- colnames(geneExp)
      colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
      
    }else{
      # print(residule)
      Pattern <- t(geneNMF$H)
      Amplitude <- geneNMF$W
      for(cl in 4:max_cl){
        geneNMF <- bignmf(geneExp, r = cl, max.iteration = 200, stop.condition = 1e-05)
        residule[cl-1,1] <- sum(abs(sum(geneExp - geneNMF$W %*% geneNMF$H)))
        dresidule[cl-2,1] <- residule[cl-2,1] - residule[cl-1,1]
        print(dresidule[1:(cl-1),])
        if(dresidule[cl-2,1] <= dresidule[cl-3,1]){
          break
        }else{
          Pattern <- t(geneNMF$H)
          Amplitude <- geneNMF$W
        }
      }
      pdf(file = options$pdf1)
      plot(residule[1:(cl-1),], type = "l", lwd = 2, xlab = "Optimal metagenes", ylab = "Residule", xaxt="n", main = "NMF statistics")
      axis(side = 1,at = 1:(cl-1)
           , labels = c(2:cl)
      )
      dev.off()
    }  
    rownames(Amplitude) <- rownames(geneExp)
    colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
    rownames(Pattern) <- colnames(geneExp)
    colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
    return(list(Amplitude = Amplitude, Pattern = Pattern, residule = residule))
  }
}
#####################################################################################
#
# Identity
#
# Calculate Index of Cell Identity (ICI)
#
#  Idan Efroni (ie10@nyu.edu)
#  Dec  14th, 2014
#####################################################################################



#####################################################################################
#
#  getMarkerList
#
#  returns a list of markers for each tissue
#
#  Arguments
#    ci - spec data structure
#    info - information threshold
#    universe - list of gene names from which to choose markers
#
# markers = getMarkerList(ci, i/2, rownames(data))
getMarkerList <- function(ci, info, universe) {
  markers=list()
  ci[[1]][ci[[1]]<0]=0
  for(i in 1:ncol(ci[[1]])) {
    ci_sub = ci[[1]][intersect(rownames(ci[[1]]),universe),i]
    ci_sub = ci_sub[ci_sub>= MIN_USEFUL_SPEC]
    
    # ci_max_sub = ci[[2]][names(ci_sub),i]
    # cum_ci <- cumsum(ci_sub[order(ci_sub, ci_max_sub, decreasing=TRUE)])
    cum_ci <- cumsum(ci_sub[order(ci_sub, decreasing=TRUE)])
    cum_ci <- cum_ci[1:head(which(cum_ci==max(cum_ci)),n=1)]
    markers[[i]] <- names(which(cum_ci<info))
  }
  names(markers)=colnames(ci[[1]])
  markers
}
#####################################################################################
#
#  getRandomBackground
#
#  returns a vector of RAND_ITER length of ICI based on random markers
#   used for significance testing
#
#  Arguments
#    cell - vector of gene expression values
#    ci_for_mark - spec scores
#    universe - list of gene names from which to choose markers
#    marker_num - number of markers to use


getRandomBackground <- function(cell, ci_for_mark, universe, marker_num) {
  
  rand_id_scores=vector()
  cell = cell[universe]
  
  for(i in 1:RAND_ITER) {
    # get marker set
    marker_set = sample(1:length(cell), marker_num)
    rand_id_scores[i] = getIdentityScore(cell, ci_for_mark, marker_set)
  }
  rand_id_scores[is.na(rand_id_scores)]=0
  rand_id_scores
}

#####################################################################################
#
#  getIdentityScore
#
#  returns ICIs for a given cell
#
#  Arguments
#    cell - vector of gene expression values
#    ci_for_mark - spec scores for the markers
#    markers - list of markers from getMarkerList

getIdentityScore <- function(cell, ci_for_mark, markers){
  mean(cell[markers]* ci_for_mark) * ((sum(cell[markers]>0))/length(markers))
}
#####################################################################################
#
#  getIdentity
#
#  returns ICIs for a data matrix
#
#  Arguments
#    data - gene expression matrix
#    ci - spec data structure
#    markers - list of markers from getMarkerList
#    returnSig - should significance be calculated
#    universe - what subset of genes should be used for randomizations
# a3 = getIdentity(data, ci, markers, returnSig=TRUE, universe=rownames(data))
getIdentity <- function(data, ci, markers, returnSig=FALSE, universe=c()) {
  
  hs_scoremat <- matrix(nrow=ncol(data), ncol=length(markers))
  
  colnames(hs_scoremat) <- names(markers)
  rownames(hs_scoremat) <- colnames(data)
  calls <- hs_scoremat
  sig <- calls
  all_markers = unlist(markers)
  
  for(cell in 1:nrow(hs_scoremat)) {
    markers_cell = data[, cell]
    for(mark in 1:length(markers)) {
      hs_scoremat[cell, mark] = getIdentityScore(data[,cell], ci[[1]][markers[[mark]],mark], markers[[mark]])
      calls[cell, mark] = sum(data[markers[[mark]],cell]>0)
      if(returnSig) {
        sig[cell,mark] <- 1-which(order(c(hs_scoremat[cell,mark], getRandomBackground(markers_cell, ci[[1]][markers[[mark]],mark], universe, length(markers[[mark]]))))==1)/(RAND_ITER+1)
      }
    }
  }
  
  hs_scoremat_norm <- hs_scoremat
  for(i in 1:nrow(hs_scoremat_norm)) { hs_scoremat_norm[i,] = hs_scoremat_norm[i,]/sum(hs_scoremat_norm[i,]) }
  hs_scoremat_norm[is.nan(hs_scoremat_norm)]=0
  if(returnSig) {
    list(hs_scoremat_norm, hs_scoremat, calls, sig, matrix(nrow=nrow(sig), ncol=ncol(sig), p.adjust(sig, "BH")))
  } else {
    list(hs_scoremat_norm, hs_scoremat)
  }
}

#####################################################################################
#
#  getICISignal
#
#  returns vector of ICI Signal for cumulative information threshold of 0.5-100 at 0.5 intervals
#
#  Arguments
#    ci - spec data structure
#    data - gene expression matrix

getICISignal <- function(ci, data) {
  ic_mn_list <- list()
  
  for(i in 1:200) {
    
    id_tmp = getIdentity(data, ci, getMarkerList(ci, i/2, rownames(data)), returnSig=FALSE)
    ic_mn_list[[i]] = id_tmp[[1]]
    cat(".")
    flush.console()
  }
  
  maxID = vector()
  for(i in 1:200) {
    maxID[i] = mean(apply(ic_mn_list[[i]],1,max))
  }
  maxID
}

#####################################################################################
#
#  getICIVariability
#
#  returns a list of vectors of ICI Variability for cumulative information threshold of 0.5-100 at 0.5 intervals (one vector per cell)
#
#   the variability at each cumulative infromartion threshold (i) is the eucledean distance between the ICI vectors for 
#     threshold i and threshold i+1
#
#  Arguments
#    ci - spec data structure
#    data - gene expression matrix

getICIVariability <- function(ci, data) {
  ic_mn_list <- list()
  for(i in 1:ncol(data)) {
    ic_mn_list[[i]] = matrix(nrow=200, ncol= ncol(ci[[1]]))
    colnames(ic_mn_list[[i]]) = colnames(ci[[1]])
  }
  
  for(i in 1:200) {
    id_tmp = getIdentity(data, ci, getMarkerList(ci, i/2, rownames(data)), returnSig=FALSE)
    for(x in 1:ncol(data)) {
      print(x)
      ic_mn_list[[x]][i,] = id_tmp[[1]][x,]
    }
    cat(".")
    flush.console()
  }
  for(i in 1:ncol(data)) {
    ic_mn_list[[i]][is.na(ic_mn_list[[i]])]=0
  }
  
  variability=list()
  for(cell in 1:ncol(data)) {
    variability[[cell]]=vector()
    for(info in 1:198) {
      variability[[cell]][info] = dist(rbind(ic_mn_list[[cell]][info,], ic_mn_list[[cell]][info+2,]))
    }
  }
  variability
}

SingleCellAnalysis <- function(){

  gene_expression <- as.matrix(read.table(options$input1, header = T, stringsAsFactors = F))
  cpu <- as.numeric(options$input7)
  if(as.numeric(options$input11) == 1){
    SpecScore <- as.matrix(read.table(options$input4, header = T, stringsAsFactors = F))
  }else{
    library(stringr)
    gene_expression_spec <- as.matrix(read.table(options$input4, header = T, stringsAsFactors = F))
    aa <- colnames(gene_expression_spec)
    bb <- unlist(sapply(1:length(aa), function(ll){
     
      return( str_replace(string = aa[ll], replacement = "", pattern= "\\.\\d+"))
    }))
    colnames(gene_expression_spec) <- bb
    SpecScore <- getAllSpec(data = gene_expression_spec, header = colnames(gene_expression_spec), cpu = cpu)
    SpecScore <- SpecScore[[1]]
  }
 
  disturb_code <- options$input5
  if(disturb_code == "1"){
    plast <- read.table(options$input6, stringsAsFactors = F, header = T)[,1]
    gene_expression <- gene_expression[setdiff(rownames(gene_expression), plast),]
  }
  
  library(Seurat)
  gene_expression_seurat <- CreateSeuratObject(gene_expression, normalization.method = "LogNormalize", do.scale = TRUE, do.center = F, display.progress = F)
  gene_expression_seurat <- FindVariableGenes(object = gene_expression_seurat, do.plot = F, display.progress = F)
  # length(genes@var.genes)
  data.use <- PrepDR(object = gene_expression_seurat, genes.use = gene_expression_seurat@var.genes, use.imputed = F, assay.type = "RNA")
  
  
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
  tsne_obj <- RunTSNE(gene_expression_seurat, dims.use = 1:ncol(Amplitude), check_duplicates = FALSE)
 
  
  gene_expression_cluster <- GetClusters(object = gene_expression_seurat)
  gene_expression_cluster[,2] <- as.numeric(gene_expression_cluster[,2])
  a1 <- GetClusters(tsne_obj)

  ci <- list()
  ci[[1]] <- SpecScore
  
  markers <- getMarkerList(ci, 20, rownames(gene_expression_seurat@raw.data))
  ICIScore <- getIdentity(data = gene_expression_seurat@raw.data, ci, markers, returnSig=FALSE, universe = rownames(gene_expression_seurat@raw.data))
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

