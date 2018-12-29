######################decomposition gene expression profile#########################################
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

# input1: input the gene expression profile
# input2: input the cpu numbers to parallelly compute
# input3: decomposition algorithms code
# output1: molecular relationships
# output2: sample relationships
# pdf1: Statistics analysis of the decomposition
option_specification = matrix(c('input1', 'i1', 2, 'character',
                                'input2', 'i2', 2, 'integer',
                                'input3', 'i3', 2, 'integer',
                                'output1', 'o1', 2, 'character',
                                'output2', 'o2', 2, 'character',
                                'pdf1', 'o3', 2, 'character'
),
byrow=TRUE, ncol=4);
# Parse options
options(warn=-1)
options = getopt(option_specification);


MatrixFactorization <- function(){
  
  if(options$input3 == 4){
    
    system(paste0("cat ", options$input1,  "> ", options$output1))
    system(paste0("cat ", options$input2, "> ", options$output2))
    
    
  }else{
    geneExp <- as.matrix(read.table(options$input1, header = T, stringsAsFactors = F))
    
    if(options$input3 == 1){
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
          cl <-  makeCluster(options$input2)
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
      
      pdf(file = options$pdf1)
      # pdf(file = "GSP_Statistics.pdf")
      plot(PCEVPlot[selectPCIdx], ylim = c(0,1), xlab = "Summary gene set patterns", ylab = "Variance Explaned (%)", main = "PCA statistics", pch = 20)
      points(PCCEVPlot[selectPCIdx], col = 2, pch = 20)
      points(alphaVec[selectPCIdx], col = 4, pch = 20)
      legend(x = length(selectPCIdx)*0.4, y = 0.6, pch = c(20,20,20), col = c(1,2,4), bty = "n",
             legend = c("Explained variance", "Cumulative explained variance", "Cronbach's alpha"))
      dev.off()
      
      ##done PCA analysis on the genes
      #save PC for further analysis
      # Pattern <- genePCA$rotation[,selectPCIdx]
      # Amplitude <- genePCA$x[,selectPCIdx]
      Amplitude <- genePCA$rotation[,selectPCIdx]
      Pattern <- genePCA$x[,selectPCIdx]
      colnames(Pattern) <- paste0("NC",1:selectGSPNum)
      colnames(Amplitude) <- paste0("NC",1:selectGSPNum)
      
      write.table(Amplitude, file = options$output1, quote = FALSE, sep = "\t")
      write.table(Pattern, file = options$output2, quote = FALSE, sep = "\t")
    }
    
    
    if(options$input3 == 2){
      # geneExp <- t(geneExp)  ##make sure samples in rows and genes in columes
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
      
      pdf(file = options$pdf1)
      # pdf(file = "GSP_Statistics.pdf")
      if(length(residule) > 1){
        plot(residule[1:(cl-1),1], xaxt = "n", yaxt = "n", 
             xlab = "Number of components", ylab = "Residue", main = "ICA statistics"
             , pch = 20, lwd = 2, type = "o")
        axis(1,1:(cl-1),(1:(cl-1))+1)
        points((cl-2), residule[cl-2,1], pch = 20, lwd = 4, col = 2)
        abline(v = cl-2, col = 2, lwd = 3, lty = 2)
      }else{
        plot(residule, xaxt = "n", yaxt = "n", 
             xlab = "Number of components", ylab = "Residue", main = "ICA statistics"
             , pch = 20, lwd = 2, type = "o")
        axis(1,1,cl-1)
      }
      
      dev.off()
      
      if(ncol(Pattern) == 2){
        
        t1 <- Sys.time()
        geneICA <- icafast(geneExp, nc = 3)
        t2 <- Sys.time()
        print(t2 - t1)
        Pattern <- geneICA$M
        Amplitude <- geneICA$S
        cl = 4  
      }
      
      ##done PCA analysis on the genes
      #save PC for further analysis
      # Pattern <- genePCA$rotation[,selectPCIdx]
      # Amplitude <- genePCA$x[,selectPCIdx]
      colnames(Pattern) <- paste0("NC",1:(cl-1))
      rownames(Pattern) <- colnames(geneExp)
      colnames(Amplitude) <- paste0("NC",1:(cl-1))
      write.table(Amplitude, file = options$output1, quote = FALSE, sep = "\t")
      write.table(Pattern, file = options$output2, quote = FALSE, sep = "\t")
    }
    
    if(options$input3 == 3){
      
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
      
      pdf(file = options$pdf1)
      # pdf(file = "GSP_Statistics.pdf")
      if(length(residule) > 1){
        plot(residule[1:(cl-1),1], xaxt = "n", yaxt = "n", 
             xlab = "Number of components", ylab = "Residue", main = "NMF statistics"
             , pch = 20, lwd = 2, type = "o")
        axis(1,1:(cl-1),(1:(cl-1))+1)
        points((cl-2), residule[cl-2,1], pch = 20, lwd = 4, col = 2)
        abline(v = cl-2, col = 2, lwd = 3, lty = 2)
      }else{
        plot(residule, xaxt = "n", yaxt = "n", 
             xlab = "Number of components", ylab = "Residue", main = "NMF statistics"
             , pch = 20, lwd = 2, type = "o")
        axis(1,1,cl-1)
      }
      dev.off()
      if(ncol(Pattern) == 2){
        
        t1 <- Sys.time()
        geneICA <- icafast(geneExp, nc = 3)
        t2 <- Sys.time()
        print(t2 - t1)
        Pattern <- geneICA$M
        Amplitude <- geneICA$S
        cl = 4  
      }
      ##done PCA analysis on the genes
      #save PC for further analysis
      # Pattern <- genePCA$rotation[,selectPCIdx]
      # Amplitude <- genePCA$x[,selectPCIdx]
      colnames(Pattern) <- paste0("NC",1:(cl-1))
      rownames(Pattern) <- colnames(geneExp)
      colnames(Amplitude) <- paste0("NC",1:(cl-1))
      rownames(Amplitude) <- rownames(geneExp)
      write.table(Amplitude, file = options$output1, quote = FALSE, sep = "\t")
      write.table(Pattern, file = options$output2, quote = FALSE, sep = "\t")
      
    }
    
  }
  
 
  
  
  
}
MatrixFactorization()
#################################################################################################
