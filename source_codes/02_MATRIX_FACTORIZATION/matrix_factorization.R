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
# input3: options for optimal metagenes
# input4: decomposition algorithms code
# input5: optimal metagenes code
# output1: molecular relationships
# output2: sample relationships
# pdf1: Statistics analysis of the decomposition
option_specification = matrix(c('input1', 'i1', 2, 'character',
                                'input2', 'i2', 2, 'integer',
                                'input3', 'i3', 2, 'character',
                                'input4', 'i4', 2, 'integer',
                                'input5', 'i5', 2, 'integer',
                                'output1', 'o1', 2, 'character',
                                'output2', 'o2', 2, 'character',
                                'pdf1', 'o3', 2, 'character'
),
byrow=TRUE, ncol=4);
# Parse options
options(warn=-1)
options = getopt(option_specification);
MatrixFactorization <- function(){
  
  if(options$input4 == 4){
  }else{
    geneExp <- as.matrix(read.table(options$input1, header = T, stringsAsFactors = F))
    
    if(options$input5 == "4"){
      if(as.numeric(options$input3) >= min(dim(geneExp))){
        cat("The specified number of metagenes is too large, we change it the all metagene!")
        options$input5 == "2"
        pdf(file = options$pdf1)
        plot(c(1,1),c(1,1),pch=20, cex=0.01, type = "o", xaxt="n",yaxt="n",xlab="",ylab="")
        text(x=1,y=1.1, labels = "The specified number of metagenes is too large")
        text(x=1,y=0.8, labels = "we change it the all metagene!")
        dev.off()
      }
    }else{
      if(as.numeric(options$input3) >= min(dim(geneExp))){
        cat("The specified number of metagenes is too large, we change it the all metagene!")
        options$input5 == "2"
      }
    }
    
    
    if(options$input4 == 1){
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
      
      if(options$input5 == "2"){
        selectGSPNum <- GSPNum
        selectPCIdx <- 1:GSPNum
        # PCEVPlot <- PCEV[selectPCIdx]
        # PCCEVPlot <- sapply(1:selectGSPNum, function(selectGSPNumIdx){sum(PCEVPlot[1:selectGSPNumIdx])})
      }else if(options$input5 == "1"){
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
      }else if(options$input5 == "4"){
        #choose PC started with specified number for further analysis based on internal consistency alpha coefficient > 0.7
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
        selectPCIdx <- union(c(1:as.numeric(options$input3)), selectPCIdx)
        PCEVPlot <- PCEV[selectPCIdx]
        PCCEVPlot <- sapply(1:selectGSPNum, function(selectGSPNumIdx){sum(PCEVPlot[1:selectGSPNumIdx])})
        
      }else{
        selectPCIdx <- 1:(as.numeric(options$input3))
        selectGSPNum <- length(selectPCIdx)
      }
      
      ########################PDF3
      if((options$input5 == "1") | (options$input5 == "4")){
        pdf(file = options$pdf1)
        # pdf(file = "GSP_Statistics.pdf")
        plot(PCEVPlot[selectPCIdx], ylim = c(0,1), xlab = "Optimal metagenes", ylab = "Variance Explaned (%)", main = "PCA statistics", pch = 20)
        points(PCCEVPlot[selectPCIdx], col = 2, pch = 20)
        points(alphaVec[selectPCIdx], col = 4, pch = 20)
        legend(x = length(selectPCIdx)*0.4, y = 0.6, pch = c(20,20,20), col = c(1,2,4), bty = "n",
               legend = c("Explained variance", "Cumulative explained variance", "Cronbach's alpha"))
        dev.off()
      }
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
      write.table(Amplitude, file = options$output1, quote = FALSE, sep = "\t")
      write.table(Pattern, file = options$output2, quote = FALSE, sep = "\t")
    }
    
    if(options$input4 == 2){
      # geneExp <- t(geneExp)  ##make sure samples in rows and genes in columes
      library("ica")
      
      if(options$input5 == "2"){
        geneICA <- icafast(geneExp, nc = min(nrow(geneExp),ncol(geneExp)))
        Pattern <- geneICA$M
        Amplitude <- geneICA$S
        rownames(Amplitude) <- rownames(geneExp)
        colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
        rownames(Pattern) <- colnames(geneExp)
        colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
      }else if(options$input5 == "1"){
        
        max_cl <- min(1000, min(dim(geneExp)))
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
            set.seed(1)
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
          pdf(file = options$pdf1)
          plot(residule[1:(cl-1),], type = "l", lwd = 2, xlab = "Optimal metagenes", ylab = "Residule", xaxt="n", main = "ICA statistics")
          axis(side = 1,at = 1:(cl-1)
               , labels = c(2:cl)
          )
          dev.off()
          rownames(Amplitude) <- rownames(geneExp)
          colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
          rownames(Pattern) <- colnames(geneExp)
          colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
          
        }  
      }else if(options$input5 == "4"){
        
        num=as.numeric(options$input3)
        max_cl <- min(num+1000, min(dim(geneExp)))
        residule <- matrix(,max_cl-num+1, 1)
        dresidule <- matrix(,max_cl-num, 1)
        set.seed(1)
        geneICA <- icafast(geneExp, nc = num)
        residule[1,1] <- mean(as.matrix(abs(geneExp - geneICA$S %*% t(geneICA$M))))
        Pattern <- geneICA$M
        Amplitude <- geneICA$S
        set.seed(1)
        geneICA <- icafast(geneExp, nc = num+1)
        residule[2,1] <- mean(as.matrix(abs(geneExp - geneICA$S %*% t(geneICA$M))))
        dresidule[1,1] <- residule[1,1] - residule[2,1]
        print(residule[1:2,])
        if(dresidule[1,1] < 0){
          pdf(file = options$pdf1)
          plot(residule[1:2,], type = "l", lwd = 2, xlab = "Optimal metagenes", ylab = "Residule", xaxt="n", main = "ICA statistics")
          axis(side = 1,at = 1:2
               , labels = c(num:(num+1)))       
          dev.off()
          rownames(Amplitude) <- rownames(geneExp)
          colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
          rownames(Pattern) <- colnames(geneExp)
          colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
          
        }else{
          # print(residule)
          Pattern <- geneICA$M
          Amplitude <- geneICA$S
          for(cl in (num+2):max_cl){
            t1 <- Sys.time()
            set.seed(1)
            geneICA <- icafast(geneExp, nc = cl)
            t2 <- Sys.time()
            print(t2 - t1)
            residule[cl-num+1,1] <- mean(as.matrix(abs(geneExp - geneICA$S %*% t(geneICA$M))))
            dresidule[cl-num,1] <- residule[cl-num,1] - residule[cl-num+1,1]
            print(dresidule[1:(cl-num+1),])
            if((dresidule[cl-num,1] <= dresidule[cl-num-1,1]) | (dresidule[cl-num,1] <= 0) ){
              break
            }else{
              Pattern <- geneICA$M
              Amplitude <- geneICA$S
            }
          }
          pdf(file = options$pdf1)
          plot(residule[1:(cl-num+1),], type = "l", lwd = 2, xlab = "Optimal metagenes", ylab = "Residule", xaxt="n", main = "ICA statistics")
          axis(side = 1,at = 1:(cl-num+1)
               , labels = c(num:cl)
          )
          dev.off()
          rownames(Amplitude) <- rownames(geneExp)
          colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
          rownames(Pattern) <- colnames(geneExp)
          colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
          
        }
        
      }else{
        set.seed(1)
        geneICA <- icafast(geneExp, nc = as.numeric(options$input3))
        Pattern <- geneICA$M
        Amplitude <- geneICA$S
        rownames(Amplitude) <- rownames(geneExp)
        colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
        rownames(Pattern) <- colnames(geneExp)
        colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
      }
      write.table(Amplitude, file = options$output1, quote = FALSE, sep = "\t")
      write.table(Pattern, file = options$output2, quote = FALSE, sep = "\t")
    }
    
    if(options$input4 == 3){
      # geneExp <- t(geneExp)  ##make sure samples in rows and genes in columes
      library(bignmf)
      
      if(options$input5 == "2"){
        geneNMF <- bignmf(geneExp, r = min(nrow(geneExp),ncol(geneExp)), max.iteration = 200, stop.condition = 1e-05)
        Pattern <- t(geneNMF$H)
        Amplitude <- geneNMF$W
        rownames(Amplitude) <- rownames(geneExp)
        colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
        rownames(Pattern) <- colnames(geneExp)
        colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
      }else if(options$input5 == "1"){
        
        max_cl <- min(1000, min(dim(geneExp)))
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
          rownames(Amplitude) <- rownames(geneExp)
          colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
          rownames(Pattern) <- colnames(geneExp)
          colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
          
        }  
      }else if(options$input5 == "4"){
        
        num=as.numeric(options$input3)
        max_cl <- min(num+1000, min(dim(geneExp)))
        residule <- matrix(,max_cl-num+1, 1)
        dresidule <- matrix(,max_cl-num, 1)
        geneNMF <- bignmf(geneExp, r = num, max.iteration = 200, stop.condition = 1e-05)
        residule[1,1] <-  sum(abs(sum(geneExp - geneNMF$W %*% geneNMF$H)))
        Pattern <- t(geneNMF$H)
        Amplitude <- geneNMF$W
        geneNMF <- bignmf(geneExp, r = num+1, max.iteration = 200, stop.condition = 1e-05)
        residule[2,1] <-  sum(abs(sum(geneExp - geneNMF$W %*% geneNMF$H)))
        dresidule[1,1] <- residule[1,1] - residule[2,1]
        print(residule[1:2,])
        if(dresidule[1,1] < 0){
          pdf(file = options$pdf1)
          plot(residule[1:2,], type = "l", lwd = 2, xlab = "Optimal metagenes", ylab = "Residule", xaxt="n", main = "NMF statistics")
          axis(side = 1,at = 1:2
               , labels = c(num:(num+1)))       
          dev.off()
          rownames(Amplitude) <- rownames(geneExp)
          colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
          rownames(Pattern) <- colnames(geneExp)
          colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
          
        }else{
          # print(residule)
          Pattern <- t(geneNMF$H)
          Amplitude <- geneNMF$W
          for(cl in (num+2):max_cl){
            geneNMF <- bignmf(geneExp, r = cl, max.iteration = 200, stop.condition = 1e-05)
            residule[cl-num+1,1] <- sum(abs(sum(geneExp - geneNMF$W %*% geneNMF$H)))
            dresidule[cl-num,1] <- residule[cl-num,1] - residule[cl-num+1,1]
            print(dresidule[1:(cl-num+1),])
            if((dresidule[cl-num,1] <= dresidule[cl-num-1,1]) | (dresidule[cl-num,1] <= 0) ){
              break
            }else{
              Pattern <- t(geneNMF$H)
              Amplitude <- geneNMF$W
            }
          }
          pdf(file = options$pdf1)
          plot(residule[1:(cl-num+1),], type = "l", lwd = 2, xlab = "Optimal metagenes", ylab = "Residule", xaxt="n", main = "NMF statistics")
          axis(side = 1,at = 1:(cl-num+1)
               , labels = c(num:cl)
          )
          dev.off()
          rownames(Amplitude) <- rownames(geneExp)
          colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
          rownames(Pattern) <- colnames(geneExp)
          colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
          
        }
        
      }else{
        geneNMF <- bignmf(geneExp, r = as.numeric(options$input3), max.iteration = 200, stop.condition = 1e-05)
        Pattern <- t(geneNMF$H)
        Amplitude <- geneNMF$W
        rownames(Amplitude) <- rownames(geneExp)
        colnames(Amplitude) <- paste0("Metagene",1:ncol(Amplitude))
        rownames(Pattern) <- colnames(geneExp)
        colnames(Pattern) <- paste0("Metagene",1:ncol(Pattern))
      }
      write.table(Amplitude, file = options$output1, quote = FALSE, sep = "\t")
      write.table(Pattern, file = options$output2, quote = FALSE, sep = "\t")
    }
  }
}
MatrixFactorization()
#################################################################################################
