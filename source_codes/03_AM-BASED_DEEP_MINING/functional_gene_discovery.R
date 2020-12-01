######################functional gene discovery#########################################
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
option_specification = matrix(c('input1', 'i1', 2, 'character',
                                'input2', 'i2', 2, 'character',
                                'input3', 'i3', 2, 'character',
                                'input4', 'i4', 2, 'character',
                                'pdf1', 'o1', 2, 'character',
                                'output1', 'o4', 2, 'character'
                      
),
byrow=TRUE, ncol=4);
# Parse options
options(warn=-1)
options = getopt(option_specification);

Ttest_MAP <- function(Amplitude = NULL, targetTrain = NULL
                      , targetTrain_nega_all = NULL
                      , targetTrain_nega = NULL){
  # Ttest
  targetTrain <- intersect(targetTrain, rownames(Amplitude))
  targetTrain_nega_all <- intersect(targetTrain_nega_all, rownames(Amplitude))
  targetTrain_nega <- intersect(targetTrain_nega, rownames(Amplitude))
  if(length(targetTrain) >= 2){
    targetPoolZscore <- matrix(NA ,ncol(Amplitude), 1)
    for (j in 1:ncol(Amplitude)) {
      dat1 <- Amplitude[targetTrain,j]
      dat2 <- Amplitude[targetTrain_nega_all,j]
      temp <- t.test(dat1, dat2)
      # targetPoolPCA[j,1] <- t.test(dat1, dat2)$statistic#T-test
      if(temp$p.value == 0){
        temp$p.value = 1.0e-100
      }
      targetPoolZscore[j,1] <- qnorm(temp$p.value/2) * if(temp$statistic > 0) -1 else 1#T-test
    }
    calCor <- apply(Amplitude[targetTrain_nega, ], 1, function(ll){
      temp <- cor.test(targetPoolZscore[,1], ll)
      if(is.na(temp$estimate)){
        return(NA)
      }else{
        return(qnorm(temp$p.value/2) * if(temp$estimate > 0) -1 else 1)
      }
    })
    #target test rank based on T-test
    tempT <- na.exclude(calCor)
    tempTOrder <- tempT[order(tempT, decreasing = TRUE)]
    TMat <- as.data.frame(cbind(names(tempTOrder), tempTOrder), stringsAsFactors = F)
    TMat[,1] <- as.character(TMat[,1])
    TMat[,2] <- as.numeric(TMat[,2])
    TMat[,2] <- (TMat[,2] - min(TMat[,2]))/(max(TMat[,2]) - min(TMat[,2]))
    rownames(TMat) <- TMat[,1]
  }else{
    TMat <- NULL
  }
  
  return(TMat)
}
FunctionalGeneDiscovery <- function(){
  
  Amplitude <- as.matrix(read.table(options$input1, header = T, stringsAsFactors = F))
  Discovery_code <- as.numeric(options$input2)
  cpu <- as.numeric(options$input3)
  if(Discovery_code == "1"){
    library(stringr)
    if(length(grep(pattern = "__cn__", x = options$input4))){  
      geneSet <-  str_replace_all(unlist(strsplit(options$input4,split="__cn__")), "-", " ")
    }
    if(length(grep(pattern = ",", x = options$input4))){  
      geneSet <-  str_replace_all(unlist(strsplit(options$input4,split=",")), "-", " ")
    }
  }else{
    geneSet <- read.table(options$input4, sep = c(",","\n"))
   if(nrow(geneSet) > 2){
     geneSet <- geneSet[,1]
   }
    
  }
    
  
    geneName <- rownames(Amplitude)
    genePoolNum <- nrow(Amplitude)
    GSPNum <- ncol(Amplitude)
    overlapTarget <- intersect(geneSet, geneName)
    targetTrain_nega_all <- setdiff(geneName, overlapTarget)
    
    TAMFResultLOOCV <- NULL
    # LOOCV
    options(stringsAsFactors=F)
    options(scipen=999)
    suppressPackageStartupMessages(library(doParallel))
    suppressPackageStartupMessages(library(foreach))
    cl <-  makeCluster(cpu)
    registerDoParallel(cl)
    #do parallel computation
    tempList <- foreach(i = 1 : length(overlapTarget)) %dopar%{
      
      ##################
      Ttest_MAP <- function(Amplitude = NULL, targetTrain = NULL
                            , targetTrain_nega_all = NULL
                            , targetTrain_nega = NULL){
        # Ttest
        targetTrain <- intersect(targetTrain, rownames(Amplitude))
        targetTrain_nega_all <- intersect(targetTrain_nega_all, rownames(Amplitude))
        targetTrain_nega <- intersect(targetTrain_nega, rownames(Amplitude))
        if(length(targetTrain) >= 2){
          targetPoolZscore <- matrix(NA ,ncol(Amplitude), 1)
          for (j in 1:ncol(Amplitude)) {
            dat1 <- Amplitude[targetTrain,j]
            dat2 <- Amplitude[targetTrain_nega_all,j]
            temp <- t.test(dat1, dat2)
            # targetPoolPCA[j,1] <- t.test(dat1, dat2)$statistic#T-test
            if(temp$p.value == 0){
              temp$p.value = 1.0e-100
            }
            targetPoolZscore[j,1] <- qnorm(temp$p.value/2) * if(temp$statistic > 0) -1 else 1#T-test
          }
          calCor <- apply(Amplitude[targetTrain_nega, ], 1, function(ll){
            temp <- cor.test(targetPoolZscore[,1], ll)
            if(is.na(temp$estimate)){
              return(NA)
            }else{
              return(qnorm(temp$p.value/2) * if(temp$estimate > 0) -1 else 1)
            }
          })
          #target test rank based on T-test
          tempT <- na.exclude(calCor)
          tempTOrder <- tempT[order(tempT, decreasing = TRUE)]
          TMat <- as.data.frame(cbind(names(tempTOrder), tempTOrder), stringsAsFactors = F)
          TMat[,1] <- as.character(TMat[,1])
          TMat[,2] <- as.numeric(TMat[,2])
          TMat[,2] <- (TMat[,2] - min(TMat[,2]))/(max(TMat[,2]) - min(TMat[,2]))
          rownames(TMat) <- TMat[,1]
        }else{
          TMat <- NULL
        }
        
        return(TMat)
      }
      ######################3
      targetTest <- overlapTarget[i]
      targetTrain <- overlapTarget[-i]
      targetTrain_nega <- setdiff(geneName, targetTrain)
      # easyMF
      TMat <- Ttest_MAP(Amplitude = Amplitude, targetTrain = targetTrain
                        , targetTrain_nega_all = targetTrain_nega_all
                        , targetTrain_nega = targetTrain_nega)
      if(length(which(TMat[,1] == targetTest)) > 0){
        Ttemp <- cbind(targetTest, TMat[which(TMat[,1] == targetTest),2], which(TMat[,1] == targetTest), "label", "easyMF")
        colnames(Ttemp) <- c("geneName", "score", "rank", "annotation", "Method")
      }else{
        Ttemp <- cbind(targetTest, NA, NA, "label", "easyMF")
        colnames(Ttemp) <- c("geneName",  "score", "rank", "annotation", "Method")
      }
      return(Ttemp = Ttemp)
    }
    names(tempList) <- c(1 : length(overlapTarget))
    stopCluster(cl)
    TAMFResultLOOCV <- do.call(rbind, tempList)
    
    AUSRCal <-  function(ratio){
      if(is.null(ratio) == F){
        rank <- 1:length(ratio)
        AUSR <- (0 + ratio[1])*(rank[1] - 0)/2
        
        for(index in 2:length(rank)){
          AUSR <- AUSR + (ratio[index-1] + ratio[index])*(rank[index] - rank[index-1])/2
        }
        return(AUSR/length(ratio))
      }else{0}
      
    }
    
    topk=1000
    AUSR <- matrix(,topk,1)
    # print(max(TAMFResultLOOCV[,4]))
    for(sortflag in 1:topk){    
      AUSR[sortflag, 1] <- length(which(as.numeric(TAMFResultLOOCV[,3]) <= sortflag))/nrow(TAMFResultLOOCV)
    }
    AUSRScore <- apply(AUSR, 2, function(cc){round(AUSRCal(cc),3)})
    # print(max(TAMFResultLOOCV[,4]))
    # pdf1
    pdf(file = options$pdf1)
    plot(c(1:topk), AUSR[, 1], type = "l", col = "red", lwd = 2, ylim = c(0,1), xlab = "Self-rank threshold", ylab = "Fraction")
    legend(x = topk*0.4, y = 0.5, legend = paste0("AUSR : ", AUSRScore), lty = c(1), col = 2
           ,border = NA, bty = "n", lwd = 2)
    dev.off()
    # print("pdf is done")
    TAMFResultUnlabel <- NULL
    ##############################
    # predict
    targetTrain <- overlapTarget
    targetTrain_nega <- setdiff(geneName, targetTrain)
    # Ttest
    TMat <- Ttest_MAP(Amplitude = Amplitude, targetTrain = targetTrain
                      , targetTrain_nega_all = targetTrain_nega_all
                      , targetTrain_nega = targetTrain_nega)
    Ttemp <- cbind(rownames(TMat), TMat[,2], 1:nrow(TMat),"Unlabel", "easyMF")
    colnames(Ttemp) <- c("geneName", "score", "rank","annotation", "Method")
    
    
    TAMFResultUnlabel <- rbind(TAMFResultUnlabel, Ttemp)
    TAMFResultUnlabel <- rbind(TAMFResultLOOCV,TAMFResultUnlabel)
    TAMFResultUnlabel <- TAMFResultUnlabel[,1:4]
    #TAMFResultUnlabel <- TAMFResultUnlabel[order(TAMFResultUnlabel[,2], decreasing = T),]
    # TAMFResultUnlabel <- cbind(TAMFResultUnlabel, geneAnnota[match(TAMFResultUnlabel[,1], geneAnnota[,1]), 2])
    # colnames(TAMFResultUnlabel)[7] <- "geneDescription"
    #TAMFResultUnlabel[,2] <- round(as.numeric(TAMFResultUnlabel[,2]), 3)
    write.table(TAMFResultUnlabel, file = options$output1, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE)
    # print("TAMFResultUnlabel is done")	
  
  
}
FunctionalGeneDiscovery()
#################################################################################################


