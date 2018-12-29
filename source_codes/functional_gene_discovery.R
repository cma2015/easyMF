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

# input1: input the amplitude matrix
# input2: functional gene discovery code, 1->New candidates priorization, 2->Biological interpretation of GWAS result
# input3: input the cpu numbers to parallelly compute
# for New candidates priorization
# input4: input the interested gene set
# output1: new candidates of the input gene set
# pdf1: AUSR plot of the gene set through the LOOCV
# for Biological interpretation of GWAS result
# input5: the SNPMapGene GWAS result
# input13: input the cpu numbers to parallelly compute
# input6: functional gene probabilities construction code, 1->default (Arabidopsis thaliana genes (TAIR10))
#                                                          2->upload the construction of functional gene probabilities
#                                                          3->calculate from functional gene annotation
# for default (Arabidopsis thaliana genes (TAIR10)), the default data are used
# for upload the construction of functional gene probabilities
# input7: input the construction of functional gene probabilities
# for calculate from functional gene annotation
# input12: input the gold standard evidence code, default is: "IDA","IMP","IPI", "IGI", "IEP", "TAS"
# input9: functional gene annotation code, 1->upload, 2->fetch from ensemblplants
# for upload
# input10: input the functional gene annotation file by personal collected
# input14: input the gene description file by personal collected
# for fetch from ensemblplants
# input11: input the specie to fetch the functional gene annotation file from ensemblplants
# output5: the functional gene probabilities construction
# output6: golden standard functional terms information
# output2: enrichment analysis of the significant genes above the significant level
# output3: prioritize the potential related genes below the significant level
# output4: filter false related genes above the significant level
# pdf2: pheatmap of the significant genes above the significant level
# pdf3: manhattan plot of the GWAS result

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
                                'input12', 'i12', 2, 'character',
                                'input13', 'i13', 2, 'character',
                                'input14', 'i14', 2, 'character',
                                'pdf1', 'o1', 2, 'character',
                                'pdf2', 'o2', 2, 'character',
                                'pdf3', 'o3', 2, 'character',
                                'output1', 'o4', 2, 'character',
                                'output2', 'o5', 2, 'character',
                                'output3', 'o6', 2, 'character',
                                'output4', 'o7', 2, 'character',
                                'output5', 'o8', 2, 'character',
                                'output6', 'o9', 2, 'character'
                      
),
byrow=TRUE, ncol=4);
# Parse options
options(warn=-1)
options = getopt(option_specification);

BuildConstruction <- function(GSPLoadings, cpu, GTOGOTermsSelect){
  GSPLoadings <- GSPLoadings
  cpu <- cpu
  geneName <- rownames(GSPLoadings)
  genePoolNum <- nrow(GSPLoadings)
  GSPNum <- ncol(GSPLoadings)
  GTOGOTermsSelect <- GTOGOTermsSelect
  domainCode <- c("P", "C", "F")
  # T1 <- Sys.time()
  TAMFResultLOOCV <- NULL
  ###############################
  # leave one out experiment
  for(codeFlag in 1:3){
    GTOGOdomain <-  GTOGOTermsSelect[which(GTOGOTermsSelect[,3] == domainCode[codeFlag]), ]
    GTOGOdomainTerm <- unique(GTOGOdomain[,2])
    for(termFlag in 1:length(GTOGOdomainTerm)){#for(termFlag in 1:length(GTOGOdomainTerm)){
      # t1 <- Sys.time()
      termgeneName <- unique(GTOGOdomain[which(GTOGOdomain[,2] == GTOGOdomainTerm[termFlag]),1])
      overlapTarget <- intersect(termgeneName, geneName)
      
      # print(paste0(GTOGOdomainTerm[termFlag], " has ", length(overlapTarget), " genes"))
      
      targetNum <- length(overlapTarget)
      unLabelNum <- genePoolNum - targetNum
      # leave one out for label gene
      targetNum <- length(overlapTarget)
      unLabelNum <- genePoolNum - targetNum
      #target pool
      targetNumPool <- 1:targetNum
      targetPool <- matrix(NA, GSPNum, 1)
      colnames(targetPool) <- c("T_test_statistic")
      #target zscore
      targetPoolZscore <- matrix(NA ,GSPNum, 1)
      colnames(targetPoolZscore) <- c("T_test")
      #test correlation
      testCor <- matrix(NA, targetNum, 1)
      colnames(testCor) <- c("T_test")
      #unLabel correlation
      unLabel <- setdiff(geneName, overlapTarget)
      unLabelCor <- matrix(NA, unLabelNum, 1)
      colnames(unLabelCor) <- c("T_test")
      #################################
      #set parallel computation parameters
      options(stringsAsFactors=F)
      options(scipen=999)
      suppressPackageStartupMessages(library(doParallel))
      suppressPackageStartupMessages(library(foreach))
      cl <-  makeCluster(cpu)
      registerDoParallel(cl)
      #do parallel computation
      tempList <- foreach(i = targetNumPool[1] : targetNumPool[targetNum]) %dopar%{
        #do leave one out test
        targetTest <- overlapTarget[i]#take one out of the pathway members to do the leave-one-out
        targetTrain <- overlapTarget[-i]#the rest as the train dataSet
        ################Do test for enrichment score###############################################
        
        for (j in 1:GSPNum) {#calculate the statistics
          dat1 <- GSPLoadings[targetTrain,j]
          dat2 <- GSPLoadings[setdiff(geneName, targetTrain),j]
          targetPool[j,1] <- t.test(dat1, dat2)$statistic#T-test
        }
        
        #calculate the Zscore
        targetPoolZscore[,1] <- ((targetPool[,1] - mean(targetPool[,1]))/sd(targetPool[,1]))#T-test Zscore
        # targetPoolZscore[,3] <- targetPoolZscore[,1] + targetPoolZscore[,2]#T-test and AD-test Zscore
        #calculate the target test correlation between the test gene's factor loadings and the target pool zscore
        testCor[i,1] <- cor(targetPoolZscore[,1], GSPLoadings[targetTest,1:GSPNum])#T-test correlation
        #calculate the correlation between the unLabel genes' factor loadings and the target pool zscore
        unLabelCor[,1] <- cor(targetPoolZscore[,1], t(GSPLoadings[unLabel, 1:GSPNum]))#T-test correlation
        #
        ##################Done test for enrichment score###############################################
        
        ####################sort the test gene and the unLabel genes####################################
        #target test rank based on T-test
        tempT <- c(testCor[i,1], unLabelCor[,1])
        names(tempT) <- c(overlapTarget[i],unLabel)
        tempTOrder <- tempT[order(tempT, decreasing = TRUE)]
        TAMFMat <- as.data.frame(cbind(names(tempTOrder), tempTOrder), stringsAsFactors = F)
        TAMFMat[,1] <- as.character(TAMFMat[,1])
        TAMFMat[,2] <- as.numeric(TAMFMat[,2])
        TAMFMat[,2] <- (TAMFMat[,2] - min(TAMFMat[,2]))/(max(TAMFMat[,2]) - min(TAMFMat[,2]))
        TAMFtemp <- cbind(targetTest, TAMFMat[which(TAMFMat[,1] == targetTest),2], GTOGOdomainTerm[termFlag], which(TAMFMat[,1] == targetTest), domainCode[codeFlag], "label")
        colnames(TAMFtemp) <- c("geneName", "score", "GO", "rank", "domain", "annotation")
        return(list(TAMFtemp = TAMFtemp))
      }
      #names(tempList) <- targetNumPool
      stopCluster(cl)
      # t2 <- Sys.time()
      # print(t2 - t1)
      for(geneFlag in 1:length(tempList)){
        TAMFResultLOOCV <- rbind(TAMFResultLOOCV, tempList[[geneFlag]]$TAMFtemp)
      }
      
    }
    
  }
  # T2 <- Sys.time()
  # print(T2 - T1)
  TAMFResultLOOCV <- TAMFResultLOOCV[order(as.numeric(TAMFResultLOOCV[,2]), decreasing = T),]
  TAMFResultUnlabel <- NULL
  ##############################
  # unlabel gene prediction
  for(codeFlag in 1:3){
    #geneNames <- unique(TAMFResultLOOCV[which(TAMFResultLOOCV[,5] == domainCode[codeFlag]),1])
    GTOGOdomain <-  GTOGOTermsSelect[which(GTOGOTermsSelect[,3] == domainCode[codeFlag]), ]
    GTOGOdomainTerm <- unique(GTOGOdomain[,2])
    for(termFlag in 1:length(GTOGOdomainTerm)){#for(termFlag in 1:length(GTOGOdomainTerm)){
      # t1 <- Sys.time()
      termgeneName <- unique(GTOGOdomain[which(GTOGOdomain[,2] == GTOGOdomainTerm[termFlag]),1])
      overlapTarget <- intersect(termgeneName, geneName)
      # print(paste0(GTOGOdomainTerm[termFlag], " has ", length(overlapTarget), " genes"))
      targetNum <- length(overlapTarget)
      unLabelNum <- genePoolNum - targetNum
      # leave one out for label gene
      targetNum <- length(overlapTarget)
      unLabelNum <- genePoolNum - targetNum
      #target pool
      targetNumPool <- 1:targetNum
      targetPool <- matrix(NA, GSPNum, 1)
      colnames(targetPool) <- c("T_test_statistic")
      #target zscore
      targetPoolZscore <- matrix(NA ,GSPNum, 1)
      colnames(targetPoolZscore) <- c("T_test")
      #test correlation
      testCor <- matrix(NA, targetNum, 1)
      colnames(testCor) <- c("T_test")
      #unLabel correlation
      unLabel <- setdiff(geneName, overlapTarget)
      unLabelCor <- matrix(NA, unLabelNum, 1)
      colnames(unLabelCor) <- c("T_test")
      #################################
      targetTrain <- overlapTarget#the rest as the train dataSet
      ################Do test for enrichment score###############################################
      # t1 <- Sys.time()
      for (j in 1:GSPNum) {#calculate the statistics
        dat1 <- GSPLoadings[targetTrain,j]
        dat2 <- GSPLoadings[setdiff(geneName, targetTrain),j]
        targetPool[j,1] <- t.test(dat1, dat2)$statistic#T-test
      }
      # t2 <- Sys.time()
      # print(t2 - t1)
      #calculate the Zscore
      targetPoolZscore[,1] <- ((targetPool[,1] - mean(targetPool[,1]))/sd(targetPool[,1]))#T-test Zscore
      #calculate the correlation between the unLabel genes' factor loadings and the target pool zscore
      unLabelCor[,1] <- cor(targetPoolZscore[,1], t(GSPLoadings[unLabel, 1:GSPNum]))#T-test correlation
      #
      ##################Done test for enrichment score###############################################
      
      ####################sort the test gene and the unLabel genes####################################
      #target test rank based on T-test
      tempT <- unLabelCor[,1]
      names(tempT) <- unLabel
      tempTOrder <- tempT[order(tempT, decreasing = TRUE)]
      TAMFMat <- as.data.frame(cbind(names(tempTOrder), tempTOrder), stringsAsFactors = F)
      TAMFMat[,1] <- as.character(TAMFMat[,1])
      TAMFMat[,2] <- as.numeric(TAMFMat[,2])
      TAMFMat[,2] <- (TAMFMat[,2] - min(TAMFMat[,2]))/(max(TAMFMat[,2]) - min(TAMFMat[,2]))
      #unlabelgene <- setdiff(geneName, overlapTarget)
      unlabelgene <- TAMFMat[,1]
      unlabelgeneMatch <- match(unlabelgene, TAMFMat[,1])
      
      TAMFUnlabelMat <- matrix(,length(unlabelgene),6)
      colnames(TAMFUnlabelMat) <- c("geneName", "score", "GO", "rank", "domain", "annotation")
      TAMFUnlabelMat <- as.data.frame(TAMFUnlabelMat)
      TAMFUnlabelMat[,1] <- unlabelgene
      TAMFUnlabelMat[,2] <- TAMFMat[unlabelgene,2]
      TAMFUnlabelMat[,4] <- unlabelgeneMatch
      TAMFUnlabelMat[,3] <- GTOGOdomainTerm[termFlag]
      TAMFUnlabelMat[,5] <- domainCode[codeFlag]
      TAMFUnlabelMat[,6] <- "unlabel"
      
      TAMFResultUnlabel <- rbind(TAMFResultUnlabel, TAMFUnlabelMat)
      # t2 <- Sys.time()
      # print(t2 - t1)
      
    }
  }
  TAMFResultUnlabel <- TAMFResultUnlabel[order(as.numeric(TAMFResultUnlabel[,2]), decreasing = T),]
  TAMFResult <- rbind(TAMFResultLOOCV, TAMFResultUnlabel)
  geneSet <- unique(TAMFResult[,3])
  TAMFProResult <- as.data.frame(TAMFResult)
  TAMFProResult[,2] <- as.numeric(TAMFProResult[,2])
  gene2GeneSet <- by(TAMFProResult$score, list(TAMFProResult$geneName, TAMFProResult$GO), function(x){return(x)})
  gene2GeneSet <- as.matrix(unlist(gene2GeneSet))
  # t2 <- Sys.time()
  # print(t2 - t1)
  return(gene2GeneSet)
  
}

FunctionalGeneDiscovery <- function(){
  
  Amplitude <- as.matrix(read.table(options$input1, header = T, stringsAsFactors = F))
  Discovery_code <- as.numeric(options$input2)
  cpu <- as.numeric(options$input3)
  if(Discovery_code == 1){
    library(stringr)
    if(length(grep(pattern = "__cn__", x = options$input4))){  
      geneSet <-  str_replace_all(unlist(strsplit(options$input4,split="__cn__")), "-", " ")
    }
    if(length(grep(pattern = ",", x = options$input4))){  
      geneSet <-  str_replace_all(unlist(strsplit(options$input4,split=",")), "-", " ")
    }
    # print(geneSet)
    geneAnnota <- as.matrix(read.delim("/galaxy/tools/TAMF/test_data/default_ara_gene_description.txt", sep = "\t", header = T, stringsAsFactors = F))
    
    geneName <- rownames(Amplitude)
    genePoolNum <- nrow(Amplitude)
    GSPNum <- ncol(Amplitude)
    
    TAMFResultLOOCV <- NULL
    ###############################
    # leave one out experiment
    overlapTarget <- intersect(geneSet, geneName)
    targetNum <- length(overlapTarget)
    unLabelNum <- genePoolNum - targetNum
    # leave one out for label gene
    targetNum <- length(overlapTarget)
    unLabelNum <- genePoolNum - targetNum
    #target pool
    targetNumPool <- 1:targetNum
    targetPool <- matrix(NA, GSPNum, 1)
    colnames(targetPool) <- c("T_test_statistic")
    #target zscore
    targetPoolZscore <- matrix(NA ,GSPNum, 1)
    colnames(targetPoolZscore) <- c("T_test")
    #test correlation
    testCor <- matrix(NA, targetNum, 1)
    colnames(testCor) <- c("T_test")
    #unLabel correlation
    unLabel <- setdiff(geneName, overlapTarget)
    unLabelCor <- matrix(NA, unLabelNum, 1)
    colnames(unLabelCor) <- c("T_test")
    #################################
    #set parallel computation parameters
    options(stringsAsFactors=F)
    options(scipen=999)
    suppressPackageStartupMessages(library(doParallel))
    suppressPackageStartupMessages(library(foreach))
    cl <-  makeCluster(cpu)
    registerDoParallel(cl)
    #do parallel computation
    tempList <- foreach(i = targetNumPool[1] : targetNumPool[targetNum]) %dopar%{
      #do leave one out test
      targetTest <- overlapTarget[i]#take one out of the gene set members to do the leave-one-out
      targetTrain <- overlapTarget[-i]#the rest as the train dataSet
      ################Do test for enrichment score###############################################
      for (j in 1:GSPNum) {#calculate the statistics
        dat1 <- Amplitude[targetTrain,j]
        dat2 <- Amplitude[setdiff(geneName, targetTrain),j]
        targetPool[j,1] <- t.test(dat1, dat2)$statistic#T-test
      }
      #calculate the Zscore
      targetPoolZscore[,1] <- ((targetPool[,1] - mean(targetPool[,1]))/sd(targetPool[,1]))#T-test Zscore
      #calculate the target test correlation between the test gene's factor loadings and the target pool zscore
      testCor[i,1] <- cor(targetPoolZscore[,1], Amplitude[targetTest,1:GSPNum])#T-test correlation
      #calculate the correlation between the unLabel genes' factor loadings and the target pool zscore
      unLabelCor[,1] <- cor(targetPoolZscore[,1], t(Amplitude[unLabel, 1:GSPNum]))#T-test correlation
      #
      ##################Done test for enrichment score###############################################
      
      ####################sort the test gene and the unLabel genes####################################
      #target test rank based on T-test
      tempT <- c(testCor[i,1], unLabelCor[,1])
      names(tempT) <- c(overlapTarget[i],unLabel)
      tempTOrder <- tempT[order(tempT, decreasing = TRUE)]
      
      TAMFMat <- as.data.frame(cbind(names(tempTOrder), tempTOrder), stringsAsFactors = F)
      TAMFMat[,1] <- as.character(TAMFMat[,1])
      TAMFMat[,2] <- as.numeric(TAMFMat[,2])
      TAMFMat[,2] <- (TAMFMat[,2] - min(TAMFMat[,2]))/(max(TAMFMat[,2]) - min(TAMFMat[,2]))
      
      TAMFtemp <- cbind(targetTest, TAMFMat[which(TAMFMat[,1] == targetTest),2], "geneSet", which(TAMFMat[,1] == targetTest), "geneSet", "label")
      colnames(TAMFtemp) <- c("geneName", "score", "GO", "rank", "domain", "annotation")
      return(list(TAMFtemp = TAMFtemp))
    }
    # names(tempList) <- targetNumPool
    stopCluster(cl)
    for(geneFlag in 1:length(tempList)){
      TAMFResultLOOCV <- rbind(TAMFResultLOOCV, tempList[[geneFlag]]$TAMFtemp)
    }
    TAMFResultLOOCV <- as.data.frame(TAMFResultLOOCV, stringsAsFactors = FALSE)
    TAMFResultLOOCV[,2] <- as.numeric(TAMFResultLOOCV[,2])
    TAMFResultLOOCV[,4] <- as.numeric(TAMFResultLOOCV[,4])
    
    
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
    AUSR <- matrix(,max(TAMFResultLOOCV[,4]),1)
    # print(max(TAMFResultLOOCV[,4]))
    for(sortflag in 1:max(TAMFResultLOOCV[,4])){    
      AUSR[sortflag, 1] <- length(which(TAMFResultLOOCV[,4] <= sortflag))/nrow(TAMFResultLOOCV)
    }
    AUSRScore <- apply(AUSR, 2, function(cc){round(AUSRCal(cc),3)})
    # print(max(TAMFResultLOOCV[,4]))
    # pdf1
    pdf(file = options$pdf1)
    plot(c(1:max(TAMFResultLOOCV[,4])), AUSR[, 1], type = "l", col = "red", lwd = 2, ylim = c(0,1), xlab = "Self-rank threshold", ylab = "Fraction")
    legend(x = max(TAMFResultLOOCV[,4])*0.4, y = 0.5, legend = paste0("AUSR : ", AUSRScore), lty = c(1), col = 2
           ,border = NA, bty = "n", lwd = 2)
    dev.off()
    # print("pdf is done")
    TAMFResultUnlabel <- NULL
    ##############################
    # unlabel gene prediction
    
    overlapTarget <- intersect(geneSet, geneName)
    targetNum <- length(overlapTarget)
    unLabelNum <- genePoolNum - targetNum
    # leave one out for label gene
    targetNum <- length(overlapTarget)
    unLabelNum <- genePoolNum - targetNum
    
    #target pool
    targetNumPool <- 1:targetNum
    targetPool <- matrix(NA, GSPNum, 1)
    colnames(targetPool) <- c("T_test_statistic")
    #target zscore
    targetPoolZscore <- matrix(NA ,GSPNum, 1)
    colnames(targetPoolZscore) <- c("T_test")
    #test correlation
    testCor <- matrix(NA, targetNum, 1)
    colnames(testCor) <- c("T_test")
    #unLabel correlation
    unLabel <- setdiff(geneName, overlapTarget)
    unLabelCor <- matrix(NA, unLabelNum, 1)
    colnames(unLabelCor) <- c("T_test")
    #################################
    targetTrain <- overlapTarget#the rest as the train dataSet
    ################Do test for enrichment score###############################################
    # t1 <- Sys.time()
    for (j in 1:GSPNum) {#calculate the statistics
      dat1 <- Amplitude[targetTrain,j]
      dat2 <- Amplitude[setdiff(geneName, targetTrain),j]
      targetPool[j,1] <- t.test(dat1, dat2)$statistic#T-test
    }
    # t2 <- Sys.time()
    # print(t2 - t1)
    #calculate the Zscore
    targetPoolZscore[,1] <- ((targetPool[,1] - mean(targetPool[,1]))/sd(targetPool[,1]))#T-test Zscore
    #calculate the correlation between the unLabel genes' factor loadings and the target pool zscore
    unLabelCor[,1] <- cor(targetPoolZscore[,1], t(Amplitude[unLabel, 1:GSPNum]))#T-test correlation
    #
    ##################Done test for enrichment score###############################################
    
    ####################sort the test gene and the unLabel genes####################################
    #target test rank based on T-test
    tempT <- unLabelCor[,1]
    names(tempT) <- unLabel
    tempTOrder <- tempT[order(tempT, decreasing = TRUE)]
    
    TAMFMat <- as.data.frame(cbind(names(tempTOrder), tempTOrder), stringsAsFactors = F)
    TAMFMat[,1] <- as.character(TAMFMat[,1])
    TAMFMat[,2] <- as.numeric(TAMFMat[,2])
    TAMFMat[,2] <- (TAMFMat[,2] - min(TAMFMat[,2]))/(max(TAMFMat[,2]) - min(TAMFMat[,2]))
    
    #unlabelgene <- setdiff(geneName, overlapTarget)
    unlabelgene <- TAMFMat[,1]
    unlabelgeneMatch <- match(unlabelgene, TAMFMat[,1])
    
    TAMFUnlabelMat <- matrix(,length(unlabelgene),6)
    colnames(TAMFUnlabelMat) <- c("geneName", "score", "GO", "rank", "domain", "annotation")
    TAMFUnlabelMat <- as.data.frame(TAMFUnlabelMat)
    TAMFUnlabelMat[,1] <- unlabelgene
    TAMFUnlabelMat[,2] <- TAMFMat[unlabelgene,2]
    TAMFUnlabelMat[,4] <- unlabelgeneMatch
    TAMFUnlabelMat[,3] <- "geneSet"
    TAMFUnlabelMat[,5] <- "geneSet"
    TAMFUnlabelMat[,6] <- "unlabel"
    
    TAMFResultUnlabel <- rbind(TAMFResultUnlabel, TAMFUnlabelMat)
    #TAMFResultUnlabel <- TAMFResultUnlabel[order(TAMFResultUnlabel[,2], decreasing = T),]
    TAMFResultUnlabel <- cbind(TAMFResultUnlabel, geneAnnota[match(TAMFResultUnlabel[,1], geneAnnota[,1]), 2])
    colnames(TAMFResultUnlabel)[7] <- "geneDescription"
    #TAMFResultUnlabel[,2] <- round(as.numeric(TAMFResultUnlabel[,2]), 3)
    write.table(TAMFResultUnlabel[, c(4,1,2,7)], file = options$output1, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE)
    # print("TAMFResultUnlabel is done")	
    
    
    
  }else{
    
    gene_pool <- as.matrix(read.table(options$input5, header = T, stringsAsFactors = F))
    sigLevel <- 10^(-as.numeric(options$input13))
    construct_code <- as.numeric(options$input6)
   
    if(construct_code == 1){# default data
      
      gene2GeneSet <- as.matrix(read.table(file = "/galaxy/tools/TAMF/test_data/default_ara_gene2geneset.txt", header = T, stringsAsFactors = F, sep = "\t"))
      geneAnnota <- as.matrix(read.delim("/galaxy/tools/TAMF/test_data/default_ara_gene_description.txt", header = T, stringsAsFactors = F, sep = "\t"))
      
    }else if(construct_code == 2){# upload gene2geneset
       
      gene2GeneSet <- as.matrix(read.table(file = options$input7, header = T, stringsAsFactors = F))
      geneAnnota <- as.matrix(read.table(options$input14, header = T, stringsAsFactors = F, sep = "\t"))
       
      }else{# calculate gene2geneset
        
        annota_code <- options$input9
        
        if(annota_code == 1){# upload the functional gene annotation
          GTOGO <- as.matrix(read.table(options$input10, header = T, stringsAsFactors = F, sep = "\t"))
          geneAnnota <- as.matrix(read.table(options$input14, header = T, stringsAsFactors = F, sep = "\t"))
          
          
        }else{# fetch functional gene annotation from ensemblplants
          plants_ensembl <- as.matrix(read.table("/galaxy/tools/TAMF/test_data/plants_ensemble.txt", header = F, stringsAsFactors = F, sep = "\t"))          
          dataSet <- plants_ensembl[as.numeric(options$input11),1]
          require(biomaRt)
          require(topGO)
          mart <- useMart(biomart = "plants_mart", dataset = dataSet, host = "plants.ensembl.org")
          GTOGO <- getBM(attributes = c("ensembl_gene_id", "go_id", "namespace_1003", "go_linkage_type", "description"), mart = mart)
          geneAnnota <- GTOGO[!duplicated(GTOGO[,1]), c(1,5)]
          GTOGO <- GTOGO[GTOGO$go_id != '', c(1,2,3,4)]
          GTOGO <- GTOGO[which((GTOGO[,3] == "biological_process") | (GTOGO[,3] == "molecular_function") | (GTOGO[,3] == "cellular_component") ),]
          colnames(GTOGO) <- c("GeneID",	"GOID",	"Domain",	"Code")
          library(stringr)
          GTOGO[,3] <- str_replace_all(string = GTOGO[,3], pattern = "biological_process", replacement = "P")
          GTOGO[,3] <- str_replace_all(string = GTOGO[,3], pattern = "molecular_function", replacement = "F")
          GTOGO[,3] <- str_replace_all(string = GTOGO[,3], pattern = "cellular_component", replacement = "C")
        }
        
        # GTOGO <- as.matrix(read.table(options$input6, header = T, stringsAsFactors = F, sep = "\t"))
        verifyCode <- c("IBA","IC","IDA","IEA","IEP","IGI","IMP","IPI","ISA","ISM","ISS","NAS","ND","RCA","TAS")
        verifyCode <-  verifyCode[as.numeric(unlist(strsplit(options$input12,split=",")))]
        
        geneOverlap <- intersect(unique(GTOGO[,1]), rownames(Amplitude))
        GTOGOCoding <- GTOGO[which(is.na(match(x = GTOGO[,1], table = geneOverlap)) == F),]
        GTOGOVerify <- GTOGOCoding[which(is.na(match(x = GTOGOCoding[,4], table = verifyCode)) == F),]
        
        GTOGOTerms <- GTOGOVerify
        WholeTerms <- matrix(NA, length(unique(GTOGOTerms[,2])), 3)
        colnames(WholeTerms) <- c("GOTerms", "GeneNum", "OntoNum")
        WholeTerms[,1] <- unique(GTOGOTerms[,2])
        rownames(WholeTerms) <- unique(GTOGOTerms[,2])
        WholeTermsGene <- list()
        for(i in 1:length(WholeTerms[,1])){
          WholeTermsGene[[WholeTerms[i]]] <- unique(GTOGOTerms[which(GTOGOTerms[,2] == WholeTerms[i,1]), 1])
          WholeTerms[i,2] <- length(WholeTermsGene[[WholeTerms[i]]])
          WholeTerms[i,3] <- unique(GTOGOTerms[which(GTOGOTerms[,2] == WholeTerms[i,1]), 3])
        }
        WholeTermsSelect <- WholeTerms[which((as.numeric(WholeTerms[,2]) >= 5) & (as.numeric(WholeTerms[,2]) <= 500)),]
        
        
        numSelect <- NULL
        for(i in 1:length(WholeTermsSelect[,1])){
          numSelect <- c(numSelect, which(GTOGOTerms[,2] == WholeTermsSelect[i,1]))
        }
        GTOGOTermsSelect <- GTOGOTerms[numSelect,]
        
        write.table(file = options$output6,paste0("Terms in Biological Process : ", length(which(WholeTermsSelect[,3] == "P"))), append = T, quote = F, row.names = F, col.names = F)
        write.table(file = options$output6,paste0("Terms in Molecular Function : ",length(which(WholeTermsSelect[,3] == "F"))), append = T, quote = F, row.names = F, col.names = F)
        write.table(file = options$output6,paste0("Terms in Cellular Component : ",length(which(WholeTermsSelect[,3] == "C"))), append = T, quote = F, row.names = F, col.names = F)
        write.table(file = options$output6,paste0("Total GO terms : ",length(unique(GTOGOTermsSelect[,2]))), append = T, quote = F, row.names = F, col.names = F)
        write.table(file = options$output6,paste0("Genes in Biological Process : ",length(unique(GTOGOTermsSelect[which(GTOGOTermsSelect[,3] == "P"),1]))), append = T, quote = F, row.names = F, col.names = F)
        write.table(file = options$output6,paste0("Genes in Molecular Function : ",length(unique(GTOGOTermsSelect[which(GTOGOTermsSelect[,3] == "F"),1]))), append = T, quote = F, row.names = F, col.names = F)
        write.table(file = options$output6,paste0("Genes in Cellular Component : ",length(unique(GTOGOTermsSelect[which(GTOGOTermsSelect[,3] == "C"),1]))), append = T, quote = F, row.names = F, col.names = F)
        write.table(file = options$output6,paste0("Total Genes : ",length(unique(GTOGOTermsSelect[,1]))), append = T, quote = F, row.names = F, col.names = F)
        
        gene2GeneSet <- BuildConstruction(GSPLoadings = Amplitude, cpu = cpu, GTOGOTermsSelect = GTOGOTermsSelect)
        write.table(x = gene2GeneSet, file = options$output5, quote = F, sep = "\t", row.names = T, col.names = T)
    }
    
    ################3
    rownames(geneAnnota) <- geneAnnota[,1]
    ################
    library(stringr)
    colnames(gene2GeneSet) <- str_replace(colnames(gene2GeneSet),"GO.","GO:")
    library(qqman)
    resMat <- matrix(,length(gene_pool[,1]),4)
    resMat[,1] <- gene_pool[,1]
    resMat[,2] <- sapply(gene_pool[,2],function(cc){sub("Chr","",cc)})
    resMat[,3] <- gene_pool[,3]
    resMat[,4] <- gene_pool[,6]
    #########manhattan plot############
    manhaData <- data.frame(CHR = as.numeric(resMat[,2] ), BP = as.numeric(resMat[,3]), P = abs(as.numeric(resMat[,4])), SNP = resMat[,1])
    pdf(file = options$pdf3)
    # pdf(file = "/home/malab12/research/Galaxy/TAMF/test_data/GWAS.pdf")
    manhattan(x = manhaData, cex = 0.5, col = c("#AFD8E7", "#040E85"), genomewideline = FALSE, suggestiveline = -log10(sigLevel))
    dev.off()
    ########################
    # gene prioritize
    positive <- unique(gene_pool[which(as.numeric(gene_pool[,6]) < sigLevel), 4])
    overlap <- intersect(rownames(gene2GeneSet), positive)
    # length(overlap)
    library(pheatmap)
    library(RColorBrewer)
    R <- cor(t(gene2GeneSet[overlap,]))
    value_break <- seq(round(min(R),1), round(max(R),1), 0.01)
    col_grade <-  colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(value_break) -1)
    pdf(file = options$pdf2)
    B <- pheatmap(R,color = col_grade,breaks = value_break, cutree_row = 2, cutree_cols = 2, border_color = NA)
    dev.off()
    # d <- dist(R,method = "euclidean")
    # hc1=hclust(d,method = "complete");
    # pdf(file = "/home/malab12/research/Galaxy/TAMF/test_data/temp2.pdf")
    # plot(hc1)
    # dev.off()
    row_cluster=cutree(B$tree_row,k=2)
    D <- list()
    D[[1]] <- names(which(row_cluster == 1))
    D[[2]] <- names(which(row_cluster == 2))
    
    # D <- rect.hclust(hc1,k = 2, border = 3)
    if(length(D[[1]]) == 1){
      highQualitySNP <- D[[2]]
    }else if(length(D[[2]]) == 1){
      highQualitySNP <- D[[1]]
    }else if((length(which(R[D[[1]],D[[1]]] > 0.5)) - length(D[[1]]))/(length(R[D[[1]],D[[1]]]) - length(D[[1]])) > (length(which(R[D[[2]],D[[2]]] > 0.5)) - length(D[[2]]))/(length(R[D[[2]],D[[2]]]) - length(D[[2]]))){
      highQualitySNP <- D[[1]]
    }else{
      highQualitySNP <- D[[2]]
    }
    # print(highQualitySNP)
    
    outgene <- setdiff(highQualitySNP, intersect(highQualitySNP, geneAnnota[,1]))
    if(length(outgene) > 0){
      temp1 <- cbind(outgene, "")
      colnames(temp1) <- colnames(geneAnnota)
      geneAnnota <- rbind(geneAnnota, temp1)
      rownames(geneAnnota) <- geneAnnota[,1]
    }
    highQualitySNP <- cbind(highQualitySNP, geneAnnota[highQualitySNP,2])
    colnames(highQualitySNP) <- c("geneName", "description")
    write.table(highQualitySNP, file = options$output4, append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
    
    ###############################################
    # gene set enrichment
    GOAnnota <- as.matrix(read.delim("/galaxy/tools/TAMF/test_data/default_GO_annota.txt", header = F, stringsAsFactors = F, sep = "\t"))
    rownames(GOAnnota) <- GOAnnota[,1]
    GOAnnotaUse <- GOAnnota[colnames(gene2GeneSet),]
    geneNameOverlap <- setdiff(intersect(rownames(gene2GeneSet), unique(gene_pool[,4])), overlap)
    options(stringsAsFactors=F)
    options(scipen=999)
    suppressPackageStartupMessages(library(doParallel))
    suppressPackageStartupMessages(library(foreach))
    cl <-  makeCluster(cpu)
    registerDoParallel(cl)
    #do parallel computation
    tempList <- foreach(j = 1 : 20) %dopar%{
    # for(j in 1:20){
      # t1 <- Sys.time()
      # print(j)
      seedflag <- (j-1)*1000 + 1
      FDRMat <- NULL
      FDRMat <- rbind(FDRMat, apply(gene2GeneSet[overlap, ], 2, sum))
      
      # options(stringsAsFactors=F)
      # options(scipen=999)
      # suppressPackageStartupMessages(library(doParallel))
      # suppressPackageStartupMessages(library(foreach))
      # cl <-  makeCluster(cpu)
      # registerDoParallel(cl)
      # #do parallel computation
      # tempList <- foreach(i = 1 : 1000) %dopar%{
      #   tempOverlap <- sample(geneNameOverlap, length(overlap))
      #   return(apply(gene2GeneSet[tempOverlap, ], 2, sum))
      # }
      # names(tempList) <- targetNumPool
      # stopCluster(cl)
      # FDRMat <- rbind(FDRMat, do.call(rbind, tempList))
      for(i in 1:1000){# for(i in 1:1000)
        set.seed(seedflag)
        tempOverlap <- sample(geneNameOverlap, length(overlap))
        FDRMat <- rbind(FDRMat, apply(gene2GeneSet[tempOverlap, ], 2, sum))
        seedflag<- seedflag+1
      }
      tmp <- sapply(1:ncol(FDRMat), function(aa){return((FDRMat[1,aa]-mean(FDRMat[2:1001,aa]))/sd(FDRMat[2:1001,aa]))})
      return(tmp)
      # t2 <- Sys.time()
      # print(t2 - t1)
    }
    stopCluster(cl)
    FDR <- do.call(rbind, tempList)
    FDR2 <- pnorm(FDR, lower.tail = FALSE)
    adjustp <- apply(FDR2, 2, function(aa){return(length(aa) * aa/rank(aa))})

    # length(which(adjustp[1,] < 0.001))
    geneSetEn <- cbind(GOAnnotaUse[adjustp[1,] < 0.01,1], adjustp[1,which(adjustp[1,] < 0.01)], GOAnnotaUse[adjustp[1,] < 0.01,2:3]) 
    colnames(geneSetEn) <- c("GO", "fdr", "function", "domain")
    geneSetEn[,2] <- format(as.numeric(geneSetEn[,2]),scientific = T, digits = 4)
    write.table(geneSetEn, file = options$output2, append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
    ##################################################
    # gene prioritize outside 
    genePrize <- highQualitySNP
    FDRout <- NULL
    FDRout <- rbind(FDRout, apply(cor(t(gene2GeneSet[overlap,]), t(gene2GeneSet[geneNameOverlap,])), 2, max))
    # t1 <- Sys.time()
    # for(i in 1:50){# for(i in 1:1000)
    #   print(i)
    #   set.seed(i)
    #   tempOverlap <- sample(geneNameOverlap, length(overlap) + 1)
    #   tempMat <- cor(t(gene2GeneSet[tempOverlap,]), t(gene2GeneSet[geneNameOverlap,]))
    #   FDRout <- rbind(FDRout, apply(tempMat, 2, function(ii){
    #     return(max(ii[which(ii < 0.999999999)]))
    #   }))
    # }
    # adjustpout <- apply(FDRout, 2, function(aa){return(length(which(aa[2:11] > aa[1]))/10)})
    # sig <- which(adjustpout < 1e-3)
    # t2 <- Sys.time()
    # print(t2 - t1)
    # t1 <- Sys.time()
    options(stringsAsFactors=F)
    options(scipen=999)
    suppressPackageStartupMessages(library(doParallel))
    suppressPackageStartupMessages(library(foreach))
    cl <-  makeCluster(cpu)
    registerDoParallel(cl)
    #do parallel computation
    tempList <- foreach(i = 1 : 1000) %dopar%{
      # for(i in 1:1000){# for(i in 1:1000)
      # print(i)
      set.seed(i)
      tempOverlap <- sample(geneNameOverlap, length(overlap) + 1)
      tempMat <- cor(t(gene2GeneSet[tempOverlap,]), t(gene2GeneSet[geneNameOverlap,]))
      return(apply(tempMat, 2, function(ii){
        return(max(ii[which(ii < 0.999999999)]))
      }))
    }
    names(tempList) <- 1 : 1000
    stopCluster(cl)
    FDRout <- rbind(FDRout, do.call(rbind, tempList))
    adjustpout <- apply(FDRout, 2, function(aa){return(length(which(aa[2:1001] > aa[1]))/1000)})
    sig <- which(adjustpout < 1e-3)
    # t2 <- Sys.time()
    # print(t2 - t1)
    sigGeneOut <- names(FDRout[1,sig])[which(FDRout[1,sig] > 0.5)]
    pvalue <- matrix(, length(sigGeneOut), 1)
    rownames(pvalue) <- sigGeneOut
    for(i in 1:length(sigGeneOut)){
      pvalue[i,1] <- min(as.numeric(gene_pool[which(gene_pool[,4] == sigGeneOut[i]),6]))
    }
    outgene <- setdiff(sigGeneOut, intersect(sigGeneOut, geneAnnota[,1]))
    if(length(outgene) > 0){
      temp1 <- cbind(outgene, "")
      colnames(temp1) <- colnames(geneAnnota)
      geneAnnota <- rbind(geneAnnota, temp1)
      rownames(geneAnnota) <- geneAnnota[,1]
    }
    sigGeneOutAnno <- cbind(geneAnnota[sigGeneOut,1],-log10(pvalue), adjustpout[sigGeneOut], FDRout[1,sigGeneOut], geneAnnota[sigGeneOut,2])
    colnames(sigGeneOutAnno) <- c("geneName", "-log10(pvalue)", "fdr", "PCC", "description")
    sigGeneOutAnno[,2] <- round(as.numeric(sigGeneOutAnno[,2]), 3)
    sigGeneOutAnno[,4] <- round(as.numeric(sigGeneOutAnno[,4]), 3)
    write.table(sigGeneOutAnno, file = options$output3, append = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
  }
  
  
  
  
}
FunctionalGeneDiscovery()
#################################################################################################


