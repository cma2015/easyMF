###################### knowledge data preparation#########################################
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

# input1: code of knowledge data from each tool: 1-> Functional gene discovery,
#                                                4->Single cell analysis, 
#                                                5->Spatial-course analysis,
#                                                6->Time-course analysis,
#                                                7->No
# input2: code of fetch source
# input3: description of the matrix source
# pdf1: Statistics analysis conducted on samples to filter low quality samples
# output1: high quality gene expression matrix for TAMF analysis
# output2: process information of the quality control
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
                                'output1', 'o1', 2, 'character',
                                'output2', 'o2', 2, 'character',
                                'output3', 'o3', 2, 'character',
                                'output4', 'o4', 2, 'character',
                                'output5', 'o5', 2, 'character'
                              
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

KnowledgePreparation <- function(){
  
  code <- as.numeric(options$input1)
  if(code == 1){
    Amplitude <- as.matrix(read.table(options$input2, header = T, stringsAsFactors = F))
    cpu <- as.numeric(options$input3)
    annota_code <- options$input4
    if(annota_code == 1){# upload the functional gene annotation
      GTOGO <- as.matrix(read.table(options$input5, header = T, stringsAsFactors = F, sep = "\t"))
      geneAnnota <- as.matrix(read.table(options$input6, header = T, stringsAsFactors = F, sep = "\t"))
    }else{# fetch functional gene annotation from ensemblplants
      plants_ensembl <- as.matrix(read.table("/galaxy/tools/TAMF/test_data/plants_ensemble.txt", header = F, stringsAsFactors = F, sep = "\t"))          
      dataSet <- plants_ensembl[as.numeric(options$input7),1]
      require(biomaRt)
      mart <- useMart(biomart = "plants_mart", dataset = dataSet, host = "plants.ensembl.org")
      GTOGO <- getBM(attributes = c("ensembl_gene_id", "go_id", "namespace_1003", "go_linkage_type", "description"), mart = mart)
      geneAnnota <- GTOGO[!duplicated(GTOGO[,1]), c(1,5)]
      write.table(geneAnnota, file = options$output3, quote = F, sep = "\t", row.names = F)
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
    verifyCode <-  verifyCode[as.numeric(unlist(strsplit(options$input8,split=",")))]
    
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
    
    write.table(file = options$output2,paste0("Terms in Biological Process : ", length(which(WholeTermsSelect[,3] == "P"))), append = T, quote = F, row.names = F, col.names = F)
    write.table(file = options$output2,paste0("Terms in Molecular Function : ",length(which(WholeTermsSelect[,3] == "F"))), append = T, quote = F, row.names = F, col.names = F)
    write.table(file = options$output2,paste0("Terms in Cellular Component : ",length(which(WholeTermsSelect[,3] == "C"))), append = T, quote = F, row.names = F, col.names = F)
    write.table(file = options$output2,paste0("Total GO terms : ",length(unique(GTOGOTermsSelect[,2]))), append = T, quote = F, row.names = F, col.names = F)
    write.table(file = options$output2,paste0("Genes in Biological Process : ",length(unique(GTOGOTermsSelect[which(GTOGOTermsSelect[,3] == "P"),1]))), append = T, quote = F, row.names = F, col.names = F)
    write.table(file = options$output2,paste0("Genes in Molecular Function : ",length(unique(GTOGOTermsSelect[which(GTOGOTermsSelect[,3] == "F"),1]))), append = T, quote = F, row.names = F, col.names = F)
    write.table(file = options$output2,paste0("Genes in Cellular Component : ",length(unique(GTOGOTermsSelect[which(GTOGOTermsSelect[,3] == "C"),1]))), append = T, quote = F, row.names = F, col.names = F)
    write.table(file = options$output2,paste0("Total Genes : ",length(unique(GTOGOTermsSelect[,1]))), append = T, quote = F, row.names = F, col.names = F)
    
    gene2GeneSet <- BuildConstruction(GSPLoadings = Amplitude, cpu = cpu, GTOGOTermsSelect = GTOGOTermsSelect)
    write.table(x = gene2GeneSet, file = options$output1, quote = F, sep = "\t", row.names = T, col.names = T)
  }
  if(code == 4){
    cpu <- as.numeric(options$input3)
    gene_expression <- as.matrix(read.table(options$input9, header = T, stringsAsFactors = F))
    aa <- colnames(gene_expression)
    bb <- unlist(sapply(1:length(aa), function(ll){
      num <- unlist(gregexpr("[.]", aa[ll]))
      if(num[1] > 0){
        return(substr(aa[ll], 1, num[length(num)] - 1))
      }
      return(aa[ll])
    }))
    colnames(gene_expression) <- bb
    source("/galaxy/tools/TAMF/test_data/spec.R")
    specscore <- getAllSpec(data = gene_expression, header = colnames(gene_expression), cpu = cpu)
    write.table(specscore[[1]], file = options$output4, quote = F, sep = "\t")
  }
  if(code == 5){
    plants_ensembl <- as.matrix(read.table("/galaxy/tools/TAMF/test_data/plants_ensemble.txt", header = F, stringsAsFactors = F, sep = "\t"))          
    dataSet <- plants_ensembl[as.numeric(options$input10),1]
    require(biomaRt)
    mart <- useMart(biomart = "plants_mart", dataset = dataSet, host = "plants.ensembl.org")
    GTOGO <- getBM(attributes = c("ensembl_gene_id", "go_id", "namespace_1003", "go_linkage_type", "description"), mart = mart)
    geneAnnota <- GTOGO[!duplicated(GTOGO[,1]), c(1,5)]
    write.table(geneAnnota, file = options$output3, quote = F, sep = "\t", row.names = F)
    prePareGO <- prePareGO(dataset = dataSet)
    save(prePareGO, file = options$output5)
  }
  if(code == 6){
    plants_ensembl <- as.matrix(read.table("/galaxy/tools/TAMF/test_data/plants_ensemble.txt", header = F, stringsAsFactors = F, sep = "\t"))          
    dataSet <- plants_ensembl[as.numeric(options$input11),1]
    require(biomaRt)
    mart <- useMart(biomart = "plants_mart", dataset = dataSet, host = "plants.ensembl.org")
    GTOGO <- getBM(attributes = c("ensembl_gene_id", "go_id", "namespace_1003", "go_linkage_type", "description"), mart = mart)
    geneAnnota <- GTOGO[!duplicated(GTOGO[,1]), c(1,5)]
    write.table(geneAnnota, file = options$output3, quote = F, sep = "\t", row.names = F)
  }

  
}
KnowledgePreparation()
#################################################################################################
