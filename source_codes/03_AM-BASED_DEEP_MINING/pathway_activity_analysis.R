######################pathway activity analysis#########################################
# Cancer classification and pathway discovery using non-negative matrix factorization
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
# input2: input the amplitude matrix
# input3: input the pattern matrix
# input4: pathway annotation
# input5: input the cpu numbers to parallelly compute
# output1: actived pathways
# pdf1: actived pathways network
option_specification = matrix(c('input1', 'i1', 2, 'character',
                                'input2', 'i2', 2, 'character',
                                'input3', 'i3', 2, 'character',
                                'input4', 'i4', 2, 'character',
                                'input5', 'i5', 2, 'character',
                                'output1', 'o1', 2, 'character'
),
byrow=TRUE, ncol=4);
# Parse options
options(warn=-1)
options = getopt(option_specification);


PathwayActiveAnalysis <- function(){
  Amplitude <- as.matrix(read.table(options$input1, sep = "\t", quote = "", header = T, stringsAsFactors = F))
  pathwayAnno <- as.matrix(read.table(options$input2, sep = "\t", quote = "", header = T, stringsAsFactors = F))
  cpu <- as.numeric(options$input3)
  
  Discovery_code <- options$input5
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
  
  ###############  Amplitude  #######################
  # pathway analysis
  pathwayAnno <- pathwayAnno[!duplicated(pathwayAnno),]
  pathwayNumber <- table(pathwayAnno[,2])
  pathway_select <- pathwayAnno[which(is.na(match(pathwayAnno[,2], names(which( (pathwayNumber >= 5))))) == F),]
  pathway <- unique(pathway_select[,2])
  pathway_name <- pathwayAnno[,2:3]
  pathway_name <- pathway_name[!duplicated(pathway_name),]
  rownames(pathway_name) <- pathway_name[,1]
  # conserve only control pathway genes amplitude matrix
  if(length(unique(intersect(pathway_select[,1], rownames(Amplitude)))) <= 1){
    information = cbind("Warning!",
                        "The pathway genes are large not enough!",
                        "The evaluation is break!")
    write.table(x = information, file = options$output1, sep = "\t", row.names = F)
  }else{
    Amplitude_control <- Amplitude[unique(intersect(pathway_select[,1], rownames(Amplitude))),]
    geneName <- rownames(Amplitude_control)
    # control pathway gene-level activity statistics
    pathStatControl <- matrix(NA, length(pathway), ncol(Amplitude_control))
    rownames(pathStatControl) <- pathway
    colnames(pathStatControl) <- colnames(Amplitude_control)
    for(i in 1:nrow(pathStatControl)){
      for(j in 1:ncol(pathStatControl)){
        flag1 <- unique(intersect(geneName, pathway_select[which(pathway_select[,2] == pathway[i]),1]))
        if((length(flag1)>=2) & (length(setdiff(geneName,flag1)) >= 2)){
          dat1 <- Amplitude_control[flag1,j]
          dat2 <- Amplitude_control[setdiff(geneName,flag1),j]
          temp <- t.test(dat1, dat2) 
          if(temp$p.value == 0){
            temp$p.value = 1.0e-100
          }
          pathStatControl[i,j] <- qnorm(temp$p.value/2) * if(temp$statistic > 0) -1 else 1
        }
      }
    }
    
    # conserve only control pathway genes amplitude matrix
    pathway_gene_set_case <- unique(intersect(geneSet, pathway_select[,1]))
    if(length(unique(intersect(pathway_gene_set_case, rownames(Amplitude)))) <= 1){
      information = cbind("Warning!",
                          "The pathway genes are large not enough!",
                          "The evaluation is break!")
      write.table(x = information, file = options$output1, sep = "\t", row.names = F)
    }else{
      Amplitude_case <- Amplitude[unique(intersect(pathway_gene_set_case, rownames(Amplitude))),]
      geneName_case <- rownames(Amplitude_case)
      # case pathway gene-level activity statisticspathway_gene_set_case <- unique(intersect(pathway_gene_set, KEGGpathway[,1]))
      pathStatCase <- matrix(NA, length(pathway), ncol(Amplitude_case))
      rownames(pathStatCase) <- pathway
      colnames(pathStatCase) <- colnames(Amplitude_case)
      
      for(i in 1:nrow(pathStatCase)){
        for(j in 1:ncol(pathStatCase)){
          flag2 <- unique(intersect(geneName_case, pathway_select[which(pathway_select[,2] == pathway[i]),1]))
          if((length(flag2)>=2) & (length(setdiff(geneName_case,flag2)) >= 2)){
            dat1 <- Amplitude_case[flag2,j]
            dat2 <- Amplitude_case[setdiff(geneName_case,flag2),j]
            temp <- t.test(dat1, dat2)
            if(temp$p.value == 0){
              temp$p.value = 1.0e-100
            }
            pathStatCase[i,j] <-  qnorm(temp$p.value/2) * if(temp$statistic > 0) -1 else 1
          }
        }
      }
      # calculate pathway activity statistics
      pathStat <- matrix(NA, nrow(pathStatCase), ncol(pathStatCase) + 3)
      rownames(pathStat) <- rownames(pathStatCase)
      colnames(pathStat) <- c("significant","annotate","p_value", colnames(pathStatCase))
      for(i in 1:nrow(pathStat)){
        pathStat[i,1] <- length(unique(intersect(pathway_select[which(pathway_select[,2] == rownames(pathStat)[i]),1], pathway_gene_set_case)))
        pathStat[i,2] <- length(unique(pathway_select[which(pathway_select[,2] == rownames(pathStat)[i]),1]))
        dat1 <- pathStatCase[i,]
        dat2 <- pathStatControl[i,]
        if(length(which(is.na(dat1) == F))){
          pathStat[i,3] <- t.test(dat1, dat2, paired = T, alternative = "greater")$p.value
          pathStat[i,4:ncol(pathStat)] <- 100*abs((dat1 - dat2)/dat1)
        }
      }
      pathStat <- pathStat[order(pathStat[,3], decreasing = F),]
      pathStat <- cbind(pathway=rownames(pathStat),pathStat)
      
      # permutation test to calculate the FDR
      options(stringsAsFactors=F)
      options(scipen=999)
      suppressPackageStartupMessages(library(doParallel))
      suppressPackageStartupMessages(library(foreach))
      cl <-  makeCluster(cpu)
      registerDoParallel(cl)
      #do parallel computation
      permutation_result <- foreach(i = 1 : 1000) %dopar%{
        # for(i in 1:1000){
        set.seed(i)
        pathway_gene_set_rand <- geneName[sample(1:length(geneName), length(geneSet))]
        # conserve only control pathway genes amplitude matrix
        pathway_gene_set_rand_case <- unique(intersect(pathway_gene_set_rand, pathway_select[,1]))
        Amplitude_rand_case <- Amplitude[pathway_gene_set_rand_case,]
        # case pathway gene-level activity statisticspathway_gene_set_case <- unique(intersect(pathway_gene_set, KEGGpathway[,1]))
        pathStatRandCase <- matrix(NA, length(pathway), ncol(Amplitude_case))
        rownames(pathStatRandCase) <- pathway
        colnames(pathStatRandCase) <- colnames(Amplitude_rand_case)
        
        for(i in 1:nrow(pathStatRandCase)){
          for(j in 1:ncol(pathStatRandCase)){
            flag2 <- unique(intersect(rownames(Amplitude_rand_case), pathway_select[which(pathway_select[,2] == rownames(pathStatRandCase)[i]),1]))
            if((length(flag2)>=2) & (length(setdiff(rownames(Amplitude_rand_case),flag2)) >= 2)){
              dat1 <- Amplitude_rand_case[flag2,j]
              dat2 <- Amplitude_rand_case[setdiff(rownames(Amplitude_rand_case),flag2),j]
              temp <- t.test(dat1, dat2)
              if(temp$p.value == 0){
                temp$p.value = 1.0e-100
              }
              pathStatRandCase[i,j] <-  qnorm(temp$p.value/2) * if(temp$statistic > 0) -1 else 1
              
            }
          }
        }
        # calculate pathway activity statistics
        pathStatRand <- matrix(NA, nrow(pathStatRandCase), ncol(pathStatRandCase) + 3)
        rownames(pathStatRand) <- rownames(pathStatRandCase)
        colnames(pathStatRand) <- c("significant","annotate","p_value", colnames(pathStatRandCase))
        for(i in 1:nrow(pathStatRand)){
          pathStatRand[i,1] <- length(unique(intersect(pathway_select[which(pathway_select[,2] == rownames(pathStatRand)[i]),1], pathway_gene_set_case)))
          pathStatRand[i,2] <- length(unique(pathway_select[which(pathway_select[,2] == rownames(pathStatRand)[i]),1]))
          dat1 <- pathStatRandCase[i,]
          dat2 <- pathStatControl[i,]
          if(length(which(is.na(dat1) == F))){
            pathStatRand[i,3] <- t.test(dat1, dat2, paired = T, alternative = "greater")$p.value
            pathStatRand[i,4:ncol(pathStatRand)] <- 100*abs((dat1 - dat2)/dat1)
          }
        }
        # pathStatRand <- pathStatRand[which(is.na(pathStatRand[,3]) == F),]
        pathStatRand <- pathStatRand[order(pathStatRand[,3], decreasing = F),]
        pathStatRand <- cbind(pathway=rownames(pathStatRand),pathStatRand)
        return(list(pathStatRand = pathStatRand))
        # permutation_result[[i]] <- list(pathStatRand = pathStatRand)
      }
      names(permutation_result) <- c(1 : 1000)
      stopCluster(cl)
      
      pathStatMat <- pathStat[,4:ncol(pathStat)]
      class(pathStatMat) <- "numeric"
      pathStat_FDR <- matrix(0, nrow(pathStat), ncol(pathStat)-3)
      rownames(pathStat_FDR) <- rownames(pathStat)
      colnames(pathStat_FDR) <- c("FDR",colnames(pathStat)[5:ncol(pathStat)])
      for(i in 1:1000){
        tempMat <- permutation_result[[i]]$pathStatRand
        tempMat <- tempMat[rownames(pathStatMat),]
        tempMat <- tempMat[,4:ncol(tempMat)]
        class(tempMat) <- "numeric"
        temp <- matrix(0, nrow(pathStat), ncol(pathStat)-3)
        temp[which(tempMat >= pathStatMat)] <- 1
        pathStat_FDR <- pathStat_FDR+ temp
        # tempMat <- tempMat/pathStatMat
        # tempMat[which(tempMat >=1)] = 1
        # tempMat[which(tempMat < 1)] = 0
        # pathStat_FDR <- pathStat_FDR + tempMat
      }
      pathStat_FDR[,1] <- 1000 -pathStat_FDR[,1]
      pathStat_FDR <- pathStat_FDR/1000
      pathStat_FDR <- cbind(pathStat[rownames(pathStat_FDR),c(1,4)],FDR=pathStat_FDR[,1]
                            ,pathStat[rownames(pathStat_FDR),c(2,3)]
                            ,pathStat_FDR[,2:ncol(pathStat_FDR)])
      pathStat_FDR <- cbind(pathStat_FDR[,1:3], Term=pathway_name[pathStat_FDR[,1],2],pathStat_FDR[,4:ncol(pathStat_FDR)])
      pathStat_FDR <- pathStat_FDR[which(is.na(pathStat[,4]) == F),]
      if(length(which(is.na(pathStat[,4]) == F)) == 0){
        information = cbind("Warning!",
                            "There is none enough evaluated pathway!")
        write.table(x = information, file = options$output1, sep = "\t", row.names = F)
      }else if(length(which(is.na(pathStat[,4]) == F)) == 1){
        pathStat_FDR <- t(as.matrix(pathStat_FDR))
        pathStat_FDR[,2] <- round(as.numeric(pathStat_FDR[,2]),3)
        write.table(x = pathStat_FDR, file = options$output1, sep = "\t", row.names = F
                    , quote = F)
      }else{
        pathStat_FDR <- pathStat_FDR[order(as.numeric(pathStat_FDR[,2])),]
        pathStat_FDR[,2] <- round(as.numeric(pathStat_FDR[,2]),3)
        write.table(x = pathStat_FDR, file = options$output1, sep = "\t", row.names = F
                    , quote = F)
      }
      
    }
    
  }
  
  
 
}
PathwayActiveAnalysis()
#################################################################################################