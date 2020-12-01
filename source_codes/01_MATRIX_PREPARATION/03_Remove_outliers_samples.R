




library("getopt")


options(stringAsfactors = F, useFancyQuotes = F)

args <- commandArgs(trailingOnly = TRUE)
option_specification <- matrix(c('input1', 'i1', 2, 'character',
                                 'input2', 'i2', 2, 'character',
                                 'input3', 'i3', 2, 'character',
                                 'output1', 'o1', 2, 'character'),
                               byrow=TRUE, ncol=4)

options = getopt(option_specification);


expressionQC <- function(geneExp = NULL, Rethreshold = 0.999, Lowthreshold=0.7){
  # geneExp
  # row -> gene
  # col -> sample
  
  log2transform=F
  if(max(geneExp) > 100){
    log2transform=T
    geneExp=log2(geneExp+1)
  }


  # average potential repeat samples
  sampleNum <- ncol(geneExp)
  sampleName = colnames(geneExp)
  #remove the repeat samples by R > 0.999
  repeatIdx <- list()
  for(sampleIdx in 1:(sampleNum-1)){
    # print(sampleIdx)
    if(length(which(unlist(repeatIdx) == sampleName[sampleIdx]))){
    }else{
      repeatInfor <- which(cor(geneExp[,sampleIdx], geneExp[,(sampleIdx+1):sampleNum]) > Rethreshold)
      if(length(repeatInfor) >0){
        cat(sampleName[sampleIdx],"\n")
        repeatIdx[[sampleName[sampleIdx]]] <- sampleName[sampleIdx + repeatInfor]
      }
    }
  }
  repeatIdx_back = repeatIdx
  repeatSampleNum <- length(repeatIdx)
  removeSampleIdx <- NULL
  averageSample <- NULL
  if(repeatSampleNum == 0){
  }else{
    repeatSampleIdx <- names(repeatIdx)
    for(Idx in 1:repeatSampleNum){
      averageSample <- cbind(averageSample, apply(as.matrix(geneExp[,c(repeatSampleIdx[Idx],repeatIdx[[Idx]])]), 1, mean))
      removeSampleIdx <- union(union(removeSampleIdx, repeatSampleIdx[Idx]),repeatIdx[[Idx]])
    } 
    geneExp <- geneExp[,-match(removeSampleIdx,sampleName)]
    new_sample_name <- colnames(geneExp)
    geneExp <- cbind(geneExp, averageSample)
    colnames(geneExp) <- c(new_sample_name,repeatSampleIdx)
  }
  
  
  # filter low-quality samples
  samplePCA <- prcomp(geneExp)
  GSPNum <- length(samplePCA$sdev)#principal component number
  PCEV <- sapply(1:GSPNum, function(PCIdx){samplePCA$sdev[PCIdx]^2/sum(samplePCA$sdev^2)})#expalined variance for each PC
  PCCEV <- sapply(1:GSPNum, function(PCIdx){sum(PCEV[1:PCIdx])})#cumulative explained variance
  SPCA <- rbind(PCEV * 100, PCCEV * 100)
  rownames(SPCA) <- c("Explained variance (%)", "Cumulative explained variance (%)")
  colnames(SPCA) <- paste0("PC", 1:length(PCEV))
 
  #calculate the correlation between the sample expression and the first PC
  sampleName <- colnames(geneExp)#undate sample name
  sampleCov <- abs(cor(geneExp, samplePCA$x[,1]))
  if(length(which(sampleCov < Lowthreshold))){
    lowQualitySample <- sampleName[which(sampleCov < Lowthreshold)]
    geneExp <- geneExp[,-match(lowQualitySample,sampleName)]
  }
  
  if(log2transform){
     geneExp = 2^geneExp - 1
  }
 
  
  return(geneExp)
  
}


rawExpMat <- as.matrix(read.table(file = options$input1, sep = "\t", quote = "", header	= T, row.names = 1))

expMat <- expressionQC(rawExpMat, options$input2, options$input3)
write.table(expMat, file = options$output1, sep = "\t", quote = F, col.names = T, row.names = T)

