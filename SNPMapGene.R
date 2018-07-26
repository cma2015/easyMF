
SNPMapGene <- function(GWASFile, AraGffFile, constructionFile, gene_poolFile){
  
  SNPMat <- as.matrix(read.table(GWASFile, sep = ",", header = T, quote = ""))
  SNPMat <- as.data.frame(SNPMat, stringsAsFactors = FALSE)
  for(col in 2:ncol(SNPMat)){
    SNPMat[,col] <- as.numeric(SNPMat[,col])
  }
  load(AraGffFile)
  gene2GeneSet <- read.table(file = constructionFile)
  gene2GeneSetRowname <- rownames(gene2GeneSet)
  gene2GeneSet <- apply(gene2GeneSet, 2, as.numeric)
  rownames(gene2GeneSet) <- gene2GeneSetRowname
  gene_pool <- NULL
  curSNP <- paste0(SNPMat[,1],"_", SNPMat[,2])
  for(j in 1:nrow(SNPMat)) {#for(j in 1:nrow(SNPMat))
    curCHR <- SNPMat[j,1]
    curPos <- as.numeric(SNPMat[j,2])
    curPvalue <- as.numeric(SNPMat[j,3])
    #curCHR <- paste0("Chr", curCHR)
    curGFF <- araGFF[which(araGFF[,1] == curCHR), ]
    geneStartVec <- curGFF[,4]
    names(geneStartVec) <- curGFF[,6]
    geneEndVec <- curGFF[,5]
    names(geneEndVec) <- curGFF[,6]
    
    curGFF[which(curGFF[,7] == "+"), 4] <- as.numeric(curGFF[which(curGFF[,7] == "+"), 4]) - 1000
    curGFF[which(curGFF[,7] == "-"), 5] <- as.numeric(curGFF[which(curGFF[,7] == "-"), 5]) + 1000
    
    curGene <- curGFF[which((as.numeric(curGFF[,4]) <= curPos) & (as.numeric(curGFF[,5]) >= curPos)), 6]
    if(length(curGene)==0){
      # next
      tempNum <- length(which(as.numeric(curGFF[,4]) <= curPos))
      dup <- curPos - as.numeric(curGFF[tempNum, 5])
      if(length(dup) == 0){
        dup <- 10000000000
      }
      if(tempNum < nrow(curGFF)){
        ddown <- as.numeric(curGFF[tempNum + 1, 4]) - curPos
      }else{
        ddown <- 10000000000
      }
      if(length(ddown) == 0){
        ddown <- 10000000000
      }
      if(dup == ddown){
        gene_pool <- rbind(gene_pool, c(curSNP[j], curCHR, curPos, curGFF[tempNum, 6], "nearest", curPvalue))
        gene_pool <- rbind(gene_pool, c(curSNP[j], curCHR, curPos, curGFF[tempNum + 1, 6], "nearest", curPvalue))
      }else if(dup > ddown){
        gene_pool <- rbind(gene_pool, c(curSNP[j], curCHR, curPos, curGFF[tempNum + 1, 6], "nearest", curPvalue))
      }else{
        gene_pool <- rbind(gene_pool, c(curSNP[j], curCHR, curPos, curGFF[tempNum, 6], "nearest", curPvalue))
      }
      
    }else{
      for(k in 1:length(curGene)){
        tmpGene <- curGene[k]
        curStart <- geneStartVec[tmpGene]
        curEnd <- geneEndVec[tmpGene]
        if((curPos >= curStart)&(curPos <= curEnd)){
          type <- "genic"
        }else{
          type <- "promoter"
        }
        gene_pool <- rbind(gene_pool, c(curSNP[j], curCHR, curPos, tmpGene, type, curPvalue))
      }
    }
    # cat(j, "\n")
  }
  colnames(gene_pool) <- c("flag", "Chr", "position", "gene", "type", "p")
  
  write.table(gene_pool, file = gene_poolFile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = T)
  
}

# SNPMapGene(GWASFile = "DTF1.csv", AraGffFile = "AraGff.RData", constructionFile = "Araconstruction.txt", gene_poolFile = "SNPMapGene.txt")
