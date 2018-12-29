###################### gene expression matrix preparation#########################################
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
# input2: gene expression name or GSM or Key Words
# pdf1: Statistics analysis conducted on samples to filter low quality samples
# output1: high quality gene expression matrix for TAMF analysis
# output2: process information of the quality control
# output3: the GSM information fetched from the Key words to filter
option_specification = matrix(c('input1', 'i1', 2, 'integer',
                                'input2', 'i2', 2, 'character',
                                'input3', 'i3', 2, 'character',
                                'pdf1', 'o1', 2, 'character',
                                'output1', 'o2', 2, 'character',
                                'output2', 'o3', 2, 'character',
                                'output3', 'o4', 2, 'character'
),
byrow=TRUE, ncol=4);
# Parse options
options(warn=-1)
options = getopt(option_specification);
getGEOSuppFiles <- function(GEO, makeDirectory=FALSE, baseDir=getwd(),
                            pattern=NULL, verbose=TRUE) {
  
  getGEO <- function(geo, makeDir, baseDir, pattern, verbose) {
    getDirListing <- function(url, verbose) {
      if (verbose) message(url)
      a <- RCurl::getURL(url)
      if (grepl("<HTML", a, ignore.case=TRUE)) {
        doc <- XML::htmlParse(a)
        links <- XML::xpathSApply(doc, "//a/@href")
        XML::free(doc)
        b <- as.character(links)[-1]
        if (verbose) message("OK")
      } else {
        tmpcon <- textConnection(a, "r")
        b <- read.table(tmpcon)
        close(tmpcon)
      }
      return(b)
    }
    
    geotype <- toupper(substr(geo, 1, 3))
    storedir <- baseDir
    fileinfo <- c()
    stub <- gsub("\\d{1,3}$", "nnn", geo, perl=TRUE)
    if (geotype == "GSM") {
      url <- sprintf(paste0("https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/",
                            "%s/suppl/"), stub, geo)
    }
    if (geotype == "GSE") {
      url <- sprintf(paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/",
                            "%s/suppl/"), stub, geo)
    }
    if (geotype == "GPL") {
      url <- sprintf(paste0("https://ftp.ncbi.nlm.nih.gov/geo/platform/%s",
                            "/%s/suppl/"), stub, geo)
    }
    fnames <- try(getDirListing(url=url, verbose=verbose), silent=TRUE)
    if (inherits(fnames, "try-error")) {
      message("No supplemental files found.")
      message("Check URL manually if in doubt")
      message(url)
      return(NULL)
    }
    
    if (!is.null(pattern)) {
      idx <- if (length(pattern) > 1 ) {
        unlist(lapply(pattern, function(p){
          grep(p, fnames)}))} else idx = grep(pattern, fnames)
          fnames <- fnames[idx]
    }
    
    if (makeDirectory) {
      suppressWarnings(dir.create(storedir <- file.path(baseDir, geo)))
    }
    for (i in fnames) {
      download.file(file.path(url, i), destfile = file.path(storedir,i),
                    mode = "wb", quiet = ifelse(verbose, FALSE, TRUE),
                    method = getOption("download.file.method.GEOquery"))
      ## fileinfo[[file.path(storedir, i)]] <-
      ## file.info(file.path(storedir, i))
      fileinfo <- c(fileinfo, file.path(storedir, i))
    }
    return(fileinfo)
  }
  if (verbose) cat("*** Downloading files ...\n")
  x <- lapply(GEO, function(x) {getGEO(geo=x, makeDir=makeDirectory,
                                       baseDir=baseDir, pattern=pattern,
                                       verbose=verbose)})
  # return(unlist(x))
}

MatrixPreparation <- function(){
  library(stringr)
  # matrix source code: 1 -> default (our default ara gene expression matrix)
  #                     2 -> user upload
  #                     3 -> fetch key words
  #                     4 -> GSM information to filter
  cpu <- as.numeric(options$input3)
  if(options$input1 == 4){
    if(length(grep(pattern = "__cn__", x = options$input2))){  
      words <-  str_replace_all(unlist(strsplit(options$input2,split="__cn__")), "-", " ")
    }
    if(length(grep(pattern = ",", x = options$input2))){  
      words <-  unlist(strsplit(str_replace_all(unlist(strsplit(options$input2,split=",")), "-", " "), ","))
    }
    library(rentrez)
    if(length(words) == 1){
      term <- words[1]
    }else if(length(words) == 2){
      term <- paste0(words[1], " AND ", words[2])
    }else{
      term <- paste0("( ", words[2])
      for(i in 3:length(words)){
        term <- paste0(term, " OR ", words[i])
      }
      term <- paste0(words[1], " AND ", term, " )")
    }    
    
    for(iter in 1:10){
      if(!inherits(try(search <- entrez_search(db="gds", term = term) , TRUE), "try-error")){
        break
      }
    }
    for(iter in 1:10){
      if(!inherits(try(search <- entrez_search(db="gds", term = term, retmax = search$count) , TRUE), "try-error")){
        break
      }
    }
    GSM_summary <- NULL
    for(i in 1:ceiling(length(search$ids)/300)){      
      for(iter in 1:10){
        if(!inherits(try( taxize_summ <- entrez_summary(db="gds", id=search$ids[((i-1)*300+1):min((300*i),length(search$ids))]) , TRUE), "try-error")){
          break
        }
      }
      for(j in 1:length(taxize_summ)){
        if(length(taxize_summ[[j]]$samples)){
          GSM_summary <- rbind(GSM_summary, cbind(taxize_summ[[j]]$samples, paste0("GPL",taxize_summ[[j]]$gpl)))
        }
      }
    }
    colnames(GSM_summary)[3] <- "Platform"
    # t2 <- Sys.time()
    # print(t2 - t1)
    write.table(GSM_summary, file = options$output3, quote = F, sep = "\t", row.names = F)
    
  }else{
    if(options$input1 == 1){
      geneExp <- as.matrix(read.table(options$input2, header = T, stringsAsFactors = F))
      # print(options$input2)
    }else if(options$input1 == 2){
      geneExp <- as.matrix(read.table(options$input2, header = T, stringsAsFactors = F))
      # print(options$input2)
      
    }else if(options$input1 == 3){
      # print(as.character(options$input2))
      gsm_str <- as.matrix(read.table(file = options$input2, header = T, stringsAsFactors = F))
      
      library(GEOquery)
      library(affy)
      
      rownames(gsm_str) <- gsm_str[,2]
      if(!dir.exists(paste0(options$input2,"_"))) dir.create(paste0(options$input2,"_"))
      setwd(paste0(options$input2,"_"))
      GPLSet <- getGEO(unique(gsm_str[,3]),destdir =".")  ##soft
      colnames(Table(GPLSet))
      mapped <- Table(GPLSet)[,c(1,2)]
      mapped[,2] <- toupper(mapped[,2])
      rownames(mapped) <- mapped[,1]
      GEO_expMat <- list()
      for(file in 1:nrow(gsm_str)){
        print(paste0(file, " in ", nrow(gsm_str)))
        if(!inherits(try(getGEOSuppFiles(gsm_str[file,2], makeDirectory = F, verbose = FALSE), TRUE), "try-error")){
          cels <- list.files(pattern = "[gz]")
          celsrm <- cels[-grep(pattern = ".CEL.gz", x = cels)]
          a1 <- file.remove(celsrm)
          cels <- cels[grep(pattern = ".CEL.gz", x = cels)]
          # library(GEOquery)
          a1 <- sapply(cels, function(ll){gunzip(ll, overwrite = T)})
          cels <- unlist(strsplit(cels,".gz"))
          for(j in 1:length(cels)){
            if(!inherits(try(ReadAffy(filenames = cels[j]),TRUE), "try-error")){
              rawdata <- ReadAffy(filenames = cels[j]) #CEL
              exprSet <- exprs(mas5(rawdata, verbose = FALSE))# r
              exprSet <- cbind(rownames(exprSet), exprSet[,1])
              expression_mat <- as.data.frame(cbind(geneID = mapped[intersect(mapped[,1], exprSet[,1]),2], Sample = exprSet[intersect(mapped[,1], exprSet[,1]),2]), stringsAsFactors = F)
              colnames(expression_mat)[2] <- gsm_str[file,1]
              expression_mat <- expression_mat[which(expression_mat[,1] != "") ,]
              if(length(cels) > 1){
                GEO_expMat[[paste0(gsm_str[file,2],"_",j)]] <- expression_mat
              }else{
                GEO_expMat[[gsm_str[file,2]]] <- expression_mat
              }
            }else{
              GEO_infor <- getGEO(gsm_str[file,2],destdir =".")##GSE
              # print(GEO_infor@dataTable@columns)
              exprSet <- GEO_infor@dataTable@table# ID
              rownames(exprSet) <- exprSet[,1]
              expression_mat <- as.data.frame(cbind(geneID = mapped[intersect(mapped[,1], exprSet[,1]),2], Sample = exprSet[intersect(mapped[,1], exprSet[,1]),2]), stringsAsFactors = F)
              colnames(expression_mat)[2] <- gsm_str[file,1]
              expression_mat <- expression_mat[which(expression_mat[,1] != "") ,]
              # print(dim(expression_mat))
              if(length(cels) > 1){
                GEO_expMat[[paste0(gsm_str[file,2],"_",j)]] <- expression_mat
              }else{
                GEO_expMat[[gsm_str[file,2]]] <- expression_mat
              }
            }
          }
          
        }else{
          print(paste0("Warning: the ", gsm_str[file,2]," can not be download automatally"))
        }
      }
      
      geo_gene_all <- NULL
      for(file in 1:length(GEO_expMat)){
        geo_gene_all <- union(geo_gene_all, GEO_expMat[[file]][,1])
      }
      expression <- matrix(0, length(geo_gene_all), length(GEO_expMat))
      rownames(expression) <- geo_gene_all
      colnames(expression) <- gsm_str[unlist(sapply(names(GEO_expMat), function(ll){
        return(unlist(strsplit(ll,"_"))[1])
      })),1]
      
      tissue <- sort(unique(gsm_str[,1]))
      # tissue <- unique(geo_str[,1])
      for(flag in 1:length(GEO_expMat)){
        expression[GEO_expMat[[flag]][,1],flag] <- as.numeric(GEO_expMat[[flag]][,2])
      }
      expression <- expression[,order(colnames(expression))]
      geneExp <- expression
      
      system(paste0("rm -r ",paste0(options$input2,"_")))
    }
    
    
    # quality control of gene expression matrix
    # filter the NA genes, retain at least 10% non-zero genes
    na_num <- apply(geneExp, 1, function(ll){
      return(length(which(is.na(ll) == T)))
    })
    geneExp <- geneExp[which(na_num == 0),]
    
    geneExp_sd <- apply(geneExp, 1, sd)
    geneExp <- geneExp[which(geneExp_sd > 0 ),]
    
    gene_expression_max <- apply(geneExp, 2, max)
    if(length(which(gene_expression_max > 100))){
      geneExp[,which(gene_expression_max > 100)] <- log2(geneExp[,which(gene_expression_max > 100)] + 1)
    }
    
    for(flag in 1:100000000){
      gene_expression_max <- apply(geneExp, 2, max)
      if(max(gene_expression_max)/min(gene_expression_max) > 2){
        geneExp <- geneExp[,-which.max(gene_expression_max)]
      }else{
        break
      }
    }
    
    nonzero_num <- apply(geneExp, 1, function(ll){
      return(length(which(ll > 0)))
    })
    geneExp <- geneExp[which(nonzero_num >= ceiling(0.1 * ncol(geneExp))),]
    
    
    # ############################Principal component analysis##############################
    # #quality control of gene expression data
    # #average the duplicated samples
    sampleNum <- ncol(geneExp)
    sampleName <- colnames(geneExp)
    geneNum <- nrow(geneExp)
    geneName <- rownames(geneExp)
    
    #remove the repeat samples by R > 0.999
    repeatIdx <- list()
    for(sampleIdx in 1:(sampleNum-1)){
      # print(sampleIdx)
      repeatInfor <- which(cor(geneExp[,sampleIdx], geneExp[,(sampleIdx+1):sampleNum]) > 0.999)
      if(length(repeatInfor) >0){
        #cat(sampleName[sampleIdx],"\n")
        repeatIdx[[sampleName[sampleIdx]]] <- sampleName[sampleIdx + repeatInfor]
      }
    }
    
    repeatSampleNum <- length(repeatIdx)
    removeSampleIdx <- NULL
    averageSample <- NULL
    if(repeatSampleNum == 0){
      write.table(file = options$output2, "There is none duplicate samples", append = T, quote = F, row.names = F, col.names = F)
    }else{
      repeatSampleIdx <- names(repeatIdx)
      for(Idx in 1:repeatSampleNum){
        write.table(file = options$output2, paste0("The sample ", repeatSampleIdx[Idx], " have duplicate sample(s): ", repeatIdx[[Idx]],", we averaged them"), append = T, quote = F, row.names = F, col.names = F)
        averageSample <- cbind(averageSample, apply(as.matrix(geneExp[,c(repeatSampleIdx[Idx],repeatIdx[[Idx]])]), 1, mean))
        removeSampleIdx <- union(union(removeSampleIdx, repeatSampleIdx[Idx]),repeatIdx[[Idx]])
      }
      write.table(file = options$output2, paste0("Finally, the following samples are removed: "), append = T, quote = F, row.names = F, col.names = F)
      write.table(file = options$output2, removeSampleIdx, append = T, quote = F, row.names = F, col.names = F)
      geneExp <- geneExp[,-match(removeSampleIdx,sampleName)]
      new_sample_name <- colnames(geneExp)
      geneExp <- cbind(geneExp, averageSample)
      colnames(geneExp) <- c(new_sample_name,repeatSampleIdx)
    }
    #done remove the repeat samples by R > 0.999
    
    #quality control gene expression data
    #remove outlier samples by the correlation of each sample expression data with
    #the first principal component score by the principal component analysis on the
    #sample with R < 0.75
    samplePCA <- prcomp(geneExp)
    GSPNum <- length(samplePCA$sdev)#principal component number
    PCEV <- sapply(1:GSPNum, function(PCIdx){samplePCA$sdev[PCIdx]^2/sum(samplePCA$sdev^2)})#expalined variance for each PC
    PCCEV <- sapply(1:GSPNum, function(PCIdx){sum(PCEV[1:PCIdx])})#cumulative explained variance
    SPCA <- rbind(PCEV * 100, PCCEV * 100)
    rownames(SPCA) <- c("Explained variance (%)", "Cumulative explained variance (%)")
    colnames(SPCA) <- paste0("PC", 1:length(PCEV))
    #write.table(t(SPCA), file = paste0(workDir, "sample-based PCA.txt"))
    
    ########################PDF1
    pdf(file = options$pdf1)
    plot(PCEV, ylim = c(0,1), xlab = "Summary gene set patterns", ylab = "Variance Explaned (%)", main = "Statistics analysis", pch = 20)
    points(PCCEV, col = 2, pch = 20)
    legend(x = 1, y = 0.6, pch = c(20,20), col = c(1,2), bty = "n",
           legend = c("Explained variance", "Cumulative explained variance"))
    dev.off()
    
    #calculate the correlation between the sample expression and the first PC
    sampleName <- colnames(geneExp)#undate sample name
    sampleCov <- abs(cor(geneExp, samplePCA$x[,1]))
    if(length(which(sampleCov < 0.75))){
      lowQualitySample <- sampleName[which(sampleCov < 0.75)]
      geneExp <- geneExp[,-match(lowQualitySample,sampleName)]
    }
    
    # save high quality gene expression matrix
    write.table(geneExp, file = options$output1, quote = FALSE, sep = "\t")
  }
  
  
}
MatrixPreparation()
#################################################################################################
