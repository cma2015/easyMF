######################clustering analysis analysis#########################################
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

# input1: input the amplitude matrix
# input2: input the automatic clustering methods
# output1: actived pathways
# pdf1: actived pathways network
option_specification = matrix(c('input1', 'i1', 2, 'character',
                                'input2', 'i2', 2, 'character',
                                'input3', 'i3', 2, 'character',
                                'output1', 'o1', 2, 'character',
                                'pdf1', 'o2', 2, 'character'
),
byrow=TRUE, ncol=4);
# Parse options
options(warn=-1)
options = getopt(option_specification);

MFCluster <- function(Pattern, methods = "mclust"){
  if(ncol(Pattern) > 30){
    num <- 2
  }else{
    num <- ncol(Pattern)
  }
  if(methods == "mclust"){
    library(mclust)
    d_clust <- Mclust(Pattern[,1:num], G=1:min(15,ncol(Pattern)))
    m.best <- dim(d_clust$z)[2]
    cat("model-based optimal number of clusters:", m.best, "\n")
    plot(d_clust, dimens = c(2,1), what = "classification", xlab = "Pattern1", ylab = "Pattern2", main = paste0(m.best, " clusters"))
    title(main = paste0("mclust :",m.best, " clusters"))
    return(d_clust$classification)
  }
  
  if(methods == "apcluster"){
    library(apcluster)
    d.apclus <- apcluster(negDistMat(r=2), Pattern[,1:num])
    cat("affinity propogation optimal number of clusters:", length(d.apclus@clusters), "\n")
    # 4
    # pheatmap(d.apclus)
    plot(d.apclus, Pattern[,1:2], xlab = "Pattern1", ylab = "Pattern2", 
         main = paste0("apcluster :",length(d.apclus@clusters), " clusters"))
    cluster <- sapply(1:length(d.apclus@clusters),function(i) {
      d.apclus@clusters[[i]][1:length(d.apclus@clusters[[i]])] <- i
      return(d.apclus@clusters[[i]])})
    return(unlist(cluster))
  }
  
  if(methods == "SSE"){
    mydata <- Pattern[,1:num]
    wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
    for (i in 2:ncol(Pattern)) wss[i] <- sum(kmeans(mydata,
                                                    centers=i)$withinss)
    oss <- order(wss, decreasing = T)
    best.cl = which(oss != 1:length(oss))[1]
    if(is.na(best.cl)) best.cl = length(oss)
    cl <- kmeans(mydata, centers = best.cl)
    plot(mydata[which(cl$cluster == 1),1], mydata[which(cl$cluster == 1),2], xlab = "Pattern1",
         ylab = "Pattern2", main = paste0("SSE :",best.cl, " clusters"), col = 1, pch = 1,
         xlim = range(mydata[,1]), ylim = range(mydata[,2]))
    for(i in 2:best.cl){
      points(mydata[which(cl$cluster == i),1], mydata[which(cl$cluster == i),2], col = i, pch = i)
    }
    return(cl$cluster)
    # fviz_cluster(cl, data = dataset)
  }
  
  if(methods == "fpc"){
    library(fpc)
    library(cluster)
    pamk.best <- pamk(Pattern[,1:num])
    cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
    aa = pam(Pattern[,1:2], pamk.best$nc)
    clusplot(pam(Pattern[,1:2], pamk.best$nc), xlab = "Pattern1",
             ylab = "Pattern2", main = paste0("fpc :",pamk.best$nc, " clusters"), sub = "")
    return(aa$clustering)
  }
  
  if(methods == "vegan"){
    
    require(vegan)
    fit <- cascadeKM(Pattern[,1:num], 1, min(20,ncol(Pattern)), iter = 1000)
    calinski.best <- as.numeric(which.max(fit$results[2,]))
    cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")
    # 5 clusters!
    plot(Pattern[which(fit$partition[,calinski.best] == 1),1], Pattern[which(fit$partition[,calinski.best] == 1),2], xlab = "Pattern1",
         ylab = "Pattern2", main = paste0("vegan :",calinski.best, " clusters"), col = 1, pch = 1,
         xlim = range(Pattern[,1]), ylim = range(Pattern[,2]))
    for(i in 2:calinski.best){
      points(Pattern[which(fit$partition[,calinski.best] == i),1], Pattern[which(fit$partition[,calinski.best] == i),2], col = i, pch = i)
    }
    return(fit$partition[,calinski.best])
    
  }
  
  if(methods == "gap"){
    library(cluster)
    aa = clusGap(Pattern[,1:num], kmeans, min(20,ncol(Pattern)), B = 100, verbose = F)
    wss = aa$Tab[,"gap"]
    oss <- order(wss, decreasing = F)
    best.cl = which(oss != 1:length(oss))[1]
    if(is.na(best.cl)) best.cl = length(oss)
    mydata <- Pattern[,1:num]
    cl <- kmeans(mydata, centers = best.cl)
    plot(mydata[which(cl$cluster == 1),1], mydata[which(cl$cluster == 1),2], xlab = "Pattern1",
         ylab = "Pattern2", main = paste0("gap :",best.cl, " clusters"), col = 1, pch = 1,
         xlim = range(mydata[,1]), ylim = range(mydata[,2]))
    for(i in 2:best.cl){
      points(mydata[which(cl$cluster == i),1], mydata[which(cl$cluster == i),2], col = i, pch = i)
    }
    return(cl$cluster)
  }
}

ClusterAnalysis <- function(){
  
  Pattern <- as.matrix(read.table(options$input1, header = T, sep = "\t", stringsAsFactors = F))
  methods <- c("mclust", "apcluster", "SSE", "fpc", "vegan", "gap")
  methods <-  methods[as.numeric(unlist(strsplit(options$input2,split=",")))]
  cluster_mat <- matrix(, nrow(Pattern), length(methods))
  rownames(cluster_mat) <- rownames(Pattern)
  colnames(cluster_mat) <- methods
  
  pdf(options$pdf1)
  par(mfrow = c(ceiling(sqrt(length(methods))), ceiling(sqrt(length(methods)))))
  for(i in 1:length(methods)){
    cluster <- MFCluster(Pattern = Pattern, methods = methods[i])
    cluster_mat[names(cluster),i] <- cluster
  }
  dev.off()
  write.table(cluster_mat, file = options$output1, quote = F, sep = "\t")
  
}
ClusterAnalysis()
#################################################################################################

