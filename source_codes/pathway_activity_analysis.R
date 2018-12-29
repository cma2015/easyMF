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
                                'output1', 'o1', 2, 'character',
                                'pdf1', 'o2', 2, 'character'
),
byrow=TRUE, ncol=4);
# Parse options
options(warn=-1)
options = getopt(option_specification);


PathwayActiveAnalysis <- function(){
  
  gene_expression <- as.matrix(read.table(options$input1, header = T, sep = "\t", stringsAsFactors = F))
  gene_expression_name <- rownames(gene_expression)
  Amplitude <- as.matrix(read.table(options$input2, sep = "\t", quote = "", header = T, stringsAsFactors = F))
  Pattern <- as.matrix(read.table(options$input3, sep = "\t", quote = "", header = T, stringsAsFactors = F))
  pathwayAnno <- as.matrix(read.table(options$input4, sep = "\t", quote = "", header = T, stringsAsFactors = F))
  cpu <- as.numeric(options$input5)
  ###############  Amplitude  #######################
  # pathway analysis
  # download data from pathwayAnno (PMN13_July; 20180702; ftp://ftp.plantcyc.org/Pathways/Data_dumps)
  m1 <- match(pathwayAnno[,3], gene_expression_name)
  m2 <- match(pathwayAnno[,4], gene_expression_name)
  m <- cbind(m1, m2)
  m4 <- apply(m, 1, function(p){
    return(p[which(p >0)[1]])
  })
  m <- cbind(m,m4)
  pathwayAnno <- cbind(pathwayAnno[which(m4 >0),], gene_expression_name[m4[which(m4 >0)]])
  colnames(pathwayAnno)[5] <- "geneName"
  pathway <- unique(pathwayAnno[,1])
  pathwayNumber <- matrix(, length(pathway),1)
  rownames(pathwayNumber) <- pathway
  for(i in 1:length(pathway)){
    pathwayNumber[i,1] <- length(which(pathwayAnno[,1] == pathway[i]))
  }
  # range(pathwayNumber)
  # length(which(pathwayNumber >= 10)) # 233
  # length(which(pathwayNumber < 10))  # 348
  pathway_select <- pathway[which( (pathwayNumber >= 10) & (pathwayNumber <= 500) )]
  pathwayAnno_Specific <- NULL
  for(ll in 1:length(pathway_select)){
    pathwayAnno_Specific <- rbind(pathwayAnno_Specific, pathwayAnno[which(pathwayAnno[,1] == pathway_select[ll]),])
  }
  
  accumul <- gene_expression %*% Pattern
  pathway <- matrix(, ncol(accumul), length(pathway_select))
  rownames(pathway) <- colnames(accumul)
  colnames(pathway) <- pathway_select
  t1 <- Sys.time()
  for (patternID in 1:ncol(accumul)) {#calculate the statistics
    for(pathwayID in 1:length(pathway_select)){
      pathway[patternID, pathwayID] <- t.test(accumul[pathwayAnno_Specific[which(pathwayAnno_Specific[,1] == pathway_select[pathwayID]),"geneName"],patternID]
                                                  ,accumul[pathwayAnno_Specific[-which(pathwayAnno_Specific[,1] == pathway_select[pathwayID]),"geneName"],patternID])$statistic
    }
  }
  pathway <- apply(pathway, 2, scale)
  t2 <- Sys.time()
  print(t2 - t1)
  
  t1 <- Sys.time()
  pathwayAnno_active <- matrix(, length(pathway_select),1)
  N <- length(unique(pathwayAnno_Specific[,"geneName"]))
  for(termID in 1:ncol(pathway)){
    annota_genes <- pathwayAnno_Specific[which(pathwayAnno_Specific[,1] == pathway_select[termID]),"geneName"]
    M <- length(annota_genes)
    p_value <- sapply(1:length(annota_genes), function(ll){
      return(cor.test(accumul[annota_genes[ll],], pathway[,termID], alternative = "greater")$p.value)
    })
    n <- length(annota_genes)
    k <- length(which(p_value < 1e-02))
    pathwayAnno_active[termID,1] <- fisher.test(data.frame(c(M-k, N-M-n+k), c(k, n-k)), alternative = "less")$p.value
  }
  t2 <- Sys.time()
  print(t2 - t1)
  # length(which(pathwayAnno_active < 1e-05))
  pathwayAnno_active_sig <- pathway_select[which(pathwayAnno_active < 1e-05)]
  
  write.table(x = pathwayAnno[match(pathwayAnno_active_sig, pathwayAnno[,1]),1:2], file = options$output1, sep = "\t", row.names = F)
  
  nodes <- NULL
  for(cl in 1:length(pathwayAnno_active_sig)){#length(cluster_result)
    nodes <- rbind(nodes,cbind(id = pathwayAnno_active_sig[cl], media.type = cl, 
                               type.label = pathwayAnno_Specific[which(pathwayAnno_Specific[,1] == pathwayAnno_active_sig[cl])[1],2]))
  }
  nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
  nodes[,2] <- as.integer(nodes[,2])
  
  links <- NULL
  t1 <- Sys.time()
  options(stringsAsFactors=F)
  options(scipen=999)
  suppressPackageStartupMessages(library(doParallel))
  suppressPackageStartupMessages(library(foreach))
  cl <-  makeCluster(cpu)
  registerDoParallel(cl)
  resMatList <- foreach(numi = 1 : (length(pathwayAnno_active_sig) - 1) ) %dopar%{#(length(genes[,1]) - 1) 
    # for(numi in 1:length(genes[,1])){
    print(numi)
    templink <- NULL
    for(numj in (numi+1):length(pathwayAnno_active_sig)){
      temp <- cor.test(pathway[,pathwayAnno_active_sig[numi]], pathway[,pathwayAnno_active_sig[numj]])
      if(temp$p.value < 0.05 ){
        templink <- rbind(templink, cbind(from = pathwayAnno_active_sig[numi], to = pathwayAnno_active_sig[numj], weight = temp$estimate))
      }
    }
    return(templink)
  }
  stopCluster(cl)
  t2 <- Sys.time()
  print(t2 - t1)
  links <- do.call(rbind, resMatList)
  rownames(links) <- NULL
  links <- as.data.frame(links)
  links[,4] <- as.numeric(links[,3])
  
  library("RColorBrewer")
  library(igraph)
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  pdf(options$pdf1)
  plot(net, layout = layout_in_circle(net))
  dev.off()
  
}
PathwayActiveAnalysis()
#################################################################################################