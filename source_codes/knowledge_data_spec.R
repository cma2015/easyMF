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
