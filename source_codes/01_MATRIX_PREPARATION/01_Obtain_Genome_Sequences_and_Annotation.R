
library("getopt")

library(rlang)
library(dplyr)
library(threadr)

options(stringAsfactors = F, useFancyQuotes = F)

args <- commandArgs(trailingOnly = TRUE)
option_specification <- matrix(c('input1', 'i1', 2, 'character',
                                 'input2', 'i2', 2, 'character',
                                 'input3', 'i3', 2, 'character',
                                 'input4', 'i4', 2, 'character',
                                 'output1', 'o1', 2, 'character',
                                 'output2', 'o2', 2, 'character'),
                               byrow=TRUE, ncol=4)

options = getopt(option_specification);

#########################
mainDic <- "/home/galaxy/tools/easyMF/01_MATRIX_PREPARATION/"

MHmakeRandomString <- function(n=1, lenght=12){
  randomString <- c(1:n) # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    lenght, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}

generateURL <- function(Type, DataType, Species){
  if(Type == "plants"){
    ftp <- paste0("ftp://ftp.ensemblgenomes.org/pub/release-", release_version, "/plants/")
    species_df <- read.delim(file = paste0(mainDic, "species_EnsemblPlants.txt"),
                             sep = '\t', header = FALSE, comment.char = "#")
  }else if(Type == "fungi"){
    ftp <- "ftp://ftp.ensemblgenomes.org/pub/release-44/fungi/"
    species_df <- read.delim(file = paste0(mainDic, "species_EnsemblFungi.txt"),
                             sep = '\t', header = FALSE, comment.char = "#")
  }else if(Type == "metazoa"){
    ftp <- "ftp://ftp.ensemblgenomes.org/pub/release-44/metazoa/"
    species_df <- read.delim(file = paste0(mainDic, "species_EnsemblMetazoa.txt"),
                             sep = '\t', header = FALSE, comment.char = "#")
  }else if(Type == "protists"){
    ftp <- "ftp://ftp.ensemblgenomes.org/pub/release-44/protists/"
    species_df <- read.delim(file = paste0(mainDic, "species_EnsemblProtists.txt"),
                             sep = '\t', header = FALSE, comment.char = "#")
  }else{
    ftp <- "ftp://ftp.ensembl.org/pub/release-97/"
    species_df <- read.delim(file = paste0(mainDic, "species_EnsemblVertebrates.txt"),
                             sep = '\t', header = FALSE, comment.char = "#")
  }
  
  species_value <- Species
  #assembly <- paste(strsplit(species_df$V3[which(species_df$V2 == Species)], " ")[[1]], collapse = "_")
  #species_value <- as.character(species_df[which(species_df$V1 == Species),]$V2)
  
  if(DataType == "Genome"){
    curFTP <- paste0(ftp, "fasta/", species_value, "/dna/")
    filenames <-  list_files_ftp(url = curFTP,
                                 credentials = "",
                                 sleep = NA,
                                 sort = FALSE,
                                 verbose = FALSE)
    url <- filenames[grep("dna_sm.toplevel.fa.gz", filenames)]
  }else if(DataType == "GTF"){
    curFTP <- paste0(ftp, "gtf/", species_value, "/")
    filenames <- list_files_ftp(url = curFTP, credentials = "", sleep = NA, sort = FALSE, verbose = FALSE)
    url <- filenames[grep(paste0(".",release_version, ".gtf.gz"),filenames)]
  }
  url
}

#########################


Type <- options$input1
release_version <- options$input2
Species <- options$input3
DataType <- options$input4


if( DataType == "Both") {
  
  for( i in c("Genome", "GTF")){
    
    url <- generateURL(Type = Type, DataType = i, Species = Species)
    
    dirname <- MHmakeRandomString()
    dir.create(paste0(mainDic, dirname), showWarnings = FALSE)
    
    if(i == "GTF"){
      outName <- paste0(paste(strsplit(Species, " ")[[1]], collapse = "_"),
                        "_", tolower(DataType), ".gtf.gz")
    }else{
      outName <- paste0(paste(strsplit(Species, " ")[[1]], collapse = "_"),
                        "_", tolower(DataType), ".fasta.gz")
    }
    
    cmd <- paste0("wget -O ", paste0(mainDic, dirname, "/"), outName, " ", url)
    system(command = cmd)
    system(command = paste0("gunzip ", paste0(mainDic,dirname, "/", outName)))
    
    outName <- gsub(".gz", "", outName)
    
    if( i == "Genome"){
      system(command = paste0("mv ", paste0(mainDic, dirname, "/", outName), " ", options$output1))
      system(command = paste0("rm -r ", paste0(mainDic, dirname)))
    } else {
      system(command = paste0("mv ", paste0(mainDic, dirname, "/", outName), " ", options$output2))
      system(command = paste0("rm -r ", paste0(mainDic, dirname)))
    }
  }
  
} else {
  url <- generateURL(Type = Type, DataType = DataType, Species = Species)
  
  dirname <- MHmakeRandomString()
  dir.create(paste0(mainDic, dirname), showWarnings = FALSE)
  
  if(DataType == "GTF"){
    outName <- paste0(paste(strsplit(Species, " ")[[1]], collapse = "_"),
                      "_", tolower(DataType), ".gtf.gz")
  }else{
    outName <- paste0(paste(strsplit(Species, " ")[[1]], collapse = "_"),
                      "_", tolower(DataType), ".fasta.gz")
  }
  cmd <- paste0("wget -O ", paste0(mainDic, dirname, "/"), outName, " ", url)
  system(command = cmd)
  system(command = paste0("gunzip ", paste0(mainDic,dirname, "/", outName)))
  
  outName <- gsub(".gz", "", outName)
  system(command = paste0("mv ", paste0(mainDic, dirname, "/", outName), " ", options$output1))
  system(command = paste0("rm -r ", paste0(mainDic, dirname)))
}

