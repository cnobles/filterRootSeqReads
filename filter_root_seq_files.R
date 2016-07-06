dependancies <- c("dplyr", "tidyr", "magrittr", "ShortRead")
null <- sapply(dependancies, function(x){
  suppressPackageStartupMessages(
    try(library(x, character.only = TRUE), silent = TRUE))
})

dependancies_present <- sapply(dependancies, function(package){
  package <- paste0("package:", package)
  logic <- package %in% search()
})

if(FALSE %in% dependancies_present){
  Unloaded_Packages <- data.frame(package=as.character(dependancies), 
                                  loaded=dependancies_present)
  write.table(Unloaded_Packages, 
              file = "Unloaded_Packages.tsv", 
              quote = FALSE, 
              row.names = FALSE)
  stop("Load required packages. Check Unloaded_Packages.tsv for missing
       dependancies.")
}else{
  remove(dependancies, dependancies_present)
  message("Required packages loaded.")
}

####    Parse directories from command line    ####
args <- commandArgs(trailingOnly = TRUE)

rootDataDir <- args[ grep("root", args) + 1 ]
branchDataDir <- args[ grep("branch", args) + 1 ]

message(paste0("rootPath: ", rootDataDir))
message(paste0("branchPath: ", branchDataDir))

####    Load required files    ####
rootDataFiles <- list.files(rootDataDir, pattern = "fastq.gz")
message(paste0("root sequencing files:\n", paste(rootDataFiles, collapse = "\n")))
rootData <- lapply(c("R1", "I1", "R2"), function(readType){
  readFastq(
    file.path(rootDataDir, grep(readType, rootDataFiles, value = TRUE)),
    withIds = TRUE)
})

branchData <- list.files(branchDataDir, pattern = "fastq.gz")
message(paste0("branch sequencing files:\n", paste(branchData, collapse = "\n")))
branchData <- readFastq(branchDataDir, withIds = TRUE)

n_root_reads <- sum(sapply(rootData, length))/3
message(paste0("Reads loaded from root: ", n_root_reads))
n_branch_reads <- length(branchData)/2
message(paste0("Reads loaded from branches: ", n_branch_reads))

####    Select only read names present in branches    ####
filteredRootData <- lapply(rootData, function(rootReads){
  rootReads[rootReads@id %in% branchData@id]
})

####    Write new files of the filtered reads    ####
system("mkdir filteredData")
lapply(1:length(filteredRootData), function(i){
  filteredData <- filteredRootData[[i]]
  readType <- c("R1", "I1", "R2")[i]
  fileName <- grep(readType, rootDataFiles, value = TRUE)
  writeFastq(filteredData, file = file.path("filteredData", fileName))
})
q()
