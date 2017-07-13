options(stringsAsFactors = FALSE)
suppressMessages(library("argparse"))
suppressMessages(library("pander"))

#' For those reviewing the code below, the following is a small style guide 
#' outlining the various formats for the code. 
#' 
#' Names with "_": objects, inlucding data.frames, GRanges, vectors, ...
#' Names in caMel format: functions or components of objects (i.e. columns 
#' within a data.frame).
#' Names with ".": arguments / options for functions

code_dir <- dirname(
  sub("--file=", "", 
      grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))

# Set up and gather command line arguments -------------------------------------
parser <- ArgumentParser(
  description = "R-based tool for filtering a list of reads form original sequencing files.")
parser$add_argument(
  "-r", "--root", nargs = "+", type = "character", default = NULL,
  help = "Files with original or root sequences (.fastq / .fastq.gz / .fasta / ...) to filter. Can select multiple files. Must be fasta or fastq format.")
parser$add_argument(
  "-b", "--branch", nargs = 1, type = "character", default = NULL,
  help = "File to filter extract filtered read ids, a subset of root sequences.")
parser$add_argument(
  "-i", "--ids", nargs = 1, type = "character", default = NULL,
  help = "Text-based file with list of read names to filter from root sequences. Header included.")
parser$add_argument(
  "-o", "--outputDir", nargs = 1, type = "character", default = "filteredData",
  help = "Output directory.")
parser$add_argument(
  "--filePattern", nargs = 1, type = "character", default = NULL,
  help = "Pattern to include in output file name. ie. filtered.[filePattern].R1.fastq.")
parser$add_argument(
  "--readTypes", nargs = "+", type = "character", default = "R1 R2 I1 I2", 
  help = "Read types (R1, R2, I1, I2) identifiers for root seqs. Default: R1 R2 I1 I2. Identifiers need to be present in same format on root sequencing files.")
parser$add_argument(
  "--readNamePattern", nargs = 1, type = "character", default = "[\\w:-]+",
  help = "Regular expression for pattern matching read names. Should not contain R1/R2/I1/I2 specific components. Default is [\\w:-]+")
parser$add_argument(
  "-c", "--cores", nargs = 1, type = "integer", default = 0,
  help = "Parallel processing option with r-parallel. Specify number of cores to use, number of cores will not be more than machine can provide or number of root files to filter, which ever is smaller.")
parser$add_argument(
  "--compress", action = "store_true", help = "Compress output with gzip.")

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Argument Conditionals
if(is.null(args$branch) & is.null(args$ids)){
  stop("Please specify file for branch or file with read ids.")}

args$readTypes <- unlist(strsplit(args$readTypes, " "))

# Print Inputs to terminal
input_table <- data.frame(
  "Variables" = paste0(names(args), " :"), 
  "Values" = sapply(1:length(args), function(i){
    paste(args[[i]], collapse = ", ")}))
input_table <- input_table[
  match(c("root :", "branch :", "ids :", "outputDir :", "filePattern :", 
          "readTypes :", "readNamePattern :", "cores :", "compress :"),
        input_table$Variables),]
pandoc.title("filterRootSeqReads Inputs")
pandoc.table(data.frame(input_table, row.names = NULL), 
             justify = c("left", "left"), 
             split.tables = Inf,
             style = "simple")

# Load additional R-packages for analysis and processing
add_packs <- c("stringr", "ShortRead")
add_packs_loaded <- suppressMessages(
  sapply(add_packs, require, character.only = TRUE))
if(!all(add_packs_loaded)){
  pandoc.table(data.frame(
    "R-Packages" = names(add_packs_loaded), 
    "Loaded" = add_packs_loaded, 
    row.names = NULL))
  stop("Check dependancies.")
}

# Load id file for filtering ---------------------------------------------------
if(!is.null(args$branch)){
  if(grepl("fastq", args$branch)){
    branchData <- readFastq(args$branch)
  }else{
    branchData <- readFasta(args$branch)
  }
  ids <- str_extract(as.character(id(branchData)), args$readNamePattern)
}else if(!is.null(args$ids)){
  ids <- read.delim(args$ids, header = TRUE)
  ids <- str_extract(as.character(ids[,1]), args$readNamePattern)
}else{
  stop("Please provided either a branch file or list of ids.")
}

## Make system folder of output ================================================
if(!dir.exists(args$outputDir)){
  system(paste0("mkdir ", args$outputDir))
  if(!dir.exists(args$outputDir)) stop("Cannont create output directory.")
}

# Filter root sequence files ---------------------------------------------------
if(args$cores == 0){
  null <- lapply(args$root, function(root, args, ids){
    rootFormat <- ifelse(grepl("fastq", root), "fastq", "fasta")
    if(rootFormat == "fastq"){
      rootReads <- readFastq(root)
    }else{
      rootReads <- readFasta(root)
    }
    pander(paste0("\nReads from ", root, " loaded, totaling ", length(rootReads), " reads."))
    rootIDs <- str_extract(as.character(id(rootReads)), args$readNamePattern)
    filteredReads <- rootReads[match(ids, rootIDs)]
    pander(paste0("\nRoot sequences filtered to ", length(filteredReads), " reads."))
    rootType <- args$readTypes[sapply(args$readTypes, grepl, x = root)]
    if(length(rootType) > 1 | length(rootType) == 0){
      stop("Root seq files do not contain read type identifiers. Change file names or identifiers.")
    }
    fileName <- paste0(
      "filteredData.", args$filePattern, ".", rootType, ".", rootFormat)
    if(args$compress) fileName <- paste0(fileName, ".gz")
    if(rootFormat == "fastq"){
      writeFastq(
        filteredReads, 
        file = file.path(args$outputDir, fileName), 
        compress = args$compress)
    }else{
      writeFasta(
        filteredReads, 
        file = file.path(args$outputDir, fileName), 
        compress = args$compress)
    }
  }, args = args, ids = ids)
}else{
  parLoaded <- suppressMessages(require(parallel))
  if(!parLoaded) stop("R-parallel not loaded, check to see if package is installed.")
  
  buster <- makeCluster(min(detectCores(), args$cores, length(args$root)))
  
  null <- parLapply(buster, args$root, function(root, args, ids){
    library(ShortRead)
    library(stringr)
    library(pander)
    rootFormat <- ifelse(grepl("fastq", root), "fastq", "fasta")
    if(rootFormat == "fastq"){
      rootReads <- readFastq(root)
    }else{
      rootReads <- readFasta(root)
    }
    pander(paste0("\nReads from ", root, " loaded, totaling ", length(rootReads), " reads."))
    rootIDs <- str_extract(as.character(id(rootReads)), args$readNamePattern)
    filteredReads <- rootReads[match(ids, rootIDs)]
    pander(paste0("\nRoot sequences filtered to ", length(filteredReads), " reads."))
    rootType <- args$readTypes[sapply(args$readTypes, grepl, x = root)]
    if(length(rootType) > 1 | length(rootType) == 0){
      stop("Root seq files do not contain read type identifiers. Change file names or identifiers.")
    }
    fileName <- paste0(
      "filteredData.", args$filePattern, ".", rootType, ".", rootFormat)
    if(args$compress) fileName <- paste0(fileName, ".gz")
    if(rootFormat == "fastq"){
      writeFastq(
        filteredReads, 
        file = file.path(args$outputDir, fileName), 
        compress = args$compress)
    }else{
      writeFasta(
        filteredReads, 
        file = file.path(args$outputDir, fileName), 
        compress = args$compress)
    }
  }, args = args, ids = ids)
  stopCluster(buster)
}
