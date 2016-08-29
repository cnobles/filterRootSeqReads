# filterRootSeqReads
Filter reads that appear in a subset of output sequence files from the original sequencing file based on read names. For example, if a subset of .fastq files have been manipulated during bioinformatic processing, then the original sequences that gave rise to the processed sequences can be obtained from the original fastq files. 

Use:
```
Rscript filterRootSeqReads.R -root path/to/original/fastq/files -branch path/to/subset/fastq/files #To initiate the script
```

The new sequence files (R1, R2, I1.fastq) will be writen into a directory called "filteredData" in the location which the script was executed. 

All matching is based on read names, no sequences are analyzed or used in the script, only loaded into memory and writen to the output file. 
