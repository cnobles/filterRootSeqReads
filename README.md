# filterRootSeqReads
Filter reads that appear in a subset of output sequence files from the original sequencing file based on read names. For example, if a subset of .fastq files have been manipulated during bioinformatic processing, then the original sequences that gave rise to the processed sequences can be obtained from the original fastq files. 

Use:
```
Rscript filterRootSeqReads.R -r path/to/rootSeqFile.fastq.gz -b path/to/subset/files.fasta
Rscript path/to/filterRootSeqReads.R -r sample.R1.fastq.gz sample.R2.fastq.gz -i read.ids.csv -o filData --compress
Rscript path/to/filterRootSeqReads.R -r R1.fasta R2.fasta I1.fasta -i ids.csv --cores 3
```

The new sequence files (R1, R2, I1.fastq) will be writen into an output directory, defaultly called "filteredData".

All matching is based on readNamePattern, no sequences are analyzed or used in the script, only loaded into memory and writen to the output file.

## Arguments

**[-r, --root]** Files with original or root sequences (.fastq / .fastq.gz / .fasta / ...) to filter. Can select multiple files. Must be fasta or fastq format.

**[-b, --branch]** File to filter extract filtered read ids, a subset of root sequences.

**[-i, --ids]** Text-based file with list of read names to filter from root sequences. Header included.

**[-o, --outputDir]** Output directory. Relative or absolute path.

**[--filePattern]** Pattern to include in output file name. ie. filtered.[filePattern].R1.fastq.

**[--readTypes]** Read types (R1, R2, I1, I2) identifiers for root seqs. Default: R1 R2 I1 I2. Identifiers need to be present in same format on root sequencing files.

**[--readNamePattern]** Regular expression for pattern matching read names. Should not contain R1/R2/I1/I2 specific components. Default is [\\w:-]+

**[-c, --cores]** Parallel processing option with r-parallel. Specify number of cores to use, number of cores will not be more than machine can provide or number of root files to filter, which ever is smaller.

**[--compress]** Compress output with gzip.

## Dependencies

* argparse
* pander
* ShortRead
* stringr
* parallel (if using multiple cores for parallel processing)
