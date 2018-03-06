## Introduction

YALFF (Yet Another Lossy FASTQ Filter) is a lossy compressor for fastq files which uses a Burrows and Wheeler Transform as an index.
This greatly reduces the amount of memory required compared to other tools such as [QUARTZ][1].

## Availability

YALFF is released under [GPLv3][2].
Because YALFF uses some BWA functions and the same format for the index, it requires the [zlib][3] dependency.

## Usage

YALFF preserves the orders of the reads of the input files so there is no need for a paired end mode.

* Smoothing a fastq file:
    
        cat file.fastq | ./yalff -d index.fa > file_smoothed.fastq
    
* Smoothing a gzipped fastq file:

        zcat file.fastq.gz | ./yalff -d index.fa | gzip > file_smoothed.fastq.gz
        
* Smoothing using an index loaded in shared memory:

        bwa shm index.fa
        zcat file.fastq.gz | ./yalff -shm index.fa | gzip > file_smoothed.fastq.gz
        bwa shm -d
        
* Other options available are:

        -k NUM	 k-mer length. [32] (max = 255 excluded)

        -m NUM	 Number of mismatches allowed per k-mer. [1]

        -c NUM	 Chunk size. The number of reads read at once on each iteration. [10000]

        -b CHAR	 Sanger threshold for a quality score to be considered. [$]

        -g CHAR	 Sanger threshold for a quality score to be considered correct independently from the dictionary. [I]

        -s NUM	 Number of bases to skip after each k-mer. A value of 0 checks all the k-mers. [0]

        -q CHAR	 Sanger value used as replacement during smoothing. [I]

        -e CHAR	 Sanger value used as an eventual replacement when a k-mer aligns badly. [j]

        -t NUM	 Number of threads available. [Hardware concurrency - 1]

        -h 	     See this help.

[1]: http://cb.csail.mit.edu/cb/quartz
[2]: http://en.wikipedia.org/wiki/GNU_General_Public_License
[3]: http://zlib.net
