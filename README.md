## Introduction

YALFF (Yet Another Lossy FASTQ Filter) is a lossy ~~compressor~~ smoother (see NOTE) for FASTQ files which uses a Burrows and Wheeler Transform as an index.
This greatly reduces the amount of memory required compared to other tools such as [QUARTZ][1]. This is because the dictionary of k-mers is stored in a linearized form, i.e. a reference sequence. The reference can be a standard reference sequence of an organism (like the hg38 for humans) or a custom sequence built from multiple k-mer counting procedures.

NOTE - the actual compression is achieved by standard lossless compressors such as gzip, bzip2 and xz. The compression ratio is increased by the smoothing procedure which basically reduces entrophy by replacing most of the quality values with a fixed value. Obviously the algorithm guarantees that the most relevant qualities for downstream analysis are kept untouched.

## Quickest solution

Version presented at [BITS Turin 2018][2] "Better quality score compression through sequence-based quality smoothing" (under review).

The fastest way to get started is by indexing a reference genome for the reads to be compressed. YALFF is 100% compliant with the index produced by **bwa**. So compliant that the whole indexing procedure and shared memory operations are handled by bwa himself with the implicit advantage of not requiring a separate construction if an index is already available.

A copy of [bwa][3] is downloaded as an external dependency (remember to clone with the '--recurse-submodules' flag).

## Most performant solution

Version to be presented at [BIOINFORMATICS 2019][16] "Indexing k-mers in linear-space for quality value compression".

It is also possible to reassemble a k-mer dictionary with the [assembler][4] developed for the [ProPhyle][5] package.
For this purpose the script folder contains the utility print_mitdb.c to print the k-mers of the dictionaries generated by [Quartz][1] and [LAVA][6].

The whole procedure can be described as follows:

1. Download the Quartz dictionary from the original links [1][7] and [2][8].
2. Clone the [LAVA repository][9], compile it and download the [hg17][10] and the associated SNP lists [dbSNPs][11], [affy][12].
3. Run LAVA on its files to construct its reference indeces.
4. Use [print_mitdb.c](scripts/print_mitdb.c) to stream the contents of the 4 dictionaries (the unswapped one for Quartz and the 3 for LAVA) to the Prophasm assembler.
5. Index the FASTA file produced by the assembler.

The commands described here can be found at the [scripts folder](./scripts).

The results produced by this index are better both in terms of overall compression and Precision and Recall than the standard reference, similarly to what happens in [GeneCodeq][13].

## Evaluation

The Precision, Recall and F-Measure are computed aligning the smoothed dataset to a reference and comparing the quality to a standard ground truth. The dataset used for comparison in our study is the Platinum genome NA12878 and relative vcf files downloaded from the Illumina ftp (which is open access and if it asks for a password just continue).

A [script](scripts/evaluation.job) is available containing the whole pipeline used for evaluation.

## Availability

~~YALFF is released under [GPL][14].~~

YALFF has been re-licenced under the more permissive [MIT](./LICENSE) licence.
Be aware that the final executable **will be GPL-ed** because of the linking at object level with bwa.

Because YALFF uses some BWA functions and the same format for the index, it requires the [zlib][15] dependency.

## Usage

YALFF preserves the orders of the reads of the input files so there is no need for a paired-end mode.

* Smoothing a fastq file:
    
        cat file.fastq | ./yalff -d index.fa > file_smoothed.fastq
    
* Smoothing a gzipped fastq file:

        zcat file.fastq.gz | ./yalff -d index.fa | gzip > file_smoothed.fastq.gz
        
* Smoothing using an index loaded in shared memory:

        bwa shm index.fa
        zcat file.fastq.gz | ./yalff -shm index.fa | gzip > file_smoothed.fastq.gz
        bwa shm -d
        
* All the options available are:

        -d STR   Reference file in fasta format.

        -k NUM	 k-mer length. [32] (max = 255 excluded)

        -m NUM	 Number of mismatches allowed per k-mer. [1]

        -c NUM	 Chunk size. The number of reads read at once on each iteration. [10000]

        -b CHAR	 Sanger threshold for a quality score to be considered. [$]

        -g CHAR	 Sanger threshold for a quality score to be considered correct independently from the dictionary. [I]

        -s NUM	 Number of bases to skip after each k-mer. A value of 0 checks all the k-mers. [0]

        -q CHAR	 Sanger value used as replacement during smoothing. [I]

        -e CHAR	 Sanger value used as an eventual replacement when a k-mer aligns badly. [j]

        -t NUM	 Number of threads available. [Hardware concurrency - 1]
        
        -sst NUM Smoothing algorithm. 0 checks all and only the k-mers considered. 1 applies a seed and extend search if a k-mer has no mismatches. [0]
        
        -shm STR Reference file loaded into shared memory.
        
        -h       See this help.

[1]: http://cb.csail.mit.edu/cb/quartz
[2]: http://bioinformatics.it/bits2018
[3]: http://bio-bwa.sourceforge.net/
[4]: https://github.com/prophyle/prophasm
[5]: https://prophyle.github.io
[6]: http://cb.csail.mit.edu/cb/lava/
[7]: http://giant.csail.mit.edu/quartz/dec200.bin.sorted.gz
[8]: http://giant.csail.mit.edu/quartz/dec200.bin.sorted.swapped.gz
[9]: https://github.com/arshajii/lava/
[10]: http://cb.csail.mit.edu/cb/lava/data/hg19.fa.gz
[11]: http://cb.csail.mit.edu/cb/lava/data/SNPs142_hg19_Common.filt.txt
[12]: http://cb.csail.mit.edu/cb/lava/data/Affymetrix_6_SNPs.txt
[13]: https://www.ncbi.nlm.nih.gov/pubmed/27354700
[14]: http://en.wikipedia.org/wiki/GNU_General_Public_License
[15]: http://zlib.net
[16]: http://www.bioinformatics.biostec.org/

