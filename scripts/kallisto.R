# Test run Kallisto
#
# > kallisto quant -i transcripts.idx -o output -b 100 reads_1.fastq.gz reads_2.fastq.gz

reads_1.fastq <- read.table("~/kallisto_windows-v0.46.1/kallisto/test/reads_1.fastq.gz", quote="\"", comment.char="")
reads_2.fastq <- read.table("~/kallisto_windows-v0.46.1/kallisto/test/reads_2.fastq.gz", quote="\"", comment.char="")
#View(reads_1.fastq)
#View(reads_2.fastq)
dim(reads_1.fastq)
dim(reads_2.fastq)

# To visualize pseudoalignments  we use --genomebam, which requires a GTF file and text file with each chrom length
# > kallisto quant -i transcripts.kidx -b 30 -o kallisto_out --genomebam --gtf transcripts.gtf.gz --chromosomes chrom.txt reads_1.fastq.gz reads_2.fastq.gz
# 
# combined read files for each cart with:
# > START copy *.gz /b C7_all.fastq.gz 

# Build kallisto index:
# > kallisto index -i C1_transcripts.idx C1_all.fastq.gz
# Run quantification alg:
# > kallisto quant -i C1_transcripts.idx -o output
########################
