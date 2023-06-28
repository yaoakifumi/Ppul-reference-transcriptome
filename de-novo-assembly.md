# De novo assembly
## Adapter trimming by Trimmomatic
remove adapter sequences and low quality reads.
adapter sequences are deposited in "adapter.fa"
```
trimmomatic PE \ 
 -phred33 \
 -trimlog 20210828trimmomatic_log.txt \
 DNBSEQ_TY1_Read1.fq.gz \
 DNBSEQ_TY1_Read2.fq.gz \
 paired_output_read1.fq.gz \
 unpaired_output_read1.fq.gz \
 paired_output_read2.fq.gz \
 unpaired_output_read2.fq.gz \
 ILLUMINACLIP:adapter.fa:2:30:10 \
 LEADING:20 \
 TRAILING:20 \
 SLIDINGWINDOW:4:15 \
 MINLEN:50
```

## Quality check by FastQC
Quality check 
```
fastqc paired_output_read1.fq.gz
fastqc paired_output_read2.fq.gz
```
## De novo assembly by Trinity


## Mapping by Trinity scripts


## Remove low-expressed contigs by Trinity scripts



