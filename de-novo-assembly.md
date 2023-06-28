# De novo assembly

Theses analysis were partly performed using NIG supercomputer system.

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

```
Trinity \
--seqType fq \
--left /<path>/<to>/paired_output_read1.fq.gz \
--right /<path>/<to>/paired_output_read2.fq.gz \
--max_memory 32G \
--CPU 12 \
--output /<path>/<to>/trinity_run
```

## Mapping by Trinity scripts
```
trinityrnaseq-v2.13.2/util/align_and_estimate_abundance.pl \
--transcripts /<path>/<to>/Trinity.fasta \
--seqType fq \
--left /<path>/<to>/paired_output_read1.fq.gz \
--right /<path>/<to>/paired_output_read2.fq.gz \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--thread_count 6 \
--prep_reference \
--output_dir  /<path>/<to>/results/align_RSEM
```

## Remove low-expressed contigs by Trinity scripts


## Extract protein-coding contigs by TransDecoder


## Integrate similar contigs by CD-HIT


## rename contigs


## 

# Quality checking

## BUSCO analysiss




