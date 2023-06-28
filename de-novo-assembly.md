Theses analysis were partly performed using [NIG supercomputer system](https://sc.ddbj.nig.ac.jp), National Institute of Genetics, Research Organaization of Information and Systemes, Japan.

# De novo assembly


## Adapter trimming by Trimmomatic
remove adapter sequences and low quality reads.
adapter sequences are deposited in "adapter.fa"

[Trimmomatic]() version 0.3.9
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
[FastQC]() version 0.11.9
```
fastqc paired_output_read1.fq.gz
fastqc paired_output_read2.fq.gz
```
## De novo assembly by Trinity

[Trinity]() version 2.13.2
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

align_and_estimate_abundance.pl
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

abundance_estimates_to_matrix.pl
```
trinityrnaseq-v2.13.2/util/abundance_estimates_to_matrix.pl \
/<path>/<to>/results/align_RSEM/RSEM.isoforms.results \
--est_method RSEM \
--gene_trans_map /<patj>/<to>/align_RSEM/Trinity.fasta.gene_trans_map 
```

## Remove low-expressed contigs by Trinity scripts
これはTPM1以下で切った時なので、後で修正する
```
trinityrnaseq-v2.13.2/util/filter_low_expr_transcripts.pl \
--transcripts /<path>/<to>/results/Trinity.fasta \
--min_expr_any 1 \
--trinity_mode \
--matrix /<path>/<to>/results/align_RSEM/RSEM.isoform.TPM.not_cross_norm  \
> /<path>/<to>/results/Ppul_lowcut.fasta 
```

## Extract protein-coding contigs by TransDecoder
[TransDecoder]() version 

## Integrate similar contigs by CD-HIT
[CD-HIT]() version

## rename contigs
[seqtk]() version


# Quality checking

## BUSCO analysis


# Functional annotations


## eggNOG-mapper
```
Run eggNOG-mapper on the web browser with default parameters
```

## KAAS

```
Run KAAS on the web browser with default parameter
```

## Reciplocal BLAST best-hit via medaka non-redundant CDSs
Pick up longest isoforms from RefSeq database
```{R}

```

Construct BLAST database
```

```

Run BLAST
```

```

Formatting BLAST results
```

```

Integrate data
```{R}

```









