These analysis were partly performed using [NIG supercomputer system](https://sc.ddbj.nig.ac.jp), National Institute of Genetics, Research Organaization of Information and Systemes, Japan.

# De novo assembly


## Adapter trimming by Trimmomatic
Remove adapter sequences and low quality reads.
Adapter sequences are deposited in "adapter.fa"

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

## Pick up renamined contigs
[seqkit]() version2.4

## rename contigs
[seqtk]() version
```

```

# Quality checking

## BUSCO analysis
[BUSCO]() version

Raw assembled contigs
```

```

Final reference transcriptome
```

```

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
### Pick up longest isoforms from RefSeq database by R
```{R}
#library
library(orthologr)
library(biomartr)

#get data from Ensembl
medaka_proteome <- biomartr::getProteome(db = "refseq", 
                    organism = "Oryzias latipes")
medaka_gff <- biomartr::getGFF(db = "refseq", 
                    organism = "Oryzias latipes")

retrieve_longest_isoforms(proteome_file = medaka_proteome, 
                          annotation_file = medaka_gff, 
                          new_file = "medaka_longest_peptide.fasta")
```

### blastp version
```

```

### Formatting BLAST results by R
```{R}
#library
library(tidyverse)
library(openxlsx)
library(janitor)

#load blast data
#Header information was attached to the tsv files before loading data.
medaka_to_ppul <- read_tsv("blastp_medaka_to_ppul.tsv")
ppul_to_medaka <- read_tsv("blastp_ppul_to_medaka.tsv")

#table curation
medaka_to_ppul <- medaka_to_ppul %>%
  select(accession, id, evalue_mtp, bit_score_mtp) 
ppul_to_medaka <- ppul_to_medaka %>%
  select(id, accession, evalue_ptm, bit_score_ptm)

medaka_to_ppul <- medaka_to_ppul %>%
  mutate(id = str_replace_all(medaka_to_ppul$id, 
                              "\\.p*\\d", ""))

#merge data
merge_blast <- medaka_to_ppul %>%
  inner_join(ppul_to_medaka, by = c("id", "accession")) %>%
  as_tibble()

#retain only results with the highest bit-score
merge_blast_nonredundant <- merge_blast %>%
  group_by(accession, id) %>%
  slice_max(bit_score_mtp, n = 1, with_ties = TRUE) %>%
  slice_max(bit_score_ptm, n = 1, with_ties = FALSE) %>%
  as_tibble()

#modify accession ID
merge_blast_nonredundant <- merge_blast_nonredundant %>%
  mutate(accession_clean = str_replace_all(merge_blast_nonredundant$accession, 
                                     "\\.\\d", ""))
#add gene name
medaka_id <- read_csv("/Users/yaoakifumi/Desktop/Research/bioinformatics/Ppul_analysis/re_assembly/annotation/ensembl_medaka_id.csv")

reciprocal_blast_annotation <- merge_blast_nonredundant %>%
  left_join(medaka_id,
            by = c("accession_clean" = "refseq_ID")) %>%
  dplyr::select(id, accession, external_gene_name, description)

#remove redundant row
reciprocal_blast_annotation <- reciprocal_blast_annotation %>%
  dplyr::distinct(accession, id, .keep_all = TRUE)

#save data
reciprocal_blast_annotation %>%
  write.xlsx("supertranscripts_reciprocal_blast_annotation.xlsx")
```

Integrate data by R
```{R}
library(tidyverse)
library(openxlsx)

load data
eggnog_data <- read.xlsx("<path>/<to>eggnog_results.xlsx")
kaas_data <- read.xlsx("<path>/<to>/KAAS_results.xlsx")
blast_data <- read.xlsx("<path>/<to>/supertranscripts_reciprocal_blast_annotation.xlsx")



```









