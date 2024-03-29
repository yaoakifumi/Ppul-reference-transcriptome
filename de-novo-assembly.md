# Parapercis pulchella reference transcriptome assembly
Reference transcriptome of the harlequin sandsmelt Parapercis pulehella was constructed from RNA-sequencing data.

These analyses were partly performed using [NIG supercomputer system](https://sc.ddbj.nig.ac.jp), National Institute of Genetics, Research Organaization of Information and Systemes, Japan.

## De novo assembly

### Adapter trimming using Trimmomatic
Remove adapter sequences and low quality reads.
Adapter sequences are deposited in "adapter.fa"

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) version 0.3.9

```shell
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

### Quality check using FastQC
Quality check 
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) version 0.11.9

```shell
fastqc paired_output_read1.fq.gz
fastqc paired_output_read2.fq.gz
```

### De novo assembly using Trinity
De novo assembly using Trinity with default parameters
[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) version 2.13.2

```shell
Trinity \
--seqType fq \
--left /<path>/<to>/paired_output_read1.fq.gz \
--right /<path>/<to>/paired_output_read2.fq.gz \
--max_memory 32G \
--CPU 12 \
--output /<path>/<to>/trinity_run
```

### Assembly statistics
[seqkit](https://bioinf.shenwei.me/seqkit/) version 2.4.0

```shell
seqkit stats -a Trinity.fasta
```

## superTranscripts construction
### mapping using salmon
[salmon](https://combine-lab.github.io/salmon/) version 1.7.0

```shell
salmon index -p 8 \
-t /<path>/<to>/Trinity.fasta \
-i /<path>/<to>/salmon_index

salmon quant -p 15 \
-i /<path>/<to>/salmon_index \
--libType A \
--dumpEq \
--hardFilter \
--skipQuant \
-1 /<path>/<to>/paired_output_read1.fq.gz \
-2 /<path>/<to>/paired_output_read2.fq.gz \
--output /<path>/<to>/salmon_raw/output
```

### Run Corset
[corset](https://github.com/Oshlack/Corset/wiki) version 1.09

```shell
corset \
-i salmon_eq_classes /<path>/<to>/salmon_raw/output/aux_info/eq_classes.txt
```

### Run Lace
[lace](https://github.com/Oshlack/Lace/wiki) version 1.00

```shell
Lace.py /<path>/<to>/Trinity.fasta \
clusters.txt \
-t -a --cores 32 \
-o /<path>/<to>/supertranscripts_lace/
```

### Rename contigs
[seqtk](https://github.com/lh3/seqtk) version 1.3

```shell
seqtk seq \
-C /<path>/<to>/supertranscripts_lace/SuperDuper.fasta  > \
/<path>/<to>/tmp.fasta

seqtk rename \
/<path>/<to>/tmp.fasta \
Pp- > \
/<path>/<to>/Ppul_supertranscripts.fasta
```

### Assembly statistics
[seqkit](https://bioinf.shenwei.me/seqkit/) version 2.4.0

```shell
seqkit stats -a Ppul_supertranscripts.fasta
```

### ORF prediction using TranDecoder
[TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki) version 5.7.0

```shell
TransDecoder.LongOrfs \
-t /<path>/<to>/Ppul_supertranscripts.fasta \
-O /<path>/<to>/transdecoder_ref_supertranscripts \
-m 100

TransDecoder.Predict \
-t /<path>/<to>/Ppul_supertranscripts.fasta \
--single_best_only \
-O /<path>/<to>/transdecoder_ref_supertranscripts
```

## Quality assessment
### BUSCO analysis
[BUSCO](https://busco.ezlab.org) version 5.1.2

#### Raw assembled contigs
```shell
busco \
-m transcriptome \
-i /<path>/<to>/Trinity.fasta \
-o busco_aftertrinity \
--auto-lineage-euk
```

#### Final superTranscripts
```shell
busco \
-m transcriptome \
-i /<path>/<to>/supertranscripts_lace/SuperDuper.fasta \
-o busco_supertranscripts_lace \
-c 10 \
--auto-lineage-euk
```

### Back-mapping using salmon
[samlon](https://combine-lab.github.io/salmon/) version 1.7.0

#### Raw assembled contigs (same as superTranscripts construction)
```
Described above
```

#### Final superTranscripts
```shell
#salmon index
salmon index -p 2 \
-t /<path>/<to>/supertranscripts_lace/SuperDuper.fasta \
-i /<path>/<to>/salmon_supertranscripts_lace/salmon_index

#mapping
salmon quant -p 15 \
-i /<path>/<to>/salmon_supertranscripts_lace/salmon_index \
--libType A --dumpEq --hardFilter --skipQuant \
-1 /<path>/<to>/paired_output_read1.fq.gz \
-2 /<path>/<to>/paired_output_read2.fq.gz \
--output /<path>/<to>/salmon_supertranscripte_lcae/output
```

### Back-mapping using STAR (for superTranscripts)
[STAR](https://github.com/alexdobin/STAR) version 2.7.8a

Multi-fasta file (SuperDuper.fasta) and annotation file (SuperDuper.gff) produced by [lace](https://github.com/Oshlack/Lace/wiki) were subjected to STAR as reference.

#### Make STAR index from SuperDuper.fasta and SuperDuper.gff
```shell
STAR \
--runMode genomeGenerate \
--genomeDir /<path>/<to>/STAR_lace_supertranscripts \
--runThreadN 8 \
--genomeFastaFiles /<path>/<to>/supertranscripts_lace/SuperDuper.fasta \
--sjdbGTFfile /<path>/<to>/supertranscripts_lace/SuperDuper.gff \
--limitGenomeGenerateRAM 64231571722
```

#### Back-mapping
```shell
STAR \
--genomeDir /<path>/<to>/STAR_lace_supertranscripts \
--runThreadN 24 \
--outFileNamePrefix /<path>/<to>/STAR_lace_supertranscripts/output/ \
--quantMode TranscriptomeSAM \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn /<path>/<to>/paired_output_read1.fq.gz \
/<path>/<to>/paired_output_read2.fq.gz \
--readFilesCommand gunzip -c
```

## Functional annotations

### eggNOG-mapper
[eggNog-mapper](http://eggnog-mapper.embl.de)

```shell
Run eggNOG-mapper on the web browser with default parameters
```

### KAAS
[KAAS-KEGG Automatic Annotation Server](https://www.genome.jp/tools/kaas/)

```shell
Run KAAS on the web browser with default parameter
```

### Reciprocal BLAST best-hit analysis via medaka non-redundant CDSs
#### Manage header name
For conducting reciprocal BLAST anaysis, header name of protein sequences (generated by TransDecoder) were modified.

```shell
sed Ppul_supertranscripts.fasta.transdecoder.pep -e 's/\.p.*//' > supertranscripts_protein.fasta
```

#### Pick up longest isoforms from RefSeq database using R
```R
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

#### blastp
[NCBI blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) version 2.13.0

```shell
#make blast database (from P. pulehcella superTranscripts)
makeblastdb \
-in /<path>/<to>/Ppul_supertranscripts.fasta.transdecoder.pep \
-out Ppul_supertranscripts_pep \
-dbtype prot \
-parse_seqids

#make blast database (from medaka non-redundant sequences)
makeblastdb \
-in /<path>/<to>/medaka_longest_peptide.fasta \
-out medaka_longest \
-dbtype prot \
-parse_seqids

#query: medaka non-redundant sequence, database: P. pulchella superTranscripts
blastp \
-query /<path>/<to>/medaka_longest_peptide.fasta \
-db /<path>/<to>/Ppul_supertranscripts_pep  \
-out /<path>/<to>/0816medaka_to_ppul.asn \
-evalue 1e-4 \
-outfmt 11 \
-max_target_seqs 1 \
-num_threads 16

#blast formatter
blast_formatter \
-archive /<path>/<to>/0816medaka_to_ppul.asn \
-outfmt 6 \
-out /<path>/<to>/0816medaka_to_ppul.tsv

#query: P. pulehcella superTranscripts, database: medaka non-redundant sequence
blastp \
-query /<path>/<to>/supertranscripts_protein.fasta \
-db /<path>/<to>/medaka_longest  \
-out /<path>/<to>/0816ppul_to_medaka.asn \
-evalue 1e-4 \
-outfmt 11 \
-max_target_seqs 1 \
-num_threads 16

#blast formatter
singularity exec /usr/local/biotools/b/blast:2.13.0--hf3cf87c_0/ \
blast_formatter \
-archive /<path>/<to>/0816ppul_to_medaka.asn \
-outfmt 6 \
-out /<path>/<to>/0816ppul_to_medaka.tsv
```

#### Formatting BLAST results using R

```R
library(tidyverse)
library(openxlsx)
library(janitor)
library(biomaRt)

#---------------------------------------------------------------
#Medaka gene annotation data 

#search ensemble database
db <- useMart("ensembl")
#view(listDatasets(db))

#get medaka dataset
medaka_data <- useDataset("olatipes_gene_ensembl", mart = db)

view(listAttributes(medaka_data))

#get data
medaka_genename <- getBM(attributes = c( 
  "refseq_peptide", "refseq_peptide_predicted", 
  "external_gene_name", "description"), 
  mart = medaka_data)

as_tibble(medaka_genename)

#make list
medaka_id <- medaka_genename %>%
  tidyr::pivot_longer(c("refseq_peptide", "refseq_peptide_predicted"), 
                      names_to = "mean", values_to = "refseq_ID") %>%
  filter(refseq_ID != "") %>%
  dplyr::select(refseq_ID, external_gene_name, description)

write_csv(medaka_id, "ensembl_medaka_id.csv")

#---------------------------------------------------------------
#integration of blast data
#Note: "ptm" means "BLAST P. pulchella to medaka" and "mtp" means "BLAST medaka to P. pulchella"

#load blast data
#Header information was attached to the tsv files before loading data.
medaka_to_ppul <- read_tsv("0816medaka_to_ppul.tsv")
ppul_to_medaka <- read_tsv("0816ppul_to_medaka.tsv")

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
medaka_id <- read_csv("/<path>/<to>/ensembl_medaka_id.csv")

reciprocal_blast_annotation <- merge_blast_nonredundant %>%
  left_join(medaka_id,
            by = c("accession_clean" = "refseq_ID")) %>%
  dplyr::select(id, accession, external_gene_name, description)

#remove redundant row
reciprocal_blast_annotation <- reciprocal_blast_annotation %>%
  dplyr::distinct(accession, id, .keep_all = TRUE)

#save the result
reciprocal_blast_annotation %>%
  write.xlsx("supertranscripts_reciprocal_blast_annotation.xlsx")
```

### anntation information integration using R

```R
library(tidyverse)
library(openxlsx)

#load data
eggnog_data <- read.xlsx("<path>/<to>eggnog_results.xlsx")
kaas_data <- read.xlsx("<path>/<to>/KAAS_results.xlsx")
blast_data <- read.xlsx("<path>/<to>/supertranscripts_reciprocal_blast_annotation.xlsx")

#data curation
eggnog_summary <- eggnog_data %>%
  select(query, Description, Preferred_name)
kaas_summary <- kaas_data %>%
  drop_na(KEGG)

#data integration
merge_tmp <- eggnog_summary %>%
  full_join(kaas_summary, by = c("query" = "id"))
merge_annotations <- merge_tmp %>%
  full_join(blast_data, by = c("query" = "id"))

#save the result
merge_annotations %>%
  write.xlsx("annotations_threemethods.xlsx")
```


## Data curation using other methods
We also conducted other methods for curation.

### CD-HIT-EST
[CD-HIT-EST](https://sites.google.com/view/cd-hit) version 4.8.1

```shell
cd-hit-est \
-i /<path>/<to>/Trinity.fasta \
-c 0.95 \
-T 15 \
-M 0 \
-o /<path>/<to>/trinity_cdhitest95.fasta
```

### EvidentialGene
[EvidentialGene](http://arthropods.eugenes.org/EvidentialGene/about/EvidentialGene_trassembly_pipe.html)
Note: Although the developer recommends multiple assembly data from different assemblers be used for EvidentialGene analysis, only the Trinity results were submitted to the analysis in this study.

```shell
perl /<path>/<to>/evigene/scripts/prot/tr2aacds4.pl \
-cdnaseq \
/<path>/<to>/Trinity.fasta \
-NCPU=10 \
-MAXMEM=32000 \
-MINCDS=100 \
-logfile
```

### ORF prediction
[TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki) version 5.7.0

CD-HIT-EST
```shell
TransDecoder.LongOrfs -t /<path>/<to>/trinity_cdhitest95.fasta \
-O /<path>/<to>/transdecoder_trinity_cdhitest95 \
-m 100

TransDecoder.Predict -t /<path>/<to>/trinity_cdhitest95.fasta \
--single_best_only \
-O /<path>/<to>/transdecoder_trinity_cdhitest95
```

EvidentialGene
```shell
TransDecoder.LongOrfs -t /<path>/<to>/okayset/Trinity.okay.tr \
-O /<path>/<to>/transdecoder_evigene \
-m 100

TransDecoder.Predict -t /<path>/<to>/okayset/Trinity.okay.tr \
--single_best_only \
-O /<path>/<to>/transdecoder_evigene
```


### Assembly statistics
[seqkit](https://bioinf.shenwei.me/seqkit/) version 2.4.0

```shell
seqkit stats -a trinity_cdhitest95.fasta.transdecoder.cds
seqkit stats -a /okayset/Trinity.okay.tr
```

### BUSCO analysis
[BUSCO](https://busco.ezlab.org) version 5.1.2

CD-HIT-EST
```shell
busco \
-m transcriptome \
-i /<path>/<to>/trinity_cdhitest95.fasta \
-o busco_cdhitest_fromtrinity \
-c 10 \
--auto-lineage-euk
```

EvidentialGene
```shell
busco \
-m transcriptome \
-i /<path>/<to>/okayset/Trinity.okay.tr \
-o busco_evidentialgene \
-c 10 \
--auto-lineage-euk
```

