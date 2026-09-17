# ============================================
# Title: Primate annotation script
# Author: Martín Santamarina García
# Date: 24-02-2026
# Description: This script retrieve annotations for primate reference genomes
# ============================================

### Load libraries
library(biomaRt)
library(karyoploteR)

### Load BioMart databases hosted by Ensembl
ensembl = useEnsembl(biomart="ensembl")
listDatasets(ensembl)

### Load BioMart databases for Macaca mulatta
ensembl = useEnsembl(biomart="ensembl", dataset="mmulatta_gene_ensembl") # Macaque genes (Mmul_10)
listAttributes(ensembl)
listAttributes(ensembl)[,1][grepl("size",listAttributes(ensembl)[,1])]

### Retrieve Gene Annotation
genes<-getBM(attributes=c('chromosome_name','start_position','end_position', 'strand', 'ensembl_gene_id','external_gene_name', "gene_biotype"),
                   # filters = 'chromosome_name', 
                   # values = '1', 
                    mart = ensembl)
genes$chr<-paste0("chr",genes$chromosome_name)
genes$gene_size<-genes$end_position-genes$start_position

dataGENES<-genes[,c('chr','start_position','end_position', 'strand', 'gene_size', 'ensembl_gene_id','external_gene_name', "gene_biotype")]
grGENES<-toGRanges(dataGENES)


### Retrieve Exon Annotation
exons<-getBM(attributes=c('chromosome_name','exon_chrom_start','exon_chrom_end', 'strand', 'ensembl_gene_id','external_gene_name', "ensembl_exon_id", "rank"),
                     #filters = 'chromosome_name', 
                     #values = '1', 
                    mart = ensembl)

exons$chr<-paste0("chr",exons$chromosome_name)
exons$exon_size<-exons$exon_chrom_end-exons$exon_chrom_start
dataEXONS<-exons[,c('chr','exon_chrom_start','exon_chrom_end', 'strand', 'exon_size', 'ensembl_gene_id','external_gene_name', "ensembl_exon_id", "rank")]
grEXONS<-toGRanges(dataEXONS)


### Load reference genome report file & select chromosomes
dataREPORT<-read.delim("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/Mmul_10_sequence_report.tsv", header = TRUE)
dataREPORT$Sequence.name

dataGENOME<-dataREPORT[dataREPORT$Sequence.name %in% paste0("chr",c(1:20,"X","Y")),c("Sequence.name","Seq.length")]
grGENOME<-toGRanges(data.frame(chr=dataGENOME$Sequence.name,start=0,end=dataGENOME$Seq.length))

# Create 10kb reference genome bins
bin_size <- 10000

grBINS <- tileGenome(
  seqlengths = setNames(dataGENOME$Seq.length, dataGENOME$Sequence.name),
  tilewidth = bin_size,
  cut.last.tile.in.chrom = TRUE
)
grBINS

dataEXONS[dataEXONS$external_gene_name %in% "TP53",]
dataEXONS[dataEXONS$external_gene_name %in% "MYC",]

### improve this file (!)
#BED<-dataEXONS[dataEXONS$external_gene_name %in% panelgenes_filt,]
#BED[,c("chr","exon_chrom_start", "exon_chrom_end", "external_gene_name", "rank")]
#BED$chr<-rename_chroms(BED$chr, from = "UCSC.style.name", to="RefSeq.seq.accession", ncbi_report = NCBI_REPORT)
#BED$gene_exon<-paste0(BED$external_gene_name, "_", BED$rank)
#BED<-BED[!duplicated(BED$gene_exon),]
#dim(BED)
#table(BED$external_gene_name)
#WRITE<-BED[c("chr","exon_chrom_start", "exon_chrom_end", "gene_exon")]
#write.table(WRITE, "~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/panel/targeted_panel_gene_exons.bed", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

