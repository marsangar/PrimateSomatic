# ============================================
# Title: Primate cosmic script
# Author: Martín Santamarina García (ms84@sanger.ac.uk)
# Description: This script implements comparative approaches between primate somatic mutation data and COSMIC mutational catalogues
# ============================================

### Load configs
source("/Users/ms84/Research/postdoc/projects/PrimateSomatic/Primate_master.r")


### Load Mutation Data
##dataSOM<-read.delim(paste0(RESDIR_VARIANTS,"/somatic_variants.txt"))
dataSOM<-read.delim("/Users/ms84/Research/postdoc/projects/PrimateSomatic/results/v0.1/Macaca_mulatta/targeted_nanoseq/variants/annotated_mutations.txt")

### Load COSMIC data
dataCOSMIC<-read.delim2(COSMIC_CENSUS_GENES_MUTATIONS)

dim(dataCOSMIC)
head(dataCOSMIC)

## Extract COSMIC genes
cosmic_genes<-unique(dataCOSMIC$GENE_SYMBOL)
cosmic_genes

### Extract macaque mutated genes
macaque_genes<-unique(dataSOM$gene)

## Get intersections
genes_common<-intersect(macaque_genes, cosmic_genes)
genes_macaque_only<-setdiff(macaque_genes, cosmic_genes)
genes_cosmic_only<-setdiff(cosmic_genes, macaque_genes)

### Visualize number of genes
gene_counts <- data.frame(
  Category = c("Common", "Macaque_only", "COSMIC_only"),
  Count = c(length(genes_common), length(genes_macaque_only), length(genes_cosmic_only))
)

ggplot(gene_counts, aes(x=Category, y=Count, fill=Category)) +
  geom_bar(stat="identity") + theme_minimal() +
  ggtitle("Comparisson of mutated genes - macaque vs COSMIC")


dataCOSMIC[1,]

# Pipeline básico para humanización de mutaciones macaque en R

# 1. Obtener lista de genes ortólogos y secuencias cDNA
library(biomaRt)
mart <- useMart("ensembl")
ensembl_macaque <- useDataset("mmulatta_gene_ensembl", mart=mart)
ensembl_human <- useDataset("hsapiens_gene_ensembl", mart=mart)

# Consulta genes ortólogos macaque a humano
orthologs <- getLDS(attributes=c("ensembl_gene_id","external_gene_name","cdna"),
                    filters="external_gene_name",
                    values=unique(dataSOM$gene),
                    mart=ensembl_macaque,
                    attributesL=c("ensembl_gene_id","external_gene_name","cdna"),
                    martL=ensembl_human)


# 2. Alinear secuencias cDNA para cada gen
# Para cada par cDNA macaque-humano alinea globalmente (sin gaps) para obtener mapeo posición a posición.
# Usa paquete Biostrings que tiene función pairwiseAlignment.
library(Biostrings)
# Ejemplo alineamiento cDNA de un gen
macaque_seq <- DNAString(orthologs$cdna[1])
human_seq <- DNAString(orthologs$cdna.1[1])
alignment <- pairwiseAlignment(macaque_seq, human_seq, type="global")


# 3. Crear el mapa posición a posición (mapeo)
Recorre el alineamiento para asignar cada posición en macaque a la posición correspondiente en humano (ignorando gaps).
Esto permite convertir coordenadas de mutaciones de macaque a humanas