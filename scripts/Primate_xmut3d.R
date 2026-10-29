############################################################################################
#### Visualize cross-species mutations on relevant (cancer) protein 3D structure models ####
############################################################################################

### Install

# Bioconductor dependencies (sequence alignment)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings", "pwalign"))

# CRAN dependencies
install.packages(c("httr2", "jsonlite", "data.table", "ggplot2", "patchwork",
                   "memoise", "r3dmol", "htmlwidgets"))

# the companion package, then this one
remotes::install_github("https://github.com/marsangar/xsmut")
remotes::install_github("https://github.com/marsangar/xsmut3d")

# optional: static PNG snapshots of the 3D view
install.packages("webshot2")


### Load libraries
library("xsmut")
library("xsmut3d")

### Rename chromosome naming between different alternative formats
rename_chroms<- function(x, from="RefSeq.seq.accession", to="Chromosome.name", ncbi_report) {
  ## This function changes the format of a chromosome vector // options are: c("Chromosome.name","Sequence.name", "RefSeq.seq.accession","GenBank.seq.accession","UCSC.style.name")
  report <- read.delim(ncbi_report, header = FALSE, sep = "\t", comment.char = "#",
                       quote = "", na.strings = "na", colClasses = "character")
  colnames(report)<-c("Sequence-Name","Sequence-Role",	"Assigned-Molecule",	"Assigned-Molecule-Location/Type",	"GenBank-Accn",	"Relationship",	"RefSeq-Accn",	"Assembly-Unit",	"Sequence-Length",	"UCSC-style-name")
  
  # Validate column exists
  if (!from %in% colnames(report) || !to %in% colnames(report)) {
    stop(paste("incorrect requested formats"))
  }
  x<-as.character(x)
  renamed<-report[match(x, report[,from]),to]
  return(renamed)
}




### Load mutation data
muts <- data.table::fread(
  system.file("extdata/example_macaque_TP53_annotmuts.tsv", package = "xsmut3d"))


### Load metadata
dataMETA<-read.xlsx("/Users/ms84/Research/postdoc/projects/PrimateSomatic/metadata/metadata_Macaca_mulatta.xlsx")

cosmic_genes<-read.delim("/Users/ms84/Research/postdoc/projects/AI/xsmut/data/cosmic/Cosmic_CancerGeneCensus_v103_GRCh38.tsv.gz")
cosmic_genes$GENE_SYMBOL

### Load mutation data
mutations <- read.delim("/Users/ms84/Research/postdoc/projects/PrimateSomatic/results/v0.2/Macaca_mulatta/targeted_nanoseq/variants/annotated_mutations.txt")
mutations$chr<-rename_chroms(mutations$chr, from = "RefSeq-Accn", to="GenBank-Accn", ncbi_report = "/Users/ms84/Research/postdoc/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/assembly_report.txt")
mutations 

gene="CBL"

for(gene in c("PIK3CA","IL7R", "TP53","NOTCH1","NOTCH3","KEAP1","CUL3","ALB","CBLB","CHD4","RB1", "PPM1D", "TP63", "RAC1", "ERBB2", "KRAS", "HRAS")){

  ### Subset to mutations on the gene
  muts<-mutations[mutations$gene==gene,]
  
### Execute the xsmut3d pipeline:
res <- xsmut3d_pipeline(
  gene        = gene,
  mutations   = muts,
  species     = "macaque",
  assembly    = "Mmul_10",
  cosmic      = "~/Research/postdoc/projects/AI/xsmut/data/cosmic/Cosmic_MutantCensus_v104_GRCh38.tsv.gz",
  census_path = "~/Research/postdoc/projects/AI/xsmut/data/cosmic/Cosmic_CancerGeneCensus_v103_GRCh38.tsv.gz",
  outdir      = paste0("/Users/ms84/Research/postdoc/projects/PrimateSomatic/results/v0.2/Macaca_mulatta/targeted_nanoseq/drivers/xsmut3d/", gene)
)
}

res                      # one-screen summary
res$view                 # rotate it (RStudio Viewer or a browser)
res$plot                 # ggsave("TP53.pdf", res$plot, width = 11, height = 9.5)
top_hits(res, 10)        # the calls worth looking at first
cat(res$confidence$text) # what you may and may not interpret
