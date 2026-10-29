######################################################################
#### Visualize cross-species mutations on relevant (cancer) genes ####
######################################################################

### Install
remotes::install_github("https://github.com/marsangar/xsmut")

### Load libraries
library("remotes")
library("openxlsx")
library("VariantAnnotation")
library("openxlsx")

source("/Users/ms84/Research/postdoc/projects/PrimateSomatic/src/genomic_toolkit.r")

library(xsmut)

devtools::document("/Users/ms84/Research/postdoc/projects/AI/xsmut")
devtools::test("/Users/ms84/Research/postdoc/projects/AI/xsmut")

### Custom functions

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


### Load metadata
dataMETA<-read.xlsx("/Users/ms84/Research/postdoc/projects/PrimateSomatic/metadata/metadata_Macaca_mulatta.xlsx")

cosmic_genes<-read.delim("/Users/ms84/Research/postdoc/projects/AI/xsmut/data/cosmic/Cosmic_CancerGeneCensus_v103_GRCh38.tsv.gz")
cosmic_genes$GENE_SYMBOL

### Load mutation data
mutations <- read.delim("/Users/ms84/Research/postdoc/projects/PrimateSomatic/results/v0.2/Macaca_mulatta/targeted_nanoseq/variants/annotated_mutations.txt")
mutations$chr<-rename_chroms(mutations$chr, from = "RefSeq-Accn", to="GenBank-Accn", ncbi_report = "/Users/ms84/Research/postdoc/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/assembly_report.txt")
mutations 

mutations[mutations$gene %in% c("NOTCH1") & mutations$sampleID %in% dataMETA[dataMETA$tissue=="Oesophagus","sample"],]
mutations[mutations$gene %in% c("NOTCH1") & mutations$sampleID %in% dataMETA[dataMETA$tissue=="Skin","sample"],]
mutations[mutations$gene %in% c("NOTCH1") & mutations$sampleID %in% dataMETA[dataMETA$tissue=="Prostate","sample"],]
mutations[mutations$gene %in% c("NOTCH1") & mutations$sampleID %in% dataMETA[dataMETA$tissue=="Bladder","sample"],]
mutations[mutations$gene %in% c("NOTCH1") & mutations$sampleID %in% dataMETA[dataMETA$tissue=="Liver","sample"],]



gene="CUL3"

unique(muts$gene)

#pdf("/Users/ms84/Research/postdoc/projects/PrimateSomatic/results/v0.2/Macaca_mulatta/targeted_nanoseq/drivers/xsmut/xsmut_report_hotspots.pdf", 12, 6)
#pdf("/Users/ms84/Research/postdoc/projects/PrimateSomatic/results/v0.2/Macaca_mulatta/targeted_nanoseq/drivers/xsmut/xsmut_report.pdf", 12, 6)

for(gene in c("TP53","NOTCH1","KEAP1","CUL3","TP63","KDM6A","IL7R","ERBB2","FAT1","ALB","CBLB","CBL", "CHD4","MTOR", "ARID2", "RB1", "ERBB4","ERBB3","NOTCH3","PPM1D","PIK3CA")){
print(gene)
  
  ### Create outdir
  outdir<-paste0("/Users/ms84/Research/postdoc/projects/PrimateSomatic/results/v0.2/Macaca_mulatta/targeted_nanoseq/drivers/xsmut/", gene)
  dir.create(file.path(outdir))
  
  ### Subset to mutations on the gene
  muts<-mutations[mutations$gene==gene,]

  ### Execute the xsmut pipeline:
  res <- xsmut::xsmut_pipeline(
    gene        = gene,
    mutations   = muts,
    aa_numbering = "human",
    highlight_colour = "red",
    highlight_style = "strong",
    hotspot_window = 1,
    hotspot_top_n = 15,
    species     = "macaque",
    assembly    = "Mmul_10",
    show_pictograms = TRUE,
    pictogram_size = 0.24,
    cosmic      = "/Users/ms84/Research/postdoc/projects/AI/xsmut/data/cosmic/Cosmic_MutantCensus_v104_GRCh38.tsv.gz",
    census_path = "/Users/ms84/Research/postdoc/projects/AI/xsmut/data/cosmic/Cosmic_CancerGeneCensus_v103_GRCh38.tsv.gz",
  
  )           
# patchwork; ggsave("TP53.pdf", res$plot, width = 12, height = 8)

  ### Display gene plot
  ggsave(file.path(outdir,  paste0(gene,".pdf")), res$gene_model_plot, width = 12, height = 8)
  ggsave(file.path(outdir,  paste0(gene,"_plot.pdf")), res$plot, width = 12, height = 8)
  ggsave(file.path(outdir,  paste0(gene,"_plot_hotspots.pdf")), res$plot_hotspots, width = 12, height = 8)
  
  ### Export mapping info
  write.xlsx(res$species_mutations, 
             file.path(outdir,  paste0(gene,"_xsmut.xlsx"))
  )
}
#dev.off()



### check results
res$gene_check         # in_census, tier, role, n_cosmic_mutations
res$model              # human canonical transcript, exons, CDS
res$cosmic_mutations   # all human COSMIC mutations for TP53
res$species_mutations  # your mutations + human_aa_pos / human_chrom / human_pos
res$plot               # patchwork; ggsave("TP53.pdf", res$plot, width = 12, height = 8)
