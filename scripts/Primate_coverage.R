# ============================================
# Title: Primate coverage script
# Author: Martín Santamarina García
# Date: 12-02-2026
# Description: This script applies assess coverage
# ============================================

### Load libraries
library(biomaRt)
library(GenomicRanges)
library(stringr)
library(dplyr)
library(gridExtra)
library("karyoploteR")
library("googledrive")
library("googlesheets4")

### Custom functions

rename_chroms<- function(x, from="RefSeq.seq.accession", to="Chromosome.name", ncbi_report) {
  ## This function changes the format of a chromosome vector // options are: c("Chromosome.name","Sequence.name", "RefSeq.seq.accession","GenBank.seq.accession","UCSC.style.name")
  report <- read.delim(ncbi_report)
  
  # Validate column exists
  if (!from %in% colnames(report) || !to %in% colnames(report)) {
    stop(paste("incorrect requested formats"))
  }
  x<-as.character(x)
  renamed<-report[match(x, report[,from]),to]
  return(renamed)
}

### Master objects

ncbi_report_macaca="~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/Mmul_10_sequence_report.tsv"



### Define palette
tissue_cols <- c(
  "Bladder" = "#F8766D",
  "Liver" = "#A3A500",
  "Oesophagus" = "#00BF7D",
  "Prostate" = "#00B0F6",
  "Skin" = "#E76BF3"
)

### Load Metadata
dataMETA<-read_sheet("https://docs.google.com/spreadsheets/d/1k3srj6Kx14Np2vBxxbhEJZQ7s6TdEdo-82_gwwgW0Jk/edit?gid=0#gid=0")
dataMETA$age<-as.numeric(dataMETA$age)
str(dataMETA)

#### Load genes in panel data.frame (Human Source)
dfPANEL<-read.delim("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/resources/panels/Macaca_mulatta/Mmul_10/Sanger_TERT-v4_TE-95148282_hg19_highstringencyfilter_buccal_gene_list.tsv")[,c("gene","target_type")]
head(dfPANEL)


ontarget<-c(0.744007743645798,
            0.751718628835458,
            0.752629749442927,
            0.756528895584605,
            0.755494333679298,
            0.580901415785987,
            0.581044655415606,
            0.567915368185581,
            0.55630132812178,
            0.606928646768616,
            0.634601820737009,
            0.620562672150147,
            0.646939940652798,
            0.564039806488363,
            0.639647882771595,
            0.579475038486336,
            0.596646294641614,
            0.64948712999598,
            0.658651282413498,
            0.621989878469722,
            0.631517489117487)

median(ontarget)
range(ontarget)


### Load Mutations
dataMUT<-read.delim()


### Load Coverage files 
dataCOV<-read.delim("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/post/MQD0001j_tds0001.cov.bed.gz", header = FALSE)
head(dataCOV)

dataCOV$chr<-rename_chroms(dataCOV$V1,from = "RefSeq.seq.accession", to="Sequence.name", ncbi_report = ncbi_report_macaca)
dataCOV$start<-dataCOV$V2
dataCOV$end<-dataCOV$V3
dataCOV$coverage<-as.integer(str_split_fixed(dataCOV$V4,";",3)[,3])

grCOV<-toGRanges(dataCOV[,c("chr","start","end","coverage")])
grCOV



######################################################################################################
### Script to plot mean coverage across samples (calculate mean coverage per 10kb genomic bins) ####
######################################################################################################

######################################################
### Load duplex coverage data for multiple samples ###
######################################################

### Directory with coverage files
cov_dir <- path.expand(
  "~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/post/"
)

### List coverage files
cov_files <- list.files(
  cov_dir,
  pattern = "*.cov.bed.gz$",
  full.names = TRUE
)

#cov_files<-cov_files[1:3] # subset for dev

### Load coverage files into a list & reformat
cov_list <- lapply(cov_files, function(f) {
  
  print(paste0("Loading:", f))
  df<-read.delim(f, header = FALSE)
  df$chr<-rename_chroms(df$V1,from = "RefSeq.seq.accession", to="Sequence.name", ncbi_report = ncbi_report_macaca)
  df$start<-df$V2
  df$end<-df$V3
  df$coverage<-as.integer(str_split_fixed(df$V4,";",3)[,3])
  df$sample <- basename(f)
  
  return(df[, c("chr", "start", "end", "coverage", "sample")])
  
})
names(cov_list) <- basename(cov_files)
str(cov_list)

### Convert to GRranges
gr_list <- lapply(cov_list, function(df) {
  GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start, end = df$end),
    coverage = df$coverage,
    sample = df$sample
  )
})

######################################
##### Compute coverage matrices ######
######################################

### Compute coverage matrix (genes × samples)
cov_genes_mat <- matrix(NA, nrow = length(grGENES), ncol = length(gr_list))
colnames(cov_genes_mat) <- names(gr_list)

for (i in seq_along(gr_list)) {
  
  hits <- findOverlaps(grGENES, gr_list[[i]])
  
  # extract values
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  cov_vals <- mcols(gr_list[[i]])$coverage[sh]
  
  # aggregate per gene
  cov_by_gene <- tapply(cov_vals, qh, mean)
  cov_genes_mat[as.integer(names(cov_by_gene)), i] <- cov_by_gene
}
cov_genes_mat


### Compute coverage matrix (exon × samples)
cov_exons_mat <- matrix(NA, nrow = length(grEXONS), ncol = length(gr_list))
colnames(cov_exons_mat) <- names(gr_list)

for (i in seq_along(gr_list)) {
  
  hits <- findOverlaps(grEXONS, gr_list[[i]])
  
  # extract values
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  cov_vals <- mcols(gr_list[[i]])$coverage[sh]
  
  # aggregate per exon
  cov_by_exon <- tapply(cov_vals, qh, mean)
  cov_exons_mat[as.integer(names(cov_by_exon)), i] <- cov_by_exon
}
cov_exons_mat
dim(cov_exons_mat)


### Compute coverage matrix (bins × samples)
cov_bins_mat <- matrix(NA, nrow = length(grBINS), ncol = length(gr_list))
colnames(cov_bins_mat) <- names(gr_list)

for (i in seq_along(gr_list)) {
  
  hits <- findOverlaps(grBINS, gr_list[[i]])
  
  # extract values
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  cov_vals <- mcols(gr_list[[i]])$coverage[sh]
  
  # aggregate per exon
  cov_by_bin <- tapply(cov_vals, qh, mean)
  cov_bins_mat[as.integer(names(cov_by_bin)), i] <- cov_by_bin
}
cov_bins_mat

###
max(cov_genes_mat[,1], na.rm = TRUE)
max(cov_exons_mat[,1], na.rm = TRUE)
max(cov_bins_mat[,1], na.rm = TRUE)

### Compute summaries across samples in the cohort (!!! double-check)
gene_median_cov <- apply(cov_genes_mat, 1, median, na.rm = TRUE)
gene_max_cov <- apply(cov_genes_mat, 1, max)
gene_min_cov <- apply(cov_genes_mat, 1, min)
gene_sd_cov <- apply(cov_genes_mat, 1, sd, na.rm = TRUE)

exon_median_cov <- apply(cov_exons_mat, 1, median, na.rm = TRUE)
exon_max_cov <- apply(cov_exons_mat, 1, max)
exon_min_cov <- apply(cov_exons_mat, 1, min)
exon_sd_cov <- apply(cov_exons_mat, 1, sd, na.rm = TRUE)

bin_median_cov <- apply(cov_bins_mat, 1, median, na.rm = TRUE)
bin_max_cov <- apply(cov_bins_mat, 1, max)
bin_min_cov <- apply(cov_bins_mat, 1, min)
bin_sd_cov <- apply(cov_bins_mat, 1, sd, na.rm = TRUE)


dataGENES$gene_median_cov<-gene_median_cov
dataGENES$gene_max_cov<-gene_max_cov
dataGENES$gene_min_cov<-gene_min_cov
dataGENES$gene_sd_cov<-gene_sd_cov

dataEXONS$exon_median_cov<-exon_median_cov
dataEXONS$exon_max_cov<-exon_max_cov
dataEXONS$exon_min_cov<-exon_min_cov
dataEXONS$exon_sd_cov<-exon_sd_cov


head(dataGENES)
head(dataEXONS)

head(dataEXONS)
head(dataEXONS[dataEXONS$external_gene_name %in% "TP53",])

### Aggregate exon information for each gene
dfAGG <- aggregate(
  exon_median_cov ~ external_gene_name,
  data = dataEXONS,
  FUN = function(x) c(
    median = median(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    n = sum(!is.na(x))
  )
)
dfAGG

dataGENES$exon_median_cov.median<-NA
dataGENES[match(dfAGG$external_gene_name, dataGENES$external_gene_name),"exon_median_cov.median"]<-dfAGG

###
median(dataGENES[dataGENES$external_gene_name %in% dfPANEL$gene,"median_exon_median_cov"], na.rm=TRUE)

median(dataEXONS[,"exon_median_cov"], na.rm = TRUE)
median(dataEXONS[dataEXONS$external_gene_name %in% dfPANEL$gene,"exon_median_cov"], na.rm = TRUE)

###################################################
### Restrict the analysis to genes in the panel ###
###################################################

### Inspect genes in panel
dfPANEL

### Subset with genes in panel & add target type
dataGENES_panel<-dataGENES[dataGENES$external_gene_name %in% dfPANEL$gene,]
dataGENES_panel$target_type<-dfPANEL[match(dataGENES_panel$external_gene_name, dfPANEL$gene),"target_type"]
grGENES_panel<-toGRanges(dataGENES_panel)

### check missing genes (!)
dfPANEL$gene[!dfPANEL$gene %in% dataGENES_panel$external_gene_name]

## Set zoom coordinates
dataGENES_panel$zoom_coordinates<-paste0(dataGENES_panel$chr,":", dataGENES_panel$start_position-1000,"-", dataGENES_panel$end_position+1000)
dataGENES_panel$chr

dataGENES_panel$exon_median_cov.median

##############################
#### Screen Coverage Tool ####
##############################

## Set plotting parameters
#plotDefaultPlotParams()
pp <- getDefaultPlotParams(plot.type=1)
pp$data1inmargin <- 5
pp$ideogramheight <- 20
#pp$topmargin <- 480

### Create Screen data.frame with target genes in chromosomes
dataSCREEN<-dataGENES_panel[dataGENES_panel$chr %in% paste0("chr",c(1:20,"X","Y")),]
#dataSCREEN<-dataSCREEN[10,] # subset for dev

head(dataSCREEN)
dim(dataSCREEN)

i=3
pdf("~/Research/postdoc/projects/PrimateSomatic/Analysis/Macaca_mulatta/targeted_nanoseq/v0.1/QC/panel_genes_screening_NEW.pdf", width = 12, height = 8)

### Iterate through the panel of genes
for(i in 1:nrow(dataSCREEN)){

  ## For each gene in the panel
  gene<-dataSCREEN$external_gene_name[i]
  print(gene)
  target_type<-dataSCREEN$target_type[i]
  print(target_type)
  gene_size<-dataSCREEN$gene_size[i]
  gene_size
  
  zoom_coordinates<-dataSCREEN$zoom_coordinates[i]
  zoom_coordinates
  
  ## Retrieve coverage per exon metrics
  dfEXONS<-mcols(grEXONS)[mcols(grEXONS)$external_gene_name %in% gene,]
  dfEXONS
  
  missing_threshold<-200
  
  exon_nb<-nrow(dfEXONS)
  
  exon_median_cov<-median(dfEXONS$median_cov, na.rm = TRUE) ### (!) Only across those profiled
  exon_min_cov<-min(dfEXONS$median_cov, na.rm = TRUE) ### improve
  exon_max_cov<-max(dfEXONS$median_cov, na.rm = TRUE) ### improve
  
  ## Plot ideogram zooming to gene coordinates
  kp<-plotKaryotype(genome = grGENOME,  main=gene, plot.type = 1, plot.params = pp, zoom = zoom_coordinates)
  kpAddCytobandLabels(kp)
  
  
  ## Report coordinates
  #mtext(zoom_coordinates, line = 0,  col = "gray")
  
  
  ## Report relevant coverage per exon metrics
  mtext(paste0("Target type: ", target_type), line = -3)
  mtext(paste0("Number of exons: ", exon_nb), line = -4)
  mtext(paste0("Median coverage: ", exon_median_cov, "X"), line = -5)
  #mtext(paste0("Median coverage: ", exon_median_cov, " (", exon_min_cov, " - ", exon_max_cov, ")"), line = -4)

  abline(h=0.8, col="darkred", lty=2 , lwd=1)

  ### Display coverage across the gene
  #kpArea(kp, chr = dataCOV$chr,  x=dataCOV$start, y=dataCOV$coverage,  ymax = 1000, col = "#007C7C")
  
  ### Plot coverage histograms (one per sample)
  n <- length(gr_list)
  
  # Color palette (one per sample)
  cols <- colorRampPalette(c("#007C7C", "#D95F02", "#7570B3"))(n)
  
  for (j in seq_along(gr_list)) {
    print("test")
    
    #r0 <- (j-1)/n
    #r1 <- j/n
    
    r0 <- (j-1)/n + 0.1
    r1 <- j/n + 0.1
    
    kpArea(
      kp,
      data = gr_list[[j]],
      y = gr_list[[j]]$coverage,
      ymax = 1000,
      r0 = r0,
      r1 = r1,
      col = cols[j],
      border = NA
    )
    kpAxis(kp, data.panel=1, ymin = 0, ymax=1000, numticks = 5, side = 1, r0 = r0, r1 = r1, cex=0.2)
  }
  
  #kpAxis(kp, data.panel=1, ymin = 0, ymax=1000, numticks = 5, side = 1)
  kpAddLabels(kp, labels="Coverage (X)", data.panel = 1, srt=90, r0 = 0.5, label.margin = 0.08)
  
  kpAddBaseNumbers(kp, tick.dist = gene_size/10, tick.len = 10, tick.col="darkred", cex=1, units = "auto", add.units="FALSE",	
                                    minor.tick.dist = gene_size/100, minor.tick.len = 5, minor.tick.col = "gray")

  ### Display gene structure
  kpLines(kp, chr=dataSCREEN$chr[i], x=c(dataSCREEN$start_position[i], dataSCREEN$end_position[i]), y=c(0.025,0.025), lwd=5, col="darkgray")
  kpRect(kp, chr = dataEXONS$chr, x0=dataEXONS$exon_chrom_start, x1=dataEXONS$exon_chrom_end, y0=0, y1=0.05, data.panel = 1, col="#D95F02")
  
  ### Legend
  legend("topright",legend = c("exon", "gene"), fill = c("#D95F02", "gray"), cex = 1)
}
dev.off()

### Assess genes
dataEXONS[dataEXONS$external_gene_name %in% "MYC",]
dataEXONS[dataEXONS$external_gene_name %in% "TP53",]


