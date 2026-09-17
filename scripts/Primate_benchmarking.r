############################################################
#### Script designed to compare variant calling methods ####
############################################################

### Load packages
library(VennDiagram)
library(grid)
library(vcfR)
library(Biostrings)
library(sigfit)
library(gridExtra)
library(ggplot2)
library(ggpubr)

#library(MutationalPatterns)

### Define functions
get_context <- function(chr, pos, ref_genome) {
  if (!(chr %in% names(ref_genome))) return(NA)
  # Extract 3bp region
  seq <- subseq(ref_genome[[chr]], pos-1, pos+1)
  as.character(seq)
}
make_spectrum <- function(vcf_file, dnaREF){
  
  vcf <- read.vcfR(vcf_file)
  
  dataVCF <- as.data.frame(vcf@fix)
  dataVCF$POS <- as.integer(dataVCF$POS)
  
  ## SNVs only
  dataVCF <- dataVCF[nchar(dataVCF$REF) == 1 & nchar(dataVCF$ALT) == 1,]
  
  ## sequence context
  dataVCF$contextREF <- apply(
    dataVCF,
    1,
    function(x)
      get_context(as.character(x["CHROM"]),
                  as.integer(x["POS"]),
                  dnaREF)
  )
  
  dataVCF$contextREF_5   <- substr(dataVCF$contextREF,1,1)
  dataVCF$contextREF_3   <- substr(dataVCF$contextREF,3,3)
  
  dataVCF$contextALT <-
    paste0(dataVCF$contextREF_5,
           dataVCF$ALT,
           dataVCF$contextREF_3)
  
  ## pyrimidine normalization
  dataVCF$contexREF_normalized <-
    ifelse(
      dataVCF$REF %in% c("A","G"),
      reverseComplement(DNAStringSet(dataVCF$contextREF)),
      dataVCF$contextREF
    )
  
  dataVCF$contextALT_normalized <-
    ifelse(
      dataVCF$REF %in% c("A","G"),
      reverseComplement(DNAStringSet(dataVCF$contextALT)),
      dataVCF$contextALT
    )
  
  dataVCF$mutation_normalized <-
    paste0(dataVCF$contexREF_normalized,
           ">",
           dataVCF$contextALT_normalized)
  
  data("cosmic_signatures_v2")
  
  mutations <- rep(0,96)
  names(mutations) <- colnames(cosmic_signatures_v2)
  
  tbl <- table(dataVCF$mutation_normalized)
  
  mutations[names(tbl)] <- as.integer(tbl)
  
  mutations
}

### Load metadata and order
dataMETA<-read.xlsx(METADATA)
dataMETA <- dataMETA[
  order(dataMETA$tissue,
        dataMETA$age,
        dataMETA$individual),
]


### Load kraken results & append to metadata
dataKRAKEN_macaca<-read.delim(paste0(KRAKEN_DIR, "KRAKEN.RESULTS_Macaca_mulatta.tsv"), header = FALSE)
colnames(dataKRAKEN_macaca)[2]<-c("kraken_macaca")
dataKRAKEN_macaca$sample<-str_split_fixed(dataKRAKEN_macaca$V1,"analysis/|_",3)[,2]
dataMETA<-dataMETA %>% dplyr::left_join(dataKRAKEN_macaca[,c("sample","kraken_macaca")], by = "sample")


### Load concordance results (pass snv & pass indel)
dataCOUNTS_snv<-read.delim("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/snv/venn_counts.tsv")
dataCOUNTS_indel<-read.delim("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/indel/venn_counts.tsv")

head(dataCOUNTS_snv)
head(dataCOUNTS_indel)

dataCOUNTS_snv$sample_id<-dataCOUNTS_snv$sample
dataCOUNTS_indel$sample_id<-dataCOUNTS_indel$sample

dataCOUNTS_snv$fraction_shared_snv<-(dataCOUNTS_snv$shared/(dataCOUNTS_snv$tdc_private+dataCOUNTS_snv$shared+dataCOUNTS_snv$tns_private))*100
dataCOUNTS_indel$fraction_shared_indel<-(dataCOUNTS_indel$shared/(dataCOUNTS_indel$tdc_private+dataCOUNTS_indel$shared+dataCOUNTS_indel$tns_private))*100

head(dataCOUNTS_snv)
head(dataCOUNTS_indel)

### Append concordance results
dataMETA<-dataMETA %>% dplyr::left_join(dataCOUNTS_snv[,c("sample_id","fraction_shared_snv")], by = "sample_id")
dataMETA<-dataMETA %>% dplyr::left_join(dataCOUNTS_indel[,c("sample_id","fraction_shared_indel")], by = "sample_id")


### Generate scatterplot for snv
p <- ggplot(
  dataMETA,
  aes(x = kraken_macaca, y = fraction_shared_snv)
) +
  geom_point(
    size = 3,
    shape = 21,
    fill = "steelblue",
    colour = "black",
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    colour = "tomato",
    fill = "tomato",
    alpha = 0.2,
    linewidth = 1
  ) +
  stat_cor(
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 5
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = "Fraction of reads assigned to macaque (kraken2)",
    y = "Fraction of concordant SNVs (called by DC & NS)",
    title = "Relationship between concordance and purity"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )
pdf(paste0("/Users/ms84/Research/postdoc/projects/PrimateSomatic/benchmarking/scatter_PASS_snv_macaque_purity.pdf"))
print(p)
dev.off()

### Generate scatterplot for indel
p <- ggplot(
  dataMETA,
  aes(x = kraken_macaca, y = fraction_shared_indel)
) +
  geom_point(
    size = 3,
    shape = 21,
    fill = "seagreen3",
    colour = "black",
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    colour = "darkorange2",
    fill = "darkorange2",
    alpha = 0.2,
    linewidth = 1
  ) +
  stat_cor(
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 5
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = "Fraction of reads assigned to macaque (kraken2)",
    y = "Fraction of concordant INDELs (called by DC & NS)",
    title = "Relationship between concordance and purity"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )
pdf(paste0("/Users/ms84/Research/postdoc/projects/PrimateSomatic/benchmarking/scatter_PASS_indel_macaque_purity.pdf"))
print(p)
dev.off()

##############################################
#### Stage 1: Variant Calling Comparisson ####
##############################################

### snv comparisson
dataCOUNTS<-read.delim("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/snv/venn_counts.tsv")
head(dataCOUNTS)

dataCOUNTS$fraction_tdc_private<-dataCOUNTS$tdc_private/(dataCOUNTS$tdc_private+dataCOUNTS$shared+dataCOUNTS$tns_private)
dataCOUNTS$fraction_shared<-dataCOUNTS$shared/(dataCOUNTS$tdc_private+dataCOUNTS$shared+dataCOUNTS$tns_private)
dataCOUNTS$fraction_tns_private<-dataCOUNTS$tns_private/(dataCOUNTS$tdc_private+dataCOUNTS$shared+dataCOUNTS$tns_private)
head(dataCOUNTS)

mean(dataCOUNTS$fraction_tdc_private)
mean(dataCOUNTS$fraction_shared)
mean(dataCOUNTS$fraction_tns_private)

range(dataCOUNTS$fraction_shared)
barplot(dataCOUNTS$fraction_shared)

pdf(paste0("/Users/ms84/Research/postdoc/projects/PrimateSomatic/benchmarking/venn_PASS_snv_macaque_increasing_purity.pdf"), 12, 24)

grid.newpage()

pushViewport(viewport(
  layout = grid.layout(7, 3)
))

for(i in seq_len(nrow(dataCOUNTS))) {
  
  vp <- viewport(
    layout.pos.row = ceiling(i / 3),
    layout.pos.col = ((i - 1) %% 3) + 1
  )
  
  pushViewport(vp)
  
  grid.draw(
    draw.pairwise.venn(
      area1 = dataCOUNTS$tdc_private[i] + dataCOUNTS$shared[i],
      area2 = dataCOUNTS$tns_private[i] + dataCOUNTS$shared[i],
      cross.area = dataCOUNTS$shared[i],
      category = c("TDC", "TNS"),
      fill = c("steelblue", "tomato"),
      alpha = c(0.5, 0.5),
      scaled = TRUE
    )
  )
  
  grid.text(
    dataCOUNTS$sample[i],
    y = unit(0.95, "npc"),
    gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  upViewport()
}
dev.off()




### indel comparisson
dataCOUNTS<-read.delim("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/indel/venn_counts.tsv")
head(dataCOUNTS)

dataCOUNTS$fraction_tdc_private<-dataCOUNTS$tdc_private/(dataCOUNTS$tdc_private+dataCOUNTS$shared+dataCOUNTS$tns_private)
dataCOUNTS$fraction_shared<-dataCOUNTS$shared/(dataCOUNTS$tdc_private+dataCOUNTS$shared+dataCOUNTS$tns_private)
dataCOUNTS$fraction_tns_private<-dataCOUNTS$tns_private/(dataCOUNTS$tdc_private+dataCOUNTS$shared+dataCOUNTS$tns_private)
head(dataCOUNTS)

mean(dataCOUNTS$fraction_tdc_private)
mean(dataCOUNTS$fraction_shared)
mean(dataCOUNTS$fraction_tns_private)

range(dataCOUNTS$fraction_shared)

barplot(dataCOUNTS$fraction_shared)



pdf(paste0("/Users/ms84/Research/postdoc/projects/PrimateSomatic/benchmarking/venn_PASS_indel_macaque_increasing_purity.pdf"), 12, 24)

grid.newpage()

pushViewport(viewport(
  layout = grid.layout(7, 3)
))

for(i in seq_len(nrow(dataCOUNTS))) {
  
  vp <- viewport(
    layout.pos.row = ceiling(i / 3),
    layout.pos.col = ((i - 1) %% 3) + 1
  )
  
  pushViewport(vp)
  
  grid.draw(
    draw.pairwise.venn(
      area1 = dataCOUNTS$tdc_private[i] + dataCOUNTS$shared[i],
      area2 = dataCOUNTS$tns_private[i] + dataCOUNTS$shared[i],
      cross.area = dataCOUNTS$shared[i],
      category = c("TDC", "TNS"),
      fill = c("steelblue", "tomato"),
      alpha = c(0.5, 0.5),
      scaled = TRUE
    )
  )
  
  grid.text(
    dataCOUNTS$sample[i],
    y = unit(0.95, "npc"),
    gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  upViewport()
}
dev.off()


### snv underpowered comparisson
dataCOUNTS<-read.delim("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/snv_underpowered/venn_counts.tsv")
head(dataCOUNTS)

dataCOUNTS$fraction_tdc_private<-dataCOUNTS$tdc_private/(dataCOUNTS$tdc_private+dataCOUNTS$shared+dataCOUNTS$tns_private)
dataCOUNTS$fraction_shared<-dataCOUNTS$shared/(dataCOUNTS$tdc_private+dataCOUNTS$shared+dataCOUNTS$tns_private)
dataCOUNTS$fraction_tns_private<-dataCOUNTS$tns_private/(dataCOUNTS$tdc_private+dataCOUNTS$shared+dataCOUNTS$tns_private)
head(dataCOUNTS)

mean(dataCOUNTS$fraction_tdc_private)
mean(dataCOUNTS$fraction_shared)
mean(dataCOUNTS$fraction_tns_private)

range(dataCOUNTS$fraction_shared)
barplot(dataCOUNTS$fraction_shared) ### 0.01320132 0.13036810

pdf(paste0("/Users/ms84/Research/postdoc/projects/PrimateSomatic/benchmarking/venn_UNDERPOWERED_snv_macaque_increasing_purity.pdf"), 12, 24)

grid.newpage()

pushViewport(viewport(
  layout = grid.layout(7, 3)
))

for(i in seq_len(nrow(dataCOUNTS))) {
  
  vp <- viewport(
    layout.pos.row = ceiling(i / 3),
    layout.pos.col = ((i - 1) %% 3) + 1
  )
  
  pushViewport(vp)
  
  grid.draw(
    draw.pairwise.venn(
      area1 = dataCOUNTS$tdc_private[i] + dataCOUNTS$shared[i],
      area2 = dataCOUNTS$tns_private[i] + dataCOUNTS$shared[i],
      cross.area = dataCOUNTS$shared[i],
      category = c("dupcaller_underpowered", "nanoseq_pass"),
      fill = c("steelblue", "tomato"),
      alpha = c(0.5, 0.5),
      scaled = TRUE
    )
  )
  
  grid.text(
    dataCOUNTS$sample[i],
    y = unit(0.95, "npc"),
    gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  upViewport()
}
dev.off()



### Plot with proper margins
pdf(
  "/Users/ms84/Research/postdoc/projects/PrimateSomatic/benchmarking/dupcaller-1.1.2_underpowered_venn_snv_unfiltered.pdf",
  width = 12,
  height = 24
)

grid.newpage()

pushViewport(
  viewport(
    layout = grid.layout(
      7,
      3,
      widths = unit(c(1, 1, 1), "null"),
      heights = unit(rep(1, 7), "null")
    )
  )
)

for (i in seq_len(nrow(dataCOUNTS))) {
  
  vp <- viewport(
    layout.pos.row = ceiling(i / 3),
    layout.pos.col = ((i - 1) %% 3) + 1,
    
    # Add space around each Venn
    x = 0.5,
    y = 0.48,
    width = 0.85,
    height = 0.80
  )
  
  pushViewport(vp)
  
  grid.draw(
    draw.pairwise.venn(
      area1 = dataCOUNTS$tdc_private[i] + dataCOUNTS$shared[i],
      area2 = dataCOUNTS$tns_private[i] + dataCOUNTS$shared[i],
      cross.area = dataCOUNTS$shared[i],
      
      category = c(
        "dupcaller_underpowered",
        "nanoseq_pass"
      ),
      
      fill = c("steelblue", "tomato"),
      alpha = c(0.5, 0.5),
      scaled = TRUE,
      
      # Fixed label positions
      cat.pos = c(180, 0),
      
      # Distance of labels from circles
      cat.dist = c(0.03, 0.03),
      
      cat.cex = 0.9
    )
  )
  
  grid.text(
    dataCOUNTS$sample[i],
    y = unit(0.95, "npc"),
    gp = gpar(
      fontsize = 10,
      fontface = "bold"
    )
  )
  
  upViewport()
}

dev.off()



### indel underpowered comparisson
dataCOUNTS<-read.delim("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/indel_underpowered/venn_counts.tsv")
head(dataCOUNTS)

dataCOUNTS$fraction_tdc_private<-dataCOUNTS$tdc_private/(dataCOUNTS$tdc_private+dataCOUNTS$shared+dataCOUNTS$tns_private)
dataCOUNTS$fraction_shared<-dataCOUNTS$shared/(dataCOUNTS$tdc_private+dataCOUNTS$shared+dataCOUNTS$tns_private)
dataCOUNTS$fraction_tns_private<-dataCOUNTS$tns_private/(dataCOUNTS$tdc_private+dataCOUNTS$shared+dataCOUNTS$tns_private)
head(dataCOUNTS)

mean(dataCOUNTS$fraction_tdc_private)
mean(dataCOUNTS$fraction_shared)
mean(dataCOUNTS$fraction_tns_private)

range(dataCOUNTS$fraction_shared) # 0.00000000 0.02670623
barplot(dataCOUNTS$fraction_shared)

pdf(paste0("/Users/ms84/Research/postdoc/projects/PrimateSomatic/benchmarking/venn_UNDERPOWERED_indel_macaque_increasing_purity.pdf"), 12, 24)

grid.newpage()

pushViewport(viewport(
  layout = grid.layout(7, 3)
))

for(i in seq_len(nrow(dataCOUNTS))) {
  
  vp <- viewport(
    layout.pos.row = ceiling(i / 3),
    layout.pos.col = ((i - 1) %% 3) + 1
  )
  
  pushViewport(vp)
  
  grid.draw(
    draw.pairwise.venn(
      area1 = dataCOUNTS$tdc_private[i] + dataCOUNTS$shared[i],
      area2 = dataCOUNTS$tns_private[i] + dataCOUNTS$shared[i],
      cross.area = dataCOUNTS$shared[i],
      category = c("TDC_underpowered", "TNS"),
      fill = c("steelblue", "tomato"),
      alpha = c(0.5, 0.5),
      scaled = TRUE
    )
  )
  
  grid.text(
    dataCOUNTS$sample[i],
    y = unit(0.95, "npc"),
    gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  upViewport()
}
dev.off()






######################################
#### Stage 2: Spectra Comparisson ####
######################################


### Section A) Compare Targeted Dupcaller vs Targeted Nanoseq

# ---- Read reference genome ----
dnaREF<-readDNAStringSet("/Users/ms84/Research/postdoc/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa")
names(dnaREF) <- str_split_fixed(names(dnaREF), " ", 2)[,1]   # clean headers


### Generate spectra plots grouped by tissue
base_dir <- "~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/snv"


## Maximum number of samples in any tissue
max_samples <- max(table(dataMETA$tissue))


pdf(
  "/Users/ms84/Research/postdoc/projects/PrimateSomatic/benchmarking/tdc_tns_comparisson_spectra_PASSsnv_unfiltered.pdf",
  width = 24,
  height = 5 * max_samples
)


for (tissue in unique(dataMETA$tissue)) {
  
  message("Plotting ", tissue)
  
  tissue_samples <- dataMETA$sample_id[
    dataMETA$tissue == tissue
  ]
  
  n_samples <- length(tissue_samples)
  
  
  ## Fixed layout based on maximum tissue size
  par(
    mfrow = c(max_samples, 3),
    mar = c(6, 6, 12, 2),
    oma = c(0, 0, 4, 0)
  )
  
  
  ## Plot samples
  for (sample in tissue_samples) {
    
    ## Get metadata
    sample_meta <- dataMETA[
      dataMETA$sample_id == sample,
    ]
    
    age <- sample_meta$age

    
    dir <- file.path(
      base_dir,
      paste0("isec_", sample)
    )
    
    
    tdc <- make_spectrum(
      file.path(dir,"0000.vcf"),
      dnaREF
    )
    
    tns <- make_spectrum(
      file.path(dir,"0001.vcf"),
      dnaREF
    )
    
    shared <- make_spectrum(
      file.path(dir,"0002.vcf"),
      dnaREF
    )
    
    
    cs_tdc <- cosine_sim(tdc, shared)
    cs_tns <- cosine_sim(tns, shared)
    
    
    plot_spectrum(
      tdc,
      name = "TDC private"
    )
    
    
    plot_spectrum(
      shared,
      name = sprintf(
        "%s\nShared",
        sample
      )
    )
    
    
    plot_spectrum(
      tns,
      name = "TNS private"
    )
  }
  
  
  ## Fill remaining rows with empty plots
  if (n_samples < max_samples) {
    
    for (i in 1:(max_samples - n_samples)) {
      plot.new()
      plot.new()
      plot.new()
    }
  }
  
  
  ## Tissue page title
  mtext(
    paste("SNV spectra -", tissue),
    side = 3,
    outer = TRUE,
    line = 1,
    cex = 2,
    font = 2
  )
  
}

dev.off()




### Section B) Compare Targeted Dupcaller Underpowered vs Targeted Dupcaller PASS

# ---- Read reference genome ----
dnaREF<-readDNAStringSet("/Users/ms84/Research/postdoc/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa")
names(dnaREF) <- str_split_fixed(names(dnaREF), " ", 2)[,1]   # clean headers


### Generate spectra plots grouped by tissue

## Maximum number of samples in any tissue
max_samples <- max(table(dataMETA$tissue))


pdf(
  "/Users/ms84/Research/postdoc/projects/PrimateSomatic/benchmarking/dupcaller-1.1.2_underpowered_spectra_snv_unfiltered.pdf",
  width = 24,
  height = 5 * max_samples
)


for (tissue in unique(dataMETA$tissue)) {
  
  message("Plotting ", tissue)
  
  tissue_samples <- dataMETA$sample_id[
    dataMETA$tissue == tissue
  ]
  
  n_samples <- length(tissue_samples)
  
  
  ## Fixed layout based on maximum tissue size
  par(
    mfrow = c(max_samples, 3),
    mar = c(6, 6, 12, 2),
    oma = c(0, 0, 4, 0)
  )
  
  
  ## Plot samples
  for (sample in tissue_samples) {
    
    ## Get metadata
    sample_meta <- dataMETA[
      dataMETA$sample_id == sample,
    ]
    
    age <- sample_meta$age
    
    
    dupcaller_underpowered <- make_spectrum(
      paste0("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.", sample,"/results.", sample, "_snv.underpowered.vcf.gz"),
      dnaREF
    )
    
    dupcaller_pass <- make_spectrum(
      paste0("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.", sample,"/results.", sample, "_snv.pass.vcf.gz"),
      dnaREF
    )
    
    nanoseq_pass <- make_spectrum(
      paste0("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/filtered/pass/", sample, ".filtered.pass.snv.vcf.gz"),
      dnaREF
    )
    
    
    ### compute cosine similarity
    #cs_underpowered <- cosine_sim(dupcaller_underpowered, nanoseq_pass)
    #cs_pass <- cosine_sim(dupcaller_pass, nanoseq_pass)
    
    
    plot_spectrum( 
      dupcaller_underpowered,
      name = "dupcaller_underpowered"
      )
    
    
    plot_spectrum( 
      dupcaller_pass, 
      name = sprintf("%s\ndupcaller_pass", sample)
    )
    
    
    plot_spectrum(
      nanoseq_pass,
      name = "nanoseq_pass"
    )
  }
  
  
  ## Fill remaining rows with empty plots
  if (n_samples < max_samples) {
    
    for (i in 1:(max_samples - n_samples)) {
      plot.new()
      plot.new()
      plot.new()
    }
  }
  
  
  ## Tissue page title
  mtext(
    paste("SNV spectra -", tissue),
    side = 3,
    outer = TRUE,
    line = 1,
    cex = 2,
    font = 2
  )
  
}

dev.off()
