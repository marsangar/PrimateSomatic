# ============================================
# Title: Primate drivers script
# Author: Martín Santamarina García (ms84@sanger.ac.uk)
# Description: This script implements dN/dS analysis to study positive selection and detect driver candidates involved in clonal expansions
# ============================================

### Load configs
source("/Users/ms84/Research/postdoc/projects/PrimateSomatic/Primate_master.r")

library(dplyr)
library(ggplot2)
library(ggrepel)


### Set samples to exclude because of contamination
dataMETA$sample
excluded_samples_tier1<-dataKRAKEN_macaca[dataKRAKEN_macaca$kraken_percent<40,"sample"]
excluded_samples_tier2<-dataKRAKEN_macaca[dataKRAKEN_macaca$kraken_percent<60,"sample"]

### Load Mutation Data
dataSOM<-read.delim(paste0(RESDIR_VARIANTS,"/somatic_variants.txt"))

### Load Genes in the panel 
dataPANEL<-read.delim(PANEL)
panelgenes<-dataPANEL$gene
panelgenes_filt<-panelgenes[!panelgenes %in% c("H3F3A", "RHOB", "NFKBIZ", "HIST1H3B", "HIST1H2BD", "HLA-A", "HLA-B", "PIM1", "CDKN2A.p16INK4a", "CDKN2A.p14arf", "FAS", "KLF5", "ERCC5", "HOXB3", "H3F3B", "MEF2B", "CEBPA", "PAK7", "FAM58A")]
  
### Load duplex coverage info
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Directory containing all gene coverage files
#dataDC<-read.delim("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.MQD0001d_tds0001/results.MQD0001d_tds0001_gene_coverage.txt", header = FALSE)
indir <- "~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/"

# Find all coverage files
files <- list.files(
  indir,
  pattern = "gene_coverage\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

# Read all files
coverage <- lapply(files, function(f){
  
  sample <- basename(dirname(f))   # e.g. results.MQD0001d_tds0001
  
  read.delim(f, header = FALSE,
             col.names = c("Gene","Coverage")) %>%
    mutate(Sample = sample)
  
}) %>%
  bind_rows()

head(coverage)

### Convert to a gene x sample matrix
cov.mat <- coverage %>%
  pivot_wider(names_from = Sample,
              values_from = Coverage) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

mat <- t(log2(cov.mat + 1))   # transpose so genes are on x-axis


pdf(paste0(RESDIR_QC,"/GeneCoverage_ComplexHeatmap.pdf"),
    width = 20,
    height = 8)

Heatmap(mat,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        column_names_rot = 90,
        column_names_gp = grid::gpar(fontsize = 5),
        row_names_gp = grid::gpar(fontsize = 10), 
        column_title = "Mean Duplex Coverage across Genes (exon-level profiling)",
        row_title = "Samples",
        )

dev.off()

## calculate the mean coverage for each gene across all samples with:
gene.mean <- coverage %>%
  group_by(Gene) %>%
  summarise(mean_coverage = mean(Coverage, na.rm = TRUE)) %>%
  arrange(desc(mean_coverage))
gene.mean.dc <- setNames(gene.mean$mean_coverage, gene.mean$Gene)

head(gene.mean.dc)


### Run dndscv on somatic callset
dndsout = dndscv(dataSOM, refdb=REFCDS, cv=NULL, dc=gene.mean.dc, gene_list=names(gene.mean.dc))
dndsout_value<-signif(dndsout$globaldnds["wall","mle"],2)

### Run dndscv on somatic callset for each sample
#list_dndsout <- lapply(
#  split(dataSOM, dataSOM$sample),
#  function(df_sample) {
#    dndscv(df_sample, refdb = REFCDS, cv = NULL)
#  }
#)

##
#dc=Duplex coverage per gene. Named Numeric Vector with values reflecting the mean duplex coverage per site per gene, and names corresponding to gene names. Use this argument only when running dNdScv on duplex sequencing data to use gene coverage in the offset of the regression model (default = NULL)

### Run dndscv on somatic callset for each sample
list_dndsout <- lapply(
  split(dataSOM, dataSOM$sample),
  function(df_sample) {
    dndscv(df_sample, refdb = REFCDS, cv = NULL)
  }
)


###
str(list_dndsout)
unique(list_dndsout$MQD0001d$sel_cv$gene_name)


## Select dnds_cv 
df<-dndsout$sel_cv[,c("gene_name", "n_syn", "n_mis", "n_non", "n_spl", "n_ind")]
dim(df)
### restrict to genes in the panel
df<-df[df$gene_name %in% dataPANEL$gene,]
dim(df)


# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = starts_with("n_"),
    names_to = "mutation_type",
    values_to = "count"
  )

# Optional: order genes by total mutation burden
gene_order <- df %>%
  mutate(total = n_syn + n_mis + n_non + n_spl + n_ind) %>%
  arrange(desc(total)) %>%
  pull(gene_name)

df_long$gene_name <- factor(df_long$gene_name, levels = gene_order)

# Rename mutation classes for prettier legend
df_long$mutation_type <- recode(
  df_long$mutation_type,
  n_syn = "Synonymous",
  n_mis = "Missense",
  n_non = "Nonsense",
  n_spl = "Splice",
  n_ind = "Indel"
)

# Stacked barplot
p<-ggplot(df_long, aes(x = gene_name, y = count, fill = mutation_type)) +
  geom_bar(stat = "identity", width = 0.8) +
  #coord_flip() +
  labs(
    title = "Mutation Landscape on the Targeted Panel",
    x = "Gene",
    y = "Mutation Count",
    fill = "Mutation Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 90, size=5, face = "bold", hjust = 1),
    legend.position = "right"
  ) +
  scale_fill_brewer(palette = "Set2")

pdf(paste0(RESDIR_DRIVERS,"/mutation_landscape.pdf"),  width = 18, height = 6)
p
dev.off()

### dndscv outputs: Table of significant genes
sel_cv = dndsout$sel_cv
signif_genes<-sel_cv[sel_cv$qglobal_cv<0.1,]
nrow(signif_genes) # 130 genes
top_signif_genes<-signif_genes[signif_genes$qglobal_cv==0,c(1:10,20)]
top_signif_genes
#write.xlsx(top_signif_genes,"/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/RG25/top_significant_genes.xlsx")





### Driver discovery in targeted sequencing data

#panel_dnds = dndscv(dataSOM, gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
#panel_dnds = dndscv(dataSOM[!dataSOM$sample %in% excluded_samples_tier1,], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
#panel_dnds = dndscv(dataSOM[!dataSOM$sample %in% excluded_samples_tier2,], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)

panel_dnds = dndscv(dataSOM, gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, dc=gene.mean.dc)


panel_dnds$genemuts
sel_cv = panel_dnds$sel_cv
print(sel_cv[sel_cv$qglobal_cv<0.1,c(1:10,19)], digits = 3)

### Plot candidates
df <- as.data.frame(panel_dnds$genemuts)

signif_genes<-sel_cv[sel_cv$qglobal_cv<0.1,]
signif_genes

# Total observed and expected nonsynonymous burden
df <- df %>%
  mutate(
    obs_nonsyn = n_mis + n_non + n_spl,
    exp_nonsyn = exp_mis + exp_non + exp_spl,
    enrichment = obs_nonsyn / exp_nonsyn,
    total_mut = n_syn + obs_nonsyn
  )

# Label strongest enrichments
df$label <- ifelse(df$enrichment > 2 & df$obs_nonsyn >= 5, df$gene_name,"")

p<-ggplot(df, aes(x = exp_nonsyn, y = obs_nonsyn)) +
  
  geom_point(
    aes(size = total_mut,
        color = enrichment),
    alpha = 0.8
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  geom_text_repel(
    aes(label = label),
    size = 3.5,
    max.overlaps = Inf
  ) +
  
  scale_x_log10() +
  scale_y_log10() +
  
  scale_color_viridis_c(
    option = "magma",
    name = "Observed / Expected"
  ) +
  
  scale_size_continuous(
    range = c(2, 10),
    name = "Total mutations"
  )  +
  
  labs(
    title = "Observed vs Expected Nonsynonymous Mutations",
    x = "Expected nonsynonymous mutations",
    y = "Observed nonsynonymous mutations"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
#pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected.pdf"),  width = 10, height = 6)
#pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected_filt_tier1.pdf"),  width = 10, height = 6)
#pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected_filt_tier2.pdf"),  width = 10, height = 6)

pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected_dc.pdf"),  width = 10, height = 6)
p
dev.off()


#tissue="Liver"
### Driver discovery in targeted sequencing data 
for(tissue in c("Bladder","Liver", "Oesophagus", "Prostate", "Skin")){
  print(paste0("driver analysis for ", tissue))

  #panel_dnds = dndscv(dataSOM[dataSOM$tissue == tissue,], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
  #panel_dnds = dndscv(dataSOM[dataSOM$tissue == tissue & !dataSOM$sample %in% excluded_samples_tier1,], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
  #panel_dnds = dndscv(dataSOM[dataSOM$tissue == tissue & !dataSOM$sample %in% excluded_samples_tier2,], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
  
  panel_dnds = dndscv(dataSOM[dataSOM$tissue == tissue,], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf,  dc=gene.mean.dc)
  
panel_dnds$genemuts
sel_cv = panel_dnds$sel_cv
print(sel_cv[sel_cv$qglobal_cv<0.1,c(1:10,19)], digits = 3)

### Plot candidates


df <- as.data.frame(panel_dnds$genemuts)

# Total observed and expected nonsynonymous burden
df <- df %>%
  mutate(
    obs_nonsyn = n_mis + n_non + n_spl,
    exp_nonsyn = exp_mis + exp_non + exp_spl,
    enrichment = obs_nonsyn / exp_nonsyn,
    total_mut = n_syn + obs_nonsyn
  )

# Label strongest enrichments
df$label <- ifelse(df$enrichment > 2 & df$obs_nonsyn >= 2, df$gene_name, "")

p<-ggplot(df, aes(x = exp_nonsyn, y = obs_nonsyn)) +
  
  geom_point(
    aes(size = total_mut,
        color = enrichment),
    alpha = 0.8
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  geom_text_repel(
    aes(label = label),
    size = 3.5,
    max.overlaps = Inf
  ) +
  
  scale_x_log10() +
  scale_y_log10() +
  
  scale_color_viridis_c(
    option = "magma",
    name = "Observed / Expected"
  ) +
  
  scale_size_continuous(
    range = c(2, 10),
    name = "Total mutations"
  )  +
  
  labs(
    title = paste0("Observed vs Expected Nonsynonymous Mutations (", tissue, ")"),
    x = "Expected nonsynonymous mutations",
    y = "Observed nonsynonymous mutations"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
#pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected_", tissue, ".pdf"),  width = 10, height = 6)
pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected_", tissue, "_dc.pdf"),  width = 10, height = 6)
print(p)
dev.off()
}


########################
### Driver plotting ####
########################

driver_plot_per_donor(dndsout$annotmuts)
driver_plot_per_donor
# source a specific function
source("~/volumes/ms84_software/team268/tools/R/driver_plot.R")
source("~/volumes/ms84_software/team268/tools/R/driver_plot_per_donor.R")


#### source every function in the shared tools directory
#list.files("~/volumes/ms84_software/team268/tools/R", full.names = TRUE, pattern = "*.R") %>%
#  purrr::walk(source)


dnds

muts <- dndsout$annotmuts |>
  dplyr::transmute(
    gene,
    sampleID,
    donor,
    seq_type,
    chr,
    pos,
    ref,
    mut,
    impact,
    duplex_vaf,
    cellfraction
  )



driver_plot
driver_plot(genes2plot=dndsout$genemuts$gene_name,
            muts = dndsout$annotmuts,
            coverage_df = gene.mean.dc)
driver_plot





# a trimmed-down real data slice (~50 driver genes, all 198 samples) extracted
# from ibd/reports/20260219_dndscv/dndscv_cache
demo_data <- readRDS("data/demo_data.rds")
list2env(demo_data, envir = environment())






### 19/10/2026

### Load Genes in the panel 
dataPANEL<-read.delim("/Users/ms84/Downloads/Sanger_TERT-v4_TE-95148282_hg19_highstringencyfilter_buccal_gene_list.tsv")
panelgenes<-dataPANEL$gene
panelgenes_filt<-panelgenes[!panelgenes %in% c("H3F3A", "RHOB", "NFKBIZ", "HIST1H3B", "HIST1H2BD", "HLA-A", "HLA-B", "PIM1", "CDKN2A.p16INK4a", "CDKN2A.p14arf", "FAS", "KLF5", "ERCC5", "HOXB3", "H3F3B", "MEF2B", "CEBPA", "PAK7", "FAM58A")]


#### Driver discovery for the ASHG annual meeting 2026 - 

dataSOM<-read.delim("/Users/ms84/Research/postdoc/projects/PrimateSomatic/results/v0.2/Macaca_mulatta/targeted_nanoseq/variants/somatic_variants.txt")
dim(dataSOM)

### without dc correction
panel_dnds = dndscv(dataSOM, gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
panel_dnds$genemuts

panel_dnds_liver = dndscv(dataSOM[dataSOM$tissue %in% "Liver",], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
panel_dnds_bladder = dndscv(dataSOM[dataSOM$tissue %in% "Bladder",], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
panel_dnds_oesophagus = dndscv(dataSOM[dataSOM$tissue %in% "Oesophagus",], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
panel_dnds_prostate = dndscv(dataSOM[dataSOM$tissue %in% "Prostate",], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
panel_dnds_skin = dndscv(dataSOM[dataSOM$tissue %in% "Skin",], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)



sel_cv = panel_dnds$sel_cv

sel_cv_liver = panel_dnds_liver$sel_cv
sel_cv_bladder = panel_dnds_bladder$sel_cv
sel_cv_oesophagus= panel_dnds_oesophagus$sel_cv
sel_cv_prostate= panel_dnds_prostate$sel_cv
sel_cv_skin = panel_dnds_skin$sel_cv






### Plot candidates
df <- as.data.frame(panel_dnds$genemuts)

signif_genes<-sel_cv[sel_cv$qglobal_cv<0.1,]
signif_genes


###
sel_cv[sel_cv$qmis_cv<0.1,]
sel_cv[sel_cv$qtrunc_cv<0.1,]
sel_cv[sel_cv$qind_cv<0.1,]
sel_cv[sel_cv$qglobal_cv<0.1,]

dim(sel_cv)



#### Set ranked objects
rank_qglobal_cv<-sel_cv[order(sel_cv$qglobal_cv, decreasing = FALSE),]
rank_qmis_cv<-sel_cv[order(sel_cv$qmis_cv, decreasing = FALSE),]
rank_qtrunc_cv<-sel_cv[order(sel_cv$qtrunc_cv, decreasing = FALSE),]
rank_qind_cv<-sel_cv[order(sel_cv$qind_cv, decreasing = FALSE),]



pdf("/Users/ms84/Research/postdoc/projects/PrimateSomatic/results/v0.2/Macaca_mulatta/targeted_nanoseq/drivers/selection_landscape_macaque.pdf")

plot("", xlim=c(0,5), ylim=c(0.5,0), main="Selection landscape in Rhesus macaque", xlab="", ylab="q value", axes=FALSE)
axis(side= 1, at=c(1,2,3,4), labels = c("global", "missense", "truncation", "indel"))
axis(side= 2)

rect(xleft = 0.95, xright = 1.05, ybottom = 1, ytop = 0, col = "black")
rect(xleft = 1.95, xright = 2.05, ybottom = 1, ytop = 0, col = "black")
rect(xleft = 2.95, xright = 3.05, ybottom = 1, ytop = 0, col = "black")
rect(xleft = 3.95, xright = 4.05, ybottom = 1, ytop = 0, col = "black")

rampalette<-colorRampPalette(c("yellow", "red"))(1000)

### plot values
points(x=rep(1, nrow(rank_qglobal_cv)), y=rank_qglobal_cv$qglobal_cv, pch=21, bg=rampalette, cex=2)
points(x=rep(2, nrow(rank_qmis_cv)), y=rank_qmis_cv$qmis_cv, pch=21, bg=rampalette, cex=2)
points(x=rep(3, nrow(rank_qtrunc_cv)), y=rank_qtrunc_cv$qtrunc_cv, pch=21, bg=rampalette, cex=2)
points(x=rep(4, nrow(rank_qind_cv)), y=rank_qind_cv$qind_cv, pch=21, bg=rampalette, cex=2)

### plot gene labels
n_labels<-5
text(rank_qglobal_cv$gene_name[1:n_labels], x=rep(1.3, n_labels), y=rank_qglobal_cv$qglobal_cv[1:n_labels], cex = 0.7)
text(rank_qmis_cv$gene_name[1:n_labels], x=rep(2.3, n_labels), y=rank_qmis_cv$qmis_cv[1:n_labels], cex = 0.7)
text(rank_qtrunc_cv$gene_name[1:n_labels], x=rep(3.3, n_labels), y=rank_qtrunc_cv$qtrunc_cv[1:n_labels], cex = 0.7)
text(rank_qind_cv$gene_name[1:n_labels], x=rep(4.3, n_labels), y=rank_qind_cv$qind_cv[1:n_labels], cex = 0.7)

dev.off()





### Explore VAF:
dataSOM[dataSOM$DUPLEX_VAF, "DUP"]




# Total observed and expected nonsynonymous burden
df <- df %>%
  mutate(
    obs_nonsyn = n_mis + n_non + n_spl,
    exp_nonsyn = exp_mis + exp_non + exp_spl,
    enrichment = obs_nonsyn / exp_nonsyn,
    total_mut = n_syn + obs_nonsyn
  )

# Label strongest enrichments
df$label <- ifelse(df$enrichment > 2 & df$obs_nonsyn >= 5, df$gene_name,"")

p<-ggplot(df, aes(x = exp_nonsyn, y = obs_nonsyn)) +
  
  geom_point(
    aes(size = total_mut,
        color = enrichment),
    alpha = 0.8
  ) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  geom_text_repel(
    aes(label = label),
    size = 3.5,
    max.overlaps = Inf
  ) +
  
  scale_x_log10() +
  scale_y_log10() +
  
  scale_color_viridis_c(
    option = "magma",
    name = "Observed / Expected"
  ) +
  
  scale_size_continuous(
    range = c(2, 10),
    name = "Total mutations"
  )  +
  
  labs(
    title = "Observed vs Expected Nonsynonymous Mutations",
    x = "Expected nonsynonymous mutations",
    y = "Observed nonsynonymous mutations"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
#pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected.pdf"),  width = 10, height = 6)
#pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected_filt_tier1.pdf"),  width = 10, height = 6)
#pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected_filt_tier2.pdf"),  width = 10, height = 6)

pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected_dc.pdf"),  width = 10, height = 6)
p
dev.off()


