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
excluded_samples_tier1<-dataKRAKEN[dataKRAKEN$kraken_percent<40,"sample"]
excluded_samples_tier2<-dataKRAKEN[dataKRAKEN$kraken_percent<60,"sample"]

### Load Mutation Data
dataSOM<-read.delim(paste0(RESDIR_VARIANTS,"/somatic_variants.txt"))

### Load Genes in the panel 
dataPANEL<-read.delim(PANEL)
panelgenes<-dataPANEL$gene
panelgenes_filt<-panelgenes[!panelgenes %in% c("H3F3A", "RHOB", "NFKBIZ", "HIST1H3B", "HIST1H2BD", "HLA-A", "HLA-B", "PIM1", "CDKN2A.p16INK4a", "CDKN2A.p14arf", "FAS", "KLF5", "ERCC5", "HOXB3", "H3F3B", "MEF2B", "CEBPA", "PAK7", "FAM58A")]
  
### Run dndscv on somatic callset
dndsout = dndscv(dataSOM, refdb=REFCDS, cv=NULL)
dndsout_value<-signif(dndsout$globaldnds["wall","mle"],2)

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
write.xlsx(top_signif_genes,"/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/RG25/top_significant_genes.xlsx")






### Driver discovery in targeted sequencing data

#panel_dnds = dndscv(dataSOM, gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
#panel_dnds = dndscv(dataSOM[!dataSOM$sample %in% excluded_samples_tier1,], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
panel_dnds = dndscv(dataSOM[!dataSOM$sample %in% excluded_samples_tier2,], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)


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
pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected_filt_tier2.pdf"),  width = 10, height = 6)
p
dev.off()


#tissue="Liver"
### Driver discovery in targeted sequencing data 
for(tissue in c("Bladder","Liver", "Oesophagus", "Prostate", "Skin")){
  print(paste0("driver analysis for ", tissue))

  #panel_dnds = dndscv(dataSOM[dataSOM$tissue == tissue,], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
  #panel_dnds = dndscv(dataSOM[dataSOM$tissue == tissue & !dataSOM$sample %in% excluded_samples_tier1,], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
  panel_dnds = dndscv(dataSOM[dataSOM$tissue == tissue & !dataSOM$sample %in% excluded_samples_tier2,], gene_list=panelgenes_filt, refdb = REFCDS, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)

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
pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected_", tissue, "_filt_tier2.pdf"),  width = 10, height = 6)
#pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected_", tissue, "_filt_tier1.pdf"),  width = 10, height = 6)
#pdf(paste0(RESDIR_DRIVERS,"/observed_vs_expected_", tissue, ".pdf"),  width = 10, height = 6)

print(p)
dev.off()
}




