### Load configs
source("/Users/ms84/Research/postdoc/projects/PrimateSomatic/Primate_master.r")

### Load vcf files
dataVCF<-load_vcf_dataset(VCF_DIR)
dataVCF$sample<-str_split_fixed(dataVCF$sample, "_", 2)[,1]

### Set vcf core fields
dataVCF$CHROM<-as.character(dataVCF$seqnames)
dataVCF$POS<-as.integer(dataVCF$start)
dataVCF$REF<-toupper(dataVCF$REF)
dataVCF$ALT<-toupper(dataVCF$ALT)


dataVCF$DPLX_NM <- sapply(dataVCF$DPLX_NM, function(v) if (length(v) > 0) v[1] else NA_integer_)

dataVCF$variant_id<-paste0(dataVCF$CHROM,"_",dataVCF$POS)
dataVCF$variant_id_complete<-paste0(dataVCF$CHROM,"_",dataVCF$POS, "_", dataVCF$REF, ">", dataVCF$ALT)


### Load metadata
dataMETA<-read.xlsx(METADATA)

### Annotate VCFs with metadata
dataVCF_annot <- dataVCF %>% dplyr::left_join(dataMETA, by = "sample")
dataVCF_annot

### Sort VCFs based on tissue
dataVCF_sorted<- dataVCF_annot %>% dplyr::arrange(tissue,sample,CHROM,POS)

### Determine if variants are shared across samples and individuals
dataVCF_sorted$shared_across_samples<-apply(dataVCF_sorted,1,function(x) ifelse(x["variant_id"] %in% dataVCF_sorted[!dataVCF_sorted$sample %in% x["sample"],"variant_id"], TRUE,FALSE))
dataVCF_sorted$shared_across_individuals<-apply(dataVCF_sorted,1,function(x) ifelse(x["variant_id"] %in% dataVCF_sorted[!dataVCF_sorted$individual %in% x["individual"],"variant_id"], TRUE,FALSE))

table(dataVCF_sorted$shared_across_samples)
table(dataVCF_sorted$shared_across_individual)

### Plot QC metrics
pdf(paste0(RESDIR_QC,"/metrics.pdf"), 12, 6)

### BAM_COV
plot_qc_metric(dataVCF_sorted, "BAM_COV", plot_type = "histogram", group_by = "all", color_by = "tissue", tissue_cols = tissue_cols)
plot_qc_metric(dataVCF_sorted, "BAM_COV", plot_type = "violin", group_by = "sample", color_by = "tissue", tissue_cols = tissue_cols)

### DUPLEX_COV
plot_qc_metric(dataVCF_sorted, "DUPLEX_COV", plot_type = "histogram", group_by = "all", color_by = "tissue", tissue_cols = tissue_cols)
plot_qc_metric(dataVCF_sorted, "DUPLEX_COV", plot_type = "violin", group_by = "sample", color_by = "tissue", tissue_cols = tissue_cols)

### BAM_VAF
plot_qc_metric(dataVCF_sorted, "BAM_VAF", plot_type = "density", group_by = "all", color_by = "tissue", tissue_cols = tissue_cols, log10_scale = TRUE)
plot_qc_metric(dataVCF_sorted, "BAM_VAF", plot_type = "density", group_by = "sample", color_by = "tissue", tissue_cols = tissue_cols, log10_scale = TRUE)

### DUPLEX_VAF
plot_qc_metric(dataVCF_sorted, "DUPLEX_VAF", plot_type = "density", group_by = "all", color_by = "tissue", tissue_cols = tissue_cols, log10_scale = TRUE)
plot_qc_metric(dataVCF_sorted, "DUPLEX_VAF", plot_type = "density", group_by = "sample", color_by = "tissue", tissue_cols = tissue_cols, log10_scale = TRUE)

dev.off()


### Perfom Somatic Filtering
dim(dataVCF)
dataSOM<-dataVCF_sorted[dataVCF_sorted$TYPE %in% c("snv", "del", "ins") &  dataVCF_sorted$shared_across_samples==FALSE & dataVCF_sorted$DPLX_NM %in% c(1, NA) & dataVCF_sorted$FILTER %in% "PASS",c("sample","CHROM", "POS", "REF", "ALT","age", "tissue", "TIMES_CALLED", "TYPE", "variant_id","variant_id_complete", "species", "shared_across_samples","shared_across_individuals", "FILTER", "DPLX_NM","DUPLEX_COV","BAM_COV", "DUPLEX_VAF", "BAM_VAF")]
dataSOM$DPLX_NM
dim(dataSOM) #
table(dataSOM$TYPE)
table(dataSOM$TYPE)


### Evaluate impact of filtering on dnds
dndsout = dndscv(dataSOM, refdb=REFCDS, cv=NULL)
dndsout_value<-signif(dndsout$globaldnds["wall","mle"],2)

list_dndsout <- lapply(
  split(dataSOM, dataSOM$sample),
  function(df_sample) {
    dndscv(df_sample, refdb = REFCDS, cv = NULL)
  }
)


globaldnds <- bind_rows(
  lapply(names(list_dndsout), function(sid) {
    res <- list_dndsout[[sid]]
    if (is.null(res$globaldnds)) return(NULL)
    res$globaldnds %>%
      as.data.frame() %>%
      mutate(sample = sid)
    
  })
)


### add tissue metadata
sample_annot <- dataSOM %>% dplyr::select(sample, tissue, age) %>% distinct() %>% dplyr::rename(sample=sample)

### join with your dN/dS table
globaldnds <- globaldnds %>% left_join(sample_annot, by = "sample")

## order samples by tissue (and age)
globaldnds <- globaldnds %>%
  arrange(tissue, age) %>%
  mutate(sample = factor(sample, levels = unique(sample)))

## Define shapes for dN/dS types (highlight wall)
shape_values <- c(
  "wmis" = 21,  # circle
  "wnon" = 24,  # triangle
  "wspl" = 22,  # square
  "wtru" = 23,  # diamond
  "wall" = 8    # keep star (no fill, but already distinct)
)

## Modified plot
pd <- position_dodge(width = 0.6)

p1 <- ggplot(globaldnds, aes(x = sample, y = mle)) +
  
  ## Error bars FIRST (background)
  geom_errorbar(
    aes(ymin = cilow, ymax = cihigh, color = tissue, group = name),
    width = 0.2,
    position = pd
  ) +
  
  ## Points SECOND (on top)
  geom_point(
    aes(fill = tissue, shape = name, group = name),
    color = "black",
    position = pd,
    size = 3,
    stroke = 0.8
  ) +
  
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  
  coord_cartesian(ylim = c(0, 3)) +
  
  scale_fill_manual(values = tissue_cols) +
  scale_color_manual(values = tissue_cols) +
  scale_shape_manual(values = shape_values) +
  
  labs(
    title = "Global dN/dS (w) across primate samples",
    subtitle = paste0("Global dN/dS: ", dndsout_value),
    x = "Sample",
    y = "w (dN/dS ratio)",
    fill = "Tissue",
    shape = "dN/dS type"
  ) +
  
  guides(color = "none") +  # remove duplicate legend
  
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16)
  )
p1
pdf(paste0(RESDIR_QC,"/globaldnds.pdf"),  width = 12, height = 7)
grid.arrange(p1)
dev.off()

### dndscv outputs: Table of significant genes
sel_cv = dndsout$sel_cv
signif_genes<-sel_cv[sel_cv$qglobal_cv<0.1,]
nrow(signif_genes) # 130 genes
top_signif_genes<-signif_genes[signif_genes$qglobal_cv==0,c(1:10,20)]
top_signif_genes
write.xlsx(top_signif_genes,"/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/RG25/top_significant_genes.xlsx")



### write Somatic Variants
write.table(dataSOM, file=paste0(RESDIR_VARIANTS,"/somatic_variants.txt"), quote = FALSE, sep="\t", row.names = FALSE)






# Load libraries
library(tidyverse)

# Example data
df <- tribble(
  ~gene_name, ~n_syn, ~n_mis, ~n_non, ~n_spl, ~n_ind,
  "ARID2", 7,13,3,1,3,
  "KMT2C", 4,19,4,0,8,
  "MGA", 6,18,2,1,5,
  "TET2", 7,14,2,0,3,
  "ERBB4",12,15,1,1,0,
  "ROBO2",9,14,0,2,1,
  "NOTCH3",5,15,1,1,4,
  "CHD4",6,13,0,2,1,
  "FAT1",6,31,0,2,6,
  "KMT2D",8,25,2,0,6,
  "KMT2E",4,14,1,0,2,
  "ZFHX3",14,19,1,0,4,
  "FAT2",8,22,0,1,5,
  "FAT4",13,26,1,0,2,
  "ATM",5,14,1,0,3,
  "MTOR",6,15,0,0,2,
  "APC",3,18,0,0,4,
  "PREX2",7,16,0,0,0,
  "MET",2,12,0,0,1,
  "KDM6A",7,11,1,3,4,
  "CBLB",2,8,1,3,1,
  "CUL3",1,7,3,1,1,
  "NF1",7,13,0,1,6
)

df<-dndsout$sel_cv[,c("gene_name", "n_syn", "n_mis", "n_non", "n_spl", "n_ind")]
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
ggplot(df_long, aes(x = gene_name, y = count, fill = mutation_type)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  labs(
    title = "Mutation Counts per Gene",
    x = "Gene",
    y = "Mutation Count",
    fill = "Mutation Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  ) +
  scale_fill_brewer(palette = "Set2")
