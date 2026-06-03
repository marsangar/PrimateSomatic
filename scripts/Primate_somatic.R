# ============================================
# Title: Primate somatic script
# Author: Martín Santamarina García (ms84@sanger.ac.uk)
# Description: This script implements qc and filtering steps to retrieve somatic variant calls
# ============================================


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

dataVCF$chr<-rename_chroms(dataVCF$CHROM, from = "RefSeq.seq.accession", to="Sequence.name", ncbi_report = NCBI_REPORT)


dataVCF$DPLX_NM <- sapply(dataVCF$DPLX_NM, function(v) if (length(v) > 0) v[1] else NA_integer_)

dataVCF$variant_id<-paste0(dataVCF$CHROM,"_",dataVCF$POS)
dataVCF$variant_id_complete<-paste0(dataVCF$CHROM,"_",dataVCF$POS, "_", dataVCF$REF, ">", dataVCF$ALT)


### Load metadata
dataMETA<-read.xlsx(METADATA)

### Load kraken results
dataKRAKEN<-read.delim(KRAKEN, header = FALSE)
colnames(dataKRAKEN)[2]<-c("kraken_percent")
dataKRAKEN$sample<-str_split_fixed(dataKRAKEN$V1,"analysis/|_",3)[,2]

### Annotate VCFs with metadata
dataVCF_annot <- dataVCF %>% dplyr::left_join(dataMETA, by = "sample")
dataVCF_annot <- dataVCF_annot %>% dplyr::left_join(dataKRAKEN[,c("sample","kraken_percent")], by = "sample")

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


### Load germline variants based according to pile-up
dataGERMLINE<-read.delim(GERMLINE)
dataGERMLINE$variant_id<-paste0(dataGERMLINE$NC_CHROM, "_", dataGERMLINE$POS)


### Perfom Somatic Filtering
dataSOM<-dataVCF_sorted[dataVCF_sorted$TYPE %in% c("snv", "del", "ins") &  dataVCF_sorted$shared_across_samples==FALSE & dataVCF_sorted$DPLX_NM %in% c(1, NA) & dataVCF_sorted$FILTER %in% "PASS" & (!dataVCF_sorted$variant_id %in% dataGERMLINE$variant_id),c("sample","CHROM", "POS", "REF", "ALT","age", "tissue", "TIMES_CALLED", "TYPE", "variant_id","variant_id_complete", "species", "shared_across_samples","shared_across_individuals", "FILTER", "DPLX_NM","DUPLEX_COV","BAM_COV", "DUPLEX_VAF", "BAM_VAF", "kraken_percent", "chr")]

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
sample_annot <- dataSOM %>% dplyr::select(sample, tissue, age, kraken_percent) %>% distinct() %>% dplyr::rename(sample=sample)

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



### Correlation between dnds & kraken percent (proxy of purity)

### Consider all variants (wall)
globaldnds_wall<-globaldnds[globaldnds$name=="wall",]

# Correlation test
cor_res <- cor.test(
  globaldnds_wall$mle,
  globaldnds_wall$kraken_percent,
  method = "pearson"
)

# Correlation label
cor_label <- paste0(
  "r = ", round(cor_res$estimate, 3),
  "\np = ", signif(cor_res$p.value, 3)
)

# Set threshold value 
threshold <- 50 ### 

# Plot
p<-ggplot(
  globaldnds_wall,
  aes(
    x = mle,
    y = kraken_percent,
    color = tissue  # column containing tissue labels
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  
  # Add sample labels
  geom_text_repel(
    aes(label = sample),   # replace with your sample-name column
    size = 3.5,
    max.overlaps = Inf
  ) +
  
  # Linear regression + 95% CI
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    fill = "gray70",
    alpha = 0.2,
    linewidth = 1
  ) + 
  
  # Horizontal threshold line
  #geom_hline(yintercept = threshold,linetype = "dashed",color = "red",linewidth = 0.3) +
  
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = cor_label,
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  
  labs(
    x = "dN/dS (MLE)",
    y = "Kraken percent macaque",
    color = "Tissue type",
    title = "Correlation between dN/dS and sample purity"
  ) +
  
  theme_classic(base_size = 14)

pdf(paste0(RESDIR_QC,"/kraken_vs_dnds.pdf"),  width = 9, height = 5)
grid.arrange(p)
dev.off()

### write Somatic Variants
write.table(dataSOM, file=paste0(RESDIR_VARIANTS,"/somatic_variants.txt"), quote = FALSE, sep="\t", row.names = FALSE)


### write Somatic Variants annotated with dndscv
write.table(dndsout$annotmuts,paste0(RESDIR_VARIANTS, "/annotated_mutations.txt"), quote = FALSE, sep="\t", row.names = FALSE)


