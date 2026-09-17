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

### Load kraken results & append to metadata
dataKRAKEN_macaca<-read.delim(paste0(KRAKEN_DIR, "KRAKEN.RESULTS_Macaca_mulatta.tsv"), header = FALSE)
dataKRAKEN_homo<-read.delim(paste0(KRAKEN_DIR, "KRAKEN.RESULTS_Homo_sapiens.tsv"), header = FALSE)
dataKRAKEN_mus<-read.delim(paste0(KRAKEN_DIR, "KRAKEN.RESULTS_Mus_musculus.tsv"), header = FALSE)
dataKRAKEN_unclassified<-read.delim(paste0(KRAKEN_DIR, "KRAKEN.RESULTS_unclassified_ok.tsv"), header = FALSE)

colnames(dataKRAKEN_macaca)[2]<-c("kraken_macaca")
colnames(dataKRAKEN_homo)[2]<-c("kraken_homo")
colnames(dataKRAKEN_mus)[2]<-c("kraken_mus")
colnames(dataKRAKEN_unclassified)[2]<-c("kraken_unclassified")

dataKRAKEN_macaca$sample<-str_split_fixed(dataKRAKEN_macaca$V1,"analysis/|_",3)[,2]
dataKRAKEN_homo$sample<-str_split_fixed(dataKRAKEN_homo$V1,"analysis/|_",3)[,2]
dataKRAKEN_mus$sample<-str_split_fixed(dataKRAKEN_mus$V1,"analysis/|_",3)[,2]
dataKRAKEN_unclassified$sample<-str_split_fixed(dataKRAKEN_unclassified$V1,"analysis/|_",3)[,2]

dataMETA<-dataMETA %>% dplyr::left_join(dataKRAKEN_macaca[,c("sample","kraken_macaca")], by = "sample")
dataMETA<-dataMETA %>% dplyr::left_join(dataKRAKEN_homo[,c("sample","kraken_homo")], by = "sample")
dataMETA<-dataMETA %>% dplyr::left_join(dataKRAKEN_mus[,c("sample","kraken_mus")], by = "sample")
dataMETA<-dataMETA %>% dplyr::left_join(dataKRAKEN_unclassified[,c("sample","kraken_unclassified")], by = "sample")
dataMETA$kraken_other<-100-(dataMETA$kraken_macaca+dataMETA$kraken_homo+dataMETA$kraken_mus+dataMETA$kraken_unclassified)
dataMETA$kraken_total<-dataMETA$kraken_macaca+dataMETA$kraken_homo+dataMETA$kraken_mus+dataMETA$kraken_unclassified+dataMETA$kraken_other

### Load clean_ns results & append to metadata
dataCLEAN_fastq<-read.delim(paste0(CLEAN_NS_DIR, "/fastq/summary_Macaca_mulatta.tsv"), header = TRUE)
dataCLEAN_bwa<-read.delim(paste0(CLEAN_NS_DIR, "/bwa/samtools_stats_summary_Macaca_mulatta.tsv"), header = TRUE)
dataCLEAN_filtered_macaca<-read.delim(paste0(CLEAN_NS_DIR, "/filtered_bams/summary_Macaca_mulatta.tsv"), header = TRUE)
dataCLEAN_filtered_homo<-read.delim(paste0(CLEAN_NS_DIR, "/filtered_bams/summary_Homo_sapiens.tsv"), header = TRUE)
colnames(dataCLEAN_filtered_macaca)<-c("sample","clean_ns_reads_Macaca_mulatta")
colnames(dataCLEAN_filtered_homo)<-c("sample","clean_ns_reads_Homo_sapiens")

dataMETA<-dataMETA %>% dplyr::left_join(dataCLEAN_fastq, by = c("sample"))
dataMETA<-dataMETA %>% dplyr::left_join(dataCLEAN_bwa, by = c("sample"))
dataMETA<-dataMETA %>% dplyr::left_join(dataCLEAN_filtered_macaca, by = c("sample"))
dataMETA<-dataMETA %>% dplyr::left_join(dataCLEAN_filtered_homo, by = c("sample"))

### Summary table
dataMETA[,c("sample", "fastq_reads", "raw_total_sequences", "reads_mapped", "reads_unmapped",  "clean_ns_reads_Macaca_mulatta", "clean_ns_reads_Homo_sapiens")]


### Get percents
dataMETA$clean_ns_reads_Macaca_mulatta_perc<-dataMETA$clean_ns_reads_Macaca_mulatta/(dataMETA$clean_ns_reads_Macaca_mulatta+dataMETA$clean_ns_reads_Homo_sapiens)
dataMETA$clean_ns_reads_Homo_sapiens_perc<-dataMETA$clean_ns_reads_Homo_sapiens/(dataMETA$clean_ns_reads_Macaca_mulatta+dataMETA$clean_ns_reads_Homo_sapiens)

df <- dataMETA[,c("sample", "reads_mapped", "clean_ns_reads_Macaca_mulatta", "clean_ns_reads_Homo_sapiens")] %>%
#df <- dataMETA[,c("sample","reads_mapped", "clean_ns_reads_Macaca_mulatta_perc", "clean_ns_reads_Homo_sapiens_perc")] %>%
  arrange(desc(reads_mapped))

plot_df <- df %>%
  pivot_longer(
    cols = c(clean_ns_reads_Macaca_mulatta, clean_ns_reads_Homo_sapiens),
    #cols = c(clean_ns_reads_Macaca_mulatta_perc, clean_ns_reads_Homo_sapiens_perc),
    names_to = "category",
    values_to = "reads"
  )

p<-ggplot(plot_df, aes(
  x = factor(sample, levels = df$sample),
  y = reads,
  fill = category
)) +
  geom_col() +
  scale_fill_manual(
    values = c(
      #clean_ns_reads_Macaca_mulatta_perc="green",
      #clean_ns_reads_Homo_sapiens_perc="red"
      clean_ns_reads_Macaca_mulatta="#00BF7D",
      clean_ns_reads_Homo_sapiens="#F8766D"
      #reads_mapped = "steelblue",
      #reads_unmapped = "firebrick",
      #clean_ns_reads = "green"
    )
  ) +
  labs(x = "Sample", y = "Reads", title = "CLEAN_NS competitive mapping from macaque targeted sequencing data") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(RESDIR_QC,"/clean_ns_competitive_barplot.pdf"), 9, 8)
p
dev.off()


### Plot KRAKEN percentages
df <- dataMETA[,c("sample", "kraken_macaca", "kraken_homo", "kraken_mus", "kraken_other", "kraken_unclassified")] %>%
  arrange(desc(kraken_macaca))


plot_df <- df %>%
  pivot_longer(
    cols = c("kraken_macaca", "kraken_homo", "kraken_mus", "kraken_other", "kraken_unclassified"),
    names_to = "category",
    values_to = "percentage"
  )

p<-ggplot(plot_df, aes(
  x = factor(sample, levels = df$sample),
  y = percentage,
  fill = category
)) +
  geom_col() +
  scale_fill_manual(
    values = c(
      kraken_macaca="#00BF7D",
      kraken_homo="#F8766D",
      kraken_mus="#00B0F6",
      kraken_other="#A3A500",
      kraken_unclassified="gray"
    )
  ) +
  labs(x = "Sample", y = "Percentage of reads", title = "KRAKEN taxonomic inference from macaque targeted sequencing data (*mutation-supporting reads)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(RESDIR_QC,"/kraken_stacked.pdf"), 9, 8)
p
dev.off()



### Annotate VCFs with metadata
dataVCF_annot <- dataVCF %>% dplyr::left_join(dataMETA, by = "sample")

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

table(dataVCF_sorted$FILTER)
prop.table(table(dataVCF_sorted$FILTER))

### Perfom Somatic Filtering
#dataSOM<-dataVCF_sorted[dataVCF_sorted$TYPE %in% c("snv", "del", "ins") &  dataVCF_sorted$shared_across_samples==FALSE & dataVCF_sorted$DPLX_NM %in% c(1, NA) & dataVCF_sorted$FILTER %in% "PASS" & (!dataVCF_sorted$variant_id %in% dataGERMLINE$variant_id),c("sample","CHROM", "POS", "REF", "ALT","age", "tissue", "TIMES_CALLED", "TYPE", "variant_id","variant_id_complete", "species", "shared_across_samples","shared_across_individuals", "FILTER", "DPLX_NM","DUPLEX_COV","BAM_COV", "DUPLEX_VAF", "BAM_VAF", "kraken_macaca", "clean_ns_reads_Homo_sapiens_perc", "chr")]
dataSOM<-dataVCF_sorted[dataVCF_sorted$TYPE %in% c("snv", "del", "ins") &  dataVCF_sorted$shared_across_samples==FALSE & dataVCF_sorted$DPLX_NM %in% c(1, NA) & (!dataVCF_sorted$variant_id %in% dataGERMLINE$variant_id),]
dataSOM<-dataVCF_sorted[dataVCF_sorted$TYPE %in% c("snv", "del", "ins") &  dataVCF_sorted$shared_across_samples==FALSE & dataVCF_sorted$DPLX_NM %in% c(1, NA) & dataVCF_sorted$FILTER %in% "PASS" & (!dataVCF_sorted$variant_id %in% dataGERMLINE$variant_id),]

dataSOM <- dataSOM %>% dplyr::select(c("sample","CHROM", "POS", "REF", "ALT"), everything())
head(dataSOM)
dataSOM$DPLX_NM
dim(dataSOM) #
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
sample_annot <- dataSOM %>% dplyr::select(sample, tissue, age, kraken_macaca, clean_ns_reads_Homo_sapiens_perc) %>% distinct() %>% dplyr::rename(sample=sample)

### join with your dN/dS table
globaldnds <- globaldnds %>% left_join(sample_annot, by = "sample")

## order samples by tissue (and age)
globaldnds <- globaldnds %>%
  arrange(tissue, age) %>%
  mutate(sample = factor(sample, levels = unique(sample)))

### check previous results
#globaldnds_old<-globaldnds





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
sample_annot <- dataSOM %>% dplyr::select(sample, tissue, age, kraken_macaca, clean_ns_reads_Homo_sapiens_perc) %>% distinct() %>% dplyr::rename(sample=sample)

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


### 
paste0(globaldnds$name,"_", globaldnds$sample)==paste0(globaldnds_old$name,"_", globaldnds_old$sample)

colnames(globaldnds_old)<-paste0(colnames(globaldnds), "_old")
globaldnds_joined<-cbind(globaldnds, globaldnds_old)
head(globaldnds_joined)
dim(globaldnds_joined)
colnames(globaldnds_joined)

globaldnds_joined_wall<-globaldnds_joined[globaldnds_joined$name=="wall",]

## sort by macaque percentage according to kraken
globaldnds_joined_wall<- globaldnds_joined_wall %>%
  arrange(kraken_macaca)  %>%
  mutate(sample = factor(sample, levels = unique(sample)))

p1 <- ggplot(globaldnds_joined_wall, aes(x = sample)) +
  
  ## Error bars old
  geom_errorbar(
    aes(ymin = cilow_old, ymax = cihigh_old, group = name),
     color = "#F8766D",
    width = 0.2,
    position = pd
  ) +
  
  ## Error bars
  geom_errorbar(
    aes(ymin = cilow, ymax = cihigh, group = name),
    color = "#00BF7D",
    width = 0.2,
    position = pd
  ) +
  
  ## Connect old -> new
  geom_segment(
    aes(
      xend = sample,
      y = mle_old,
      yend = mle,
      group = name
    ),
    color="black",
    position = pd,
    linewidth = 1,
  ) +
  
  ## Points old
  geom_point(
    aes( group = name, y = mle_old),
    shape= 21,
    fill = "#F8766D",
    position = pd,
    size = 5,
    stroke = 0.8
  ) +
  
  ## Points
  geom_point(
    aes(group = name, y=mle),
    shape=21,
    fill ="#00BF7D",
    position = pd,
    size = 5,
    stroke = 0.8
  ) +
  
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  
  coord_cartesian(ylim = c(0, 3)) +
  
  scale_fill_manual(values = tissue_cols) +
  scale_color_manual(values = tissue_cols) +
  scale_shape_manual(values = shape_values) +
  
  labs(
    title = "Decontamination of macaque targeted sequencing",
    subtitle = "Global dN/dS: 0.85 (original) > 1.2 (decontaminated) // sorted by increasing purity (KRAKEN)",
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
pdf(paste0(RESDIR_QC,"/globaldnds_decontamination.pdf"),  width = 12, height = 7)
grid.arrange(p1)
dev.off()

globaldnds_joined_wall[globaldnds_joined_wall$sample=="MQD0006h",]
globaldnds_joined_wall[globaldnds_joined_wall$sample=="MQD0004i",]












### Correlation between dnds & kraken macaca (proxy of purity)

### Consider all variants (wall)
globaldnds_wall<-globaldnds[globaldnds$name=="wall",]

# Correlation test
cor_res <- cor.test(
  globaldnds_wall$mle,
  globaldnds_wall$kraken_macaca,
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
    y = kraken_macaca,
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


### Correlation between dnds & clean_ns

### Consider all variants (wall)
globaldnds_wall<-globaldnds[globaldnds$name=="wall",]

# Correlation test
cor_res <- cor.test(
  globaldnds_wall$mle,
  globaldnds_wall$clean_ns_reads_Homo_sapiens_perc,
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
    y = clean_ns_reads_Homo_sapiens_perc,
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
    y = "clean_ns_reads_Homo_sapiens_perc",
    color = "Tissue type"
  ) +
  
  theme_classic(base_size = 14)

pdf(paste0(RESDIR_QC,"/clean_ns_vs_dnds.pdf"),  width = 9, height = 5)
grid.arrange(p)
dev.off()


### Consider all variants (wall)
globaldnds_wall<-globaldnds[globaldnds$name=="wall",]

# Correlation test
cor_res <- cor.test(
  globaldnds_wall$mle,
  globaldnds_wall$kraken_macaca,
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
    y = kraken_macaca,
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

### Correlation between clean_ns & kraken percent

## A) Macaca_mulatta

# Correlation test
cor_res <- cor.test(
  dataMETA$clean_ns_reads_Macaca_mulatta_perc,
  dataMETA$kraken_macaca,
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
p_macaca<-ggplot(
  dataMETA,
  aes(
    x = clean_ns_reads_Macaca_mulatta_perc,
    y = kraken_macaca,
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
    x = "Macaque percentage (clean_ns)",
    y = "Macaque percentage (kraken)",
    color = "Tissue type",
    title = "KRAKEN vs CLEAN_NS // macaque targeted sequencing"
  ) +
  
  theme_classic(base_size = 14)
pdf(paste0(RESDIR_QC,"/kraken_vs_clean_ns_Macaca.pdf"),  width = 9, height = 6)
p
dev.off()


## B) Homo_sapiens

# Correlation test
cor_res <- cor.test(
  dataMETA$clean_ns_reads_Homo_sapiens_perc,
  dataMETA$kraken_homo,
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
p_homo<-ggplot(
  dataMETA,
  aes(
    x = clean_ns_reads_Homo_sapiens_perc,
    y = kraken_homo,
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
    x = "Human percentage (clean_ns)",
    y = "Human percentage (kraken)",
    color = "Tissue type",
    title = "KRAKEN vs CLEAN_NS // macaque targeted sequencing"
  ) +
  
  theme_classic(base_size = 14)


pdf(paste0(RESDIR_QC,"/kraken_vs_clean_ns_competitive.pdf"),  width = 9, height = 8)
grid.arrange(p_macaca,p_homo)
dev.off()


head(dataSOM)

### write Somatic Variants
write.table(dataSOM[,c("sample","CHROM", "POS", "REF", "ALT","FILTER", "TYPE", "TRI", "TIMES_CALLED", "DUPLEX_VAF", "DUPLEX_COV", "chr", "variant_id", "variant_id_complete",  "sample_id", "age", "tissue"),],
            file=paste0(RESDIR_VARIANTS,"/somatic_variants.txt"), quote = FALSE, sep="\t", row.names = FALSE)


### write Somatic Variants annotated with dndscv
write.table(dndsout$annotmuts,paste0(RESDIR_VARIANTS, "/annotated_mutations.txt"), quote = FALSE, sep="\t", row.names = FALSE)

### Save RDS object
#save.image("/Users/ms84/Research/postdoc/projects/PrimateSomatic/objects/Primate_somatic.RData")
