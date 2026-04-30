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
pdf(paste0(RESULT_DIR,"/qc/metrics.pdf"), 12, 6)

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
dataSOM<-dataVCF_sorted[dataVCF_sorted$TYPE %in% c("snv", "del", "ins") &  dataVCF_sorted$SHARED_ACROSS_SAMPLES==FALSE & dataVCF_sorted$DPLX_NM %in% c(1, NA) & dataVCF_sorted$FILTER %in% "PASS",c("SAMPLE","CHROM", "POS", "REF", "ALT","AGE", "TISSUE", "TIMES_CALLED", "TYPE", "varID", "SHARED_ACROSS_SAMPLES", "FILTER", "DPLX_NM","DUPLEX_COV","BAM_COV", "DUPLEX_VAF", "BAM_VAF")]
dataSOM$DPLX_NM
dim(dataSOM) #
table(dataSOM$TYPE)
table(dataSOM$TYPE)


### Evaluate impact of filtering on dnds
dndsout = dndscv(dataSOM, refdb=REFCDS, cv=NULL)
dndsout_value<-signif(dndsout$globaldnds["wall","mle"],2)

list_dndsout <- lapply(
  split(dataSOM, dataSOM$SAMPLE),
  function(df_sample) {
    dndscv(df_sample, refdb = REFCDS, cv = NULL)
  }
)

library(dplyr)

globaldnds <- bind_rows(
  lapply(names(list_dndsout), function(sid) {
    res <- list_dndsout[[sid]]
    if (is.null(res$globaldnds)) return(NULL)
    res$globaldnds %>%
      as.data.frame() %>%
      mutate(sample = sid)
    
  })
)


globaldnds$name <- factor(globaldnds$name, levels = c("wmis","wnon","wspl","wtru","wall"))

### Plot dnds across samples
p<-ggplot(globaldnds, aes(x = sample, y = mle, color = name)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = cilow, ymax = 6),width = 0.2,position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = dndsout_value, linetype = "dashed", color = "red") +
  
  scale_y_continuous(breaks = seq(0, 6, by = 0.5))  +
  #scale_color_manual(values = c(1:9)) +
  labs(
    title = paste0("DPLX_NM=", i),
    subtitle = paste0("Global dN/dS=", dndsout_value),
    x = "Sample",
    y = "w (dN/dS ratio)",
    color = "dN/dS"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16)
  )
p
print(paste0("Adding plot:", i))
plotlist_dnds[[i]]<-p
}

plotlist_dnds<-plotlist
pdf("~/Research/postdoc/projects/PrimateSomatic/06_plots/prelim/params/paramsDPLX_NM_dnds.pdf", width = 10, height = 18)
plotlist_dnds[[1]] / plotlist_dnds[[2]] / plotlist_dnds[[3]]
dev.off()




### write Somatic Variants
write.table(dataSOM, file=paste0(RESULT_DIR,"/somatic.txt"), quote = FALSE, sep="\t", row.names = FALSE)







