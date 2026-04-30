# ============================================
# Genomic Toolkit Script
# --------------------------------------------
# Author: Martín Santamarina García
# Contact: ms3242@cam.ac.uk
# Date: 24-04-2026
#
# Description:
# Toolkit with functions for genomic analysis


### Load libraries
library(ggplot2)
library(VariantAnnotation)
library(dplyr)


### Load all vcfs from a directory
load_vcf_dataset <- function(vcf_dir) {
  
  vcf_files <- list.files(vcf_dir, pattern = "\\.vcf(\\.gz)?$", full.names = TRUE)
  
  print(paste("Found", length(vcf_files), "VCF files"))
  
  vcf_list <- lapply(vcf_files, function(file) {
    
    vcf <- readVcf(file)
    print(paste("Loaded:", basename(file)))
    
    # Core variant info
    df <- as.data.frame(rowRanges(vcf))
    df$REF <- as.character(ref(vcf))
    df$ALT <- sapply(alt(vcf), function(x) as.character(x[1]))
    
    # Extract INFO field
    info_df <- as.data.frame(info(vcf))
    
    # Combine everything
    df <- cbind(df, info_df)
    
    # Add sample name
    df$sample <- tools::file_path_sans_ext(basename(file))
    
    return(df)
  })
  
  vcf_df <- bind_rows(vcf_list)
  
  return(vcf_df)
}


### Plot qc metric
plot_qc_metric <- function(
    df,
    metric,
    plot_type = c("density", "histogram", "boxplot", "violin"),
    group_by = c("all", "sample", "tissue"),
    color_by = c("none", "tissue", "sample"),
    log10_scale = FALSE,
    bins = 100,
    title = NULL,
    tissue_cols = NULL,
    sample_cols = NULL,
    sample_order = NULL,
    auto_las = TRUE,
    max_labels_horizontal = 8,
    short_labels = FALSE,
    show_stats = c("none", "mean", "median", "all"),
    stats_color = "black",
    stats_size = 0.4
) {
  
  plot_type <- match.arg(plot_type)
  group_by  <- match.arg(group_by)
  color_by  <- match.arg(color_by)
  show_stats <- match.arg(show_stats)
  
  stopifnot(metric %in% colnames(df))
  
  library(ggplot2)
  
  # ---------------------------
  # OPTIONAL LABEL SHORTENING
  # ---------------------------
  if (short_labels && "sample" %in% colnames(df)) {
    df$sample <- gsub("MQD", "M", df$sample)
  }
  
  # ---------------------------
  # LOG SCALE SAFETY
  # ---------------------------
  if (log10_scale) {
    df <- df[df[[metric]] > 0, ]
  }
  
  # ---------------------------
  # ORDERING
  # ---------------------------
  if ("sample" %in% colnames(df)) {
    if (is.null(sample_order)) {
      sample_order <- unique(df$sample)
    }
    df$sample <- factor(df$sample, levels = sample_order)
  }
  
  if ("tissue" %in% colnames(df)) {
    df$tissue <- factor(df$tissue, levels = unique(df$tissue))
  }
  
  # ---------------------------
  # TITLE
  # ---------------------------
  if (is.null(title)) {
    title <- paste(plot_type, "of", metric)
  }
  
  # ---------------------------
  # BASE AESTHETICS
  # ---------------------------
  if (plot_type %in% c("boxplot", "violin")) {
    
    if (group_by == "all") {
      stop("boxplot/violin require group_by = 'sample' or 'tissue'")
    }
    
    x_var <- if (group_by == "sample") "sample" else "tissue"
    
    p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[metric]]))
    
  } else {
    
    p <- ggplot(df, aes(x = .data[[metric]]))
  }
  
  # ---------------------------
  # GEOMETRIES
  # ---------------------------
  if (plot_type == "density") {
    
    p <- p + geom_density(alpha = 0.5)
    
  } else if (plot_type == "histogram") {
    
    p <- p + geom_histogram(bins = bins, alpha = 0.6)
    
  } else if (plot_type == "boxplot") {
    
    p <- p + geom_boxplot(outlier.alpha = 0.3)
    
  } else if (plot_type == "violin") {
    
    p <- p +
      geom_violin(alpha = 0.6, trim = FALSE) +
      geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.7)
  }
  
  # ---------------------------
  # COLOR MAPPING
  # ---------------------------
  if (color_by == "tissue") {
    p <- p + aes(fill = tissue)
  } else if (color_by == "sample") {
    p <- p + aes(fill = sample)
  }
  
  # ---------------------------
  # STATS (NEW)
  # ---------------------------
  if (show_stats != "none" && is.numeric(df[[metric]])) {
    
    # GLOBAL STATS (density / histogram)
    if (plot_type %in% c("density", "histogram")) {
      
      if (show_stats %in% c("mean", "all")) {
        p <- p + geom_vline(
          xintercept = mean(df[[metric]], na.rm = TRUE),
          linetype = "dashed",
          color = stats_color
        )
      }
      
      if (show_stats %in% c("median", "all")) {
        p <- p + geom_vline(
          xintercept = median(df[[metric]], na.rm = TRUE),
          linetype = "solid",
          color = stats_color
        )
      }
    }
    
    # GROUPED STATS (boxplot / violin)
    if (plot_type %in% c("boxplot", "violin")) {
      
      if (show_stats %in% c("mean", "all")) {
        p <- p + stat_summary(fun = mean, geom = "point",
                              shape = 23, size = 2,
                              fill = stats_color)
      }
      
      if (show_stats %in% c("median", "all")) {
        p <- p + stat_summary(fun = median, geom = "point",
                              shape = 95, size = 4,
                              color = stats_color)
      }
    }
  }
  
  # ---------------------------
  # STYLING
  # ---------------------------
  p <- p +
    theme_classic() +
    labs(
      title = title,
      x = if (plot_type %in% c("boxplot", "violin")) group_by else metric,
      y = if (plot_type %in% c("boxplot", "violin")) metric else "Count",
      fill = color_by
    )
  
  # ---------------------------
  # AUTO LAS
  # ---------------------------
  if (auto_las && group_by != "all") {
    
    n_levels <- length(unique(df[[group_by]]))
    
    axis_angle <- if (n_levels > max_labels_horizontal) 90 else 0
    
    p <- p + theme(
      axis.text.x = element_text(angle = axis_angle, hjust = 1),
      axis.text.y = element_text(angle = 0)
    )
  }
  
  # ---------------------------
  # LOG SCALE
  # ---------------------------
  if (log10_scale) {
    if (plot_type %in% c("boxplot", "violin")) {
      p <- p + scale_y_log10()
    } else {
      p <- p + scale_x_log10()
    }
  }
  
  # ---------------------------
  # COLOR SCALES
  # ---------------------------
  if (color_by == "tissue" && !is.null(tissue_cols)) {
    p <- p + scale_fill_manual(values = tissue_cols)
  }
  
  if (color_by == "sample" && !is.null(sample_cols)) {
    p <- p + scale_fill_manual(values = sample_cols)
  }
  
  # ---------------------------
  # FACETING
  # ---------------------------
  if (plot_type %in% c("density", "histogram")) {
    
    if (group_by == "sample") {
      p <- p + facet_wrap(~sample, scales = "free_y")
    }
    
    if (group_by == "tissue") {
      p <- p + facet_wrap(~tissue, scales = "free_y")
    }
  }
  
  return(p)
}
