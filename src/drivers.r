# -----------------------------
# Filter variants in driver genes
# -----------------------------
filter_driver_genes <- function(df, driver_genes) {
  df %>%
    dplyr::filter(gene %in% driver_genes)
}

# -----------------------------
# Count mutations per gene
# -----------------------------
count_mutations_per_gene <- function(df) {
  df %>%
    dplyr::count(gene, sort = TRUE)
}
