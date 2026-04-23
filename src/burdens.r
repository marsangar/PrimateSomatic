# -----------------------------
# Calculate mutation burden per sample
# -----------------------------
calculate_burden <- function(df, sample_col = "sampleID") {
  df %>%
    dplyr::group_by(.data[[sample_col]]) %>%
    dplyr::summarise(n_mutations = n(), .groups = "drop")
}

# -----------------------------
# Normalize burden by callable bases
# -----------------------------
normalize_burden <- function(burden_df, callable_df) {
  burden_df %>%
    dplyr::left_join(callable_df, by = "sampleID") %>%
    dplyr::mutate(
      burden_per_mb = n_mutations / (callable_bases / 1e6)
    )
}

# -----------------------------
# Add metadata (age, tissue)
# -----------------------------
annotate_burden <- function(burden_df, metadata) {
  dplyr::left_join(burden_df, metadata, by = "sampleID")
}