library(sigfit)
data("cosmic_signatures_v2")

# -----------------------------
# Create trinucleotide matrix (SBS96)
# -----------------------------
build_sbs96_matrix <- function(vcf_df) {
  table(vcf_df$context96) |> as.matrix()
}

# -----------------------------
# Normalize spectra
# -----------------------------
normalize_spectra <- function(mat) {
  sweep(mat, 2, colSums(mat), "/")
}

# -----------------------------
# Plot simple spectrum
# -----------------------------
plot_spectrum_simple <- function(mat, palette) {
  df <- as.data.frame(mat)
  df$context <- rownames(df)
  
  ggplot2::ggplot(df, ggplot2::aes(x = context, y = V1, fill = context)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
}
