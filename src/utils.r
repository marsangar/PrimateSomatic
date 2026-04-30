# -----------------------------
# Create directory safely
# -----------------------------
make_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

# -----------------------------
# Save table (csv + rds)
# -----------------------------
save_table <- function(df, path, name) {
  make_dir(path)
  write.csv(df, file.path(path, paste0(name, ".csv")), row.names = FALSE)
  saveRDS(df, file.path(path, paste0(name, ".rds")))
}

# -----------------------------
# Save plot
# -----------------------------
save_plot <- function(plot, path, name, width = 8, height = 6) {
  make_dir(path)
  ggsave(filename = file.path(path, paste0(name, ".pdf")),
         plot = plot, width = width, height = height)
}

# -----------------------------
# Log message with timestamp
# -----------------------------
log_message <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg))
}

