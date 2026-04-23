#### config/paths.R (AUTO-DETECTION).
#### This script decides whether you’re on HPC or local

library(yaml)

paths_config <- read_yaml("config/paths.yaml")

# --- Detect environment automatically ---
detect_env <- function() {
  sysname <- Sys.info()[["sysname"]]
  nodename <- Sys.info()[["nodename"]]
  
  # Heuristics (customize to your HPC)
  if (grepl("login|compute|node|slurm", nodename) ||
      grepl("Linux", sysname) && grepl("lustre", getwd())) {
    return("hpc")
  } else {
    return("local")
  }
}

ENV <- detect_env()

# Override if explicitly set
if (!is.null(paths_config$default)) {
  ENV <- paths_config$default
}

message("Running in environment: ", ENV)

# --- Load paths ---
PATHS <- paths_config[[ENV]]

BASE_DIR    <- PATHS$base_dir
DATA_DIR    <- PATHS$data_dir
RESULTS_DIR <- PATHS$results_dir
SCRATCH_DIR <- PATHS$scratch_dir

