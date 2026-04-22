# ============================================
# PrimateSomatic Master Script
# Author: Martín Santamarina García
# Contact: ms3242@cam.ac.uk
# Date: 2026-04-21
# Project: Comparative analyses of somatic mutational processes in primates across lifespans
#
# Description:
# Main orchestration script for duplex sequencing analysis:
# - Germline / somatic filtering
# - dN/dS estimation
# - Driver screening
# - Mutational signatures
# - Mutation burden analysis
#
# Species:
# - Macaca mulatta
# - Pan troglodytes
#
# Technologies:
# - NanoSeq
# - Targeted NanoSeq
#
# Output:
# - Versioned results: v0.1,v0.2,v0.3... (exploratory) // v1.0, v1.1, v1.2.. (publication) 
#
# Reproducibility:
# - renv-managed environment
# ============================================

#### Clean environment
rm(list = ls())
gc()

#### Load package manager
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::activate()

#### Load core packages
library(dplyr)
library(data.table)
library(ggplot2)
library(openxlsx)
library(yaml)

#### Load project configuration
paths_config <- yaml::read_yaml("config/paths.yaml")

project_config <- yaml::read_yaml("config/project.yaml")
species_config <- yaml::read_yaml("config/species.yaml")
tech_config    <- yaml::read_yaml("config/technology.yaml")

source("config/palette.R")

#### Load internal functions
source("src/utils.R")

#### Load external toolkit

#### Define global parameters
SPECIES    <- "Macaca_mulatta"        # or "Pan_troglodytes"
REFERENCE <- "Mmul_10"
TECH       <- "nanoseq"        # or "targeted_nanoseq"
VERSION    <- "v0.1"       # v0.1 / v0.2 / v0.3
RUN_MODE   <- "full"           # full / test / debug
