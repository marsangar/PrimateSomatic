# ============================================
# PrimateSomatic — Master Analysis Script
# --------------------------------------------
# Author: Martín Santamarina García
# Contact: ms3242@cam.ac.uk
# Date: 2026-04-21
#
# Project:
# Comparative analysis of somatic mutational processes 
# across primate species and tissues
#
# Description:
# Central orchestrator for duplex sequencing analysis:
#   • Germline / somatic filtering
#   • dN/dS estimation
#   • Driver mutation screening
#   • Mutational signature analysis
#   • Mutation burden estimation
#
# Species:
#   • Macaca mulatta (Rhesus macaque)
#   • Pan troglodytes (Chimpanzee)
#
# Technologies:
#   • NanoSeq
#   • Targeted NanoSeq
#
# Outputs:
#   • Versioned results:
#       - v0.x → exploratory / preliminary
#       - v1.x → manuscript-ready
#
# Reproducibility:
#   • renv-managed environment
#   • YAML-based configuration
#   • HPC/local portability via paths.yaml
# ============================================

#### Clean environment
rm(list = ls())
gc()
renv::restore()

#### Load core libraries
library(data.table)
library(dplyr)
library(dndscv)
library(ggplot2)
library(gridExtra)
library(openxlsx)
library(yaml)
library(stringr)


#### Load configuration
paths_config <- yaml::read_yaml("config/paths.yaml")
project_config <- yaml::read_yaml("config/project.yaml")
species_config <- yaml::read_yaml("config/species.yaml")
tech_config    <- yaml::read_yaml("config/technology.yaml")

source("config/paths.r")
source("config/palette.R")

#### Load internal modules
source("src/utils.r")
source("src/genomic_toolkit.r")

#### Define global parameters
VERSION    <- "v0.1"       # v0.1 / v0.2 / v0.3
SPECIES    <- "Macaca_mulatta"        # "Pan_troglodytes"
TECH       <- "targeted_nanoseq"        # "nanoseq"
VCF_DIR    <- "~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/filtered/"
COV_DIR    <- NA
METADATA   <- "metadata/metadata_Macaca_mulatta.xlsx"
GERMLINE   <- "~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/GERMLINE/pileup_genotype_COMBINED.txt"
RESULT_DIR <- file.path("results",VERSION, SPECIES, TECH)
REF  <- "/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa"
REFCDS<-"resources/RefCDS/RefCDS_Macaca_mulatta.Mmul_10.rda"


#### Derived paths
DATASET_ID <- paste(SPECIES, TECH, VERSION, sep = "_")
log_message(paste("Running:", DATASET_ID))

#### Create RESULT_DIR (if not existing)
make_dir(RESULT_DIR)

if(TECH=="targeted_nanoseq"){
  make_dir(file.path(RESULT_DIR, "qc"))
  make_dir(file.path(RESULT_DIR, "variants/filtered"))
  make_dir(file.path(RESULT_DIR, "variants/annotated"))
  make_dir(file.path(RESULT_DIR, "burdens"))
  make_dir(file.path(RESULT_DIR, "spectra"))
  make_dir(file.path(RESULT_DIR, "drivers/dnds"))
  make_dir(file.path(RESULT_DIR, "drivers/tier1"))
  make_dir(file.path(RESULT_DIR, "drivers/tier1"))
}

if(TECH=="nanoseq"){
  make_dir(file.path(RESULT_DIR, "qc"))
  make_dir(file.path(RESULT_DIR, "variants/filtered"))
  make_dir(file.path(RESULT_DIR, "variants/annotated"))
  make_dir(file.path(RESULT_DIR, "burdens/"))
  make_dir(file.path(RESULT_DIR, "spectra/"))
  make_dir(file.path(RESULT_DIR, "signatures/"))
}

