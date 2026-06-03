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
library(ggrepel)
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
COV_DIR    <- "~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/post/"
METADATA   <- "metadata/metadata_Macaca_mulatta.xlsx"
GERMLINE   <- "~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/GERMLINE/pileup_genotype_COMBINED.txt"
RESDIR <- file.path("results",VERSION, SPECIES, TECH)
REF  <- "/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa"
NCBI_REPORT  <- "~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/Mmul_10_sequence_report.tsv"
REFCDS<-"resources/RefCDS/RefCDS_Macaca_mulatta.Mmul_10.rda"
PANEL<-"~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/resources/panels/Macaca_mulatta/Mmul_10/Sanger_TERT-v4_TE-95148282_hg19_highstringencyfilter_buccal_gene_list.tsv"
KRAKEN<-"~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/KRAKEN.RESULTS.tsv"


#### Derived paths
DATASET_ID <- paste(SPECIES, TECH, VERSION, sep = "_")
log_message(paste("Running:", DATASET_ID))

#### Create RESDIR (if not existing)
make_dir(RESDIR)

RESDIR_QC<-paste0(RESDIR,"/", "qc")
RESDIR_VARIANTS<-paste0(RESDIR,"/", "variants")
RESDIR_BURDENS<-paste0(RESDIR,"/", "burdens")
RESDIR_SPECTRA<-paste0(RESDIR,"/", "spectra")
RESDIR_SIGNATURES<-paste0(RESDIR,"/", "signatures")
RESDIR_DRIVERS<-paste0(RESDIR,"/", "drivers")
RESDIR_GENES<-paste0(RESDIR,"/", "genes")



if(TECH=="targeted_nanoseq"){
  make_dir(file.path(RESDIR_QC))
  make_dir(file.path(RESDIR_VARIANTS))
  make_dir(file.path(RESDIR_SPECTRA))
  make_dir(file.path(RESDIR_DRIVERS))
  make_dir(file.path(RESDIR_GENES))
}

if(TECH=="nanoseq"){
  make_dir(file.path(RESDIR_QC))
  make_dir(file.path(RESDIR_VARIANTS))
  make_dir(file.path(RESDIR_BURDENS))
  make_dir(file.path(RESDIR_SPECTRA))
  make_dir(file.path(RESDIR_SIGNATURES))
}

