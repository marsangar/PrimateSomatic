# ============================================
# Title: Primate Signatures script
# Author: Martín Santamarina García
# Date: 20-11-2025
# Description: This script is designed to retrieve mutational spectra from filtered VCF files
# ============================================

# Note: This code is only reliable is the density of somatic mutations is low (potential issues with clustered variants when retrieving context)

# ---- Install Libraries ----
#devtools::install_github("kgori/sigfit", build_vignettes = TRUE,
#                         build_opts = c("--no-resave-data", "--no-manual"))

# ---- Load Libraries ----
library("Biostrings")
library("GenomicRanges")
library("dplyr")
library("sigfit")
library("stringr")
library("vcfR")


# ---- Define Functions ---- 
get_context <- function(chr, pos, ref_genome) {
  if (!(chr %in% names(ref_genome))) return(NA)
  # Extract 3bp region
  seq <- subseq(ref_genome[[chr]], pos-1, pos+1)
  as.character(seq)
}

# ---- Read reference genome ----
dnaREF<-readDNAStringSet("/Users/ms84/Research/postdoc/reference_genomes/Pan_troglodytes/NHGRI_mPanTro3-v2.1_pri/reference_files/genome.fa")
names(dnaREF) <- str_split_fixed(names(dnaREF), " ", 2)[,1]   # clean headers

# ---- Read VCF ----
#vcf<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/pan_troglodytes/CHPD0002/CHPD0002b_ds0002/CHPD0002b_ds0002.vcf.gz")

vcf<-read.vcfR("/Users/ms84/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/benchmarking/latest/isec_CHPD0002b_ds0002_dcns/0000.vcf")
vcf<-read.vcfR("/Users/ms84/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/benchmarking/latest/isec_CHPD0002b_ds0002_dcns/0001.vcf")
vcf<-read.vcfR("/Users/ms84/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/benchmarking/latest/isec_CHPD0002b_ds0002_dcns/0002.vcf")
vcf<-read.vcfR("/Users/ms84/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/benchmarking/latest/isec_CHPD0002b_ds0002_dcns/0003.vcf")

vcf<-read.vcfR("/Users/ms84/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/benchmarking/latest/isec_CHPD0002f_ds0001_dcns/0000.vcf")
vcf<-read.vcfR("/Users/ms84/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/benchmarking/latest/isec_CHPD0002f_ds0001_dcns/0001.vcf")
vcf<-read.vcfR("/Users/ms84/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/benchmarking/latest/isec_CHPD0002f_ds0001_dcns/0002.vcf")
vcf<-read.vcfR("/Users/ms84/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/benchmarking/latest/isec_CHPD0002f_ds0001_dcns/0003.vcf")

dataVCF<-as.data.frame(vcf@fix)
dataVCF$POS<-as.integer(dataVCF$POS)

### Work only with SNVs
#dataVCF<-dataVCF[grepl("TYPE=snv", dataVCF$INFO),]
dataVCF<-dataVCF[nchar(dataVCF$REF)==1 & nchar(dataVCF$ALT)==1,]
dim(dataVCF)

# ---- Extract upstream & downstream sequence context from reference
dataVCF$contextREF<-apply(dataVCF, 1, function(x) get_context(as.character(x["CHROM"]), as.integer(x["POS"]), dnaREF))
dataVCF$contextREF_5   = substr(dataVCF$contextREF, 1, 1)
dataVCF$contextREF_pos = substr(dataVCF$contextREF, 2, 2)
dataVCF$contextREF_3   = substr(dataVCF$contextREF, 3, 3)

dataVCF$contextALT<-paste0(dataVCF$contextREF_5,dataVCF$ALT,dataVCF$contextREF_3)

### --- Normalize to pyrimidine context ----------------------------------
dataVCF$contexREF_normalized<-ifelse(dataVCF$REF %in% c("A","G"), reverseComplement(DNAStringSet(dataVCF$contextREF)), dataVCF$contextREF)
dataVCF$contextALT_normalized<-ifelse(dataVCF$REF %in% c("A","G"), reverseComplement(DNAStringSet(dataVCF$contextALT)), dataVCF$contextALT)

dataVCF$mutation_normalized<-paste0(dataVCF$contexREF_normalized,">",dataVCF$contextALT_normalized)
table(dataVCF$mutation_normalized)


### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations<-rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations)<-colnames(cosmic_signatures_v2)
mutations

### Add mutation counts
tbl<-as.integer(table(dataVCF$mutation_normalized))
names(tbl)<-names(table(dataVCF$mutation_normalized))
  
mutations[names(tbl)]<-tbl


### Plot raw spectra
#pdf("/Users/ms84/Research/postdoc/projects/PrimateSomatic/06_plots/prelim/benchmarking/privateDupcaller_CHPD0002b_ds0002.pdf", 24,8)
#pdf("/Users/ms84/Research/postdoc/projects/PrimateSomatic/06_plots/prelim/benchmarking/privateNanoseq_CHPD0002b_ds0002.pdf", 24,8)
#pdf("/Users/ms84/Research/postdoc/projects/PrimateSomatic/06_plots/prelim/benchmarking/shared_CHPD0002b_ds0002.pdf", 24,8)

#pdf("/Users/ms84/Research/postdoc/projects/PrimateSomatic/06_plots/prelim/benchmarking/privateDupcaller_CHPD0002f_ds0001.pdf", 24,8)
#pdf("/Users/ms84/Research/postdoc/projects/PrimateSomatic/06_plots/prelim/benchmarking/privateNanoseq_CHPD0002f_ds0001.pdf", 24,8)
#pdf("/Users/ms84/Research/postdoc/projects/PrimateSomatic/06_plots/prelim/benchmarking/shared_CHPD0002f_ds0001.pdf", 24,8)
plot_spectrum(mutations, name = "SNV Spectra")
dev.off()

### Compare the spectra of shared calls
plot_spectrum(mutations_CHPD0002b_ds0002, name = "SNV Spectra") ### n=814 SNVs
plot_spectrum(mutations_CHPD0002f_ds0001, name = "SNV Spectra") ### n=792 SNVs
cosine_sim(mutations_CHPD0002b_ds0002, mutations_CHPD0002f_ds0001) # 0.66


### Compare the spectra of shared calls
priv_CHPD0002b_ds0002<-mutations
shared_CHPD0002b_ds0002<-mutations

cosine_sim(priv_CHPD0002b_ds0002, shared_CHPD0002b_ds0002) # 0.89 Cosine similarity

priv_CHPD0002f_ds0001<-mutations
shared_CHPD0002f_ds0001<-mutations

cosine_sim(priv_CHPD0002f_ds0001, shared_CHPD0002f_ds0001) # 0.97 Cosine similarity

privNanoSeq_CHPD0002f_ds0001<-mutations

cosine_sim(privNanoSeq_CHPD0002f_ds0001, shared_CHPD0002f_ds0001) # 0.97 Cosine similarity



### SIGFIT TEST ###
data("cosmic_signatures_v2")
set.seed(1)
probs <- c(0.4, 0.3, 0.2, 0.1) %*% as.matrix(cosmic_signatures_v2[c(1, 3, 7, 11), ])
mutations <- matrix(rmultinom(1, 20000, probs), nrow = 1)
colnames(mutations) <- colnames(cosmic_signatures_v2)

plot_spectrum(mutations, name = "Simulated counts")





# Example: generate a gradient from light pink to purple
library(RColorBrewer)
library(sigfit)
# SBS 96 mutation types
my_sbs_colors <- colorRampPalette(c("#FBB4C4", "#D65A9A", "#7A0177"))(96)

# Then pass it to plot_spectrum
plot_spectrum(spectra = mutations, colors = my_sbs_colors)




