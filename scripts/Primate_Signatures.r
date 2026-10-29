#########################################################
#### Script designed to implement signature analysis ####
#########################################################

### Load packages
library(sigfit)
library(Biostrings)
library(vcfR)

### Load functions

### 
get_context <- function(chr, pos, ref_genome) {
  if (!(chr %in% names(ref_genome))) return(NA)
  # Extract 3bp region
  seq <- subseq(ref_genome[[chr]], pos-1, pos+1)
  as.character(seq)
}

#------------------------------------------------------------
# make_spectrum_vcf
#
# Generate a 96-channel SBS mutational spectrum from mutations listed on a VCF 
#
# Args:
#   vcf_file : Character string giving the path to a VCF file containing mutations
#   dnaREF : DNAStringSet object used by get_context() to extract the trinucleotide sequence context.
#
# Returns:
#   A named integer vector of length 96 containing the counts of each pyrimidine-normalized SBS mutation type (COSMIC ordering).
#------------------------------------------------------------
make_spectrum_vcf <- function(vcf_file, dnaREF){
  
  ## Read VCF file
  vcf <- read.vcfR(vcf_file)
  
  dataMUT <- as.data.frame(vcf@fix)
  dataMUT$POS <- as.integer(dataMUT$POS)
  
  ## SNVs only
  dataMUT <- dataMUT[nchar(dataMUT$REF) == 1 & nchar(dataMUT$ALT) == 1,]
  
  ## sequence context
  dataMUT$contextREF <- apply(
    dataMUT,
    1,
    function(x)
      get_context(as.character(x["CHROM"]),
                  as.integer(x["POS"]),
                  dnaREF)
  )
  
  dataMUT$contextREF_5   <- substr(dataMUT$contextREF,1,1)
  dataMUT$contextREF_3   <- substr(dataMUT$contextREF,3,3)
  
  dataMUT$contextALT <-
    paste0(dataMUT$contextREF_5,
           dataMUT$ALT,
           dataMUT$contextREF_3)
  
  ## pyrimidine normalization
  dataMUT$contexREF_normalized <-
    ifelse(
      dataMUT$REF %in% c("A","G"),
      reverseComplement(DNAStringSet(dataMUT$contextREF)),
      dataMUT$contextREF
    )
  
  dataMUT$contextALT_normalized <-
    ifelse(
      dataMUT$REF %in% c("A","G"),
      reverseComplement(DNAStringSet(dataMUT$contextALT)),
      dataMUT$contextALT
    )
  
  dataMUT$mutation_normalized <-
    paste0(dataMUT$contexREF_normalized,
           ">",
           dataMUT$contextALT_normalized)
  
  data("cosmic_signatures_v2")
  
  mutations <- rep(0,96)
  names(mutations) <- colnames(cosmic_signatures_v2)
  
  tbl <- table(dataMUT$mutation_normalized)
  
  mutations[names(tbl)] <- as.integer(tbl)
  
  mutations
}


#------------------------------------------------------------
# make_spectrum_df
#
# Generate a 96-channel SBS mutational spectrum from mutations listed on a data.frame
#
# Args:
#   df : data.frame object containing mutations. mus have "CHROM", "POS", "REF", "ALT" columns.
#   dnaREF : DNAStringSet object used by get_context() to extract the trinucleotide sequence context.
#
# Returns:
#   A named integer vector of length 96 containing the counts of each pyrimidine-normalized SBS mutation type (COSMIC ordering).
#------------------------------------------------------------
make_spectrum_df <- function(df, dnaREF){
  
  dataMUT<-df
  
  ## SNVs only
  dataMUT <- dataMUT[nchar(dataMUT$REF) == 1 & nchar(dataMUT$ALT) == 1,]
  
  ## sequence context
  dataMUT$contextREF <- apply(
    dataMUT,
    1,
    function(x)
      get_context(as.character(x["CHROM"]),
                  as.integer(x["POS"]),
                  dnaREF)
  )
  
  dataMUT$contextREF_5   <- substr(dataMUT$contextREF,1,1)
  dataMUT$contextREF_3   <- substr(dataMUT$contextREF,3,3)
  
  dataMUT$contextALT <-
    paste0(dataMUT$contextREF_5,
           dataMUT$ALT,
           dataMUT$contextREF_3)
  
  ## pyrimidine normalization
  dataMUT$contexREF_normalized <-
    ifelse(
      dataMUT$REF %in% c("A","G"),
      reverseComplement(DNAStringSet(dataMUT$contextREF)),
      dataMUT$contextREF
    )
  
  dataMUT$contextALT_normalized <-
    ifelse(
      dataMUT$REF %in% c("A","G"),
      reverseComplement(DNAStringSet(dataMUT$contextALT)),
      dataMUT$contextALT
    )
  
  dataMUT$mutation_normalized <-
    paste0(dataMUT$contexREF_normalized,
           ">",
           dataMUT$contextALT_normalized)
  
  data("cosmic_signatures_v2")
  
  mutations <- rep(0,96)
  names(mutations) <- colnames(cosmic_signatures_v2)
  
  tbl <- table(dataMUT$mutation_normalized)
  
  mutations[names(tbl)] <- as.integer(tbl)
  
  mutations
}












# ---- Read reference genome ----
dnaREF<-readDNAStringSet("/Users/ms84/Research/postdoc/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa")
names(dnaREF) <- str_split_fixed(names(dnaREF), " ", 2)[,1]   # clean headers


# Calculate trinucleotide frequencies and mutational opportunities
freqs = colSums(trinucleotideFrequency(dnaREF))
freqs.pyr = sapply(which(substr(names(freqs), 2, 2) %in% c("C", "T")), function(i) {
  rcomp = rev.comp(names(freqs)[i])
  freqs[i] + freqs[rcomp]
})
opportunities_macaque = normalise(freqs.pyr[substr(sigfit:::mut_types(), 1, 3)])


#### Load COSMIC signatures
data("cosmic_signatures_v2")
cosmic_signatures_v2
#data("cosmic_signatures_v3")
#cosmic_signatures_v3

# Plot COSMIC signature 1
barplot(cosmic_signatures_v2[1,])

### Normalize to macaque genome opportunities
cosmic_macaque <- convert_signatures(
  signatures = cosmic_signatures_v2,
  opportunities_from = "human-genome",
  opportunities_to = opportunities_macaque
)

### Normalize to panel opportunities
cosmic_panel <- convert_signatures(
  signatures = cosmic_signatures_v2,
  opportunities_from = "human-genome",
  opportunities_to = panel_opp
)



### Load metadata and order
dataMETA<-read.xlsx(METADATA)
dataMETA <- dataMETA[
  order(dataMETA$tissue,
        dataMETA$age,
        dataMETA$individual),
]


########################################
#### Fitting signatures to COSMIC ######
########################################

base_dir <- "~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/snv"


## Generate panel for reconstruction
pdf(
  file.path(RESDIR_SIGNATURES, "reconstruction_cosmic_macaque.pdf"),
  width = 25,
  height = 25
)
par(mar = c(10, 8, 8, 4))

for (tissue in unique(dataMETA$tissue)) {
  
  message("Plotting ", tissue)
  
  tissue_samples <- dataMETA$sample_id[
    dataMETA$tissue == tissue
  ]
  
  for (sample in tissue_samples) {
    
    sample_meta <- dataMETA[dataMETA$sample_id == sample, ]
    age <- sample_meta$age
    
    ## Retrieve the VCF
    vcf <- file.path(
      base_dir,
      paste0("isec_", sample),
      "0002.vcf"
    )
    
    ## Build mutational spectrum
    mut <- make_spectrum(vcf, dnaREF)
    
    mutations <- matrix(mut, nrow = 1)
    colnames(mutations) <- names(mut)
    rownames(mutations) <- sample
    
  
    ## Fit signatures to cosmic
    mcmc_samples_fit <- fit_signatures(
      counts = mutations,
      signatures = cosmic_macaque,
      #iter = 2000,
      iter = 5000,
      warmup = 1000,
      chains = 4
    )
    
    ## Plot reconstruction
    plot_reconstruction(mcmc_samples = mcmc_samples_fit)
    
    ## Add page title
    mtext(
      sprintf(
        "%s | %s | Age %.1f years",
        tissue,
        sample,
        age
      ),
      side = 3,
      outer = TRUE,
      line = -2,
      cex = 1.4,
      font = 2
    )
  }
}

dev.off()



### Generate panel for exposures
pdf(
  file.path(RESDIR_SIGNATURES, "exposures_cosmic_macaque.pdf"),
  width = 20,
  height = 10
)
par(mar = c(10, 8, 8, 4))

for (tissue in unique(dataMETA$tissue)) {
  
  message("Plotting ", tissue)
  
  tissue_samples <- dataMETA$sample_id[
    dataMETA$tissue == tissue
  ]
  
  for (sample in tissue_samples) {
    
    sample_meta <- dataMETA[dataMETA$sample_id == sample, ]
    age <- sample_meta$age
    
    ## Retrieve the VCF
    vcf <- file.path(
      base_dir,
      paste0("isec_", sample),
      "0002.vcf"
    )
    
    ## Build mutational spectrum
    mut <- make_spectrum(vcf, dnaREF)
    
    mutations <- matrix(mut, nrow = 1)
    colnames(mutations) <- names(mut)
    rownames(mutations) <- sample
    
    
    ## Fit signatures to cosmic
    mcmc_samples_fit <- fit_signatures(
      counts = mutations,
      signatures = cosmic_macaque,
      #iter = 2000,
      iter = 5000,
      warmup = 1000,
      chains = 4
    )
    
    ## Plot exposures
    plot_exposures(mcmc_samples = mcmc_samples_fit)
    
    ## Add page title
    mtext(
      sprintf(
        "%s | %s | Age %.1f years",
        tissue,
        sample,
        age
      ),
      side = 3,
      outer = TRUE,
      line = -2,
      cex = 1.4,
      font = 2
    )
  }
}

dev.off()












########################################
#### Exploratory Signature Analysis ####
########################################

### Load somatic mutation dataset
dataSOM<-read.delim(paste0(RESDIR_VARIANTS, "/somatic_variants.txt"))
dim(dataSOM)

### Display raw mutation spectra for different tissues
mut_bladder <- make_spectrum_df(dataSOM[dataSOM$tissue %in% "Bladder",], dnaREF)
mut_liver <- make_spectrum_df(dataSOM[dataSOM$tissue %in% "Liver",], dnaREF)
mut_oesophagus <- make_spectrum_df(dataSOM[dataSOM$tissue %in% "Oesophagus",], dnaREF)
mut_prostate <- make_spectrum_df(dataSOM[dataSOM$tissue %in% "Prostate",], dnaREF)
mut_skin <- make_spectrum_df(dataSOM[dataSOM$tissue %in% "Skin",], dnaREF)

###
pdf(paste0(RESDIR_SIGNATURES,"/spectrum_tissues.pdf"), 20, 10)
par(mar = c(10, 8, 8, 4) )
plot_spectrum(mut_bladder, name = "bladder")
plot_spectrum(mut_liver, name = "liver")
plot_spectrum(mut_oesophagus, name = "oesophagus")
plot_spectrum(mut_prostate, name = "prostate")
plot_spectrum(mut_skin, name = "skin")
dev.off()



###
mutations<-matrix(mut, nrow = 1)
colnames(mutations)<-colnames(mut)


# Fit signatures
mcmc_samples_fit <- fit_signatures(
  counts = mutations,
  signatures = cosmic_signatures_v2,
  iter = 2000,
  warmup = 1000,
  chains = 4
)

### Retrieving signature exposures
exposures <- retrieve_pars(mcmc_samples_fit, 
                           par = "exposures", 
                           hpd_prob = 0.90)
names(exposures)


### Visualisation
pdf(paste0(RESDIR_SIGNATURES,"/spectrum.pdf"), 20, 10)
par(mar = c(10, 8, 8, 4) )
plot_spectrum(mutations)
dev.off()

pdf(paste0(RESDIR_SIGNATURES,"/exposures.pdf"), 20, 10)
par(mar = c(10, 8, 8, 4) )
plot_exposures(mcmc_samples = mcmc_samples_fit)
dev.off()

pdf(paste0(RESDIR_SIGNATURES,"/reconstruction.pdf"), 25, 25)
par(mar = c(10, 8, 8, 4) )
plot_reconstruction(mcmc_samples = mcmc_samples_fit)
dev.off()



####
# Normalise COSMIC signatures to genome-independent representation
sigs.fixed = convert_signatures(cosmic_signatures_v3,
                                opportunities_from="human-genome")






