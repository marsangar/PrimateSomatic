#########################################################
#### Script designed to implement signature analysis ####
#########################################################

### Load packages
library(sigfit)
library(Biostrings)

'''''
make_spectrum <- function(vcf_file, dnaREF){
  
  vcf <- read.vcfR(vcf_file)
  
  dataVCF <- as.data.frame(vcf@fix)
  dataVCF$POS <- as.integer(dataVCF$POS)
  
  ## SNVs only
  dataVCF <- dataVCF[nchar(dataVCF$REF) == 1 & nchar(dataVCF$ALT) == 1,]
  
  ## sequence context
  dataVCF$contextREF <- apply(
    dataVCF,
    1,
    function(x)
      get_context(as.character(x["CHROM"]),
                  as.integer(x["POS"]),
                  dnaREF)
  )
  
  dataVCF$contextREF_5   <- substr(dataVCF$contextREF,1,1)
  dataVCF$contextREF_3   <- substr(dataVCF$contextREF,3,3)
  
  dataVCF$contextALT <-
    paste0(dataVCF$contextREF_5,
           dataVCF$ALT,
           dataVCF$contextREF_3)
  
  ## pyrimidine normalization
  dataVCF$contexREF_normalized <-
    ifelse(
      dataVCF$REF %in% c("A","G"),
      reverseComplement(DNAStringSet(dataVCF$contextREF)),
      dataVCF$contextREF
    )
  
  dataVCF$contextALT_normalized <-
    ifelse(
      dataVCF$REF %in% c("A","G"),
      reverseComplement(DNAStringSet(dataVCF$contextALT)),
      dataVCF$contextALT
    )
  
  dataVCF$mutation_normalized <-
    paste0(dataVCF$contexREF_normalized,
           ">",
           dataVCF$contextALT_normalized)
  
  data("cosmic_signatures_v2")
  
  mutations <- rep(0,96)
  names(mutations) <- colnames(cosmic_signatures_v2)
  
  tbl <- table(dataVCF$mutation_normalized)
  
  mutations[names(tbl)] <- as.integer(tbl)
  
  mutations
}
'''''

### Load functions
get_context <- function(chr, pos, ref_genome) {
  if (!(chr %in% names(ref_genome))) return(NA)
  # Extract 3bp region
  seq <- subseq(ref_genome[[chr]], pos-1, pos+1)
  as.character(seq)
}


#------------------------------------------------------------
# make_spectrum
#
# Generate a 96-channel SBS mutational spectrum from either:
#   1) a VCF file, or
#   2) a data.frame containing SNVs.
#
# Args:
#   x      : Character string giving the path to a VCF file, or a
#            data.frame with columns CHROM, POS, REF and ALT.
#   dnaREF : BSgenome/DNA reference object used by get_context() to
#            extract the trinucleotide sequence context.
#
# Returns:
#   A named integer vector of length 96 containing the counts of each
#   pyrimidine-normalized SBS mutation type (COSMIC ordering).
#------------------------------------------------------------

make_spectrum <- function(x, dnaREF){
  
  ## Read input
  if (is.character(x) && length(x) == 1) {
    
    ## Input is a VCF file
    vcf <- read.vcfR(x)
    dataVCF <- as.data.frame(vcf@fix)
    
  } else if (is.data.frame(x)) {
    
    ## Input is already a data.frame of variants
    dataVCF <- x
    
  } else {
    
    stop("Input must be either a VCF file path or a data.frame.")
    
  }
  
  ## Check that the required columns are present
  required_cols <- c("CHROM", "POS", "REF", "ALT")
  
  if (!all(required_cols %in% colnames(dataVCF))) {
    stop(
      "Input data.frame must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }
  
  dataVCF$POS <- as.integer(dataVCF$POS)
  
  ## Keep SNVs only
  dataVCF <- dataVCF[
    nchar(dataVCF$REF) == 1 &
      nchar(dataVCF$ALT) == 1,
  ]
  
  ## Retrieve trinucleotide sequence context
  dataVCF$contextREF <- mapply(
    get_context,
    chrom = as.character(dataVCF$CHROM),
    pos = dataVCF$POS,
    MoreArgs = list(dnaREF = dnaREF)
  )
  
  ## Build reference and alternate trinucleotide contexts
  dataVCF$contextREF_5 <- substr(dataVCF$contextREF, 1, 1)
  dataVCF$contextREF_3 <- substr(dataVCF$contextREF, 3, 3)
  
  dataVCF$contextALT <- paste0(
    dataVCF$contextREF_5,
    dataVCF$ALT,
    dataVCF$contextREF_3
  )
  
  ## Normalize mutations to the pyrimidine representation
  idx <- dataVCF$REF %in% c("A", "G")
  
  dataVCF$contexREF_normalized <- dataVCF$contextREF
  dataVCF$contexREF_normalized[idx] <-
    as.character(reverseComplement(
      DNAStringSet(dataVCF$contextREF[idx])
    ))
  
  dataVCF$contextALT_normalized <- dataVCF$contextALT
  dataVCF$contextALT_normalized[idx] <-
    as.character(reverseComplement(
      DNAStringSet(dataVCF$contextALT[idx])
    ))
  
  dataVCF$mutation_normalized <- paste0(
    dataVCF$contexREF_normalized,
    ">",
    dataVCF$contextALT_normalized
  )
  
  ## Load COSMIC SBS96 mutation classes
  data("cosmic_signatures_v2")
  
  ## Initialize the spectrum
  mutations <- integer(96)
  names(mutations) <- colnames(cosmic_signatures_v2)
  
  ## Count mutations belonging to each SBS96 class
  tbl <- table(dataVCF$mutation_normalized)
  mutations[names(tbl)] <- as.integer(tbl)
  
  return(mutations)
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












### Load somatic mutations
mut <- make_spectrum("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/snv/isec_MQD0004k_tds0002/0002.vcf",
                     dnaREF)
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



