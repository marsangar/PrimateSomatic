# ============================================
# Title: Primate dndscv script
# Author: Martín Santamarina García
# Date: 10-09-2025
# Description: This script applies dndscv to primate somatic dataset
# ============================================

# ---- Install Libraries ----
#install.packages("devtools")
#install.packages("vcfR")

#library(devtools); install_github("im3sanger/dndscv")

# ---- Load Libraries ----
library("dndscv")
library("vcfR")
library("dplyr")
library("ggplot2")
library("karyoploteR")
library("gridExtra")
library("openxlsx")
library("patchwork")
library("stringr")

mem.maxVSize()

test_palette <- c(
  a       = "#A6CEE3",  # soft blue
  b       = "#B2DF8A",  # soft green
  c       = "#FDBF6F",  # soft orange
  d       = "#CAB2D6",  # soft lavender
  e       = "#FFFF99",  # pale yellow
  f       = "#FF9DA7",  # soft pink
  g       = "#CCEBC5",  # mint green
  h       = "#D9D9D9",  # light gray
  i       = "#B3CDE3",  # muted light blue
  j       = "#FBB4AE",  # pastel rose
  k       = "#DECBE4",  # pastel violet
  l       = "#FED9A6"   # warm light orange
)


# Optional: define consistent colors for tissues
tissue_cols <- c(
  "Bladder" = "#F8766D",
  "Liver" = "#A3A500",
  "Oesophagus" = "#00BF7D",
  "Prostate" = "#00B0F6",
  "Skin" = "#E76BF3"
)


# ---- Custom functions ----

rename_chromosomes_Mmul_10<-function(chrom){
  ### This function renames macaque chromosome naming from UCSC.style.name to RefSeq.seq.accession
  # input: a vector with UCSC.style.name chromosome/sequence identifiers
  # output: a vector with RefSeq.seq.accession chromosome/sequence identifiers
  
  report<-read.delim("/Users/ms84/Research/postdoc/projects/PrimateSomatic/02_data/databases/macaca_mulatta/Mmul_10_sequence_report.tsv")
  output<-report[match(chrom, report$UCSC.style.name),"RefSeq.seq.accession"]
  return(output)
}


### Rename BIOMART cds annotation to match our reference fasta convention
#dataBIOMART<-read.delim("/Users/ms84/Research/postdoc/projects/PrimateSomatic/02_data/databases/macaca_mulatta/BioMart_Macaca_mulatta_Mmul_10.txt")
#table(dataBIOMART$Chromosome.scaffold.name)
#dataBIOMART$Chromosome.scaffold.name<-rename_chromosomes_Mmul_10(paste0("chr", dataBIOMART$Chromosome.scaffold.name))
#table(dataBIOMART$Chromosome.scaffold.name)
#dataBIOMART_renamed<-dataBIOMART
#write.table(dataBIOMART_renamed, "/Users/ms84/Research/postdoc/projects/PrimateSomatic/02_data/databases/macaca_mulatta/BioMart_Macaca_mulatta_Mmul_10_renamed.txt", sep = "\t", quote = FALSE)


#### Using buildref
#path_cds_table = "/Users/ms84/Research/postdoc/projects/PrimateSomatic/02_data/databases/macaca_mulatta/BioMart_Macaca_mulatta_Mmul_10_renamed.txt"
#path_genome_fasta = "/Users/ms84/Research/postdoc/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa"
#buildref(cdsfile=path_cds_table, genomefile=path_genome_fasta, outfile = "RefCDS_Macaca_mulatta.Mmul_10.rda", excludechrs="MT")
# Warning message: 21591 unique gene IDs (column 1) found. 17105 unique gene names (column 2) found. Consider combining gene names and gene IDs or replacing gene names by gene IDs to avoid losing genes (see useids argument in ? buildref)


### Load mGAP database
#dataMGAP<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/02_data/databases/macaca_mulatta/mGap.Rhesus_macaque.v3.0.vcf.gz")


### Define sample vector
macaca_samples<-c("MQD0001d",
                  "MQD0001j",
                  "MQD0002e",
                  "MQD0002f",
                  "MQD0002h",
                  "MQD0002l",
                  "MQD0003d",
                  "MQD0003e",
                  "MQD0003g",
                  "MQD0003j",
                  "MQD0004d",
                  "MQD0004e",
                  "MQD0004i",
                  "MQD0004k",
                  "MQD0006e",
                  "MQD0006f",
                  "MQD0006h",
                  "MQD0006l")

# Conflictive samples: MQD0001h_tds0001, MQD0006j_tds0001, MQD0002j_tds0001

### Load Variant Calls
vcf_MQD0001d<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0001d_tds0001.filtered.vcf.gz")
vcf_MQD0001j<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0001j_tds0001.filtered.vcf.gz")
vcf_MQD0002e<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0002e_tds0001.filtered.vcf.gz")
vcf_MQD0002f<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0002f_tds0001.filtered.vcf.gz")
vcf_MQD0002h<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0002h_tds0001.filtered.vcf.gz")
vcf_MQD0002l<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0002l_tds0001.filtered.vcf.gz")
vcf_MQD0003d<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0003d_tds0001.filtered.vcf.gz")
vcf_MQD0003e<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0003e_tds0001.filtered.vcf.gz")
vcf_MQD0003g<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0003g_tds0001.filtered.vcf.gz")
vcf_MQD0003j<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0003j_tds0001.filtered.vcf.gz")
vcf_MQD0004d<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0004d_tds0001.filtered.vcf.gz")
vcf_MQD0004e<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0004e_tds0002.filtered.vcf.gz")
vcf_MQD0004i<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0004i_tds0001.filtered.vcf.gz")
vcf_MQD0004k<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0004k_tds0002.filtered.vcf.gz")
vcf_MQD0006e<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0006e_tds0001.filtered.vcf.gz")
vcf_MQD0006f<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0006f_tds0001.filtered.vcf.gz")
vcf_MQD0006h<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0006h_tds0001.filtered.vcf.gz")
vcf_MQD0006l<-read.vcfR("/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/export/MQD0006l_tds0001.filtered.vcf.gz")

dataVCF_MQD0001d<-as.data.frame(vcf_MQD0001d@fix)
dataVCF_MQD0001j<-as.data.frame(vcf_MQD0001j@fix)
dataVCF_MQD0002e<-as.data.frame(vcf_MQD0002e@fix)
dataVCF_MQD0002f<-as.data.frame(vcf_MQD0002f@fix)
dataVCF_MQD0002h<-as.data.frame(vcf_MQD0002h@fix)
dataVCF_MQD0002l<-as.data.frame(vcf_MQD0002l@fix)
dataVCF_MQD0003d<-as.data.frame(vcf_MQD0003d@fix)
dataVCF_MQD0003e<-as.data.frame(vcf_MQD0003e@fix)
dataVCF_MQD0003g<-as.data.frame(vcf_MQD0003g@fix)
dataVCF_MQD0003j<-as.data.frame(vcf_MQD0003j@fix)
dataVCF_MQD0004d<-as.data.frame(vcf_MQD0004d@fix)
dataVCF_MQD0004e<-as.data.frame(vcf_MQD0004e@fix)
dataVCF_MQD0004i<-as.data.frame(vcf_MQD0004i@fix)
dataVCF_MQD0004k<-as.data.frame(vcf_MQD0004k@fix)
dataVCF_MQD0006e<-as.data.frame(vcf_MQD0006e@fix)
dataVCF_MQD0006f<-as.data.frame(vcf_MQD0006f@fix)
dataVCF_MQD0006h<-as.data.frame(vcf_MQD0006h@fix)
dataVCF_MQD0006l<-as.data.frame(vcf_MQD0006l@fix)

dataVCF_MQD0001d$sampleID<-"MQD0001d"
dataVCF_MQD0001j$sampleID<-"MQD0001j"
dataVCF_MQD0002e$sampleID<-"MQD0002e"
dataVCF_MQD0002f$sampleID<-"MQD0002f"
dataVCF_MQD0002h$sampleID<-"MQD0002h"
dataVCF_MQD0002l$sampleID<-"MQD0002l"
dataVCF_MQD0003d$sampleID<-"MQD0003d"
dataVCF_MQD0003e$sampleID<-"MQD0003e"
dataVCF_MQD0003g$sampleID<-"MQD0003g"
dataVCF_MQD0003j$sampleID<-"MQD0003j"
dataVCF_MQD0004d$sampleID<-"MQD0004d"
dataVCF_MQD0004e$sampleID<-"MQD0004e"
dataVCF_MQD0004i$sampleID<-"MQD0004i"
dataVCF_MQD0004k$sampleID<-"MQD0004k"
dataVCF_MQD0006e$sampleID<-"MQD0006e"
dataVCF_MQD0006f$sampleID<-"MQD0006f"
dataVCF_MQD0006h$sampleID<-"MQD0006h"
dataVCF_MQD0006l$sampleID<-"MQD0006l"

dataVCF_MQD0001d$age<-15.37
dataVCF_MQD0001j$age<-15.37
dataVCF_MQD0002e$age<-5.41
dataVCF_MQD0002f$age<-5.41
dataVCF_MQD0002h$age<-5.41
dataVCF_MQD0002l$age<-5.41
dataVCF_MQD0003d$age<-17.74
dataVCF_MQD0003e$age<-17.74
dataVCF_MQD0003g$age<-17.74
dataVCF_MQD0003j$age<-17.74
dataVCF_MQD0004d$age<-19.53
dataVCF_MQD0004e$age<-19.53
dataVCF_MQD0004i$age<-19.53
dataVCF_MQD0004k$age<-19.53
dataVCF_MQD0006e$age<-16.18
dataVCF_MQD0006f$age<-16.18
dataVCF_MQD0006h$age<-16.18
dataVCF_MQD0006l$age<-16.18
  
dataVCF_MQD0001d$tissue<-"Liver"
dataVCF_MQD0001j$tissue<-"Bladder"
dataVCF_MQD0002e$tissue<-"Liver"
dataVCF_MQD0002f$tissue<-"Oesophagus"
dataVCF_MQD0002h$tissue<-"Prostate"
dataVCF_MQD0002l$tissue<-"Bladder"
dataVCF_MQD0003d$tissue<-"Liver"
dataVCF_MQD0003e$tissue<-"Oesophagus"
dataVCF_MQD0003g$tissue<-"Prostate"
dataVCF_MQD0003j$tissue<-"Bladder"
dataVCF_MQD0004d$tissue<-"Liver"
dataVCF_MQD0004e$tissue<-"Oesophagus"
dataVCF_MQD0004i$tissue<-"Skin"
dataVCF_MQD0004k$tissue<-"Bladder"
dataVCF_MQD0006e$tissue<-"Liver"
dataVCF_MQD0006f$tissue<-"Oesophagus"
dataVCF_MQD0006h$tissue<-"Prostate"
dataVCF_MQD0006l$tissue<-"Bladder"

dataVCF<-rbind(dataVCF_MQD0001d,
               dataVCF_MQD0001j,
               dataVCF_MQD0002e,
               dataVCF_MQD0002f,
               dataVCF_MQD0002h,
               dataVCF_MQD0002l,
               dataVCF_MQD0003d,
               dataVCF_MQD0003e,
               dataVCF_MQD0003g,
               dataVCF_MQD0003j,
               dataVCF_MQD0004d,
               dataVCF_MQD0004e,
               dataVCF_MQD0004i,
               dataVCF_MQD0004k,
               dataVCF_MQD0006e,
               dataVCF_MQD0006f,
               dataVCF_MQD0006h,
               dataVCF_MQD0006l)

str(dataVCF)
dim(dataVCF)
head(dataVCF)

table(dataVCF$CHROM)
dataVCF$INFO[2]
dataVCF[715,]

### Convert to uppercase and define variant id
dataVCF$CHROM
dataVCF$POS<-as.integer(dataVCF$POS)
#dataVCF$REF<-dataVCF$REF
#dataVCF$ALT<-dataVCF$ALT
dataVCF$REF<-toupper(dataVCF$REF)
dataVCF$ALT<-toupper(dataVCF$ALT)

#dataVCF$variantID<-paste0(dataVCF$CHROM,"_",dataVCF$POS, "_", dataVCF$REF, ">", dataVCF$ALT)
dataVCF$variantID<-paste0(dataVCF$CHROM,"_",dataVCF$POS)



### Set relevant info as independent columns
dataVCF$SHARED_ACROSS_SAMPLES<-apply(dataVCF,1,function(x) ifelse(x["variantID"] %in% dataVCF[!dataVCF$sampleID %in% x["sampleID"],"variantID"], TRUE,FALSE))
table(dataVCF$SHARED_ACROSS_SAMPLES)
dataVCF[dataVCF$SHARED_ACROSS_SAMPLES=="FALSE",]
dataVCF[dataVCF$variantID=="NC_041754.1_83599401",c("variantID","sampleID","SHARED_ACROSS_SAMPLES")]
dataVCF[dataVCF$variantID=="NC_041769.1_27457812",c("variantID", "sampleID","SHARED_ACROSS_SAMPLES")]
dataVCF[dataVCF$variantID=="NC_041755.1_194394201",c("variantID", "sampleID","SHARED_ACROSS_SAMPLES")]
dataVCF[dataVCF$variantID=="NC_041755.1_108560386",c("variantID", "sampleID","SHARED_ACROSS_SAMPLES")]

dataVCF$TIMES_CALLED<-sub(".*TIMES_CALLED=([^;]+).*", "\\1", dataVCF$INFO)
dataVCF$TYPE<-sub(".*TYPE=([^;]+).*", "\\1", dataVCF$INFO)
dataVCF$DUPLEX_COV<-as.integer(sub(".*DUPLEX_COV=([^;]+).*", "\\1", dataVCF$INFO))
dataVCF$BAM_COV<-as.integer(sub(".*BAM_COV=([^;]+).*", "\\1", dataVCF$INFO))
dataVCF$DUPLEX_VAF<-as.numeric(sub(".*DUPLEX_VAF=([^;]+).*", "\\1", dataVCF$INFO))
dataVCF$BAM_VAF<-as.numeric(sub(".*BAM_VAF=([^;]+).*", "\\1", dataVCF$INFO))
dataVCF$DPLX_NM<-as.numeric(sub(".*DPLX_NM=([^;]+).*", "\\1", dataVCF$INFO))


hist(sort(table(dataVCF$variantID)), breaks = 100)
dim(dataVCF)
prop.table(table(dataVCF$RECURRENT))*100 #
prop.table(table(dataVCF$TIMES_CALLED))*100 # 
prop.table(table(dataVCF$TYPE))*100 #
prop.table(table(dataVCF$DPLX_NM))*100

median(dataVCF$DUPLEX_COV, na.rm = TRUE)
quantile(dataVCF$DUPLEX_COV, 0.9, na.rm=TRUE)
quantile(dataVCF$DUPLEX_COV, 0.1, na.rm=TRUE)


### Filtering with MASK file
#library(data.table)
#dataMASK<-fread("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/maskdir/output/MASKS/SNP+NOISE.NV2.SeenNV1.Macaca_mulatta.bed.gz")


### FILTERING: Keep only PASS
dataGERMLINE<-read.delim("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/GERMLINE/pileup_genotype_COMBINED.txt")
dim(dataGERMLINE)
dataGERMLINE$SAMPLE<-str_split_fixed(dataGERMLINE$SAMPLE, "_", 2)[,1]
#dataGERMLINE$variantID<-paste0(dataGERMLINE$NC_CHROM, "_", dataGERMLINE$POS, "_", dataGERMLINE$REF, ">", dataGERMLINE$ALT)
dataGERMLINE$variantID<-paste0(dataGERMLINE$NC_CHROM, "_", dataGERMLINE$POS)

dataVCF[1,]
dataVCF[100,]
dataVCF[200,]
dataVCF[300,]
dataVCF[400,]
dataVCF[500,]

head(dataGERMLINE)
dataGERMLINE[100,]
dataGERMLINE[200,]
dataGERMLINE[300,]
dataGERMLINE[400,]
dataGERMLINE[500,]

dataGERMLINE[dataGERMLINE$variantID %in% "NC_041754.1_222607540",]
dataGERMLINE[dataGERMLINE$variantID %in% "NC_041756.1_47075696",]
dataGERMLINE[dataGERMLINE$variantID %in% "NC_041756.1_47075696",]

### check VAF distribution
for(i in macaca_samples){
  hist(dataGERMLINE[dataGERMLINE$SAMPLE %in% i, "VAF"], main = i)
}
hist(dataGERMLINE$VAF)

dim(dataVCF)
prop.table(table(dataVCF$TYPE))

### Filtering Variants
### (sapply(macaca_samples, function())
dim(dataVCF)
#dataVCF<-dataVCF[dataVCF$variantID %in% dataGERMLINE$variantID,] # filter to exclude 
#dataVCF<-dataVCF[dataVCF$FILTER %in% "PASS",] # filter to keep only pass variants
#dataVCF<-dataVCF[dataVCF$DPLX_NM %in% "1",] # filter to to keep duplex nm == 1 
#dataVCF<-dataVCF[dataVCF$DUPLEX_COV < quantile(dataVCF$DUPLEX_COV, 0.9, na.rm=TRUE) & dataVCF$DUPLEX_COV > quantile(dataVCF$DUPLEX_COV, 0.1, na.rm=TRUE),] # filter to exclude coverage outlayers
dim(dataVCF)
dataVCF$DUPLEX_COV
prop.table(table(dataVCF$TYPE))

table(dataVCF$TYPE)
table(is.na(dataVCF[dataVCF$TYPE %in% c("snv"),"DPLX_NM"]))
table(is.na(dataVCF[dataVCF$TYPE %in% c("dnv"),"DPLX_NM"]))
table(is.na(dataVCF[dataVCF$TYPE %in% c("del"),"DPLX_NM"]))
table(is.na(dataVCF[dataVCF$TYPE %in% c("ins"),"DPLX_NM"]))

### Explore impact of DPLX_NM filter
prop.table(table(dataVCF[dataVCF$DPLX_NM %in% 1,c("TYPE")]))
prop.table(table(dataVCF[dataVCF$DPLX_NM %in% 2,c("TYPE")]))
prop.table(table(dataVCF[dataVCF$DPLX_NM %in% 3,c("TYPE")]))

dataVCF[dataVCF$sampleID %in% "MQD0002h" & dataVCF$DPLX_NM %in% 3,][1,] # macaque
dataVCF[dataVCF$sampleID %in% "MQD0002h" & dataVCF$DPLX_NM %in% 3,][2,] # ?
dataVCF[dataVCF$sampleID %in% "MQD0002h" & dataVCF$DPLX_NM %in% 3,][3,] # ?

############################################
#### *** Select mutation data.frame *** ####
############################################
dim(dataVCF)
dataMUT<-dataVCF[dataVCF$TYPE %in% c("snv", "del", "ins") &  dataVCF$SHARED_ACROSS_SAMPLES==FALSE & dataVCF$DPLX_NM %in% c(1, NA) & dataVCF$FILTER %in% "PASS",c("sampleID","CHROM", "POS", "REF", "ALT","age", "tissue", "TIMES_CALLED", "TYPE", "variantID", "SHARED_ACROSS_SAMPLES", "FILTER", "DPLX_NM","DUPLEX_COV","BAM_COV", "DUPLEX_VAF", "BAM_VAF")]
dataMUT$DPLX_NM
dim(dataMUT) #
table(dataMUT$TYPE)





###########################################
#### *** Assess impact of filters  *** ####
###########################################

### Assess the impact of each filter step on each sample
for(i in macaca_samples){
  print(nrow(dataVCF[dataVCF$TYPE %in% c("snv", "del", "ins") & dataVCF$sampleID==i,]))
  print(nrow(dataVCF[dataVCF$SHARED_ACROSS_SAMPLES==FALSE & dataVCF$sampleID==i,]))
  print(nrow(dataVCF[dataVCF$DPLX_NM %in% c(1, NA) & dataVCF$sampleID==i,]))
  print(nrow(dataVCF[dataVCF$FILTER=="PASS" & dataVCF$sampleID==i,]))
  ### With filters applied:
  print(nrow(dataMUT[dataMUT$sampleID==i,]))
}


### Evaluate impact of filtering on dnds
dndsout = dndscv(dataMUT, refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)

dndscv(dataMUT[dataMUT$TYPE,], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)


dndsout = dndscv(dataMUT[dataMUT$DPLX_NM %in% 1,], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout = dndscv(dataMUT[dataMUT$DPLX_NM %in% 2,], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout = dndscv(dataMUT[dataMUT$DPLX_NM %in% 3,], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout = dndscv(dataMUT[is.na(dataMUT$DPLX_NM),], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)

dndsout = dndscv(dataMUT[!dataMUT$FILTER=="PASS",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)

dndsout = dndscv(dataMUT[dataMUT$RECURRENT==TRUE,], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout = dndscv(dataMUT[dataMUT$RECURRENT==FALSE,], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout = dndscv(dataMUT[dataMUT$RECURRENT==FALSE & dataMUT$DPLX_NM %in% 1,], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)

dndsout$nbreg$theta
dndsout$globaldnds

p6+p6

### Evaluate impact of filtering on correlation with age
params<-c(3,2,1)

plotlist_age<-vector("list", length(params))
names(plotlist_age)<-params

for(i in params){
  print(paste0("Testing param: ", i))

  
  ### Test DPLX_NM
  dataMUT<-dataVCF[dataVCF$SHARED_ACROSS_SAMPLES==FALSE & dataVCF$DPLX_NM %in% i & dataVCF$FILTER %in% "PASS",c("sampleID","CHROM", "POS", "REF", "ALT","age", "tissue", "TIMES_CALLED", "TYPE", "variantID", "SHARED_ACROSS_SAMPLES", "FILTER", "DPLX_NM","DUPLEX_COV","BAM_COV", "DUPLEX_VAF", "BAM_VAF")]
  print(nrow(dataMUT))
  
  # Step 1: count mutations per sample 
  mutation_summary <- dataMUT %>%
  group_by(sampleID,tissue, age) %>%
  summarise(n_mutations = n(), .groups = "drop")

# Step 2: fit linear model
fit <- lm(n_mutations ~ age, data = mutation_summary)

# Step 3: get model stats9
summary_fit <- summary(fit)
intercept <- coef(fit)[1]
slope <- coef(fit)[2]
ci <- confint(fit)  # confidence intervals
cor_val <- cor(mutation_summary$age, mutation_summary$n_mutations, method = "pearson")
p_val <- summary_fit$coefficients[2,4]

# 4. Scatter plot, color by tissue, shape by individual
p<-ggplot(mutation_summary, aes(x = age, y = n_mutations, color = tissue, label=sampleID)) +
  geom_point(size = 5) +
  geom_text(size = 3, vjust = -1, col="black") +
  geom_smooth(method = "lm", se = TRUE, inherit.aes = FALSE,
              aes(x = age, y = n_mutations), color = "black", linetype = "dashed") +
  #coord_cartesian(ylim = c(0, 2200)) +   # Zoom in to 0–0.1 range
  coord_cartesian(ylim = c(0, 600)) +   # Zoom in to 0–0.1 range
  scale_y_continuous(breaks = seq(0, max(mutation_summary$n_mutations, na.rm=TRUE), by = 200))  +
  scale_x_continuous(breaks = seq(0, max(mutation_summary$age, na.rm=TRUE), by = 1))  +
  labs(
    title = paste0("DPLX_NM=", i),
    subtitle = paste0("Number of variants: ", sum(mutation_summary$n_mutations)  , " // Pearson cor=", signif(cor_val,2), " ; pval=", signif(p_val,2)),
    x = "Age",
    y = "Number of Variants",
    color = "Tissue",
    shape = "Sample"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )
print(paste0("Adding plot:", i))
plotlist_age[[i]]<-p
}

pdf("~/Research/postdoc/projects/PrimateSomatic/06_plots/prelim/params/paramsDPLX_NM_age.pdf", width = 10, height = 18)
plotlist_age[[1]] / plotlist_age[[2]] / plotlist_age[[3]]
dev.off()



### Evaluate impact of filtering on dn/ds
params<-c(3,2,1)

plotlist_dnds<-vector("list", length(params))
names(plotlist_dnds)<-params

for(i in params){
  print(paste0("Testing param: ", i))
  
  
  ### Test DPLX_NM
  dataMUT<-dataVCF[dataVCF$SHARED_ACROSS_SAMPLES==FALSE & dataVCF$DPLX_NM %in% i & dataVCF$FILTER %in% "PASS",c("sampleID","CHROM", "POS", "REF", "ALT","age", "tissue", "TIMES_CALLED", "TYPE", "variantID", "SHARED_ACROSS_SAMPLES", "FILTER", "DPLX_NM","DUPLEX_COV","BAM_COV", "DUPLEX_VAF", "BAM_VAF")]
  print(nrow(dataMUT))
  
  ### Compute overall dndsout
  dndsout = dndscv(dataMUT, refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout$globaldnds
  
  ### Compute sample specific dnds
  dndsout_MQD0001d<-dndscv(dataMUT[dataMUT$sampleID=="MQD0001d",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0001j<-dndscv(dataMUT[dataMUT$sampleID=="MQD0001j",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  #dndsout_MQD0002e<-dndscv(dataMUT[dataMUT$sampleID=="MQD0002e",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0002f<-dndscv(dataMUT[dataMUT$sampleID=="MQD0002f",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0002h<-dndscv(dataMUT[dataMUT$sampleID=="MQD0002h",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0002l<-dndscv(dataMUT[dataMUT$sampleID=="MQD0002l",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0003d<-dndscv(dataMUT[dataMUT$sampleID=="MQD0003d",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0003e<-dndscv(dataMUT[dataMUT$sampleID=="MQD0003e",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0003g<-dndscv(dataMUT[dataMUT$sampleID=="MQD0003g",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0003j<-dndscv(dataMUT[dataMUT$sampleID=="MQD0003j",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0004d<-dndscv(dataMUT[dataMUT$sampleID=="MQD0004d",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0004e<-dndscv(dataMUT[dataMUT$sampleID=="MQD0004e",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0004i<-dndscv(dataMUT[dataMUT$sampleID=="MQD0004i",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0004k<-dndscv(dataMUT[dataMUT$sampleID=="MQD0004k",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0006e<-dndscv(dataMUT[dataMUT$sampleID=="MQD0006e",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0006f<-dndscv(dataMUT[dataMUT$sampleID=="MQD0006f",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0006h<-dndscv(dataMUT[dataMUT$sampleID=="MQD0006h",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  dndsout_MQD0006l<-dndscv(dataMUT[dataMUT$sampleID=="MQD0006l",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
  
list_dndsout<-list(dndsout_MQD0001d,
                   dndsout_MQD0001j,
                   #dndsout_MQD0002e,
                   dndsout_MQD0002f,
                   dndsout_MQD0002h,
                   dndsout_MQD0002l,
                   dndsout_MQD0003d,
                   dndsout_MQD0003e,
                   dndsout_MQD0003g,
                   dndsout_MQD0003j,
                   dndsout_MQD0004d,
                   dndsout_MQD0004e,
                   dndsout_MQD0004i,
                   dndsout_MQD0004k,
                   dndsout_MQD0006e,
                   dndsout_MQD0006f,
                   dndsout_MQD0006h,
                   dndsout_MQD0006l)


names(list_dndsout)<-c("MQD0001d",
                       "MQD0001j",
                       #"MQD0002e",
                       "MQD0002f",
                       "MQD0002h",
                       "MQD0002l",
                       "MQD0003d",
                       "MQD0003e",
                       "MQD0003g",
                       "MQD0003j",
                       "MQD0004d",
                       "MQD0004e",
                       "MQD0004i",
                       "MQD0004k",
                       "MQD0006e",
                       "MQD0006f",
                       "MQD0006h",
                       "MQD0006l")


names(list_dndsout)


# Combine all your results into one data frame
globaldnds <- bind_rows(
  dndsout_MQD0001d$globaldnds |> as.data.frame() |> mutate(sample = "MQD0001d"),
  dndsout_MQD0001j$globaldnds |> as.data.frame() |> mutate(sample = "MQD0001j"),
  #dndsout_MQD0002e$globaldnds |> as.data.frame() |> mutate(sample = "MQD0002e"),
  dndsout_MQD0002f$globaldnds |> as.data.frame() |> mutate(sample = "MQD0002f"),
  dndsout_MQD0002h$globaldnds |> as.data.frame() |> mutate(sample = "MQD0002h"),
  dndsout_MQD0002l$globaldnds |> as.data.frame() |> mutate(sample = "MQD0002l"),
  dndsout_MQD0003d$globaldnds |> as.data.frame() |> mutate(sample = "MQD0003d"),
  dndsout_MQD0003e$globaldnds |> as.data.frame() |> mutate(sample = "MQD0003e"),
  dndsout_MQD0003g$globaldnds |> as.data.frame() |> mutate(sample = "MQD0003g"),
  dndsout_MQD0003j$globaldnds |> as.data.frame() |> mutate(sample = "MQD0003j"),
  dndsout_MQD0004d$globaldnds |> as.data.frame() |> mutate(sample = "MQD0004d"),
  dndsout_MQD0004e$globaldnds |> as.data.frame() |> mutate(sample = "MQD0004e"),
  dndsout_MQD0004i$globaldnds |> as.data.frame() |> mutate(sample = "MQD0004i"),
  dndsout_MQD0004k$globaldnds |> as.data.frame() |> mutate(sample = "MQD0004k"),
  dndsout_MQD0006e$globaldnds |> as.data.frame() |> mutate(sample = "MQD0006e"),
  dndsout_MQD0006f$globaldnds |> as.data.frame() |> mutate(sample = "MQD0006f"),
  dndsout_MQD0006h$globaldnds |> as.data.frame() |> mutate(sample = "MQD0006h"),
  dndsout_MQD0006l$globaldnds |> as.data.frame() |> mutate(sample = "MQD0006l"),
)

globaldnds$name <- factor(globaldnds$name, 
                          levels = c("wmis","wnon","wspl","wtru","wall"))

### dnds value
dndsout_value<-signif(dndsout$globaldnds["wall","mle"],2)

### Plot dnds across samples
p<-ggplot(globaldnds, aes(x = sample, y = mle, color = name)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = cilow, ymax = 6),width = 0.2,position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = dndsout_value, linetype = "dashed", color = "red") +
  
  scale_y_continuous(breaks = seq(0, 6, by = 0.5))  +
  #scale_color_manual(values = c(1:9)) +
  labs(
    title = paste0("DPLX_NM=", i),
    subtitle = paste0("Global dN/dS=", dndsout_value),
    x = "Sample",
    y = "w (dN/dS ratio)",
    color = "dN/dS"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16)
  )
p
print(paste0("Adding plot:", i))
plotlist_dnds[[i]]<-p
}

plotlist_dnds<-plotlist
pdf("~/Research/postdoc/projects/PrimateSomatic/06_plots/prelim/params/paramsDPLX_NM_dnds.pdf", width = 10, height = 18)
plotlist_dnds[[1]] / plotlist_dnds[[2]] / plotlist_dnds[[3]]
dev.off()


########################
### Quality Control  ###
########################

library(ggplot2)
library(scales)

#display ggplot2 default hex color codes from 1 to 8
for(i in 1:8){
  print(hue_pal()(i))
}

# Optional: define consistent colors for tissues
tissue_cols <- c(
  "Bladder" = "#F8766D",
  "Liver" = "#A3A500",
  "Oesophagus" = "#00BF7D",
  "Prostate" = "#00B0F6",
  "Skin" = "#E76BF3"
)


### Order by tissue type
dataMUT_tissue <- dataMUT %>%
  mutate(sampleID = factor(sampleID,
                      levels = unique(sampleID[order(tissue)]))
  )

#### Explore Coverage at Site
p3<-ggplot(dataMUT_tissue, aes(x = sampleID, y = DUPLEX_COV, fill = tissue)) +
    geom_violin(trim = FALSE, alpha = 0.4) +
    geom_boxplot(width = 0.1, outlier.size = 0.5) +
    #facet_wrap(~ tissue, scales = "free_x") +
    labs(title = "Coverage at Candidate Somatic Variant Sites", x = "Sample", y = "Duplex Coverage") +
    theme_minimal(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
pdf("~/Research/postdoc/projects/PrimateSomatic/Analysis/Macaca_mulatta/targeted_nanoseq/v0.1/QC/macaca_targeted_nanoseq_DUPLEX_COV.pdf", width = 12, height = 6)
grid.arrange(p3)
dev.off()

p4<-ggplot(dataMUT_tissue, aes(x = sampleID, y = BAM_COV, fill = tissue)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  labs(title = "Coverage at Candidate Somatic Variant Sites", x = "Sample", y = "BAM Coverage") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf("~/Research/postdoc/projects/PrimateSomatic/Analysis/Macaca_mulatta/targeted_nanoseq/v0.1/QC/macaca_targeted_nanoseq_BAM_COV.pdf", width = 12, height = 6)
grid.arrange(p4)
dev.off()


#### Explore DUPLEX VAF

# Create all plots first and store in a list
plots_DUPLEX_VAF <- list()

for (i in levels(dataMUT_tissue$sampleID)) {
  message("Processing sample: ", i)
  
  df <- dataMUT_tissue[dataMUT_tissue$sampleID == i, ]
  
  p1 <- ggplot(df, aes(x = DUPLEX_VAF, fill = tissue)) +
    geom_density(alpha = 0.5, adjust = 1) +
    coord_cartesian(xlim = c(0, 0.1)) +
    scale_fill_manual(values = tissue_cols, drop = FALSE) +
    labs(
      title = i,
      x = "DUPLEX VAF",
      y = "Density",
      subtitle = paste0(
        unique(df$tissue), " // ",
        nrow(df), " variants"
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "none"  # redundant since each plot is one tissue
    )
  
  plots_DUPLEX_VAF[[i]] <- p1
}

# Write to PDF with 18 plots per page (3 x 6)
pdf("~/Research/postdoc/projects/PrimateSomatic/Analysis/Macaca_mulatta/targeted_nanoseq/v0.1/QC/macaca_targeted_nanoseq_DUPLEX_VAF.pdf", width = 12, height = 6)
grid.arrange(grobs = plots_DUPLEX_VAF,  nrow = 3, ncol=6)
dev.off()



load("/Users/ms84/Research/postdoc/projects/PrimateSomatic/RefCDS_Macaca_mulatta.Mmul_10.rda")
str(RefCDS)


# Step 1: count mutations per sample
mutation_summary <- dataMUT %>%
  group_by(sampleID,tissue, age) %>%
  summarise(n_mutations = n(), .groups = "drop")

# Step 2: fit linear model
fit <- lm(n_mutations ~ age, data = mutation_summary)

# Step 3: get model stats9
summary_fit <- summary(fit)
intercept <- coef(fit)[1]
slope <- coef(fit)[2]
ci <- confint(fit)  # confidence intervals
cor_val <- cor(mutation_summary$age, mutation_summary$n_mutations, method = "pearson")
p_val <- summary_fit$coefficients[2,4]


# 4. Scatter plot, color by tissue, shape by individual
p6<-ggplot(mutation_summary, aes(x = age, y = n_mutations, color = tissue, label=sampleID)) +
  geom_point(size = 5) +
  geom_text(size = 3, vjust = -1, col="black") +
  geom_smooth(method = "lm", se = TRUE, inherit.aes = FALSE,
              aes(x = age, y = n_mutations), color = "black", linetype = "dashed") +
  #coord_cartesian(ylim = c(0, 2200)) +   # Zoom in to 0–0.1 range
  coord_cartesian(ylim = c(0, 700)) +   # Zoom in to 0–0.1 range
  scale_y_continuous(breaks = seq(0, max(mutation_summary$n_mutations, na.rm=TRUE), by = 200))  +
  scale_x_continuous(breaks = seq(0, max(mutation_summary$age, na.rm=TRUE), by = 1))  +
  scale_color_manual(values = tissue_cols) +
  labs(
    title = "Correlation between number of variants & age",
    subtitle = paste0("Pearson's correlation=", signif(cor_val,2), " // p-val=", signif(p_val,2)),
    x = "Age",
    y = "Number of Variants",
    color = "Tissue",
    shape = "Sample"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

pdf("~/Research/postdoc/projects/PrimateSomatic/Analysis/Macaca_mulatta/targeted_nanoseq/v0.1/QC/macaca_targeted_nanoseq_Age-Variants.pdf", width = 10, height = 6)
grid.arrange(p6)
dev.off()


##########################
#### DNDS_CV ANALYSIS ####
##########################

dndsout_MQD0001d<-dndscv(dataMUT[dataMUT$sampleID=="MQD0001d",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0001j<-dndscv(dataMUT[dataMUT$sampleID=="MQD0001j",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0002e<-dndscv(dataMUT[dataMUT$sampleID=="MQD0002e",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0002f<-dndscv(dataMUT[dataMUT$sampleID=="MQD0002f",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0002h<-dndscv(dataMUT[dataMUT$sampleID=="MQD0002h",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0002l<-dndscv(dataMUT[dataMUT$sampleID=="MQD0002l",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0003d<-dndscv(dataMUT[dataMUT$sampleID=="MQD0003d",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0003e<-dndscv(dataMUT[dataMUT$sampleID=="MQD0003e",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0003g<-dndscv(dataMUT[dataMUT$sampleID=="MQD0003g",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0003j<-dndscv(dataMUT[dataMUT$sampleID=="MQD0003j",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0004d<-dndscv(dataMUT[dataMUT$sampleID=="MQD0004d",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0004e<-dndscv(dataMUT[dataMUT$sampleID=="MQD0004e",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0004i<-dndscv(dataMUT[dataMUT$sampleID=="MQD0004i",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0004k<-dndscv(dataMUT[dataMUT$sampleID=="MQD0004k",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0006e<-dndscv(dataMUT[dataMUT$sampleID=="MQD0006e",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0006f<-dndscv(dataMUT[dataMUT$sampleID=="MQD0006f",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0006h<-dndscv(dataMUT[dataMUT$sampleID=="MQD0006h",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout_MQD0006l<-dndscv(dataMUT[dataMUT$sampleID=="MQD0006l",], refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)

list_dndsout<-list(dndsout_MQD0001d,
                   dndsout_MQD0001j,
                   dndsout_MQD0002e,
                   dndsout_MQD0002f,
                   dndsout_MQD0002h,
                   dndsout_MQD0002l,
                   dndsout_MQD0003d,
                   dndsout_MQD0003e,
                   dndsout_MQD0003g,
                   dndsout_MQD0003j,
                   dndsout_MQD0004d,
                   dndsout_MQD0004e,
                   dndsout_MQD0004i,
                   dndsout_MQD0004k,
                   dndsout_MQD0006e,
                   dndsout_MQD0006f,
                   dndsout_MQD0006h,
                   dndsout_MQD0006l)

names(list_dndsout)<-c("MQD0001d",
                          "MQD0001j",
                          "MQD0002e",
                          "MQD0002f",
                          "MQD0002h",
                          "MQD0002l",
                          "MQD0003d",
                          "MQD0003e",
                          "MQD0003g",
                          "MQD0003j",
                          "MQD0004d",
                          "MQD0004e",
                          "MQD0004i",
                          "MQD0004k",
                          "MQD0006e",
                          "MQD0006f",
                          "MQD0006h",
                          "MQD0006l")


names(list_dndsout)


outdir="/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/MITS_2026/dndscv/"
for (sample in names(list_dndsout)) {
  
  dnds_obj <- list_dndsout[[sample]]
  
  # Choose the appropriate results table
  if (!is.null(dnds_obj$sel_cv)) {
    out_table <- dnds_obj$sel_cv
  } else {
    out_table <- dnds_obj$sel
  }
  
  write.xlsx(
    out_table,
    file = paste0(outdir,"dndsout_", sample, ".xlsx"),
    rowNames = FALSE
  )
}

### Evaluation of the that values
dndsout_MQD0001d$nbreg$theta
dndsout_MQD0001j$nbreg$theta
dndsout_MQD0002e$nbreg$theta
dndsout_MQD0002f$nbreg$theta
dndsout_MQD0002h$nbreg$theta
dndsout_MQD0002l$nbreg$theta
dndsout_MQD0003d$nbreg$theta
dndsout_MQD0003e$nbreg$theta
dndsout_MQD0003g$nbreg$theta
dndsout_MQD0003j$nbreg$theta
dndsout_MQD0004d$nbreg$theta
dndsout_MQD0004e$nbreg$theta
dndsout_MQD0004i$nbreg$theta
dndsout_MQD0004k$nbreg$theta
dndsout_MQD0006e$nbreg$theta
dndsout_MQD0006f$nbreg$theta
dndsout_MQD0006h$nbreg$theta
dndsout_MQD0006l$nbreg$theta

### Export results
OUTDIR="/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/dnds"
for(i in names(list_dndsout)[7:length(list_dndsout)]){
  write.xlsx(list_dndsout[[i]], paste0(OUTDIR, "/dndsout_",i,".xlsx"))
}


dndsout = dndscv(dataMUT, refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)


library(dplyr)
library(ggplot2)

# Combine all your results into one data frame
globaldnds <- bind_rows(
  dndsout_MQD0001d$globaldnds |> as.data.frame() |> mutate(sample = "MQD0001d"),
  dndsout_MQD0001j$globaldnds |> as.data.frame() |> mutate(sample = "MQD0001j"),
  dndsout_MQD0002e$globaldnds |> as.data.frame() |> mutate(sample = "MQD0002e"),
  dndsout_MQD0002f$globaldnds |> as.data.frame() |> mutate(sample = "MQD0002f"),
  dndsout_MQD0002h$globaldnds |> as.data.frame() |> mutate(sample = "MQD0002h"),
  dndsout_MQD0002l$globaldnds |> as.data.frame() |> mutate(sample = "MQD0002l"),
  dndsout_MQD0003d$globaldnds |> as.data.frame() |> mutate(sample = "MQD0003d"),
  dndsout_MQD0003e$globaldnds |> as.data.frame() |> mutate(sample = "MQD0003e"),
  dndsout_MQD0003g$globaldnds |> as.data.frame() |> mutate(sample = "MQD0003g"),
  dndsout_MQD0003j$globaldnds |> as.data.frame() |> mutate(sample = "MQD0003j"),
  dndsout_MQD0004d$globaldnds |> as.data.frame() |> mutate(sample = "MQD0004d"),
  dndsout_MQD0004e$globaldnds |> as.data.frame() |> mutate(sample = "MQD0004e"),
  dndsout_MQD0004i$globaldnds |> as.data.frame() |> mutate(sample = "MQD0004i"),
  dndsout_MQD0004k$globaldnds |> as.data.frame() |> mutate(sample = "MQD0004k"),
  dndsout_MQD0006e$globaldnds |> as.data.frame() |> mutate(sample = "MQD0006e"),
  dndsout_MQD0006f$globaldnds |> as.data.frame() |> mutate(sample = "MQD0006f"),
  dndsout_MQD0006h$globaldnds |> as.data.frame() |> mutate(sample = "MQD0006h"),
  dndsout_MQD0006l$globaldnds |> as.data.frame() |> mutate(sample = "MQD0006l"),
)

globaldnds$name <- factor(globaldnds$name, 
                          levels = c("wmis","wnon","wspl","wtru","wall"))

globaldnds[grepl("wall", rownames(globaldnds)),"mle"]

### dnds value
dndsout_value<-signif(dndsout$globaldnds["wall","mle"],2)


#### Prepare for dn/ds plots

### add tissue metadata
sample_annot <- dataMUT %>% select(sampleID, tissue, age) %>% distinct() %>% dplyr::rename(sample=sampleID)

### join with your dN/dS table
globaldnds <- globaldnds %>% left_join(sample_annot, by = "sample")

## order samples by tissue (and age)
globaldnds <- globaldnds %>%
  arrange(tissue, age) %>%
  mutate(sample = factor(sample, levels = unique(sample)))

## Define shapes for dN/dS types (highlight wall)
shape_values <- c(
  "wmis" = 21,  # circle
  "wnon" = 24,  # triangle
  "wspl" = 22,  # square
  "wtru" = 23,  # diamond
  "wall" = 8    # keep star (no fill, but already distinct)
)

## Modified plot
pd <- position_dodge(width = 0.6)

p6 <- ggplot(globaldnds, aes(x = sample, y = mle)) +
  
  ## Error bars FIRST (background)
  geom_errorbar(
    aes(ymin = cilow, ymax = cihigh, color = tissue, group = name),
    width = 0.2,
    position = pd
  ) +
  
  ## Points SECOND (on top)
  geom_point(
    aes(fill = tissue, shape = name, group = name),
    color = "black",
    position = pd,
    size = 3,
    stroke = 0.8
  ) +
  
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  
  coord_cartesian(ylim = c(0, 3)) +
  
  scale_fill_manual(values = tissue_cols) +
  scale_color_manual(values = tissue_cols) +
  scale_shape_manual(values = shape_values) +
  
  labs(
    title = "Global dN/dS (w) across primate samples",
    subtitle = paste0("Global dN/dS: ", dndsout_value),
    x = "Sample",
    y = "w (dN/dS ratio)",
    fill = "Tissue",
    shape = "dN/dS type"
  ) +
  
  guides(color = "none") +  # remove duplicate legend
  
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16)
  )
p6
pdf("~/Research/postdoc/projects/PrimateSomatic/Analysis/Macaca_mulatta/targeted_nanoseq/v0.1/QC/macaca_targeted_nanoseq_dNdS-samples.filtered_MITS2026.pdf", width = 12, height = 7)
grid.arrange(p6)
dev.off()

p6<-ggplot(globaldnds, aes(x = sample, y = mle, color = name)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  #geom_hline(yintercept = dndsout_value, linetype = "dashed", color = "red") +
  
  coord_cartesian(ylim = c(0, 3)) +
  #scale_y_continuous(breaks = seq(0, 2, by = 0.2))  +
  #scale_color_manual(values = sigfit_palette) +
  labs(
    title = "Global dN/dS (w) across primate samples",
    subtitle = paste0("Global dN/dS: ", dndsout_value),
    x = "Sample",
    y = "w (dN/dS ratio)",
    color = "dN/dS"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16)
  )
p6
pdf("~/Research/postdoc/projects/PrimateSomatic/Analysis/Macaca_mulatta/targeted_nanoseq/v0.1/QC/macaca_targeted_nanoseq_dNdS-samples.filtered.pdf", width = 12, height = 7)
grid.arrange(p6)
dev.off()




### dndscv outputs: Table of significant genes
sel_cv = dndsout$sel_cv
signif_genes<-sel_cv[sel_cv$qglobal_cv<0.1,]
nrow(signif_genes) # 130 genes
top_signif_genes<-signif_genes[signif_genes$qglobal_cv==0,c(1:10,20)]
top_signif_genes
write.xlsx(top_signif_genes,"/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/RG25/top_significant_genes.xlsx")


print(sel_cv[sel_cv$qglobal_cv<0.1,c(1:10,19)], digits = 3)

print(head(sel_cv), digits = 3)

sel_cv_MQD0001d<-dndsout_MQD0001d$sel_cv
sel_cv_MQD0001j<-dndsout_MQD0001j$sel_cv
sel_cv_MQD0002e<-dndsout_MQD0002e$sel_cv
sel_cv_MQD0002f<-dndsout_MQD0002f$sel_cv
sel_cv_MQD0002h<-dndsout_MQD0002h$sel_cv
sel_cv_MQD0002l<-dndsout_MQD0002l$sel_cv
sel_cv_MQD0003d<-dndsout_MQD0003d$sel_cv
sel_cv_MQD0003e<-dndsout_MQD0003e$sel_cv
sel_cv_MQD0003g<-dndsout_MQD0003g$sel_cv
sel_cv_MQD0003j<-dndsout_MQD0003j$sel_cv
sel_cv_MQD0004d<-dndsout_MQD0004d$sel_cv
sel_cv_MQD0004e<-dndsout_MQD0004e$sel_cv
sel_cv_MQD0004i<-dndsout_MQD0004i$sel_cv
sel_cv_MQD0004k<-dndsout_MQD0004k$sel_cv
sel_cv_MQD0006e<-dndsout_MQD0006e$sel_cv
sel_cv_MQD0006f<-dndsout_MQD0006f$sel_cv
sel_cv_MQD0006h<-dndsout_MQD0006h$sel_cv
sel_cv_MQD0006l<-dndsout_MQD0006l$sel_cv


signif_genes_MQD0001d = sel_cv_MQD0001d[sel_cv_MQD0001d$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0001j = sel_cv_MQD0001j[sel_cv_MQD0001j$qglobal_cv<0.1, c("gene_name","qglobal_cv")]

signif_genes_MQD0002e = sel_cv_MQD0002e[sel_cv_MQD0002e$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0002f = sel_cv_MQD0002f[sel_cv_MQD0002f$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0002h = sel_cv_MQD0002h[sel_cv_MQD0002h$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0002l = sel_cv_MQD0002l[sel_cv_MQD0002h$qglobal_cv<0.1, c("gene_name","qglobal_cv")]

signif_genes_MQD0003d = sel_cv_MQD0003d[sel_cv_MQD0003d$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0003e = sel_cv_MQD0003e[sel_cv_MQD0003e$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0003g = sel_cv_MQD0003g[sel_cv_MQD0003g$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0003j = sel_cv_MQD0003j[sel_cv_MQD0003j$qglobal_cv<0.1, c("gene_name","qglobal_cv")]

signif_genes_MQD0004d = sel_cv_MQD0004d[sel_cv_MQD0004d$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0004e = sel_cv_MQD0004e[sel_cv_MQD0004e$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0004i = sel_cv_MQD0004i[sel_cv_MQD0004i$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0004k = sel_cv_MQD0004k[sel_cv_MQD0004k$qglobal_cv<0.1, c("gene_name","qglobal_cv")]

signif_genes_MQD0006e = sel_cv_MQD0006e[sel_cv_MQD0006e$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0006f = sel_cv_MQD0006f[sel_cv_MQD0006f$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0006h = sel_cv_MQD0006h[sel_cv_MQD0006h$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
signif_genes_MQD0006l = sel_cv_MQD0006l[sel_cv_MQD0006l$qglobal_cv<0.1, c("gene_name","qglobal_cv")]



### Print significant genes
signif_genes = sel_cv[sel_cv$qglobal_cv<0.05, c("gene_name","qglobal_cv")]
nrow(signif_genes) # 126 genes
rownames(signif_genes) = NULL
print(signif_genes)


### Annotated
dim(dndsout$annotmuts)
head(dndsout$annotmuts)
prop.table(table(dndsout$annotmuts$impact))



dndsout = dndscv(dataMUT, refdb="RefCDS_Macaca_mulatta.Mmul_10.rda", cv=NULL)
dndsout
str(dndsout)

dndsout$sel_cv

### Global test of selection
print(dndsout$globaldnds)



#GRdataMUT<-toGRanges(dataMUT[,c("CHROM","POS", "POS")])

GRdataMUT<-toGRanges(dndscv[,c("CHROM","POS", "POS")])

dim(dataMUT)

### Create dataANNOTMUTS
dataANNOTMUTS<-dndsout$annotmuts

dataVCF$mutID<-paste0(dataVCF$CHROM,"_",dataVCF$POS,"_",dataVCF$REF,"_",dataVCF$ALT, "_", dataVCF$sampleID)
dataANNOTMUTS$mutID<-paste0(dataANNOTMUTS$chr,"_",dataANNOTMUTS$pos,"_",dataANNOTMUTS$ref,"_",dataANNOTMUTS$mut, "_", dataANNOTMUTS$sampleID)

dataANNOTMUTS$DUPLEX_COV<-dataVCF[match(dataANNOTMUTS$mutID, dataVCF$mutID),"DUPLEX_COV"]
dataANNOTMUTS$DUPLEX_VAF<-dataVCF[match(dataANNOTMUTS$mutID, dataVCF$mutID),"DUPLEX_VAF"]
dataANNOTMUTS$tissue<-dataVCF[match(dataANNOTMUTS$sampleID, dataVCF$sampleID),"tissue"]
dataANNOTMUTS$age<-dataVCF[match(dataANNOTMUTS$sampleID, dataVCF$sampleID),"age"]
dataANNOTMUTS
head(dataANNOTMUTS)

# Create a frequency table
impact_counts <- table(dataANNOTMUTS$impact)

pdf("~/Research/postdoc/projects/PrimateSomatic/06_plots/prelim/macaca_targeted_nanoseq_circlechart.pdf", width = 7, height = 7)

# Create a pie chart
pie(
  impact_counts,
  main = "Impact of Candidate Somatic Variants",
  col = test_palette,
  labels = paste(names(impact_counts), "\n", impact_counts)
)

# Optional: add a legend
legend(
  "bottomright",
  legend = names(impact_counts),
  fill = test_palette,
  cex = 0.8
)
dev.off()


table(dataANNOTMUTS$impact)
#top_impact<-c("Essential_Splice", "Nonsense", "Stop_loss")
top_impact<-c("Missense","Essential_Splice", "Nonsense", "Stop_loss","Missense", "no-SNV")
top_impact

dataIMPACT<-dataANNOTMUTS[dataANNOTMUTS$impact %in% top_impact,]
dim(dataIMPACT)


### Add COSMIC - Cancer Gene Census information
dataCOSMIC_cgc<-read.delim("~/volumes/ms84_lustre/scratch126/casm/teams/team294/projects/cseg/resources/databases/cosmic/cancer_gene_census/Cosmic_CancerGeneCensus_v103_GRCh38.tsv")
head(dataCOSMIC_cgc)

dataIMPACT$COSMIC_gene_id<-dataCOSMIC_cgc[match(dataIMPACT$gene, dataCOSMIC_cgc$GENE_SYMBOL),"COSMIC_GENE_ID"]
dataIMPACT$COSMIC_name<-dataCOSMIC_cgc[match(dataIMPACT$gene, dataCOSMIC_cgc$GENE_SYMBOL),"NAME"]
dataIMPACT$COSMIC_tumour_types_somatic<-dataCOSMIC_cgc[match(dataIMPACT$gene, dataCOSMIC_cgc$GENE_SYMBOL),"TUMOUR_TYPES_SOMATIC"]
dataIMPACT$COSMIC_tumour_types_germline<-dataCOSMIC_cgc[match(dataIMPACT$gene, dataCOSMIC_cgc$GENE_SYMBOL),"TUMOUR_TYPES_GERMLINE"]
dataIMPACT$COSMIC_tissue_type<-dataCOSMIC_cgc[match(dataIMPACT$gene, dataCOSMIC_cgc$GENE_SYMBOL),"TISSUE_TYPE"]
dataIMPACT$COSMIC_role_in_cancer<-dataCOSMIC_cgc[match(dataIMPACT$gene, dataCOSMIC_cgc$GENE_SYMBOL),"ROLE_IN_CANCER"]
dataIMPACT$COSMIC_mutation_types<-dataCOSMIC_cgc[match(dataIMPACT$gene, dataCOSMIC_cgc$GENE_SYMBOL),"MUTATION_TYPES"]

### Reorganize the data.frames with dplyr
dataIMPACTsorted<-dataIMPACT %>%
  mutate(age=as.numeric(age), tissue=factor(tissue, levels=c("Liver","Bladder", "Oesophagus", "Prostate", "Skin")), impact = factor(impact, levels = c("Essential_Splice", "Nonsense", "Stop_loss","Missense", "no-SNV"))) %>%
  arrange(tissue, impact, age) 

head(dataIMPACTsorted)


dataWRITE<-dataIMPACTsorted[,c("sampleID","age","tissue","gene","impact","DUPLEX_VAF", "DUPLEX_COV","chr","pos","ntchange","aachange","codonsub", "COSMIC_gene_id", "COSMIC_name",  "COSMIC_tumour_types_somatic", "COSMIC_tumour_types_germline", "COSMIC_tissue_type", "COSMIC_role_in_cancer", "COSMIC_mutation_types")]
dim(dataWRITE)
dim(dataWRITE)
write.xlsx(dataWRITE, "~/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/MITS_2026/somatic_mutations_targeted_nanoseq_macaca.xlsx")

### Set duplex vaf threshold based on quantiles
hist(dataIMPACTsorted$DUPLEX_VAF, breaks = 100)
duplex_vaf_thres<-quantile(dataIMPACTsorted$DUPLEX_VAF, probs=seq(0, 1, 0.1))["0%"]
duplex_vaf_thres

### Review COSMIC somatic tumour terminology
Liver_tumours_types_somatic<-c("liver","hepatoblastoma", "cholangiocarcinoma", "hepatocellular", "hepatocellular carcinoma", "hepatocellular carcinoma", "hepatic adenoma","biliary tract", "biliary tract carcinoma", "gallbladder carcinoma")
Bladder_tumours_types_somatic<-c("bladder", "bladder cancer", "bladder carcinoma", "urothelial cancer", "urothelial cell carcinoma")
Oesophagus_tumours_types_somatic<-c("oesophagus", "oesophagus cancer", "oesophageal SCC", "oesophageal squamous cell carcinoma")
Prostate_tumours_types_somatic<-c("prostate", "prostate cancer", "prostate carcinoma", "prostae adenocarcinoma")
Skin_tumours_types_somatic<-c("skin", "skin cancer", "skin basal cell", "skin basal cell carcinoma", "basal cell carcinoma", "skin SCC", "skin squamous cell", "skin squamous cell carcinoma-burn scar related", "squamous cell carcinoma", "cutaneous melanoma", "desmoplastic melanoma", "malignant melanoma", "melanoma", "mucosal melanoma", "uveal melanoma", "Spitzoid tumour", "melanocytic nevus")


### Select candidate drivers for each tissue type
candidate_drivers_Liver<-dataIMPACTsorted[dataIMPACTsorted$tissue=="Liver" & grepl(paste(Liver_tumours_types_somatic,collapse="|"), dataIMPACTsorted$COSMIC_tumour_types_somatic, ignore.case = TRUE) &
                                                       ((dataIMPACTsorted$impact=="Essential_Splice" & grepl("S", dataIMPACTsorted$COSMIC_mutation_types))|(dataIMPACTsorted$impact=="Nonsense" & grepl("N", dataIMPACTsorted$COSMIC_mutation_types))|(dataIMPACTsorted$impact=="Missense" & grepl("Mis", dataIMPACTsorted$COSMIC_mutation_types))) &
                                                        dataIMPACTsorted$DUPLEX_VAF >= duplex_vaf_thres,] 
candidate_drivers_Bladder<-dataIMPACTsorted[dataIMPACTsorted$tissue=="Bladder" & grepl(paste(Bladder_tumours_types_somatic,collapse="|"), dataIMPACTsorted$COSMIC_tumour_types_somatic, ignore.case = TRUE) &
                                                      ((dataIMPACTsorted$impact=="Essential_Splice" & grepl("S", dataIMPACTsorted$COSMIC_mutation_types))|(dataIMPACTsorted$impact=="Nonsense" & grepl("N", dataIMPACTsorted$COSMIC_mutation_types))|(dataIMPACTsorted$impact=="Missense" & grepl("Mis", dataIMPACTsorted$COSMIC_mutation_types))) &
                                                        dataIMPACTsorted$DUPLEX_VAF >= duplex_vaf_thres,]
candidate_drivers_Oesophagus<-dataIMPACTsorted[dataIMPACTsorted$tissue=="Oesophagus" & grepl(paste(Oesophagus_tumours_types_somatic,collapse="|"), dataIMPACTsorted$COSMIC_tumour_types_somatic, ignore.case = TRUE) &
                                                      #((dataIMPACTsorted$impact=="Essential_Splice" & grepl("S", dataIMPACTsorted$COSMIC_mutation_types))|(dataIMPACTsorted$impact=="Nonsense" & grepl("N", dataIMPACTsorted$COSMIC_mutation_types))|(dataIMPACTsorted$impact=="Missense" & grepl("Mis", dataIMPACTsorted$COSMIC_mutation_types))) &
                                                      dataIMPACTsorted$DUPLEX_VAF >= duplex_vaf_thres,]
candidate_drivers_Prostate<-dataIMPACTsorted[dataIMPACTsorted$tissue=="Prostate" & grepl(paste(Prostate_tumours_types_somatic,collapse="|"), dataIMPACTsorted$COSMIC_tumour_types_somatic, ignore.case = TRUE) &
                                                      ((dataIMPACTsorted$impact=="Essential_Splice" & grepl("S", dataIMPACTsorted$COSMIC_mutation_types))|(dataIMPACTsorted$impact=="Nonsense" & grepl("N", dataIMPACTsorted$COSMIC_mutation_types))|(dataIMPACTsorted$impact=="Missense" & grepl("Mis", dataIMPACTsorted$COSMIC_mutation_types))) &
                                                      dataIMPACTsorted$DUPLEX_VAF >= duplex_vaf_thres,]
candidate_drivers_Skin<-dataIMPACTsorted[dataIMPACTsorted$tissue=="Skin" & grepl(paste(Skin_tumours_types_somatic,collapse="|"), dataIMPACTsorted$COSMIC_tumour_types_somatic, ignore.case = TRUE) &
                                                      ((dataIMPACTsorted$impact=="Essential_Splice" & grepl("S", dataIMPACTsorted$COSMIC_mutation_types))|(dataIMPACTsorted$impact=="Nonsense" & grepl("N", dataIMPACTsorted$COSMIC_mutation_types))|(dataIMPACTsorted$impact=="Missense" & grepl("Mis", dataIMPACTsorted$COSMIC_mutation_types))) &
                                                      dataIMPACTsorted$DUPLEX_VAF >= duplex_vaf_thres,]


candidate_drivers_Liver[,c("gene","impact")]
candidate_drivers_Bladder[,c("gene","impact")]
candidate_drivers_Oesophagus[,c("gene","impact")]
candidate_drivers_Prostate[,c("gene","impact")]
candidate_drivers_Skin[,c("gene","impact")]

relevant_fields<-c("sampleID","age","tissue","gene","impact","DUPLEX_VAF", "DUPLEX_COV","chr","pos","ntchange","aachange","codonsub", "COSMIC_tumour_types_somatic",  "COSMIC_role_in_cancer", "COSMIC_mutation_types")
write.xlsx(candidate_drivers_Liver[,relevant_fields], "~/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/MITS_2026/candidate_drivers_Liver.xlsx")
write.xlsx(candidate_drivers_Bladder[,relevant_fields], "~/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/MITS_2026/candidate_drivers_Bladder.xlsx")
write.xlsx(candidate_drivers_Oesophagus[,relevant_fields], "~/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/MITS_2026/candidate_drivers_Oesophagus.xlsx")
write.xlsx(candidate_drivers_Prostate[,relevant_fields], "~/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/MITS_2026/candidate_drivers_Prostate.xlsx")
write.xlsx(candidate_drivers_Skin[,relevant_fields], "~/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/MITS_2026/candidate_drivers_Skin.xlsx")


### Load packages
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

table(dataIMPACTsorted$impact)

## Prepare the data
df <- dataIMPACTsorted %>%
  select(sampleID, tissue, gene, impact, age) %>%
  distinct()

### Keep only significant genes
signif_genes_ranked<-signif_genes$gene_name

#dim(df)
#df<-df[df$gene %in% signif_genes$gene_name,]
#dim(df)

## Define impact severity ranking
impact_rank <- c(
  "None"= 0,
  "Missense"= 1,
  "no-SNV" = 2,
  "Essential_Splice" = 3,
  "Nonsense" = 4,
  "Stop_loss" = 4
)

## Define colours for impact types
impact_colors <- c(
  "None"              = "white",
  "Missense"          = "#FDB462",
  "Essential_Splice"  = "#377EB8",
  "Nonsense"          = "#E41A1C",
  "Stop_loss"         = "#984EA3",
  "no-SNV"= "darkgray"
)

impact_colors <- c(
  "None"              = "white",     # keep neutral
  "Missense"          = "#80B1D3",   # soft blue
  "Essential_Splice"  = "#1F78B4",   # strong blue
  "Nonsense"          = "#00441B",   # deep green (high severity)
  "Stop_loss"         = "#081D58",   # dark navy (high severity)
  "no-SNV"            = "#4D4D4D"    # dark grey (neutral but visible)
)

impact_colors <- c(
  "None"              = "#F7F7F7",  # almost white
  "Missense"          = "#F4A3B4",  # soft eosin pink
  "Essential_Splice"  = "#D65A9A",  # magenta (transition)
  "Nonsense"          = "#7B3294",  # hematoxylin purple
  "Stop_loss"         = "#4A1486",  # deep purple (most severe)
  "no-SNV"            = "darkgray"   # Dark grey
)

## collapse multiple variants per sample-gene (! Warning)
df_collapsed <- df %>%
  mutate(
  impact=as.character(impact),
  impact_rank = impact_rank[impact]) %>%
  group_by(tissue, age, sampleID, gene) %>%
  summarise(
    impact = impact[which.max(impact_rank)],
    .groups = "drop"
  )


## Pivot to wide format:
heatmap_df <- df_collapsed %>%
  pivot_wider(
    names_from = gene,
    values_from = impact,
    values_fill = list(impact = "None")
  )

## Create matrix 
mat <- heatmap_df %>%
  select(-sampleID, -tissue, -age) %>%
  as.matrix()

rownames(mat) <- heatmap_df$tissue

## Ensure genes are in desired order
#ranking_genes <- signif_genes$gene_name
ranking_genes <- signif_genes$gene_name[1:100]

## Keep only genes present in matrix (avoids errors)
ranking_genes <- ranking_genes[ranking_genes%in% colnames(mat)]

## Reorder columns
mat <- mat[, ranking_genes, drop = FALSE]
dim(mat)


### Map colors to each sample
age_bar_cols <- tissue_cols[heatmap_df$tissue]
### Create row annotation (tissue)
age_annotation <- rowAnnotation(
  Age = anno_barplot(
    heatmap_df$age,
    gp = gpar(fill = age_bar_cols, col = NA),
    border = FALSE
  )
)

### Plot heatmap
p1<-Heatmap(
  mat,
  name = "Impact",
  col = impact_colors,
  na_col = "white",
  right_annotation = age_annotation,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 5),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_title = "Samples",
  column_title = "Genes",
  border=TRUE,
  rect_gp=gpar(col = "gray"),
  row_split = heatmap_df$tissue
)
pdf("~/Research/postdoc/projects/PrimateSomatic/Analysis/Macaca_mulatta/targeted_nanoseq/v0.1/drivers/driver_screening_heatmap_MITS2026.pdf", width = 16, height = 5)
p1
dev.off()



save.image(file = "~/Research/postdoc/projects/PrimateSomatic/03_scripts/r/environments/PrimateSomatic_dndscv.Rd")






















#### (IMPROVE THIS !!!)
dataANNOTMUTS_impact<-dataANNOTMUTS_impact[!duplicated(dataANNOTMUTS_impact$gene),]
dim(dataANNOTMUTS_impact)

GRdataMUT_impact<-toGRanges(dataANNOTMUTS_impact[,c("chr","pos","pos","gene")])
GRdataMUT_impact


library("openxlsx")
table(dataANNOTMUTS_impact$impact)

dataWRITE<-dataANNOTMUTS_impact[,c("gene","sampleID","age","tissue","impact","DUPLEX_VAF","chr","pos","ntchange","aachange","codonsub")]
duplex_vaf_thres<-0.001162
dataWRITE<-dataWRITE[dataWRITE$DUPLEX_VAF>duplex_vaf_thres,]
dim(dataWRITE)
dataWRITE$gene
#write.xlsx(dataANNOTMUTS_impact, "~/Research/postdoc/projects/PrimateSomatic/06_plots/prelim/macaca_targeted_ANNOTMUTS_impact.filtered.xlsx")
write.xlsx(dataWRITE, "~/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/RG25/macaca_targeted_ANNOTMUTS_impact.filtered.xlsx")
dim(dataWRITE)












###############################
#### KaryoploteR Analysis #####
###############################

### Load reference genome report file & select chromosomes
dataREPORT<-read.delim("~/Research/postdoc/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/Mmul_10_sequence_report.tsv", header = TRUE)
dataREPORT[1:22,]
chromosomes<-dataREPORT[dataREPORT$UCSC.style.name %in% paste0("chr", c(1:20,"X","Y")),"RefSeq.seq.accession"]

dataSIZE<-read.delim("~/Research/postdoc/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/Mmul_10_sequence_report.bed", header = FALSE)
head(dataSIZE)

### Load Gene annotation from ENSEMBL
dataGENE<-read.delim("~/Downloads/Mmul_10_ensGene.txt")
dataGENE$chrom<-dataREPORT[match(dataGENE$chrom, dataREPORT$UCSC.style.name), "RefSeq.seq.accession"]
dataGENE<-dataGENE[dataGENE$chrom %in% chromosomes,]
GRdataGENE<-toGRanges(dataGENE[,c("chrom", "txStart", "txEnd")])

### Load target panel distribution
dataPANEL<-read.delim("~/Research/postdoc/tools/Targeted_NanoSeq/panel/targeted_panel_Mmul_10.sorted.merged_500.bed", header = FALSE)
names(dataPANEL)<-c("chr","start","end")
GRdataPANEL<-toGRanges(dataPANEL)

### Set mutation calls
GRdataMUT<-toGRanges(dataMUT[,c("CHROM","POS","POS","TYPE")])

### Set Marker genes
head(dataWRITE)
GRdataWRITE<-toGRanges(dataWRITE[,c("chr","pos","pos","gene")])


## Set plot parameters
pp <- getDefaultPlotParams(plot.type=2)
pp$ideogramheight <- 100
pp$topmargin <- 480


## Macaca mulatta // rheMac10
pdf("~/Research/postdoc/projects/PrimateSomatic/06_plots/prelim/RG25/macaca_targeted_nanoseq_karyoploter.pdf", width = 7, height = 8)
kp_Mmul_10<-plotKaryotype(genome = GRdataSIZE,  main="Macaca mulatta // Mmul_10", plot.type = 2, plot.params = pp, chromosomes = chromosomes)
#kpPlotRegions(kp_Mmul_10, GRdataPANEL, col ="firebrick2", data.panel = 1, r1 = 2)

kpPlotRainfall(kp_Mmul_10, GRdataMUT, col = "firebrick2", data.panel = 1)

kpPlotMarkers(kp_Mmul_10, data=GRdataWRITE, labels=GRdataWRITE$gene, cex=0.5, adjust.label.position = TRUE, text.orientation = "horizontal", y=1.5)
kpPoints(kp_Mmul_10, GRdataWRITE, y = 0, cex=1.2, pch=18 , bg="black", data.panel = 1)
#kpPoints(kp_Mmul_10, GRdataWRITE, y = 0, cex=1, pch=25 , bg="black", data.panel = 2)

legend("bottomright",
       legend = c("Candidate somatic variants", "Genes with candidate driver mutations"),
       pch = c(20, 18),
       col = c("firebrick2", "black"),
       pt.cex = c(1.5, 1.5),
       cex = 0.9,
       bty = "n",
       title = "Genomic Annotations")
dev.off()

#kpPlotDensity(kp_Mmul_10, GRdataGENE, col="lightgray", window.size = 1000000, data.panel = 1, alpha=0.5)
#kpPlotDensity(kp_Mmul_10, GRdataMUT, col = "firebrick2", data.panel = 2)

