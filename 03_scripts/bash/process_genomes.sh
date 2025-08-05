#!/bin/bash


# DOWNLOADING AND PREPROCESSING REFERENCE GENOMES FOR CROSS-SPECIES II


# Function to download and process a reference genome;
# needs variables $SPECIES, $ASSEMBLY, $URL




process_genome() {

    echo -e "\n\n-----------------------------"
    echo $SPECIES
    echo -e "-----------------------------\n"
    mkdir -p /lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/$SPECIES/$ASSEMBLY/reference_files
    cd /lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/$SPECIES/$ASSEMBLY/reference_files

    wget $URL -O genome.fa.gz
    gunzip -f genome.fa.gz

    module load pcap-core
    module load samtools-1.19
    bsub -M 100000 -R "select[mem>100000] rusage[mem=100000]" -o $PWD/log.%J "echo bwa-mem2 index; bwa-mem2 index genome.fa; echo; echo samtools dict; samtools dict -a $ASSEMBLY -s $SPECIES genome.fa > genome.fa.dict; echo; echo ref_freqs.R; Rscript /software/team294/ms84/nanoseq/Set_reference_genomes/ref_freqs.R genome.fa"

}





# # rhesus_macaque
# #------------------------------
#SPECIES=macaca_mulatta
#ASSEMBLY=Mmul_10
#URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/339/765/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.fna.gz 
#process_genome
# #------------------------------

# # Komodo dragon
# #------------------------------
#SPECIES=varanus_komodoensis
#ASSEMBLY=ASM479886v1
#URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/798/865/GCF_004798865.1_ASM479886v1/GCF_004798865.1_ASM479886v1_genomic.fna.gz
#process_genome
# #------------------------------

# Indian jumping ant
#------------------------------
#SPECIES=harpegnathos_saltator
#ASSEMBLY=Hsal_v8.6
#URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/227/715/GCF_003227715.2_Hsal_v8.6/GCF_003227715.2_Hsal_v8.6_genomic.fna.gz
#process_genome
# #------------------------------

# # Domestic dog
# #------------------------------
#SPECIES=canis_lupus_familiaris
#ASSEMBLY=UU_Cfam_GSD_1.0
#URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/100/685/GCF_011100685.1_UU_Cfam_GSD_1.0/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna.gz
#process_genome
# #------------------------------

# Gray squirrel
#------------------------------
#SPECIES=sciurus_carolinensis
#ASSEMBLY=mSciCar1.2
#URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/686/445/GCF_902686445.1_mSciCar1.2/GCF_902686445.1_mSciCar1.2_genomic.fna.gz
#process_genome
#------------------------------

# # Fruit fly
# #------------------------------
#SPECIES=drosophila_melanogaster
#ASSEMBLY=Release_6_plus_ISO1_MT
#URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
#process_genome
# #------------------------------

# # Sperm whale
# #------------------------------
#SPECIES=physeter_macrocephalus
#ASSEMBLY=ASM283717v5
#URL=ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/837/175/GCF_002837175.3_ASM283717v5/GCF_002837175.3_ASM283717v5_genomic.fna.gz
#process_genome
# #------------------------------

# # Desert locust
# #------------------------------
#SPECIES=schistocerca_gregaria
#ASSEMBLY=iqSchGreg1.2
#URL=ftp.ncbi.nlm.nih.gov/genomes/all/GCF/023/897/955/GCF_023897955.1_iqSchGreg1.2/GCF_023897955.1_iqSchGreg1.2_genomic.fna.gz
#process_genome
# #------------------------------

# # HUMBOLDT PENGUIN
# #------------------------------
#SPECIES=spheniscus_humboldti
#ASSEMBLY=bSphHub1.pri
#URL=ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/474/245/GCA_027474245.1_bSphHub1.pri/GCA_027474245.1_bSphHub1.pri_genomic.fna.gz
#process_genome
# #------------------------------

# # HYACINTH MACAW
# #------------------------------
#SPECIES=anodorhynchus_hyacinthinus
#ASSEMBLY=ASM993644v2
#URL=ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/936/445/GCA_009936445.2_ASM993644v2/GCA_009936445.2_ASM993644v2_genomic.fna.gz
#process_genome
# #------------------------------

# # FERRET
# #------------------------------
#SPECIES=mustela_putorius_furo
#ASSEMBLY=ASM1176430v1.1
#URL=ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/764/305/GCF_011764305.1_ASM1176430v1.1/GCF_011764305.1_ASM1176430v1.1_genomic.fna.gz
#process_genome
# #------------------------------

# # KING COBRA
# #------------------------------
#SPECIES=ophiophagus_hannah
#ASSEMBLY=OphHan1.0
#URL=ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/516/915/GCA_000516915.1_OphHan1.0/GCA_000516915.1_OphHan1.0_genomic.fna.gz
#process_genome
# #------------------------------

# # CHIMPANZEE
# #------------------------------
#SPECIES=pan_troglodytes
#ASSEMBLY=NHGRI_mPanTro3-v2.1_pri
#URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz
#process_genome
# #------------------------------

# # JAPANESE QUAIL
# #------------------------------
#SPECIES=coturnix_japonica
#ASSEMBLY=Coturnix_japonica_2.1
#URL=ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/577/835/GCF_001577835.2_Coturnix_japonica_2.1/GCF_001577835.2_Coturnix_japonica_2.1_genomic.fna.gz
#process_genome
# #------------------------------

# # OSTRICH
# #------------------------------
#SPECIES=struthio_camelus_australis
#ASSEMBLY=ASM69896v1
#URL=ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/698/965/GCF_000698965.1_ASM69896v1/GCF_000698965.1_ASM69896v1_genomic.fna.gz
#process_genome
# #------------------------------

# # HYDRACTINIA SYMBIOLONGICARPUS
# #------------------------------
#SPECIES=hydractinia_symbiolongicarpus
#ASSEMBLY=HSymV2.1
#URL=ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/227/915/GCF_029227915.1_HSymV2.1/GCF_029227915.1_HSymV2.1_genomic.fna.gz
#process_genome
# #------------------------------


# # Common marmoset
# #------------------------------
#SPECIES=callitrhix_jacchus
#ASSEMBLY=mCalJa1.2.pat.X
#URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/100/555/GCF_011100555.1_mCalJa1.2.pat.X/GCF_011100555.1_mCalJa1.2.pat.X_genomic.fna.gz
#process_genome
# #------------------------------

# # Olive baboon
# #------------------------------
#SPECIES=papio_anubis
#ASSEMBLY=Panubis1.1
#URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/728/515/GCA_008728515.2_Panubis1.1/GCA_008728515.2_Panubis1.1_genomic.fna.gz
#process_genome
# #------------------------------


# # English oak
# #------------------------------
SPECIES=quercus_robur
ASSEMBLY=dhQueRobu3.1
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/932/294/415/GCF_932294415.1_dhQueRobu3.1/GCF_932294415.1_dhQueRobu3.1_genomic.gff.gz
process_genome
# #------------------------------


# # Thale cress
# #------------------------------
SPECIES=arabidopsis_thaliana
ASSEMBLY=TAIR10.1
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz
process_genome
# #------------------------------
