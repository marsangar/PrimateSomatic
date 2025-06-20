# Comparative analyses of somatic mutational processes in primates across lifespans
This repository hosts data, scripts, and documentation related to the primate somatic mutation (PMS) project

## Project Overview
The vast diversity in lifespans among organisms provides a remarkable natural experiment in which to explore the
evolutionary innovations that have shaped the extensive variation of this phenotype. Non-human primates in
particular represent a critical taxon for understanding the evolution of human lifespan due to their close phylogenetic
proximity and broad range of lifespans. Here, we propose to characterize the somatic mutational landscape of aging
across multiple tissues in diverse primate species spanning million years of evolution.


![dataBPMR_combined](06_plots/prelim/dataBPMR_combined.png)


## Objectives

- To create a pan-primate somatic mutational atlas
- To characterize somatic clonal evolution across primates
- To identify genetic determinants of primate somatic evolution


## Primate Species

### Foundational species
- Common marmosets (_Callitrhix jacchus_) // calJac4 assembly
- Reshus macaques (_Macaca mulatta_) // rheMac10 assembly
- Olive baboons (_Papio anubis_) // papAnu4 assembly
- Humans (_Homo sapiens_) // hg38 assembly

### Additional species
- Chimpanzee (_Pan troglodytes_) // panTro6 assembly
- Western lowland gorilla (_Gorilla gorilla gorilla_) // gorGor6 assembly
- White-naped mangabey (_Cercocebus lunulatus_) // PGDP_CerLun assembly
- Squirrel monkeys (_Saimiri boliviensis boliviensis_) // saiBol1 assembly
- Pied tamarins (_Saguinus bicolor_) // PGDP_SagBic assembly
- Slender loris (_Loris lydekkeranius_) // PGDP_LorLyd assembly
- Ring-tailed lemur (_Lemur catta_) // mLemCat1.pri assembly

## Target tissue types
- Cardiac muscle (cardiomyocytes)
- Skeletal muscle (skeletal myofibers)
- Peripheral blood (bood cells)
- Liver (hepatocytes)
- Kidney (renal endothelial cells)
- Colon (epithelial cells and infiltrating lymphocytes)
- Skin (epithelial cells)
- Esophagus (epithelial cells)
- Bladder (epithelial cells)
- Cerebral cortex (pyramidal neurons)

## Requirements
- R (version 4.5.0)
- Bioconductor (version 3.21)

See renv.lock for comprehensive R package list

## Clone the Repository

```bash
git clone https://github.com/marsangar/PrimateSomatic.git
cd PrimateSomatic
```

## Contact
Martín Santamarina García

- Research Associate - Comparative Somatic Evolutionary Genomics / Department of Genetics / University of Cambridge
- Contingent Worker - Cancer predisposition and ageing / CASM / Wellcome Sanger Institute

ms3242@cam.ac.uk / ms84@sanger.ac.uk

## Acknowledgements
We thank the collaborators, funding agencies, and facilities that support this project:

NIH Funding Reference: 1R01AG087974 / 00011897
CUFS project code: PCAG/521 
 
- *Wellcome Sanger Institute*
- *University of Berkeley*
- *University of Cambridge*

- *Medical Research Council Centre for Macaques (CFM), UK*
- *Biomedical Primate Research Centre (BPRC), Netherlands*
- *Southwest National Primate Research Center (SNPRC), USA*
- *Zoological Society of London (ZSL), UK*
- *Coriell Institute for Medical Research (CIMR), USA*




