#### Martín Santamarina García
#### 28/03/2026

####  Transfer files from FARM to CSD3 

### Reference genomes
rsync -avh --progress \
/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Varanus_komodoensis/ASM479886v1/reference_files/* \
ms3242@login.hpc.cam.ac.uk:/home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/reference_genomes/Varanus_komodoensis/ASM479886v1/reference_files

rsync -avh --progress \
/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Aldabrachelys_gigantea/AldGig_1.0/reference_files/* \
ms3242@login.hpc.cam.ac.uk:/home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/reference_genomes/Aldabrachelys_gigantea/AldGig_1.0/reference_files

rsync -avh --progress \
/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Ophiophagus_hannah/Ophiophagus_hannah_FG1930_v3.1/reference_files/* \
ms3242@login.hpc.cam.ac.uk:/home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/reference_genomes/Ophiophagus_hannah/Ophiophagus_hannah_FG1930_v3.1/reference_files
