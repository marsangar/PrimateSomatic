#### Martín Santamarina García
#### 08/07/2025

### This script is used to keep track of relevant shell commands while working at FARM

### Farm connection
ssh ms84@farm22-head1

### Checking quota
lfs quota -h /lustre/scratch125

### LSF commands
bjobs
bjobs -p
bjobs -o estart_time
bjobs -o start_time



########################
### NanoSeq Cleaning ###
########################
#WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/nanoseq/grandparent
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/Set1
cd $WORKDIR

#### Save CRAM files with sample name
mkdir -p CRAM
cut -f1,2 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/add_rb/${ID}*.cram CRAM/${SAMPLE}.cram  ; done 
cut -f1,2 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/add_rb/${ID}*.cram.crai CRAM/${SAMPLE}.cram.crai  ; done 

#### Save NEAT CRAM files with sample name
mkdir -p NEAT_CRAM
cut -f1,2 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/dedup/${ID}*.neat.cram NEAT_CRAM/${SAMPLE}.neat.cram  ; done 
cut -f1,2 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/dedup/${ID}*.neat.cram.crai NEAT_CRAM/${SAMPLE}.neat.cram.crai  ; done 


### Clean work dir
find work/ -type f ! -name ‘*.command.out’ ! -name ‘*.command.err’ -delete
##################





### Copy files with rsync
rsync -avh --dry-run /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/
rsync -avh /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/

# dry-run
rsync -avh --dry-run /lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/ /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/
# actual copy
rsync -avh /lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/ /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/
sent 641.90G bytes  received 8.64K bytes  445.30M bytes/sec
total size is 641.75G  speedup is 1.00

### Major Filesystems:

### NFS
/nfs/users/nfs_m/ms84 # HOME
/nfs/casm/team294rr/ms84

### WAREHOUSE
/warehouse/casm/team294rr/ms84

### LUSTRE
/lustre/scratch125/casm/teams/team294/projects/cseg

/lustre/scratch125/casm/staging/team294/ms84/nanoseq

### SOFTWARE
/software/team294/ms84

### iRODS
module add ISG/IRODS/1.0
iinit

ils 

##################################
### Confluence iRODS commands ####
##################################
module load ISG/IRODS/1.0
iinit

# Search data by Study Name 
imeta qu -z seq -d study like '%Primate%' and target = 1 and manual_qc = 1


#  Downloading data from iRODS
cibackup get -h


### 
module add ISG/IRODS/1.0


Luego
iinit

Si tienes acceso, puedes ver lo que hay de Pacbio con:
ils /seq/pacbio/r84093_20250409_090424/1_D01

Si quieres descargarlo, con:
iget -r /seq/pacbio/r84093_20250409_090424/1_D01


#### Import data from IRODs
module load dataImportExport

#exportData.pl -n -st -l live -s 7061STDY13446126 -o /lustre/scratch125/casm/teams/team294/projects/Apert_Syndrome_Oxford_Sperm/Output/irods_data -study 7061 -t uc

exportData.pl -n -st -l live -s 6416STDY14759778 -o /lustre/scratch125/casm/staging/team294/ms84/PrimateSomatic -study 6416 -t uc
exportData.pl -n -st -l live -s 6416STDY14759779 -o /lustre/scratch125/casm/staging/team294/ms84/PrimateSomatic -study 6416 -t uc
exportData.pl -n -st -l live -s 6416STDY14759777 -o /lustre/scratch125/casm/staging/team294/ms84/PrimateSomatic -study 6416 -t uc



6416STDY14759778
14:02
6416STDY14759779
14:03
6416STDY14759777

### Change directory permissions
chmod -R 775 /lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes
chmod -R 775 /lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/Ovis_aries/Oar_rambouillet_v1.0/reference_files

chmod -R g+w /nfs/casm/team294rr/cseg

chmod 755 /nfs/casm/team294rr/cseg/projects/PrimateSomatic


### Useful links

#### Mount drives ####
- sshfs_config.zip
- cgp_mount.sh

### Mount Drives
/Users/rs30/cgp_mount.sh /Users/rs30/sshfs_config/network_home 
/Users/rs30/cgp_mount.sh /Users/rs30/sshfs_config/lustre

/Users/ms84/cgp_mount.sh /Users/ms84/sshfs_config/network_home 
/Users/ms84/cgp_mount.sh /Users/ms84/sshfs_config/lustre




### Hannah's Zoom recording (QC, signature analysis, etc.)
https://sanger.zoom.us/rec/share/PSIt2rgRbG3AKwfBdCDq3ggLQZawzPytnAAINrt23MmAdDYhM_-8M5OgXGvfUW6V.POhwQZaP2XWmGVVU



### Adrian Baez Cross-species NanoSeq


### Run Nanoseq on chimpanzee data
cd /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/nanoseq/chimpanzee
bsub <  /nfs/casm/team294rr/cseg/projects/PrimateSomatic/03_scripts/bash/runNanoseq_chimpanzee.sh

### Run MASK pipeline

### Step 1
cd /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/nanoseq/chimpanzee
/software/team294/ms84/nanoseq/SNPmask_pipeline_nextflow/1_CoverageHist.sh /nfs/casm/team294rr/cseg/projects/PrimateSomatic/03_scripts/configs/nanoseq/chimpanzee/SPECIES_INFO.txt

### Step 2
cd /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/nanoseq/chimpanzee
/software/team294/ms84/nanoseq/SNPmask_pipeline_nextflow/2_VariantCalling.sh /nfs/casm/team294rr/cseg/projects/PrimateSomatic/03_scripts/configs/nanoseq/chimpanzee/SPECIES_INFO.txt

### Step 3
cd /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/nanoseq/chimpanzee
/software/team294/ms84/nanoseq/SNPmask_pipeline_nextflow/3_bsub_command.sh /nfs/casm/team294rr/cseg/projects/PrimateSomatic/03_scripts/configs/nanoseq/chimpanzee/SPECIES_INFO.txt

### Step 4
cd /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/nanoseq/chimpanzee
/software/team294/ms84/nanoseq/SNPmask_pipeline_nextflow/4_MaskMerging.sh /nfs/casm/team294rr/cseg/projects/PrimateSomatic/03_scripts/configs/nanoseq/chimpanzee/SPECIES_INFO.txt



#### Save NEAT CRAM files with sample name
cut -f1,2 -d "," input/samples_Set1.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/dedup/${ID}*.neat.cram NEAT_CRAM/${SAMPLE}.neat.cram  ; done 

#### Save CRAM files with sample name
cut -f1,2 -d "," input/samples_Set1.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/dedup/${ID}*.cram CRAM/${SAMPLE}.cram  ; done 



### find work/ -type f ! -name ‘*.command.out’ ! -name ‘*.command.err’ -delete
 

### Save NEAT CRAM files with sample name (Pan Troglodytes Nanoseq)
cut -f1,2 -d "," /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/nanoseq/set1/input/samples_set1.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/dedup/${ID}*.neat.cram NEAT_CRAM/${SAMPLE}.neat.cram  ; done 
 
mv   work/*/*/dedup/6416STDY14759777.neat.cram NEAT_CRAM/CHPD0002b_ds0001.neat.cram


### Save CRAM files with sample name (Pan Troglodytes Nanoseq)
cut -f1,2 -d "," /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/nanoseq/set1/input/samples_set1.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v | while read -r SAMPLE ID ; do mv work/*/*/add_rb/${ID}*.cram CRAM/${SAMPLE}.cram  ; done

mv   work/*/*/add_rb/6416STDY14759777.cram NEAT_CRAM/CHPD0002b_ds0001.cram


###################################
###### Software Installation ######
###################################

### Install kraken2

## Move to software directory
cd /software/team294/cseg/

## clone kraken2 repository
git clone https://github.com/DerrickWood/kraken2.git

## Move to kraken2 directory
cd kraken2
## Build kraken2
KRAKEN2_DIR=/software/team294/cseg/kraken2
./install_kraken2.sh $KRAKEN2_DIR

### To make things easier for you, you may want to copy/symlink the following
#files into a directory in your PATH:
#  /software/team294/cseg/kraken2/kraken2
#  /software/team294/cseg/kraken2/kraken2-build
#  /software/team294/cseg/kraken2/kraken2-inspect

### Add /software/team294/cseg/bin to your PATH
cp $KRAKEN2_DIR/kraken2{,-build,-inspect} /software/team294/cseg/bin

