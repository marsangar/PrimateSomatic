#### Martín Santamarina García
#### 19/08/2025

#### Run Nanoseq Targeted on Squirrel (TEST)
cd /lustre/scratch125/casm/teams/team294/projects/cseg/SquirrelSentinel/Sciurus_carolinensis/nanoseq_targeted/Test
bsub < /nfs/casm/team294rr/cseg/projects/SquirrelSentinel/Sciurus_carolinensis/nanoseq_targeted/launch/runNanoseqTargeted_Test.sh 


### Sort, merged and set +/- 500 bp buffer 
module load bedtools2-2.29.0/
sort -k1,1 -k2,2n targeted_panel_mSciCar1.2.bed > targeted_panel_mSciCar1.2.sorted.bed
bedtools merge -i targeted_panel_mSciCar1.2.sorted.bed > targeted_panel_mSciCar1.2.sorted.merged.bed
cat targeted_panel_mSciCar1.2.sorted.merged.bed |awk 'OFS="\t" {print $1, $2-500, $3+500}' > targeted_panel_mSciCar1.2.sorted.merged_500.bed


### Prepare exclude panel
module load bedtools2-2.29.0/
BED_a=/lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/Sciurus_carolinensis/mSciCar1.2/reference_files/mSciCar1.2_sequence_report.bed
BED_b=/lustre/scratch125/casm/teams/team294/projects/cseg/SquirrelSentinel/Sciurus_carolinensis/nanoseq_targeted/panel/targeted_panel_mSciCar1.2.sorted.merged_500.bed
OUTPUT=/lustre/scratch125/casm/teams/team294/projects/cseg/SquirrelSentinel/Sciurus_carolinensis/nanoseq_targeted/panel/exclude_panel_mSciCar1.2.sorted.merged_500.bed
bedtools subtract -a $BED_a -b $BED_b > $OUTPUT
bgzip $OUTPUT
tabix -p bed $OUTPUT.gz

#### Run Nanoseq Targeted on Squirrel (TEST)
cd /lustre/scratch125/casm/teams/team294/projects/cseg/SquirrelSentinel/Sciurus_carolinensis/nanoseq_targeted/Test
bsub < /lustre/scratch125/casm/teams/team294/projects/cseg/SquirrelSentinel/Sciurus_carolinensis/nanoseq_targeted/Test/runNano.bsub
Job <487287> is submitted to queue <week>.

#### Run Nanoseq Targeted on Squirrel (SET1)
module add ISG/IRODS/1.0
iinit
cd /lustre/scratch125/casm/teams/team294/projects/cseg/SquirrelSentinel/Sciurus_carolinensis/nanoseq_targeted/Set1
bsub < /lustre/scratch125/casm/teams/team294/projects/cseg/SquirrelSentinel/Sciurus_carolinensis/nanoseq_targeted/Set1/runNano.bsub



#### Stage Squirrel data

#### Import data from IRODs
module load dataImportExport

exportData.pl -n -st -l live -s 7949STDY15350561 -o /lustre/scratch126/casm/staging/team294_cseg/users/ms84/Sciurus_carolinensis/export -study 7949 -t uc
exportData.pl -n -st -l live -s 7949STDY15350573 -o /lustre/scratch126/casm/staging/team294_cseg/users/ms84/Sciurus_carolinensis/export -study 7949 -t uc
exportData.pl -n -st -l live -s 7949STDY15350583 -o /lustre/scratch126/casm/staging/team294_cseg/users/ms84/Sciurus_carolinensis/export -study 7949 -t uc

############################
### DupCaller Conversion ###
############################

### Create symbolic links to the CRAM files in the working directory
ln -s /lustre/scratch126/casm/staging/team294_cseg/users/ms84/Sciurus_carolinensis/export/7949/7949STDY15350561/raw_lane_bams/7949STDY15350561_50131_7_plx2.ua.cram 7949STDY15350561_50131_7_plx2.ua.cram
ln -s /lustre/scratch126/casm/staging/team294_cseg/users/ms84/Sciurus_carolinensis/export/7949/7949STDY15350573/raw_lane_bams/7949STDY15350573_50131_8_plx3.ua.cram 7949STDY15350573_50131_8_plx3.ua.cram
ln -s /lustre/scratch126/casm/staging/team294_cseg/users/ms84/Sciurus_carolinensis/export/7949/7949STDY15350583/raw_lane_bams/7949STDY15350583_50131_8_plx6.ua.cram 7949STDY15350583_50131_8_plx6.ua.cram


### Apply irods2dupcaller conversion
module load samtools-1.19.2
IRODS2DUPCALLER=/software/team294/cseg/DupCaller_preprocessing/semi/irods2dupcaller.sh

WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/SquirrelSentinel/Sciurus_carolinensis/dupcaller/CRAM
cd $WORKDIR
mkdir -p logs
for SAMPLE in 7949STDY15350561_50131_7_plx2.ua 7949STDY15350573_50131_8_plx3.ua 7949STDY15350583_50131_8_plx6.ua
do
 bsub5000 -e logs/$SAMPLE.err -o logs/$SAMPLE.err "bash $IRODS2DUPCALLER $SAMPLE.cram"
done


### Transfer from FARM to CSD3
rsync -avh --progress \
dupcaller-7949STDY15350583_50131_8_plx6.ua.bam* \
dupcaller-7949STDY15350573_50131_8_plx3.ua.bam* \
dupcaller-7949STDY15350561_50131_7_plx2.ua.bam* \
ms3242@login.hpc.cam.ac.uk:/rds/user/ms3242/hpc-work/SquirrelSentinel/dupcaller

md5sum dupcaller-*.bam* > dupcaller.md5




### 26/03/2026

### Run nanoseq on pilot squirrel data (same matched-normal) just to get proper remapping
module load ISG/IRODS/1.0
iinit
cd /lustre/scratch126/casm/teams/team294/projects/cseg/SquirrelSentinel/Sciurus_carolinensis/nanoseq/pilot
bsub < runNext.sh

#### Save CRAM files with sample name
cut -f1,2 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/dedup/${ID}*.cram CRAM/${SAMPLE}.cram  ; done 
cut -f1,2 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/dedup/${ID}*.cram.crai CRAM/${SAMPLE}.cram.crai  ; done 

### Apply nanoseq2dupcaller conversion (FULL)
module load samtools-1.19.2
NANOSEQ2DUPCALLER=/software/team294/cseg/DupCaller_preprocessing/full/nanoseq2dupcaller.sh

WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/SquirrelSentinel/Sciurus_carolinensis/nanoseq/pilot/CRAM
cd $WORKDIR
mkdir -p logs
for SAMPLE in EGSD0002f_lo0001 EGSD0013e_lo0001 EGSD0013d_lo0005
do
bsub80000 -e logs/$SAMPLE.err -o logs/$SAMPLE.err "bash $NANOSEQ2DUPCALLER $SAMPLE.cram"
done

### Get md5sum of the converted files
md5sum dupcaller-*.bam* > dupcaller.md5


### Transfer from FARM to CSD3

rsync -avh --progress \
dupcaller-EGSD0002f_lo0001.bam* \
dupcaller-EGSD0013e_lo0001.bam* \
dupcaller-EGSD0013d_lo0005.bam* \
ms3242@login.hpc.cam.ac.uk:/home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/projects/SquirrelProject/Sciurus_carolinensis/dupcaller/2_alignment