#######################################
### 02/10/2025 - Martín Santamarina ###
#######################################

######################################
#### Dupcaller on Pan_troglodytes ####
######################################

### install DupCaller (DONE)
# git clone --branch dev https://github.com/YuheCheng62/DupCaller.git
# cd DupCaller
# pip install .

### Activate the environment
module load conda
conda activate /nfs/users/nfs_m/ms84/.conda/envs/Dupcaller // 10/10/2025
module load samtools-1.19.2

### Create a reference index for use with DupCaller (DONE)
# Index the reference genome.
# REFDIR=/lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/Pan_troglodytes/NHGRI_mPanTro3-v2.1_pri/reference_files
# cd $REFDIR
# bsub80000 -e kk_dupcaller_index.err -o kk_dupcaller_index.out DupCaller.py index -f genome.fa

### Set Working directory
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/dupcaller/runs/latest


### Apply nanoseq2dupcaller conversion script to the CRAM files (semi-processed > irods2dupcaller.sh, full-processed > nanoseq2dupcaller.sh)
for SAMPLE in CHPD0002b_ds0001 CHPD0002b_ds0002 CHPD0002f_ds0001
do
CRAM="$SAMPLE".cram
bsub80000 -q long -e nanoseq2dupcaller_$SAMPLE.err -o nanoseq2dupcaller_$SAMPLE.err bash nanoseq2dupcaller.sh  $CRAM
done
Job <836905> is submitted to queue <long>.
Job <836906> is submitted to queue <long>.
Job <836907> is submitted to queue <long>.


### Dupcaller Variant Calling // rescue
sample="CHPD0002f_ds0001"
bsub80000 -q long -n 16 -R "span[ptile=16]" -e rescue_$sample.dupcallerV2.err -o rescue_$sample.dupcallerV2.err DupCaller.py call \
                   -b /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/dupcaller/set1/dupcaller-CHPD0002f_ds0001.bam \
                   -f /lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/Pan_troglodytes/NHGRI_mPanTro3-v2.1_pri/reference_files/genome.fa \
                   -r NC_072398.2 NC_086015.1 NC_072401.2 NC_072402.2 NC_072403.2 NC_072404.2 NC_072405.2 NC_072406.2 NC_072407.2 NC_072408.2 NC_072409.2 NC_072410.2 NC_072411.2 NC_072412.2 NC_072413.2 NC_072414.2 NC_072415.2 NC_072416.2 NC_072417.2 NC_072418.2 NC_072419.2 NC_072420.2 NC_086016.1 NC_072421.2 NC_072422.2    \
                   -p 16 \
                   -n /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/dupcaller/set1/dupcaller-CHPD0002b_ds0001.bam \
                   -m /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/nanoseq/chimpanzee/output/MASKS/SNP+NOISE.NV2.SeenNV1.chimpanzee.bed.gz \
                   -o rescue_DupCallerV2.$sample \
                   -ax 50 \
                   -tt 8 \
                   -tr 8 \
                   -d 15 \
                   -nm 4 \
                   --rescue TRUE


sample="CHPD0002b_ds0002"
bsub80000 -q long -n 16 -R "span[ptile=16]" -e rescue_$sample.dupcallerV2.err -o rescue_$sample.dupcallerV2.err DupCaller.py call \
                   -b /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/dupcaller/set1/dupcaller-CHPD0002b_ds0002.bam \
                   -f /lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/Pan_troglodytes/NHGRI_mPanTro3-v2.1_pri/reference_files/genome.fa \
                   -r NC_072398.2 NC_086015.1 NC_072401.2 NC_072402.2 NC_072403.2 NC_072404.2 NC_072405.2 NC_072406.2 NC_072407.2 NC_072408.2 NC_072409.2 NC_072410.2 NC_072411.2 NC_072412.2 NC_072413.2 NC_072414.2 NC_072415.2 NC_072416.2 NC_072417.2 NC_072418.2 NC_072419.2 NC_072420.2 NC_086016.1 NC_072421.2 NC_072422.2    \
                   -p 16 \
                   -n /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/dupcaller/set1/dupcaller-CHPD0002b_ds0001.bam \
                   -m /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/nanoseq/chimpanzee/output/MASKS/SNP+NOISE.NV2.SeenNV1.chimpanzee.bed.gz \
                   -o rescue_DupCallerV2.$sample \
                   -ax 50 \
                   -tt 8 \
                   -tr 8 \
                   -d 15 \
                   -nm 4 \
                   --rescue TRUE


Job <872414> is submitted to queue <long>.
Job <872415> is submitted to queue <long>.


### Compare NanoSeq & DupCaller 
module load bcftools-1.9

SAMPLE=CHPD0002f_ds0001
DUPCALLER_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/dupcaller/set1/DupCallerV2."$SAMPLE"/DupCallerV2."$SAMPLE"_snv.vcf.gz
NANOSEQ_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/nanoseq/set1/output/outNextflow/NanoSeq/"$SAMPLE"/"$SAMPLE".vcf.gz
bsub20000 -e bcftools_isec.err -o bcftools_isec.err bcftools isec $DUPCALLER_VCF $NANOSEQ_VCF -p isec_"$SAMPLE"_DC_NS


SAMPLE=CHPD0002b_ds0002
DUPCALLER_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/dupcaller/set1/DupCallerV2."$SAMPLE"/DupCallerV2."$SAMPLE"_snv.vcf.gz
NANOSEQ_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/nanoseq/set1/output/outNextflow/NanoSeq/"$SAMPLE"/"$SAMPLE".vcf.gz
bsub20000 -e bcftools_isec.err -o bcftools_isec.err bcftools isec $DUPCALLER_VCF $NANOSEQ_VCF -p isec_"$SAMPLE"_DC_NS


for SAMPLE in CHPD0002f_ds0001 CHPD0002b_ds0002
do
echo "Evaluating "$SAMPLE" executions:"
echo "Dupcaller SNVs"
zcat /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/dupcaller/set1/DupCallerV2."$SAMPLE"/DupCallerV2."$SAMPLE"_snv.vcf.gz |grep -v "^#" |wc -l
echo "Dupcaller SNVs shared with NanoSeq"
cat /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/benchmarking/isec_"$SAMPLE"_DC_NS/0002.vcf |grep -v "^#" |wc -l
echo
done




####  Dupcaller reference genome index for Mus_musculus
REFDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Mus_musculus/GRCm39/reference_files
cd $REFDIR
bsub80000 -e kk_dupcaller_index.err -o kk_dupcaller_index.out DupCaller.py index -f genome.fa


#### Dupcaller genome index for Drosophila_melanogaster
REFDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Drosophila_melanogaster/Release_6_plus_ISO1_MT/reference_files
cd $REFDIR
bsub80000 -e kk_dupcaller_index.err -o kk_dupcaller_index.out DupCaller.py index -f genome.fa


#### Dupcaller genome index for Ovis aries
REFDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Ovis_aries/Oar_rambouillet_v1.0/reference_files
cd $REFDIR
bsub80000 -e kk_dupcaller_index.err -o kk_dupcaller_index.out DupCaller.py index -f genome.fa


