#### Martín Santamarina García
#### 


########################################
#### Run Nanoseq on Pan troglodytes ####
########################################

module load ISG/IRODS/1.0
iinit
cd /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/nanoseq/set1
bsub <  runNext.sh



#### Chimpanzee Nanoseq variant calling results:
/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/nanoseq/chimpanzee/output/outNextflow/NanoSeq/CHPD0002b_ds0002/CHPD0002b_ds0002.vcf.gz # liver
/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/nanoseq/chimpanzee/output/outNextflow/NanoSeq/CHPD0002f_ds0001/CHPD0002f_ds0001.vcf.gz # intestine

###################################
#### Nanoseq Targeted Analysis ####
###################################

### Turn panel
cd /Users/ms84/Research/postdoc/projects/PrimateSomatic/02_data/databases/homo_sapiens
cat Sanger_TERT-v4_TE-95148282_hg19_highstringencyfilter_buccal_gene_list.tsv |awk 'OFS="\t" {print "chr"$2,$3,$4}' |tail -n +2 > Sanger_TERT-v4_TE-95148282_hg19_highstringencyfilter_buccal_gene_list.bed

### UCSC liftover & change chromosome names

### Sort, merged and set +/- 500 bp buffer 
module load bedtools2-2.29.0/
sort -k1,1 -k2,2n targeted_panel_Mmul_10.bed > targeted_panel_Mmul_10.sorted.bed
bedtools merge -i targeted_panel_Mmul_10.sorted.bed > targeted_panel_Mmul_10.sorted.merged.bed
cat targeted_panel_Mmul_10.sorted.merged.bed |awk 'OFS="\t" {print $1, $2-500, $3+500}' > targeted_panel_Mmul_10.sorted.merged_500.bed

### Prepare exclude panel
module load bedtools2-2.29.0/
BED_a=/lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/Mmul_10_sequence_report.bed
BED_b=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/panel/targeted_panel_Mmul_10.sorted.merged_500.bed
OUTPUT=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/panel/exclude_panel_Mmul_10.sorted.merged_500.bed
bedtools subtract -a $BED_a -b $BED_b > $OUTPUT
bgzip $OUTPUT
tabix -p bed $OUTPUT.gz


#### Run Nanoseq Targeted on Reshus (Set1)
module add ISG/IRODS/1.0
iinit
cd /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1
bsub < /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1/runNano.bsub
Job <698982> is submitted to queue <week>.


#### Run Nanoseq Targeted on Reshus (Set1_recover)
module add ISG/IRODS/1.0
iinit
cd /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1_recover
bsub < /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1_recover/runNano.bsub
Job <700145> is submitted to queue <week>.

#### Run Nanoseq Targeted on Reshus (Set1_failed)
module add ISG/IRODS/1.0
iinit
cd /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1_failed
bsub < /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1_failed/runNano.bsub
Job <700145> is submitted to queue <week>.

#### Run Nanoseq Targeted on Reshus (Set1_recover_1)
module add ISG/IRODS/1.0
iinit
cd /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1_recover_1
bsub < /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1_recover_1/runNano.bsub


### Download Nanoseq Targeted output
WORKDIR=/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/Set1
cd $WORKDIR
#for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0004d_tds0001 MQD0004i_tds0001 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001
#for SAMPLE in MQD0003g_tds0001 MQD0003j_tds0001
for SAMPLE in MQD0002h_tds0001  MQD0002l_tds0001  MQD0003d_tds0001  MQD0003e_tds0001 MQD0004e_tds0002  MQD0004k_tds0002  MQD0006l_tds0001
do
mkdir -p $WORKDIR/"$SAMPLE"

### VCF files
scp ms84@farm22-head1:/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1_recover_1/outNextflow/NanoSeq/$SAMPLE/$SAMPLE.vcf.gz ./$SAMPLE

### Trinuc profiles
scp ms84@farm22-head1:/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1_recover_1/outNextflow/NanoSeq/$SAMPLE/post/$SAMPLE.trinuc-profiles.pdf ./$SAMPLE
done




WORKDIR=/Users/ms84/Research/postdoc/projects/PrimateSomatic/05_results/prelim/macaca_mulatta/nanoseq_targeted/Set1
cd $WORKDIR

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0004d_tds0001 MQD0004i_tds0001 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001
do
mkdir -p $WORKDIR/"$SAMPLE"

### VCF files
scp ms84@farm22-head1:/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1/outNextflow/NanoSeq/$SAMPLE/$SAMPLE.vcf.gz ./$SAMPLE

### Trinuc profiles
scp ms84@farm22-head1:/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1/outNextflow/NanoSeq/$SAMPLE/post/$SAMPLE.trinuc-profiles.pdf ./$SAMPLE
done


for SAMPLE in MQD0003g_tds0001 MQD0003j_tds0001
do
mkdir -p $WORKDIR/"$SAMPLE"

### VCF files
scp ms84@farm22-head1:/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1_recover/outNextflow/NanoSeq/$SAMPLE/$SAMPLE.vcf.gz ./$SAMPLE

### Trinuc profiles
scp ms84@farm22-head1:/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1_recover/outNextflow/NanoSeq/$SAMPLE/post/$SAMPLE.trinuc-profiles.pdf ./$SAMPLE
done


for SAMPLE in MQD0002h_tds0001  MQD0002l_tds0001  MQD0003d_tds0001  MQD0003e_tds0001 MQD0004e_tds0002  MQD0004k_tds0002  MQD0006l_tds0001
do
mkdir -p $WORKDIR/"$SAMPLE"

### VCF files
scp ms84@farm22-head1:/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1_recover_1/outNextflow/NanoSeq/$SAMPLE/$SAMPLE.vcf.gz ./$SAMPLE

### Trinuc profiles
scp ms84@farm22-head1:/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1_recover_1/outNextflow/NanoSeq/$SAMPLE/post/$SAMPLE.trinuc-profiles.pdf ./$SAMPLE
done





#### Filter Nanoseq variant calls with MASK
module load bcftools-1.9

### Create a CHR_MAP
cd /lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/Macaca_mulatta/Mmul_10/reference_files
cat Mmul_10_sequence_report.tsv |awk 'OFS="\t" {print $4, $9}' |head -n 23|  tail -n +2 > chr_map.txt

### Rename chromosomes in the MASK file
VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/masks/Macaca_mulatta/Mmul_10/mGAP/mGap.Rhesus_macaque.v3.0.vcf.gz
CHR_MAP=/lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/chr_map.txt
bsub20000 -e bcftools_annotate.err -o bcftools_annotate.err bcftools annotate $VCF --rename-chrs $CHR_MAP -O z -o mGap.Rhesus_macaque.v3.0_renamed.vcf.gz

### Index MASK file
MASK=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/masks/Macaca_mulatta/Mmul_10/mGAP/mGap.Rhesus_macaque.v3.0_renamed.vcf.gz
bsub20000 -e bcftools_index.err -o bcftools_index.err bcftools index $MASK


module load bcftools-1.9
cd /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/export


### Count variants on each VCF
#for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001  
#do
#VCF=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/export/$SAMPLE.vcf.gz
#echo $SAMPLE
#zcat $VCF | grep -v "^#" | wc -l
#zcat $VCF | grep -v "^#" | grep "PASS" | wc -l
#zcat $VCF | grep -v "^#" | grep "PASS" | grep "TYPE=snv" | wc -l
#zcat $VCF | grep -v "^#" | grep "PASS" | grep "TYPE=dnv" | wc -l
#zcat $VCF | grep -v "^#" | grep "PASS" | grep "TYPE=mnv" | wc -l
#zcat $VCF | grep -v "^#" | grep "PASS" | grep "TYPE=ins" | wc -l
#zcat $VCF | grep -v "^#" | grep "PASS" | grep "TYPE=del" | wc -l
#done

### Count variants on each filtered VCF
#for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001  
#do
#FILTERED_VCF=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/export/$SAMPLE.filtered.vcf.gz
#zcat $FILTERED_VCF | grep -v "^#" | wc -l
#done

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001  
do
VCF=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/export/$SAMPLE.vcf.gz
echo $SAMPLE
zcat $VCF | grep -v "^#" | wc -l
done

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001  
do
VCF=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/export/$SAMPLE.vcf.gz
MASK=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/masks/Macaca_mulatta/Mmul_10/mGAP/mGap.Rhesus_macaque.v3.0_renamed.vcf.gz
OUTPUT=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/export/$SAMPLE.filtered.vcf.gz
bsub20000 -q small -e bcftools_view_$SAMPLE.err -o bcftools_view_$SAMPLE.err bcftools view $VCF -T ^$MASK -O z -o $OUTPUT
done

#bsub20000 -e bcftools_isec_$SAMPLE.err -o bcftools_isec_$SAMPLE.err bcftools isec $VCF $MASK -o $OUTPUT -O z
#bsub20000 -e bcftools_isec_$SAMPLE.err -o bcftools_isec_$SAMPLE.err bcftools isec $VCF $MASK -p dir


for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001  
do
mv $SAMPLE.vcf_filtered.gz $SAMPLE_filtered.vcf.gz
done


#### Perform QC analysis on targeted nanoseq data
WORKDIR=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/QC
cd $WORKDIR

module load bedtools2-2.29.0/

for SAMPLE in MQD0001d_tds0001
do
BED=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/panel/targeted_panel_Mmul_10.sorted.merged_500.bed.fixed
BAM=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/Set1/NEAT_CRAM/"$SAMPLE".neat.cram
GENOME=/lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/Mmul_10_chrom_sizes.bed
OUTPUT=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/QC/"$SAMPLE"_panel_coverage.txt
bsub20000 -e bedtools_"$SAMPLE".err -o bedtools_"$SAMPLE".err bedtools coverage -a $BED -b $BAM -mean -g $GENOME -sorted
done



module load bedtools2-2.29.0/
bedtools sort -i /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/panel/targeted_panel_Mmul_10.sorted.merged_500.bed \
-g /lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/Mmul_10_chrom_sizes.bed \
> /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/panel/targeted_panel_Mmul_10.sorted.merged_500.bed.fixed




#########################
##### export module #####
#########################

### Coverage Evaluations
cd /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/export/post
cp ../../Set1/outNextflow/NanoSeq/*/post/*.cov.bed.gz* .
cp ../../Set1_recover/outNextflow/NanoSeq/*/post/*.cov.bed.gz* .
cp ../../Set1_recover_1/outNextflow/NanoSeq/*/post/*.cov.bed.gz* .


WORKDIR=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/export/post
cd $WORKDIR
DX_COVERAGE=/software/team294/ms84/scripts/perl/dx_coverage_per_captured_region.pl
BED=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/panel/targeted_panel_Mmul_10_withgenes.bed
#COVERAGE=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/export/post/MQD0001d_tds0001.cov.bed.gz
COVERAGE=*.cov.bed.gz
module load labProjectAdmin 
module load bedtools2-2.29.0/
perl $DX_COVERAGE $BED $COVERAGE


#DX_DROPOUTS=/software/team294/ms84/scripts/r/dx_comparisons2detect_dropouts.R
#Rscript $DX_DROPOUTS
 
scp ms84@farm22-head2:/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/export/post/coverage_info.persample.tsv .




################################
#### NSDROP QUALITY CONTROL ####
################################

WORKDIR=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/QC/nsdrop
cd $WORKDIR

### Create manifest for nsdrop
for SAMPLE in MQD0001j MQD0001d MQD0006f MQD0006e MQD0004i MQD0004d MQD0002f MQD0002e MQD0006h MQD0003j MQD0003g MQD0002h MQD0004k MQD0004e MQD0003e MQD0003d MQD0002l MQD0006l
do
echo "$SAMPLE"  /lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/export/post/"$SAMPLE".cov.bed.gz
done > PrimateSomatic_manifest_nsdrop.tsv

### run nsdrop on rhesus macaque targeted data
module load nsdrop/0.3.0 
BED=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/panel/targeted_panel_Mmul_10.sorted.merged_500.bed
MANIFEST=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/QC/nsdrop/PrimateSomatic_manifest_nsdrop.tsv
nsdrop -o out/ -b $BED $MANIFEST



#######################################################
#### Target NanoSeq Cross Sample Germline Filtering ###
#######################################################

### Get the CRAMS
CRAMDIR=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/CRAM
cd $CRAMDIR
cut -f1,2 -d "," ../Set1_recover/input/samples_Set1_recover.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v | while read -r SAMPLE ID ; do mv ../Set1_recover/work/*/*/add_rb/${ID}*.cram ${SAMPLE}.cram  ; done
cut -f1,2 -d "," ../Set1_recover/input/samples_Set1_recover.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v | while read -r SAMPLE ID ; do mv ../Set1_recover/work/*/*/add_rb/${ID}*.cram.crai ${SAMPLE}.cram.crai  ; done

cut -f1,2 -d "," ../Set1_recover_1/input/samples_Set1_recover_1.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v | while read -r SAMPLE ID ; do mv ../Set1_recover_1/work/*/*/add_rb/${ID}*.cram ${SAMPLE}.cram  ; done
cut -f1,2 -d "," ../Set1_recover_1/input/samples_Set1_recover_1.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v | while read -r SAMPLE ID ; do mv ../Set1_recover_1/work/*/*/add_rb/${ID}*.cram.crai ${SAMPLE}.cram.crai  ; done


### Load environment
module load conda
conda activate /nfs/users/nfs_m/ms84/.conda/envs/Dupcaller

WORKDIR=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/GERMLINE
cd $WORKDIR

### Genotype positions across samples
BAM=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/CRAM/MQD0002h_tds0001.cram
BED=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/GERMLINE/dataMUT.bed
OUTPUT=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/GERMLINE/test_pileup_genotype.txt
bsub20000 -e pileup.err -o pileup.out "python /software/team294/ms84/genomic_toolkit/scripts/pileup_genotype.py $BAM $BED > $OUTPUT"



SAMPLE=MQD0001d_tds0001
for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
BAM=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/CRAM/"$SAMPLE".cram
BED=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/panel/targeted_panel_Mmul_10_withgenes.bed
OUTPUT=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/GERMLINE/pileup_genotype_"$SAMPLE".txt
bsub20000 -e pileup_"$SAMPLE".err -o pileup_"$SAMPLE".out "python /software/team294/ms84/genomic_toolkit/scripts/pileup_genotype.py $BAM $BED > $OUTPUT"
done

### Combine into a single file
OUTPUT="pileup_genotype_COMBINED.txt"
> "$OUTPUT" 
for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
FILE=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/GERMLINE/pileup_genotype_"$SAMPLE".txt
awk -v sample="$SAMPLE" 'BEGIN { OFS="\t" } { print $0, sample }' "$FILE" >> "$OUTPUT"
done

### Combine into a single file
OUTPUT="pileup_genotype_COMBINED.txt"
> "$OUTPUT"
# --- Set header ---
FIRST_FILE="/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/GERMLINE/pileup_genotype_MQD0001d_tds0001.txt"
head -n 1 "$FIRST_FILE" | awk 'BEGIN{OFS="\t"}{print $0, "SAMPLE"}' >> "$OUTPUT"

# --- Process every sample
for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
    FILE=/lustre/scratch125/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/nanoseq_targeted/GERMLINE/pileup_genotype_"$SAMPLE".txt

    # Añadir contenido (saltando la primera línea del header) + columna SAMPLE
    tail -n +2 "$FILE" | awk -v sample="$SAMPLE" 'BEGIN{OFS="\t"}{print $0, sample}' >> "$OUTPUT"
done





#### Inspect Variant Calling Filters
for SAMPLE in PD*; do echo $SAMPLE; zcat "$SAMPLE"/"$SAMPLE".vcf.gz|grep -v "^#" | awk '{print $7}'|sort |uniq -c;echo ""; done

### Summarize percentages
zcat PD*/PD*.vcf.gz | grep -v "^#" | awk '
{
    count[$7]++
    total++
}
END {
    for (f in count) {
        printf "%s\t%d\t%.2f%%\n", f, count[f], 100*count[f]/total
    }
}'



#################################################
#### Generate MASK file for targeted nanoseq ####
#################################################

### Set workdir
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/maskdir
cd $WORKDIR

### Create directory structure
mkdir -p work/00/NEAT_CRAM/dedup
mkdir -p input

### Place (symlink) normal crams in ./work/00/NEAT_CRAM/dedup
### Place config samples.csv in ./input
### Place config SPECIES_INFO.txt in .


module load samtools-1.19.2/python-3.11.6
for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
samtools index $SAMPLE.neat.cram
done



## RENAME NEAT_CRAMS
tail -n +2 /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/maskdir/input/samples.csv | while IFS=',' read -r ID STDY _; do rename "s/^${ID}/${STDY}/" ${ID}.neat.cram*; done



#################################
### MASK generation pipeline ####
#################################

### Define input and working directory
INPUT=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/maskdir/SPECIES_INFO.txt
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/maskdir
cd $WORKDIR

### Launch Step 1 (Coverage profiling)
bash /software/team294/cseg/SNPmask_pipeline_nextflow/1_CoverageHist.sh $INPUT

### Launch Step 2 (Variant calling)
bash /software/team294/cseg/SNPmask_pipeline_nextflow/2_VariantCalling.sh $INPUT

### Launch Step 3 (Mask Files Building Step) // Once 1 and 2 are completed
bash /software/team294/cseg/SNPmask_pipeline_nextflow/3_bsub_command.sh $INPUT

### Launch Step 4 (Merge Mask files) // Once 1, 2 and 3 are completed
bash /software/team294/cseg/SNPmask_pipeline_nextflow/4_MaskMerging.sh $INPUT


### 29032026

#######################
#### KRAKEN ANALYS ####
#######################
module load samtools-1.19

# Some aliases:
alias bsub40000='bsub -M 40000 -R "select[mem>40000] rusage[mem=40000]" '
alias bsub2000='bsub -M 2000 -R "select[mem>2000] rusage[mem=2000]" '

#### Set working directory
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken
cd $WORKDIR

# Get reads supporting mutations 
mkdir -p /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/supporting_mutations
mkdir -p /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/supporting_mutations/logs

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
echo $SAMPLE
STDERR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/supporting_mutations/logs/$SAMPLE.err
VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/$SAMPLE.filtered.vcf.gz
NEAT_CRAM=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/NEAT_CRAM/$SAMPLE.neat.cram
OUTPUT=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/supporting_mutations/$SAMPLE.fq
bsub2000 -e $STDERR -o $STDERR "perl /software/team294/cseg/bin/rbs_supporting_mutations.pl $VCF $NEAT_CRAM > $OUTPUT"
done



# Run Kraken analyis
mkdir -p /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/analysis
mkdir -p /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/analysis/logs

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do 
echo "Processing $SAMPLE"
STDERR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/analysis/logs/$SAMPLE.err
KRAKEN_DB=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/databases/KRAKEN_DB/buccals
UNCLASSIFIED_OUT=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/analysis/$SAMPLE.UNCLASS.kraken
CLASSIFIED_OUT=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/analysis/$SAMPLE.CLASS.kraken
REPORT=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/analysis/$SAMPLE.REPORT.kraken
FASTQ=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/supporting_mutations/$SAMPLE.fq
OUTPUT=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken/analysis/$SAMPLE.MAIN.kraken
echo "Supporting reads are: $FASTQ"
echo "Running Kraken"
bsub40000 -e $STDERR -o $STDERR "kraken2 --db $KRAKEN_DB --unclassified-out $UNCLASSIFIED_OUT -classified-out $CLASSIFIED_OUT --report $REPORT $FASTQ > $OUTPUT"
done

# Summarise Kraken results: (Update according to database)
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/kraken
cd $WORKDIR
for a in ./analysis/*REPORT.kraken; do echo $a `egrep "(Bos taurus)|(Sus scrofa)|(Meleagris gallopavo)|(Gallus gallus)|(Oryctolagus cuniculus)|(Mus musculus)|(Equus caballus)|(Ovis aries)|(Macaca mulatta)|(Homo sapiens)" $a|sort -k2,2n| tail -1`; done|tr ' ' '\t' | sort -k2,2n -r > KRAKEN.RESULTS.tsv














