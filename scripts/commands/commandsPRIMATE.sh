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

for a in ./analysis/*REPORT.kraken; do echo $a `egrep "(Macaca mulatta)" $a|sort -k2,2n| tail -1`; done|tr ' ' '\t' | sort -k2,2n -r > KRAKEN.RESULTS_Macaca_mulatta.tsv
for a in ./analysis/*REPORT.kraken; do echo $a `egrep "(Homo sapiens)" $a|sort -k2,2n| tail -1`; done|tr ' ' '\t' | sort -k2,2n -r > KRAKEN.RESULTS_Homo_sapiens.tsv
for a in ./analysis/*REPORT.kraken; do echo $a `egrep "(Mus musculus)" $a|sort -k2,2n| tail -1`; done|tr ' ' '\t' | sort -k2,2n -r > KRAKEN.RESULTS_Mus_musculus.tsv
for a in ./analysis/*REPORT.kraken; do echo $a `egrep "(Mus musculus)" $a|sort -k2,2n| tail -1`; done|tr ' ' '\t' | sort -k2,2n -r > KRAKEN.RESULTS_Mus_musculus.tsv

for a in ./analysis/*REPORT.kraken; do echo $a `egrep "(unclassified)" $a|sort -k2,2n| tail -1`; done|tr ' ' '\t' | sort -k2,2n -r | cut -f 1,2 > KRAKEN.RESULTS_unclassified.tsv


##############################
###### clean_ns Pipeline #####
##############################

See https://gitlab.internal.sanger.ac.uk/yi1/clean_ns

# 1- Prepare and cd to a dedicated staging area where symlinks to filtered bams are stored under ./contaminated/ under the name {sampleID}.bam (i.e. contaminated/PD69593a_tds0001.bam)
done

# 2- Prepare sample_list_preproc.txt under the root of the staging area. This is a text file with each line a PD ID (i.e. PD69593a_tds0001).
done

# 3 - cd to the staging area. Extract fastqs: bash {path to the repo}/01_write_fastq.sh
cd /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/clean_ns
bash /software/team294/cseg/clean_ns/src/01_write_fastq.sh

# 4 - Map to a chimeric reference genome:
cd /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/clean_ns
bash /software/team294/cseg/clean_ns/src/02_map_to_chimeric.sh /lustre/scratch126/casm/teams/team294/projects/cseg/resources/hybrid_genomes/macaca_human/macaca_human.fa

# 5 - Extract reads mapped to macaca contigs:
cd /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/clean_ns
bash /software/team294/cseg/clean_ns/src/03_extract_macaca_reads.sh 

## *also check reads mapped to human contigs
bash /software/team294/cseg/clean_ns/src/03_extract_human_reads_short.sh 

# 6 Prepare random reads in bundles:
cd /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/clean_ns
bash /software/team294/cseg/clean_ns/src/04_dedupe_bam.sh


#############################
#### clean_ns evaluation ####
#############################

05/06/2026

### Check number of reads in fastq:
cd /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/clean_ns/fastq
bsub count.sh


### Run samtools stats on each BAM
cd /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/clean_ns/bwa
module load samtools-1.19.2
for BAM in *.bam
do
echo $BAM
bsub5000 -e log/samtools_stats_$BAM.err -o log/samtools_stats_$BAM.err "samtools stats $BAM > ${BAM%.bam}.stats"
done

### Get summary:
{
    echo -e "sample\traw_total_sequences\treads_mapped\treads_unmapped"

    for f in *.stats; do
        sample=${f%%_tds*}

        raw=$(awk -F'\t' '$2=="raw total sequences:" {print $3}' "$f")
        mapped=$(awk -F'\t' '$2=="reads mapped:" {print $3}' "$f")
        unmapped=$(awk -F'\t' '$2=="reads unmapped:" {print $3}' "$f")

        echo -e "${sample}\t${raw}\t${mapped}\t${unmapped}"
    done
} > samtools_stats_summary_Macaca_mulatta.tsv



### Check number of macaque derived reads:
cd /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/clean_ns/filtered_bams

REPORT=summary_Macaca_mulatta.tsv
printf '%s\t%s\n' "sample" "clean_ns_reads" > "$REPORT"

for SAMPLE_ID in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
 do 
SAMPLE=${SAMPLE_ID%%_*}
 echo "screening $SAMPLE"
READS=$(wc -l < "${SAMPLE_ID}_macaca_reads.txt")
printf '%s\t%s\n' "$SAMPLE" "$READS" >> "$REPORT"
done

### Check number of human derived reads:
cd /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/clean_ns/filtered_bams

REPORT=summary_Homo_sapiens.tsv
printf '%s\t%s\n' "sample" "clean_ns_reads" > "$REPORT"

for SAMPLE_ID in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
 do 
SAMPLE=${SAMPLE_ID%%_*}
 echo "screening $SAMPLE"
READS=$(wc -l < "${SAMPLE_ID}_human_reads.txt")
printf '%s\t%s\n' "$SAMPLE" "$READS" >> "$REPORT"
done







#### Picard test on chimpanzee dupcaller-compliant bams
module load gatk-4.5.0.0/python-3.12.0
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/dupcaller/runs/latest
cd $WORKDIR

for SAMPLE in dupcaller-CHPD0002f_ds0001  dupcaller-CHPD0002b_ds0002
do
bsub40000 -e picard.${SAMPLE}.err -o picard.${SAMPLE}.err gatk MarkDuplicates -I ${SAMPLE}.bam -O ${SAMPLE}.mkdped.bam -M ${SAMPLE}.mkdp_metrics.txt \
    --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$" \
    --DUPLEX_UMI --TAGGING_POLICY OpticalOnly --BARCODE_TAG DB
done





##################################
#### Filter VCFs to PASS only ####
##################################
module load bcftools-1.9

workdir=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/filtered/pass
cd $workdir 
for sample in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
bcftools view -f PASS ../$sample.filtered.vcf.gz -Oz -o $sample.filtered.pass.vcf.gz
done



####################################
### Run Decontamination Pipeline ###
####################################
workdir=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/filtered/pass/decontamination
cd $workdir
mkdir -p $workdir/logs

module load bwa-0.7.17 
module load samtools-1.19.2/python-3.11.6 

for sample in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
 do echo $sample
 bsub20000 -e $workdir/logs/$sample.err -o $workdir/logs/$sample.err -n 2 -R "span[ptile=2]" perl /software/team294/cseg/decontamination-pipeline/rbs_supporting_mutations.v2_macaca.pl  \
   /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/filtered/pass/$sample.filtered.pass.vcf.gz \
   /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/NEAT_CRAM/$sample.neat.cram \
   $sample
 done


#######################################################
#### nanoseq2dupcaller conversion of macaque CRAMs ####
#######################################################

### 22/06/2026
module load samtools-1.19.2
nanoseq2dupcaller=/software/team294/cseg/DupCaller_preprocessing/full/nanoseq2dupcaller.sh
workdir=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/CRAM
cd $workdir

### Apply nanoseq2dupcaller conversion script to the CRAM files
for sample in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
cram="$sample".cram
bsub80000 -q long -e nanoseq2dupcaller_$sample.err -o nanoseq2dupcaller_$sample.err bash $nanoseq2dupcaller  $cram
done

for sample in MQD0001h_tds0001 MQD0002j_tds0001 MQD0006j_tds0001
do
cram="$sample".cram
bsub80000 -q long -e nanoseq2dupcaller_$sample.err -o nanoseq2dupcaller_$sample.err bash $nanoseq2dupcaller  $cram
done


#########################################
#### pileup_genotype across samples  ####
#########################################

module load python-3.12.0/perl-5.38.0 
module load conda
conda activate /nfs/users/nfs_m/ms84/.conda/envs/Dupcaller


### VAF 0.3
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.3
mkdir -p $WORKDIR
cd $WORKDIR

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
BAM=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/CRAM/"$SAMPLE".cram
BED=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/panel/targeted_panel_Mmul_10_withgenes.bed
REF=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa
VAF=0.3
OUTPUT=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.3/pileup_genotype_"$SAMPLE".txt


bsub20000 -e pileup_"$SAMPLE".err -o pileup_"$SAMPLE".out "python /software/team294/cseg/genotyping/pileup_genotype_v0.1.py \
  --bam $BAM \
  --bed $BED \
  --ref $REF \
  --vaf $VAF > $OUTPUT" 
done

### VAF 0.2
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.2
mkdir -p $WORKDIR
cd $WORKDIR

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
BAM=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/CRAM/"$SAMPLE".cram
BED=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/panel/targeted_panel_Mmul_10_withgenes.bed
REF=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa
VAF=0.2
OUTPUT=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.2/pileup_genotype_"$SAMPLE".txt


bsub20000 -e pileup_"$SAMPLE".err -o pileup_"$SAMPLE".out "python /software/team294/cseg/genotyping/pileup_genotype_v0.1.py \
  --bam $BAM \
  --bed $BED \
  --ref $REF \
  --vaf $VAF > $OUTPUT" 
done

### VAF 0.1
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.1
mkdir -p $WORKDIR
cd $WORKDIR

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
BAM=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/CRAM/"$SAMPLE".cram
BED=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/panel/targeted_panel_Mmul_10_withgenes.bed
REF=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa
VAF=0.1
OUTPUT=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.1/pileup_genotype_"$SAMPLE".txt

bsub20000 -e pileup_"$SAMPLE".err -o pileup_"$SAMPLE".out "python /software/team294/cseg/genotyping/pileup_genotype_v0.1.py \
  --bam $BAM \
  --bed $BED \
  --ref $REF \
  --vaf $VAF > $OUTPUT" 
done


### Combine into a single file

### VAF 0.1
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.1
cd $WORKDIR

OUTPUT="pileup_genotype_COMBINED.txt"
> "$OUTPUT"
# --- Set header ---
FIRST_FILE="/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.1/pileup_genotype_MQD0001d_tds0001.txt"
head -n 1 "$FIRST_FILE" | awk 'BEGIN{OFS="\t"}{print $0, "SAMPLE"}' >> "$OUTPUT"

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
    FILE=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.1/pileup_genotype_"$SAMPLE".txt

    # Añadir contenido (saltando la primera línea del header) + columna SAMPLE
    tail -n +2 "$FILE" | awk -v sample="$SAMPLE" 'BEGIN{OFS="\t"}{print $0, sample}' >> "$OUTPUT"
done

### VAF 0.2
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.2
cd $WORKDIR

OUTPUT="pileup_genotype_COMBINED.txt"
> "$OUTPUT"
# --- Set header ---
FIRST_FILE="/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.2/pileup_genotype_MQD0001d_tds0001.txt"
head -n 1 "$FIRST_FILE" | awk 'BEGIN{OFS="\t"}{print $0, "SAMPLE"}' >> "$OUTPUT"

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
    FILE=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.2/pileup_genotype_"$SAMPLE".txt

    # Añadir contenido (saltando la primera línea del header) + columna SAMPLE
    tail -n +2 "$FILE" | awk -v sample="$SAMPLE" 'BEGIN{OFS="\t"}{print $0, sample}' >> "$OUTPUT"
done

### VAF 0.3
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.3
cd $WORKDIR

OUTPUT="pileup_genotype_COMBINED.txt"
> "$OUTPUT"
# --- Set header ---
FIRST_FILE="/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.3/pileup_genotype_MQD0001d_tds0001.txt"
head -n 1 "$FIRST_FILE" | awk 'BEGIN{OFS="\t"}{print $0, "SAMPLE"}' >> "$OUTPUT"

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
    FILE=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/germline_pileup/vaf0.3/pileup_genotype_"$SAMPLE".txt

    # Añadir contenido (saltando la primera línea del header) + columna SAMPLE
    tail -n +2 "$FILE" | awk -v sample="$SAMPLE" 'BEGIN{OFS="\t"}{print $0, sample}' >> "$OUTPUT"
done



#######################################################
#### Run DupCaller on Macaque targeted sequencing  ####
#######################################################

### Short Tandem Repeat File

### Download STR annotation from UCSC
Mmul_10_str.bed  (https://genome.ucsc.edu/cgi-bin/hgTables)

### Copy to annotation directory
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/annotation
cd $WORKDIR

### Retrieve relevant columns
cat Mmul_10_str.bed |grep -v "^#" |awk 'OFS="\t" {print $6,$7,$8,$11,$2, $10}'  > str_regions_tmp.bed

### Change chromosome format from ucsc-like to refseq
awk '
BEGIN{FS=OFS="\t"}
FNR==NR {
    if ($0 ~ /^[^#]/) {
        ucsc=$10
        ref=$7

        gsub(/\r/, "", ucsc)
        gsub(/\r/, "", ref)

        if (ucsc != "" && ref != "" && ucsc != "na")
            map[ucsc]=ref
    }
    next
}
{
    if ($1 in map)
        $1=map[$1]
    print
}
' ../reference_files/GCF_003339765.1_Mmul_10_assembly_report.txt str_regions_tmp.bed > str_regions.bed


### Sort, compress and index
module load bedtools2-2.29.0/
sortBed -i str_regions.bed -faidx ../reference_files/genome.fa.fai > str_regions_sorted.bed
bgzip str_regions_sorted.bed
tabix str_regions_sorted.bed.gz


### ms84
### 01/07/2026

##########################################
#### Run DupCaller Singularity Image #####
########################################## 

### Index the reference genome
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files
cd $WORKDIR

module load singularity
DUPCALLER_IMAGE=/software/CASM/singularity/dupcaller/dupcaller-1.1.2-main-300626.sif

bsub80000 -e dupcaller_index.err -o dupcaller_index.out singularity exec \
  -B /lustre/scratch126/casm/teams/team294/projects/cseg \
  $DUPCALLER_IMAGE \
  DupCaller.py index \
  -f /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa \
  -s /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/annotation/str_regions_sorted.bed.gz



### Run dupcaller in targeted mode without matched normal
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1
cd $WORKDIR

module load singularity
DUPCALLER_IMAGE=/software/CASM/singularity/dupcaller/dupcaller-1.1.2-main-300626.sif
CONFIG=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/input/samples.csv

### Dupcaller call
tail -n +2 "$CONFIG" | while IFS=',' read -r SAMPLE DUPLEX NORMAL
do
bsub80000 -q long -n 16 -R "span[ptile=16]" -e dupcaller_call_$SAMPLE.err -o dupcaller_call_$SAMPLE.err singularity exec \
    -B /lustre/scratch126/casm/teams/team294/projects/cseg \
    $DUPCALLER_IMAGE \
    DupCaller.py call \
    -b /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/CRAM/dupcaller-$SAMPLE.bam \
    -f /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa \
    -r NC_041754.1 NC_041755.1 NC_041756.1 NC_041757.1 NC_041758.1 NC_041759.1 NC_041760.1 NC_041761.1 NC_041762.1 NC_041763.1 NC_041764.1 NC_041765.1 NC_041766.1 NC_041767.1 NC_041768.1 NC_041769.1 NC_041770.1 NC_041771.1 NC_041772.1 NC_041773.1 NC_041774.1 NC_027914.1 \
    -R /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/panel/targeted_panel_Mmul_10.sorted.merged_500.bed.gz \
    -g /lustre/scratch126/casm/teams/team294/projects/cseg/resources/masks/Macaca_mulatta/Mmul_10/mGAP/mGap.Rhesus_macaque.v3.0_renamed.vcf.gz \
    -o results.$SAMPLE \
    -p 16 \
    -maf 0.1 \
    -ax 10 \
    -d 25
done

### Save pass snvs
module load bcftools-1.9
CONFIG=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/input/samples.csv
tail -n +2 "$CONFIG" | while IFS=',' read -r SAMPLE DUPLEX NORMAL
do
VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.$SAMPLE/results.${SAMPLE}_snv.vcf
PASS_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.$SAMPLE/results.${SAMPLE}_snv.pass.vcf
bcftools view -f PASS $VCF |bcftools sort -o $PASS_VCF
bgzip $PASS_VCF
bcftools index $PASS_VCF.gz
done

### Save underpowered snvs
module load bcftools-1.9
CONFIG=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/input/samples.csv
tail -n +2 "$CONFIG" | while IFS=',' read -r SAMPLE DUPLEX NORMAL
do
VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.$SAMPLE/results.${SAMPLE}_snv.vcf
UNDERPOWERED_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.$SAMPLE/results.${SAMPLE}_snv.underpowered.vcf
bcftools view -f underpowered $VCF |bcftools sort -o $UNDERPOWERED_VCF
bgzip $UNDERPOWERED_VCF
bcftools index $UNDERPOWERED_VCF.gz
done



### Save pass indels
module load bcftools-1.9
CONFIG=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/input/samples.csv
tail -n +2 "$CONFIG" | while IFS=',' read -r SAMPLE DUPLEX NORMAL
do
VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.$SAMPLE/results.${SAMPLE}_indel.vcf
PASS_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.$SAMPLE/results.${SAMPLE}_indel.pass.vcf
bcftools view -f PASS $VCF |bcftools sort -o $PASS_VCF
bgzip $PASS_VCF
bcftools index $PASS_VCF.gz
done

### Save underpowered indels
module load bcftools-1.9
CONFIG=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/input/samples.csv
tail -n +2 "$CONFIG" | while IFS=',' read -r SAMPLE DUPLEX NORMAL
do
VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.$SAMPLE/results.${SAMPLE}_indel.vcf
UNDERPOWERED_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.$SAMPLE/results.${SAMPLE}_indel.underpowered.vcf
bcftools view -f underpowered $VCF |bcftools sort -o $UNDERPOWERED_VCF
bgzip $UNDERPOWERED_VCF
bcftools index $UNDERPOWERED_VCF.gz
done



### Dupcaller estimate
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1
cd $WORKDIR

module load singularity
DUPCALLER_IMAGE=/software/CASM/singularity/dupcaller/dupcaller-1.1.2-main-300626.sif
CONFIG=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/input/samples.csv
GENEBED=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/panel/targeted_panel_gene_exons.sorted.bed.gz
### Dupcaller estimate
tail -n +2 "$CONFIG" | while IFS=',' read -r SAMPLE DUPLEX NORMAL
do
bsub5000 -q normal -e dupcaller_estimate_$SAMPLE.err -o dupcaller_estimate_$SAMPLE.err singularity exec \
    -B /lustre/scratch126/casm/teams/team294/projects/cseg \
    $DUPCALLER_IMAGE \
    DupCaller.py estimate \
    -i results.$SAMPLE \
    -f /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa \
    -r NC_041754.1 NC_041755.1 NC_041756.1 NC_041757.1 NC_041758.1 NC_041759.1 NC_041760.1 NC_041761.1 NC_041762.1 NC_041763.1 NC_041764.1 NC_041765.1 NC_041766.1 NC_041767.1 NC_041768.1 NC_041769.1 NC_041770.1 NC_041771.1 NC_041772.1 NC_041773.1 NC_041774.1 NC_027914.1 \
    -gb $GENEBED
done


### Dupcaller summarize
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1
cd $WORKDIR

module load singularity
DUPCALLER_IMAGE=/software/CASM/singularity/dupcaller/dupcaller-1.1.2-main-300626.sif
bsub5000 -q normal -e dupcaller_summarize.err -o dupcaller_summarize.err singularity exec \
    -B /lustre/scratch126/casm/teams/team294/projects/cseg \
    $DUPCALLER_IMAGE \
    DupCaller.py summarize \
    -i results.MQD0001d_tds0001 results.MQD0001h_tds0001 results.MQD0001j_tds0001 results.MQD0002e_tds0001 results.MQD0002f_tds0001 results.MQD0002h_tds0001 results.MQD0002j_tds0001 results.MQD0002l_tds0001 results.MQD0003d_tds0001 results.MQD0003e_tds0001 results.MQD0003g_tds0001 results.MQD0003j_tds0001 results.MQD0004d_tds0001 results.MQD0004e_tds0002 results.MQD0004i_tds0001 results.MQD0004k_tds0002 results.MQD0006e_tds0001 results.MQD0006f_tds0001 results.MQD0006h_tds0001 results.MQD0006j_tds0001 results.MQD0006l_tds0001 \
    -o cohort_summary.txt




### SigProfilerPlotting
module load python-3.11.6/perl-5.38.0

export PATH=/nfs/users/nfs_m/ms84/.local/lib/python3.11/site-packages:$PATH
MATRIX_PATH=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/cohort_summary_SBS96_corrected.txt
OUTPUT_PATH=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/sigprofilerplotting
PROJECT=macaque
PLOT_TYPE=96
SigProfilerPlotting plotSBS $MATRIX_PATH $OUTPUT_PATH $PROJECT $PLOT_TYPE



### Plot with SigProfiler
python3
import sigProfilerPlotting as sigPlt
 sigPlt.plotSBS(matrix_path, output_path, project, plot_type, percentage=False)





### Run dupcaller in targeted mode without matched normal (params)
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/params
cd $WORKDIR

module load singularity
DUPCALLER_IMAGE=/software/CASM/singularity/dupcaller/dupcaller-1.1.2-main-300626.sif
CONFIG=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/params/input/samples.csv

### Dupcaller Variant Calling
tail -n +2 "$CONFIG" | while IFS=',' read -r SAMPLE DUPLEX NORMAL
do
bsub80000 -q long -n 16 -R "span[ptile=16]" -e dupcaller_call_$SAMPLE.err -o dupcaller_call_$SAMPLE.err singularity exec \
    -B /lustre/scratch126/casm/teams/team294/projects/cseg \
    $DUPCALLER_IMAGE \
    DupCaller.py call \
    -b /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/CRAM/dupcaller-$SAMPLE.bam \
    -f /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa \
    -r NC_041754.1 NC_041755.1 NC_041756.1 NC_041757.1 NC_041758.1 NC_041759.1 NC_041760.1 NC_041761.1 NC_041762.1 NC_041763.1 NC_041764.1 NC_041765.1 NC_041766.1 NC_041767.1 NC_041768.1 NC_041769.1 NC_041770.1 NC_041771.1 NC_041772.1 NC_041773.1 NC_041774.1 NC_027914.1 \
    -R /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/panel/targeted_panel_Mmul_10.sorted.merged_500.bed.gz \
    -g /lustre/scratch126/casm/teams/team294/projects/cseg/resources/masks/Macaca_mulatta/Mmul_10/mGAP/mGap.Rhesus_macaque.v3.0_renamed.vcf.gz \
    -o results.$SAMPLE \
    -p 16 \
    -maf 0.1 \
    -ax 50 \
    -d 25
done










### Compare results:
for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
echo
echo "Analyzing sample $SAMPLE"
echo
echo "dupcaller ax 10:"
cat /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.$SAMPLE/results."$SAMPLE"_snv.vcf |grep -v "^#" |awk 'OFS="\t"{if($7=="PASS") print $0}' |wc -l
echo "dupcaller ax 50:"
cat /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/params/results.$SAMPLE/results."$SAMPLE"_snv.vcf |grep -v "^#" |awk 'OFS="\t"{if($7=="PASS") print $0}' |wc -l
echo
echo
done



OUT=counts_dupcaller.tsv
echo -e "sample\tax10\tax50" > $OUT

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do

ax10=$(grep -v "^#" /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.$SAMPLE/results.${SAMPLE}_snv.vcf \
  | awk '$7=="PASS"' | wc -l)

ax50=$(grep -v "^#" /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/params/results.$SAMPLE/results.${SAMPLE}_snv.vcf \
  | awk '$7=="PASS"' | wc -l)

echo -e "${SAMPLE}\t${ax10}\t${ax50}" >> $OUT

done


OUT=counts_TDC_TNS.tsv
echo -e "sample\tTDC\tTNS" > $OUT

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do

TDC=$(grep -v "^#" /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.$SAMPLE/results.${SAMPLE}_snv.vcf \
  | awk '$7=="PASS"' | wc -l)

TNS=$(zcat /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/${SAMPLE}.vcf.gz \
  | awk '$0 !~ /^#/ && $7=="PASS" && $0 ~ /TYPE=snv/' \
  | wc -l)

echo -e "${SAMPLE}\t${TDC}\t${TNS}" >> $OUT

done


14/07/2026







#######################################################
### Compare targeted_dupcaller and targeted_nanoseq ###
#######################################################

### run comparisson for snv (PASS vs PASS)
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/snv
mkdir -p $WORKDIR
cd $WORKDIR

### Run bcftools isec
module load bcftools-1.9

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
DUPCALLER_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.${SAMPLE}/results.${SAMPLE}_snv.pass.vcf.gz
NANOSEQ_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/filtered/pass/${SAMPLE}.filtered.pass.snv.vcf.gz
bcftools isec $DUPCALLER_VCF $NANOSEQ_VCF -p isec_${SAMPLE}
done


### run comparisson for indel  (PASS vs PASS)
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/indel
mkdir -p $WORKDIR
cd $WORKDIR

### Run bcftools isec
module load bcftools-1.9

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
DUPCALLER_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.${SAMPLE}/results.${SAMPLE}_indel.pass.vcf.gz
NANOSEQ_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/filtered/pass/${SAMPLE}.filtered.pass.indel.vcf.gz
bcftools isec $DUPCALLER_VCF $NANOSEQ_VCF -p isec_${SAMPLE}
done

### creata summary table sorting by increasing purity levels
echo -e "sample\ttdc_private\ttns_private\tshared" > venn_counts.tsv
for SAMPLE in MQD0002f_tds0001 MQD0004i_tds0001 MQD0006l_tds0001 MQD0002l_tds0001 MQD0004d_tds0001 MQD0003e_tds0001 MQD0002h_tds0001 MQD0001j_tds0001 MQD0003j_tds0001 MQD0006f_tds0001 MQD0003g_tds0001 MQD0006e_tds0001 MQD0002e_tds0001 MQD0004e_tds0002 MQD0004k_tds0002 MQD0003d_tds0001 MQD0001d_tds0001 MQD0006h_tds0001
do
    DIR="isec_${SAMPLE}"
    tdc=$(grep -vc '^#' "$DIR/0000.vcf")
    tns=$(grep -vc '^#' "$DIR/0001.vcf")
    shared=$(grep -vc '^#' "$DIR/0002.vcf")
    echo -e "${SAMPLE}\t${tdc}\t${tns}\t${shared}" >> venn_counts.tsv
done








### run comparisson for snv (underpowered vs PASS)
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/snv_underpowered
mkdir -p $WORKDIR
cd $WORKDIR

### Run bcftools isec
module load bcftools-1.9

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
DUPCALLER_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.${SAMPLE}/results.${SAMPLE}_snv.underpowered.vcf.gz
NANOSEQ_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/filtered/pass/${SAMPLE}.filtered.pass.snv.vcf.gz
bcftools isec $DUPCALLER_VCF $NANOSEQ_VCF -p isec_${SAMPLE}
done

### creata summary table sorting by increasing purity levels
echo -e "sample\ttdc_private\ttns_private\tshared" > venn_counts.tsv
for SAMPLE in MQD0002f_tds0001 MQD0004i_tds0001 MQD0006l_tds0001 MQD0002l_tds0001 MQD0004d_tds0001 MQD0003e_tds0001 MQD0002h_tds0001 MQD0001j_tds0001 MQD0003j_tds0001 MQD0006f_tds0001 MQD0003g_tds0001 MQD0006e_tds0001 MQD0002e_tds0001 MQD0004e_tds0002 MQD0004k_tds0002 MQD0003d_tds0001 MQD0001d_tds0001 MQD0006h_tds0001
do
    DIR="isec_${SAMPLE}"
    tdc=$(grep -vc '^#' "$DIR/0000.vcf")
    tns=$(grep -vc '^#' "$DIR/0001.vcf")
    shared=$(grep -vc '^#' "$DIR/0002.vcf")
    echo -e "${SAMPLE}\t${tdc}\t${tns}\t${shared}" >> venn_counts.tsv
done


### run comparisson for indel  (underpowered vs PASS)
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/indel_underpowered
mkdir -p $WORKDIR
cd $WORKDIR

### Run bcftools isec
module load bcftools-1.9

for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
DUPCALLER_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_dupcaller/set1/results.${SAMPLE}/results.${SAMPLE}_indel.underpowered.vcf.gz
NANOSEQ_VCF=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/filtered/pass/${SAMPLE}.filtered.pass.indel.vcf.gz
bcftools isec $DUPCALLER_VCF $NANOSEQ_VCF -p isec_${SAMPLE}
done

### creata summary table sorting by increasing purity levels
echo -e "sample\ttdc_private\ttns_private\tshared" > venn_counts.tsv
for SAMPLE in MQD0002f_tds0001 MQD0004i_tds0001 MQD0006l_tds0001 MQD0002l_tds0001 MQD0004d_tds0001 MQD0003e_tds0001 MQD0002h_tds0001 MQD0001j_tds0001 MQD0003j_tds0001 MQD0006f_tds0001 MQD0003g_tds0001 MQD0006e_tds0001 MQD0002e_tds0001 MQD0004e_tds0002 MQD0004k_tds0002 MQD0003d_tds0001 MQD0001d_tds0001 MQD0006h_tds0001
do
    DIR="isec_${SAMPLE}"
    tdc=$(grep -vc '^#' "$DIR/0000.vcf")
    tns=$(grep -vc '^#' "$DIR/0001.vcf")
    shared=$(grep -vc '^#' "$DIR/0002.vcf")
    echo -e "${SAMPLE}\t${tdc}\t${tns}\t${shared}" >> venn_counts.tsv
done









### Inspect results (burdens)
for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
echo "Analyzing sample: $SAMPLE"
cat /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/dupcaller/Set1/rescue_dupcaller.$SAMPLE/rescue_dupcaller."$SAMPLE"_sbs_burden.txt |grep "Corrected burden" |head -n 1
echo
done

### Inspect results (snv calls)
for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
echo
echo "Analyzing sample $SAMPLE"
echo "dupcaller:"
cat /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/dupcaller/Set1/rescue_dupcaller.$SAMPLE/rescue_dupcaller."$SAMPLE"_snv.vcf |grep -v "^#" |awk 'OFS="\t"{if($7=="PASS") print $0}' |wc -l
echo "nanoseq:"
zcat /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/$SAMPLE.vcf.gz  |grep -v "^#" |awk 'OFS="\t"{if($7=="PASS") print $0}'|grep "TYPE=snv"| wc -l
echo
echo
done

### Inspect results (indel calls)
for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
echo
echo "Analyzing sample $SAMPLE"
echo "dupcaller:"
cat /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/dupcaller/Set1/rescue_dupcaller.$SAMPLE/rescue_dupcaller."$SAMPLE"_indel.vcf |grep -v "^#" |awk 'OFS="\t"{if($7=="PASS") print $0}' |wc -l
echo "nanoseq:"
zcat /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/$SAMPLE.vcf.gz  |grep -v "^#" |awk 'OFS="\t"{if($7=="PASS") print $0}'|grep -E "TYPE=ins|TYPE=del"| wc -l
echo
echo
done

### Inspect results (indel calls)
for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
zcat /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/VCF/$SAMPLE.vcf.gz  |grep -v "^#" |awk 'OFS="\t"{if($7=="PASS") print $0}'|grep -E "TYPE=ins|TYPE=del"
done



for SAMPLE in MQD0001d_tds0001 MQD0001j_tds0001 MQD0002e_tds0001 MQD0002f_tds0001 MQD0002h_tds0001 MQD0002l_tds0001 MQD0003d_tds0001 MQD0003e_tds0001 MQD0003g_tds0001 MQD0003j_tds0001 MQD0004d_tds0001 MQD0004e_tds0002 MQD0004i_tds0001 MQD0004k_tds0002 MQD0006e_tds0001 MQD0006f_tds0001 MQD0006h_tds0001 MQD0006l_tds0001
do
cat /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/benchmarking/tdc_tns_comparisson/indel/isec_$SAMPLE/0001.vcf |grep -v "^#" |cut -f 4 |sort |uniq -c
done





### 07/09/2026

###################################
#### Upgraded Mask Generation #####
###################################

################################################
### Mask generation for chimpanzee (nanoseq) ###
################################################

### Load modules
module load singularity
module load bedtools2-2.29.0/

### Set maskdir directory
OUTDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/nanoseq/set1/maskdir
mkdir -p $OUTDIR
cd $OUTDIR

### Step one: Counts per-position allele counts against a reference genome into a SQLite database with mmrate-count
bsub80000 -q normal -n 16 -R "span[ptile=16]" -e mmrate-count.err -o mmrate-count.err singularity run \
  -B /lustre/scratch126/casm/teams/team294/projects/cseg \
  /software/CASM/singularity/nanoseq-mask/mmrate-count.sif \
  /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Pan_troglodytes/NHGRI_mPanTro3-v2.1_pri/reference_files/genome.fa \
  /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/nanoseq/set1/NEAT_CRAM/CHPD0002b_ds0001.neat.cram \
  /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/nanoseq/set1/maskdir\
  CHPD0002b_ds0001

### Step two: Generates genome-wide noise/SNP exclusion masks from per-sample allele-count databases with mmrate-query 
bsub5000 -q normal -e mmrate-query.err -o mmrate-query.err singularity run \
  -B /lustre/scratch126/casm/teams/team294/projects/cseg \
    /software/CASM/singularity/nanoseq-mask/mmrate-query.sif masks CHPD0002b_ds0001.counts.db \
  --config /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/nanoseq/set1/MASK/config.json \
  --output-dir /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Pan_troglodytes/nanoseq/set1/maskdir

### Collapse into contiguous regions
sort -k1,1 -k2,2n noise.bed | bedtools merge > mask.bed



######################################################
### Mask generation for macaque (targeted_nanoseq) ###
######################################################

### Load modules
module load singularity
module load bedtools2-2.29.0/

### Set maskdir directory
REF=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Macaca_mulatta/Mmul_10/reference_files/genome.fa
OUTDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/maskdir
mkdir -p $OUTDIR
cd $OUTDIR

### Step one: Counts per-position allele counts against a reference genome into a SQLite database with mmrate-count
tail -n +2 /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/cohort/config/samples.csv | while IFS=',' read -r SAMPLE DUPLEX NORMAL
do
BAM=/lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/NEAT_CRAM/"$SAMPLE".neat.cram

bsub80000 -q normal -n 16 -R "span[ptile=16]" -e mmrate-count_"$SAMPLE".err -o mmrate-count_"$SAMPLE".err singularity run \
  -B /lustre/scratch126/casm/teams/team294/projects/cseg \
  /software/CASM/singularity/nanoseq-mask/mmrate-count.sif \
  $REF \
  $BAM \
  $OUTDIR \
  $SAMPLE
done

### (!) Two samples failed due to issues with NEAT_CRAM
module load samtools-1.19.2/python-3.11.6
samtools quickcheck -v /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/NEAT_CRAM/MQD0002l_tds0001.neat.cram
samtools quickcheck -v /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/NEAT_CRAM/MQD0003d_tds0001.neat.cram

MQD0002l_tds0001.neat.cram was missing EOF block when one should be present.
MQD0003d_tds0001.neat.cram was missing EOF block when one should be present.
##########################

### Step two: Generates genome-wide noise/SNP exclusion masks from per-sample allele-count databases with mmrate-query 
bsub5000 -q normal -e mmrate-query.err -o mmrate-query.err singularity run \
  -B /lustre/scratch126/casm/teams/team294/projects/cseg \
    /software/CASM/singularity/nanoseq-mask/mmrate-query.sif masks *.counts.db \
  --config /lustre/scratch126/casm/teams/team294/projects/cseg/resources/masks/config/mask_parameters.json \
  --output-dir /lustre/scratch126/casm/teams/team294/projects/cseg/PrimateSomatic/Macaca_mulatta/targeted_nanoseq/maskdir

### Collapse into contiguous regions
sort -k1,1 -k2,2n noise.bed | bedtools merge > mask.bed
