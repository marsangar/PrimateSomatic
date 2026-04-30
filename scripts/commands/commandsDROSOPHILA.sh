#### Martín Santamarina García
#### 23/01/2026

#### Drosophila commands
cat Drosophila_sequenced.txt |grep "Adult" |awk 'OFS=","{print $11, $13, $13}' > samples_adult.txt

cat Drosophila_sequenced.txt |grep -E 'embryo|instar|Pupa' |awk 'OFS=","{print $11, $13, $13}' > samples_preadult.txt

cat Drosophila_sequenced.txt |grep "Grand" |awk 'OFS=","{print $11, $13, $13}' > samples_grandparent.txt

### Run NanoSeq
cd 
module add ISG/IRODS/1.0
iinit
bsub < runNext.sh


### Ree movall but 
find work/ -type f ! -name ‘*.command.out’ ! -name ‘*.command.err’ -delete




### Transfer results

#### Save CRAM files with sample name
cut -f1,2 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/add_rb/${ID}*.cram CRAM/${SAMPLE}.cram  ; done 
cut -f1,2 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/add_rb/${ID}*.cram.crai CRAM/${SAMPLE}.cram.crai  ; done 

#### Save NEAT CRAM files with sample name
cut -f1,2 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/dedup/${ID}*.neat.cram NEAT_CRAM/${SAMPLE}.neat.cram  ; done 
cut -f1,2 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2 | grep "id" -v  | while read -r SAMPLE ID ; do mv work/*/*/dedup/${ID}*.neat.cram.crai NEAT_CRAM/${SAMPLE}.neat.cram.crai  ; done 



#### DupCaller Analysis on Drosophila
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/dupcaller/CRAM
cd $WORKDIR

for SAMPLE in FFLYD0136b_ds0001 FFLYD0166b_ds0001 FFLYD0219b_ds0001
do
cp /lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/nanoseq/preadult/CRAM/${SAMPLE}.cram $WORKDIR/
cp /lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/nanoseq/preadult/CRAM/${SAMPLE}.cram.crai $WORKDIR/
cp /lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/nanoseq/adult/CRAM/${SAMPLE}.cram $WORKDIR/
cp /lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/nanoseq/adult/CRAM/${SAMPLE}.cram.crai $WORKDIR/
done

### Apply nanoseq2dupcaller conversion (ARCHIVE)
module load samtools-1.19.2
NANOSEQ2DUPCALLER=/software/team294/cseg/DupCaller_preprocessing/full/archive/nanoseq2dupcaller.sh 

for SAMPLE in FFLYD0136b_ds0001 FFLYD0166b_ds0001 FFLYD0219b_ds0001
do
CRAM="$SAMPLE".cram
bsub80000 -q long -e logs/nanoseq2dupcaller_$SAMPLE.err -o logs/nanoseq2dupcaller_$SAMPLE.err bash $NANOSEQ2DUPCALLER $CRAM
done



### Apply nanoseq2dupcaller conversion (NEW)
module load samtools-1.19.2
NANOSEQ2DUPCALLER=/software/team294/cseg/DupCaller_preprocessing/full/nanoseq2dupcaller.sh 

for SAMPLE in FFLYD0136b_ds0001 FFLYD0166b_ds0001 FFLYD0219b_ds0001
do
CRAM="$SAMPLE".cram
bsub80000 -q long -e logs/nanoseq2dupcaller_$SAMPLE.err -o logs/nanoseq2dupcaller_$SAMPLE.err bash $NANOSEQ2DUPCALLER $CRAM
done

for SAMPLE in FFLYD0136b_ds0001 FFLYD0166b_ds0001 FFLYD0219b_ds0001
do
cp $SAMPLE.nano2dupcaller.tmp/dupcaller-$SAMPLE.bam.bai .
done



### Index the output
module load samtools-1.19.2
for SAMPLE in FFLYD0136b_ds0001 FFLYD0166b_ds0001 FFLYD0219b_ds0001
do
bsub20000 -q small -e samtools_index_$SAMPLE.err -o samtools_index_$SAMPLE.err samtools index dupcaller-$SAMPLE.bam
done



#######################
#### Run DupCaller ####
#######################

### Activate the environment
module load conda
conda activate /software/team294/cseg/conda_envs/dupcaller
module load samtools-1.19.2

### Go to working directory
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/dupcaller/CRAM_new
cd $WORKDIR

## Dupcaller with matched-normal
for SAMPLE in FFLYD0136b_ds0001 FFLYD0166b_ds0001 FFLYD0219b_ds0001
do
bsub80000 -q long -n 16 -R "span[ptile=16]" -e dupcaller_$SAMPLE.err -o dupcaller_$SAMPLE.err DupCaller.py call \
                   --bam /lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/dupcaller/CRAM_new/dupcaller-$SAMPLE.bam \
                   --reference /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Drosophila_melanogaster/Release_6_plus_ISO1_MT/reference_files/genome.fa \
                   --regions NC_004354.4 NT_033779.5 NT_033778.4 NT_037436.4 NT_033777.3 NT_033777.3 NC_024512.1   \
                   --normalBam /lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/dupcaller/CRAM_new/dupcaller-$SAMPLE.bam \
                   --noise /lustre/scratch126/casm/teams/team294/projects/cseg/resources/masks/Drosophila_melanogaster/Release_6_plus_ISO1_MT/SNP+NOISE.NV2.SeenNV1.Drosophila_melanogaster.bed.gz \
                   --output dupcaller.$SAMPLE \
                   --threads 16 \
                   --rescue TRUE
done

## Dupcaller without matched-normal (maxAF)
for SAMPLE in FFLYD0136b_ds0001 FFLYD0166b_ds0001 FFLYD0219b_ds0001
do
bsub80000 -q long -n 16 -R "span[ptile=16]" -e  maxAF.dupcaller_$SAMPLE.err -o  maxAF.dupcaller_$SAMPLE.err DupCaller.py call \
                   --bam /lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/dupcaller/CRAM/dupcaller-$SAMPLE.bam \
                   --reference /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Drosophila_melanogaster/Release_6_plus_ISO1_MT/reference_files/genome.fa \
                   --regions NC_004354.4 NT_033779.5 NT_033778.4 NT_037436.4 NT_033777.3 NT_033777.3 NC_024512.1   \
                   --noise /lustre/scratch126/casm/teams/team294/projects/cseg/resources/masks/Drosophila_melanogaster/Release_6_plus_ISO1_MT/SNP+NOISE.NV2.SeenNV1.Drosophila_melanogaster.bed.gz \
                   --output maxAF.dupcaller.$SAMPLE \
                   --threads 16 \
                   --maxAF 0.1 \
                   --rescue TRUE
done



### Assess number of variants calls
for SAMPLE in FFLYD0136b_ds0001 FFLYD0166b_ds0001 FFLYD0219b_ds0001
do
echo "Evaluating "$SAMPLE" executions:"
echo "Dupcaller SNVs"
cat /lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/dupcaller/CRAM_new/dupcaller.$SAMPLE/dupcaller."$SAMPLE"_snv.vcf |grep -v "^#" |grep "PASS" |wc -l
cat /lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/dupcaller/CRAM_new/dupcaller.$SAMPLE/dupcaller."$SAMPLE"_indel.vcf |grep -v "^#" |grep "PASS" |wc -l
done


#### Estimate DupCaller Burdens
for SAMPLE in FFLYD0136b_ds0001 FFLYD0166b_ds0001 FFLYD0219b_ds0001
do
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/dupcaller/CRAM_new/"$SAMPLE"
cd $WORKDIR
bsub5000 -q normal -e burdens_$SAMPLE.err -o burdens_$SAMPLE.err DupCaller.py estimate \
        -i dupcaller.$SAMPLE \
        -f /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Drosophila_melanogaster/Release_6_plus_ISO1_MT/reference_files/genome.fa \
        -r NC_004354.4 NT_033779.5 NT_033778.4 NT_037436.4 NT_033777.3 NT_033777.3 NC_024512.1
done




#######################
#### SAMTOOLS STATS ###
#######################

### Set working directory
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/nanoseq/preadult/irods_data
#WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/nanoseq/adult/irods_data
cd $WORKDIR

### move CRAM files to the working directory
mv ../output/irods_data/6416/*/*/*.cram .



### Generate samtools stats for each CRAM file
module load samtools-1.19.2
for SAMPLE in *.cram
do
echo $SAMPLE
bsub5000 -e logs/samtools_stats_$SAMPLE.err -o logs/samtools_stats_$SAMPLE.err "samtools stats $SAMPLE > ${SAMPLE%.cram}.stats"
done







### Plot samtools stats for each CRAM file

cd /lustre/scratch126/casm/teams/team294/projects/cseg/FlySoma/Drosophila_melanogaster/irods_data/stats

module load samtools-1.19.2
for SAMPLE in *.stats
do
echo $SAMPLE
plot-bamstats -p ../plots/${SAMPLE%.stats} $SAMPLE
done



###########################################
#### Remove NanoSeq intermediate files ####
###########################################

### helpful commands
find work/ -type f ! -name ‘*.command.out’ ! -name ‘*.command.err’ -delete