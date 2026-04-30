### 16022026
### Martín 


### Activate the environment
module load conda
conda activate /nfs/users/nfs_m/ms84/.conda/envs/Dupcaller
module load samtools-1.19.2

### Index the Ovis aries reference genome.
REFDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Ovis_aries/Oar_rambouillet_v1.0/reference_files/
cd $REFDIR
bsub80000 -e kk_dupcaller_index.err -o kk_dupcaller_index.out DupCaller.py index -f genome.fa



for SAMPLE in CHPD0002b_ds0001 CHPD0002b_ds0002 CHPD0002f_ds0001
do
CRAM="$SAMPLE".cram
bsub80000 -q long -e nanoseq2dupcaller_$SAMPLE.err -o nanoseq2dupcaller_$SAMPLE.err bash nanoseq2dupcaller.sh  $CRAM
done

Job <145613> is submitted to queue <long>.



###### Apply nanoseq2dupcaller conversion to Ovis aries samples
IRODS2DUPCALLER=/software/team294/cseg/DupCaller/DupCaller_preprocessing/semi/irods2dupcaller.sh
WORKDIR=/lustre/scratch127/casm/projects/mutographs/ac36/sheep/staged_bams_Nanoseq
module load samtools-1.19.2

for SAMPLE in OARID0076b_ds0001
do
SAMPLEDIR="$WORKDIR"/3466/"$SAMPLE"/mapped_sample
cd $SAMPLEDIR
bsub40000 -q long -e irods2dupcaller_$SAMPLE.err -o irods2dupcaller_$SAMPLE.out bash $IRODS2DUPCALLER  $SAMPLE.sample.dupmarked.bam
done


module load samtools-1.19.2
for SAMPLE in OARID0076b_ds0001
do
bsub20000 -q small -e samtools_index_$SAMPLE.err -o samtools_index_$SAMPLE.err samtools index dupcaller-$SAMPLE.sample.dupmarked.bam
done


### Apply nanoseq2dupcaller conversion to Ovis aries samples using a sample list
while read SAMPLE; do
bsub80000 -q long -e nanoseq2dupcaller_$SAMPLE.err -o nanoseq2dupcaller_$SAMPLE.out bash $NANOSEQ2DUPCALLER  $SAMPLE.sample.dupmarked.bam
done < stageBam_sample_list.txt





### Dupcaller Variant Calling // rescue
sample="OARID0076b_ds0001"
bsub80000 -q long -n 16 -R "span[ptile=16]" -e $sample.dupcaller.err -o $sample.dupcaller.err DupCaller.py call \
                   -b /lustre/scratch127/casm/projects/mutographs/ac36/sheep/staged_bams_Nanoseq/3466/OARID0076b_ds0001/mapped_sample/test/dupcaller-OARID0076b_ds0001.sample.dupmarked.bam \
                   -f /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Ovis_aries/Oar_rambouillet_v1.0/reference_files/genome.fa \
                   -r 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 X \
                   -p 16 \
                   -n /lustre/scratch127/casm/projects/mutographs/ac36/sheep/staged_bams_WGS/3446/OARID0076b/mapped_sample/OARID0076b.sample.dupmarked.bam \
                   -m /lustre/scratch124/casm/references/pipeline_ref/Ovis_aries/Oar_rambouillet_v1.0/botseq/SNP+NOISE.NV2.SeenNV1.Sheep.v2.bed.gz \
                   -o dupcaller.$sample \
                   -ax 50 \
                   -tt 8 \
                   -tr 8 \
                   -d 15 \
                   -nm 4 \
                   --rescue TRUE


### DupCaller Burden Estimation
for SAMPLE in OARID0076b_ds0001
do
WORKDIR=/lustre/scratch127/casm/projects/mutographs/ac36/sheep/staged_bams_Nanoseq/3466/"$SAMPLE"/mapped_sample/test/
cd $WORKDIR
bsub5000 -q normal -e burdens_$SAMPLE.err -o burdens_$SAMPLE.err DupCaller.py estimate \
        -i dupcaller.$SAMPLE \
        -f /lustre/scratch126/casm/teams/team294/projects/cseg/resources/reference_genomes/Ovis_aries/Oar_rambouillet_v1.0/reference_files/genome.fa \
        -r 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 X
done




