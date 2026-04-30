#### Martín Santamarina García
#### 16/04/2026


### Run nanoseq short-mode
module load ISG/IRODS/1.0
iinit
cd /lustre/scratch126/casm/teams/team294/projects/cseg/ReptileSomatic/Varanus_komodoensis/nanoseq/set1
bsub < runNext.sh

cd /lustre/scratch126/casm/teams/team294/projects/cseg/ReptileSomatic/Ophiophagus_hannah/nanoseq/set1
bsub < runNext.sh

cd /lustre/scratch126/casm/teams/team294/projects/cseg/ReptileSomatic/Aldabrachelys_gigantea/nanoseq/set1
bsub < runNext.sh


### Save CRAM files with sample name
for SPECIES in Varanus_komodoensis Ophiophagus_hannah Aldabrachelys_gigantea
do
echo "Processing $SPECIES"
cd /lustre/scratch126/casm/teams/team294/projects/cseg/ReptileSomatic/${SPECIES}/nanoseq/set1

#### Save CRAM files
mkdir -p CRAM
mkdir -p CRAM/duplex
mkdir -p CRAM/normal

## Duplex samples 
cut -f1,2,3 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2,3 | grep "id" -v  | while read -r ID DUPLEX NORMAL ; do mv work/*/*/add_rb/${DUPLEX}*.cram CRAM/duplex/${ID}.cram  ; done 
cut -f1,2,3 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2,3 | grep "id" -v  | while read -r ID DUPLEX NORMAL; do mv work/*/*/add_rb/${DUPLEX}*.cram.crai CRAM/duplex/${ID}.cram.crai  ; done 

## Normal samples
cut -f1,2,3 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2,3 | grep "id" -v  | while read -r ID DUPLEX NORMAL ; do mv work/*/*/add_rb/${NORMAL}*.cram CRAM/normal/${ID}_normal.cram  ; done 
cut -f1,2,3 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2,3 | grep "id" -v  | while read -r ID DUPLEX NORMAL; do mv work/*/*/add_rb/${NORMAL}*.cram.crai CRAM/normal/${ID}_normal.cram.crai  ; done 

#### Save NEAT_CRAM files
mkdir -p NEAT_CRAM
mkdir -p NEAT_CRAM/duplex
mkdir -p NEAT_CRAM/normal

## Duplex samples 
cut -f1,2,3 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2,3 | grep "id" -v  | while read -r ID DUPLEX NORMAL ; do mv work/*/*/dedup/${DUPLEX}*.neat.cram NEAT_CRAM/duplex/${ID}.neat.cram  ; done 
cut -f1,2,3 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2,3 | grep "id" -v  | while read -r ID DUPLEX NORMAL; do mv work/*/*/dedup/${DUPLEX}*.neat.cram.crai NEAT_CRAM/duplex/${ID}.neat.cram.crai  ; done 

## Normal samples
cut -f1,2,3 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2,3 | grep "id" -v  | while read -r ID DUPLEX NORMAL ; do mv work/*/*/dedup/${NORMAL}*.neat.cram NEAT_CRAM/normal/${ID}_normal.neat.cram  ; done 
cut -f1,2,3 -d "," input/samples.csv | sed "s/[,,:]/\t/g" | cut -f1,2,3 | grep "id" -v  | while read -r ID DUPLEX NORMAL; do mv work/*/*/dedup/${NORMAL}*.neat.cram.crai NEAT_CRAM/normal/${ID}_normal.neat.cram.crai  ; done 

done


### Clean work dir after saving CRAM and NEAT_CRAM files
#find work/ -type f ! -name ‘*.command.out’ ! -name ‘*.command.err’ -delete
find /lustre/scratch126/casm/teams/team294/projects/cseg/ReptileSomatic/Varanus_komodoensis/nanoseq/set1/work/ -type f ! -name ‘*.command.out’ ! -name ‘*.command.err’ -delete
find /lustre/scratch126/casm/teams/team294/projects/cseg/ReptileSomatic/Ophiophagus_hannah/nanoseq/set1/work/ -type f ! -name ‘*.command.out’ ! -name ‘*.command.err’ -delete
find /lustre/scratch126/casm/teams/team294/projects/cseg/ReptileSomatic/Aldabrachelys_gigantea/nanoseq/set1/work/ -type f ! -name ‘*.command.out’ ! -name ‘*.command.err’ -delete


### Apply nanoseq2dupcaller conversion (FULL)
module load samtools-1.19.2
NANOSEQ2DUPCALLER=/software/team294/cseg/DupCaller_preprocessing/full/nanoseq2dupcaller.sh

## Convert Varanus_komodoensis
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/ReptileSomatic/Varanus_komodoensis/nanoseq/set1/CRAM/duplex
cd $WORKDIR
mkdir -p logs
for SAMPLE in KDRGD0001b_ds0001 KDRGD0001b_ds0002
do
bsub80000 -e logs/$SAMPLE.err -o logs/$SAMPLE.err "bash $NANOSEQ2DUPCALLER $SAMPLE.cram"
done

## Convert Ophiophagus_hannah
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/ReptileSomatic/Ophiophagus_hannah/nanoseq/set1/CRAM/duplex
cd $WORKDIR
mkdir -p logs
for SAMPLE in KCOBD0001b_ds0001 KCOBD0001b_ds0002 KCOBD0001b_ds0003
do
bsub80000 -e logs/$SAMPLE.err -o logs/$SAMPLE.err "bash $NANOSEQ2DUPCALLER $SAMPLE.cram"
done

## Convert Aldabrachelys_gigantea
WORKDIR=/lustre/scratch126/casm/teams/team294/projects/cseg/ReptileSomatic/Aldabrachelys_gigantea/nanoseq/set1/CRAM/duplex
cd $WORKDIR
mkdir -p logs
for SAMPLE in AGD0001b_ds0001 AGD0002b_ds0001 AGD0003b_ds0001 AGD0004b_ds0001 AGD0005b_ds0001 AGD0006b_ds0001 AGD0007b_ds0001 AGD0008b_ds0001 AGD0009b_ds0001 AGD0010b_ds0001
do
bsub80000 -e logs/$SAMPLE.err -o logs/$SAMPLE.err "bash $NANOSEQ2DUPCALLER $SAMPLE.cram"
done
