#!/bin/bash
#BSUB -q basement
#BSUB -J NF_chimp_test
#BSUB -n 1
#BSUB -R "select[mem>5000] rusage[mem=5000]" -M5000
#BSUB -o log.nf.%J

# USER OPTIONS ---------------------------------------------------------------------------#
# For additional variant calling options and default values, see beginning of the script:
# /nfs/dog_n_devil/adrian/software/NanoSeq_Nextflow/NanoSeq_main.nf
#-----------------------------------------------------------------------------------------#
# Input file paths
INPUT=/nfs/casm/team294rr/cseg/projects/PrimateSomatic/03_scripts/configs/nanoseq/chimpanzee_config.csv
REF=/lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/pan_troglodytes/NHGRI_mPanTro3-v2.1_pri/reference_files/genome.fa
TRINUC=/lustre/scratch125/casm/teams/team294/projects/cseg/reference_genomes/pan_troglodytes/NHGRI_mPanTro3-v2.1_pri/reference_files/ref_freqs.txt
#MASK=/nfs/casm/team294rr/cseg/projects/PrimateSomatic/03_scripts/configs/nanoseq/chimpanzee_mask_file.bed.gz
MASK=/nfs/casm/team294rr/cseg/projects/PrimateSomatic/03_scripts/configs/nanoseq/empty_mask_file.bed.gz


# Output directory
OUTPUT=$PWD/output
mkdir -p $OUTPUT

# Remap option: if REMAP="--remap", CRAMs are remapped against reference genome
REMAP="--remap"
#REMAP=""

# Short-run mode: if SHORT="--short", pipeline stops after NANOSEQ_EFFI step
#SHORT="--short"
SHORT="--short"

# Keep-temp option: if KEEP_TEMP="--keep_temp", `work` directory will NOT be deleted after success
#                   (this option is not needed if already specifying SHORT="--short")
#KEEP_TEMP="--keep_temp"
KEEP_TEMP=""
#-----------------------------------------------------------------------------------------#
echo "Testing Stage 1"

# Load modules
module purge
#module add R                 # R/4.1.0 on farm5
module add singularity       # singularity/3.9.0 on farm5
module load HGI/softpack/groups/hgi/ubeR
module add samtools-1.19.2   # (farm22)
module add nextflow-23.10.0  # (farm22)
export NXF_VER=22.04.5
NANOSEQ=/nfs/dog_n_devil/adrian/software/NanoSeq_Nextflow/NanoSeq_main.nf
echo "Testing Stage 2"

# Check mandatory parameters
if [ -z "$INPUT" ] || [ -z "$OUTPUT" ] || [ -z "$REF" ] || [ -z "$MASK" ] || [ -z "$TRINUC" ]; then
    echo -e "\nERROR! Some of the following variables are missing: INPUT, OUTPUT, REF, MASK, TRINUC\n"
    exit
fi

# Parse study ID
STUDY=`tail -n+2 $INPUT | cut -f2 -d, | head -c4`
echo -e "\nStudy: $STUDY\n"

if [ "$SHORT" = "--short" ] || [ "$KEEP_TEMP" = "--keep_temp" ]; then
    echo -e "NOTE: Temporary files will be preserved after success ('work' directory)\n"
else
    echo -e "NOTE: Temporary files will be deleted after success ('work' directory)\n"
fi

# Launch pipeline
nextflow run $NANOSEQ \
    --jobs 200 -qs 20000 -profile lsf_singularity -resume $KEEP_TEMP $SHORT $REMAP \
    --sample_sheet $INPUT --outDir $OUTPUT --study $STUDY \
    --ref $REF --noise_bed $MASK --post_triNuc $TRINUC \


# TO LAUNCH (FARM22):
#    module load IRODS/1.0; iinit
#    bsub < runNanoseq_komodo.bsub
