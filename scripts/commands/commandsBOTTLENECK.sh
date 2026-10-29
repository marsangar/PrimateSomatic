####
#### 

#################################################################
#### Tutorial: download files from a remote server with lftp ####
#################################################################

### lftp is a command-line tool for transferring files between your computer and remote servers, particularly via SFTP, FTP, or FTPS. It is practical for HPC servers and large transfers because it allows for resuming downloads and performing robust transfers.

## Example: Connecting to the CGS ftp server via Linux command line

### Step 1. Request a username and password from your project manager.
Username: xxxxxxxxxx
Password: xxxxxxxxxx

### Step 2. Connect to the remote server using the following commands substituting your user name and password for “<user>” and “<password>” in the final line.
lftp
set ssl:verify-certificate no
set ftp:ssl-force true
set ftp:ssl-protect-data true
open -u <user>,<password> 131.111.189.67:991

### Step 3. Download datasets using lftp commands

## Explore and navigate local directory
lpwd # show local directory
lcd # change local directory

# Download file:
get file.txt

# Upload file:
put file.txt

# Download multiple files:
mget -c "*.bam" /local/path # important: -c for continuous transfer if connection breaks

# Upload multiple files:
mput -c "*.bam" /local/path # important: -c for continuous transfer if connection breaks

# Download a directory and its content
mirror -c sequencing_data  # important: -c for continuous transfer if connection breaks

# To exit the lftp session
exit





### Your files should be available to view in the home directory.


### 06/05/2026 ###

### Run dupcaller on a array of samples
sbatch --array=1-$(wc -l < ../config/samples.txt) slurm_dupcaller_pipeline




#### 28/09/2026

### Explore dupcaller results
CONFIG=/home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/projects/BottleneckOptimization/Canis_lupus_familiaris/nanoseq/NovaSeq_CRUK_data/config/samples.txt
WORKDIR=/home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/projects/BottleneckOptimization/Canis_lupus_familiaris/nanoseq/NovaSeq_CRUK_data/dupcaller
cd $WORKDIR




echo "Mutation burden"
#echo "Number of PASS SNVs:"
#echo "Number of PASS INDELs:"

while read SAMPLE
do
cat /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/projects/BottleneckOptimization/Canis_lupus_familiaris/nanoseq/NovaSeq_CRUK_data/dupcaller/"$SAMPLE"/"$SAMPLE"_sbs_burden.txt |awk 'FNR == 5'
#cat /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/projects/BottleneckOptimization/Canis_lupus_familiaris/nanoseq/NovaSeq_CRUK_data/dupcaller/"$SAMPLE"/"$SAMPLE"_snv.vcf |grep -v "#" |grep "PASS" |wc -l
#cat /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/projects/BottleneckOptimization/Canis_lupus_familiaris/nanoseq/NovaSeq_CRUK_data/dupcaller/"$SAMPLE"/"$SAMPLE"_indel.vcf |grep -v "#"|grep "PASS" |wc -l
done < $CONFIG


###################################
### Summarize dupcaller results ###
###################################

### ms84
### Important: double-check if there are changes dupcaller output files format - code requires to be updated if so

#### Specify config and working directory
CONFIG=/home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/projects/BottleneckOptimization/Canis_lupus_familiaris/nanoseq/NovaSeq_CRUK_data/config/samples.txt
WORKDIR=/home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/projects/BottleneckOptimization/Canis_lupus_familiaris/nanoseq/NovaSeq_CRUK_data/dupcaller

cd $WORKDIR


### Generate summary table for stats across samples 
OUTPUT="${WORKDIR}/summary_stats.txt"

# Header
printf "Sample\tNumber_Read_Families\tNumber_Passfilter_Reads\tNumber_Effective_Read_Families\tEffective_Coverage\tUnmasked_Coverage\tEffective_Indel_Coverage\tUnmasked_Indel_Coverage\tPer_Read_Family_Coverage\tPassfilter_Duplication_Rate\tEfficiency\n" \
    > "$OUTPUT"
while IFS= read -r SAMPLE; do

    # Skip empty lines and comments
    [[ -z "$SAMPLE" || "$SAMPLE" =~ ^# ]] && continue

    STATS_FILE="$WORKDIR/$SAMPLE/${SAMPLE}_stats.txt"

    awk -v sample="$SAMPLE" '
        BEGIN { OFS="\t" }

        {
            # Everything after the metric name is the value
            value = $NF
            values[NR] = value
        }

        END {
            print sample,
                  values[1],
                  values[2],
                  values[3],
                  values[4],
                  values[5],
                  values[6],
                  values[7],
                  values[8],
                  values[9],
                  values[10]
        }
    ' "$STATS_FILE" >> "$OUTPUT"

done < "$CONFIG"




### Generate summary table for sbs burden across samples 
OUTPUT="${WORKDIR}/summary_sbs_burden.txt"

# Header
printf "Sample\tUncorrected_burden\tUncorrected_burden_95_lower\tUncorrected_burden_95_upper\tUncorrected_mutation_number\tCorrected_burden\tCorrected_burden_95_lower\tCorrected_burden_95_upper\tmutation_number_per_genome\tgenome_coverage\tUnmasked_burden\tUnmasked_burden_95_lower\tUnmasked_burden_95_upper\n" \
    > "$OUTPUT"
while IFS= read -r SAMPLE; do

    # Skip empty lines and comments
    [[ -z "$SAMPLE" || "$SAMPLE" =~ ^# ]] && continue

    SBS_BURDEN_FILE="$WORKDIR/$SAMPLE/${SAMPLE}_sbs_burden.txt"

    awk -v sample="$SAMPLE" '
        BEGIN { OFS="\t" }

        {
            # Everything after the metric name is the value
            value = $NF
            values[NR] = value
        }

        END {
            print sample,
                  values[1],
                  values[2],
                  values[3],
                  values[4],
                  values[5],
                  values[6],
                  values[7],
                  values[8],
                  values[9],
                  values[10],
                  values[11],
                  values[12]


        }
    ' "$SBS_BURDEN_FILE" >> "$OUTPUT"

done < "$CONFIG"




### Generate summary table for indel burden across samples 
OUTPUT="${WORKDIR}/summary_indel_burden.txt"

# Header
printf "Sample\tIndel_burden\tIndel_burden_95_lower\tIndel_burden_95_upper\tIndel_number\tUnmasked_Indel_burden\tUnmasked_Indel_burden_95_lower\tUnmasked_Indel_burden_95_upper\n" \
    > "$OUTPUT"
while IFS= read -r SAMPLE; do

    # Skip empty lines and comments
    [[ -z "$SAMPLE" || "$SAMPLE" =~ ^# ]] && continue

    INDEL_BURDEN_FILE="$WORKDIR/$SAMPLE/${SAMPLE}_indel_burden.txt"

    awk -v sample="$SAMPLE" '
        BEGIN { OFS="\t" }

        {
            # Everything after the metric name is the value
            value = $NF
            values[NR] = value
        }

        END {
            print sample,
                  values[1],
                  values[2],
                  values[3],
                  values[4],
                  values[5],
                  values[6],
                  values[7]


        }
    ' "$INDEL_BURDEN_FILE" >> "$OUTPUT"

done < "$CONFIG"
















### Generate summary table
while read SAMPLE
do
    BURDEN=$(awk 'FNR == 5' "$WORKDIR/$SAMPLE/${SAMPLE}_sbs_burden.txt")
    SNVs=$(grep -v "#" "$WORKDIR/$SAMPLE/${SAMPLE}_snv.vcf" | grep "PASS" | wc -l)
    INDELs=$(grep -v "#" "$WORKDIR/$SAMPLE/${SAMPLE}_indel.vcf" | grep "PASS" | wc -l)
    DUPLICATION=$(awk 'FNR == 9' "$WORKDIR/$SAMPLE/${SAMPLE}_stats.txt")
    EFFICIENCY=$(awk 'FNR == 10' "$WORKDIR/$SAMPLE/${SAMPLE}_stats.txt")

    echo -e "$SAMPLE\t$BURDEN\t$SNVs\t$INDELs\t$DUPLICATION\t$EFFICIENCY"
done < $CONFIG


### Show specific result
while read SAMPLE
do
    BURDEN=$(awk 'FNR == 5' "$WORKDIR/$SAMPLE/${SAMPLE}_sbs_burden.txt")
    SNVs=$(grep -v "#" "$WORKDIR/$SAMPLE/${SAMPLE}_snv.vcf" | grep "PASS" | wc -l)
    INDELs=$(grep -v "#" "$WORKDIR/$SAMPLE/${SAMPLE}_indel.vcf" | grep "PASS" | wc -l)
    DUPLICATION=$(awk 'FNR == 9' "$WORKDIR/$SAMPLE/${SAMPLE}_stats.txt")
    EFFICIENCY=$(awk 'FNR == 10' "$WORKDIR/$SAMPLE/${SAMPLE}_stats.txt")
done < $CONFIG

for SAMP