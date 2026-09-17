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


#