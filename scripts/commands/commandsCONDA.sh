#######################################
### 30/04/2026 - Martín Santamarina ###
#######################################

##########################################
#### Miniconda3 installation for CSD3 ####
##########################################

### Download miniconda3
cd /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/miniconda3

### Make Conda usable
# file: software/miniconda3/etc/profile.d/conda.sh already exists
source /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/miniconda3/etc/profile.d/conda.sh

### Configure Conda to store envs in your shared path
#nano /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/miniconda3/.condarc
#envs_dirs:
#  - /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/conda_envs
#
#pkgs_dirs:
#  - /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/conda_pkgs

mkdir -p /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/conda_envs
mkdir -p /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/conda_pkgs

### Check conda info
conda info



#####################################
#### Conda environments for CSD3 ####
#####################################

cd /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/conda_envs

### remove conda evnvironment
#conda remove -n dupcaller_v1.0.4_cb2f352 --all

#### dupcaller conda environments ####

mkdir -p /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/conda_envs/dupcaller
cd /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/conda_envs/dupcaller


### conda environment for dupcaller_cb2f352
conda create -p /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/conda_envs/dupcaller/dupcaller_cb2f352 python=3.12 bioconda::tabix
cd dupcaller_cb2f352/

conda activate /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/conda_envs/dupcaller/dupcaller_cb2f352
pip install git+https://github.com/AlexandrovLab/DupCaller.git@cb2f352

### conda environment for dupcaller_2984779
conda create -p /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/conda_envs/dupcaller/dupcaller_2984779 python=3.12 bioconda::tabix
cd dupcaller_2984779/

conda activate /home/ms3242/rds/rds-dog-de-novo-mXW0tnsCK5M/software/conda_envs/dupcaller/dupcaller_2984779
pip install git+https://github.com/AlexandrovLab/DupCaller.git@2984779