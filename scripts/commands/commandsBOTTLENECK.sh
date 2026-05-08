##################
### 06/05/2026 ###
##################

### Run dupcaller on a array of samples
sbatch --array=1-$(wc -l < ../config/samples.txt) slurm_dupcaller_pipeline
