#!/bin/bash
#SBATCH --account=PDS0325
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=itermag-test
#SBATCH --output=itermag-test-%j.out

# Ensure we are in the directory where the job was submitted
#cd $SLURM_SUBMIT_DIR

# Activate the correct conda environment
source /fs/ess/PDS0325/bioinformatic_tools/miniforge3/bin/activate iterMAG-test

# Run the pipeline
itermag run \
    -1 data/reads/2022_IncS_FSU.v72.S25_60d_r4.subsampled1mil.1.fq.gz \
    -2 data/reads/2022_IncS_FSU.v72.S25_60d_r4.subsampled1mil.2.fq.gz \
    -o test_run \
    -t 40
