#!/bin/bash
#SBATCH --account=PDS0325
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --job-name=itermag-hifi
#SBATCH --output=itermag-hifi-%j.out

set -euo pipefail

# Submit from the repository root:
# sbatch tests/submit_hifi_test.sh /path/sample.hifi.fastq.gz /path/output [genome-directory]
cd "$SLURM_SUBMIT_DIR"
if (( $# < 2 || $# > 3 )); then
    echo "Usage: sbatch tests/submit_hifi_test.sh READS OUTPUT [GENOMES]" >&2
    exit 2
fi
source /fs/ess/PDS0325/bioinformatic_tools/miniforge3/bin/activate iterMAG-test

extra=()
if (( $# == 3 )); then
    extra+=(--genomes "$3")
fi
itermag run --read-type hifi --long-reads "$1" --output "$2" \
    --threads "$SLURM_CPUS_PER_TASK" --max_iterations 5 \
    --myloasm-env /fs/ess/PDS0325/bioinformatic_tools/miniforge3/envs/myloasm-0.7.0 \
    "${extra[@]}"
