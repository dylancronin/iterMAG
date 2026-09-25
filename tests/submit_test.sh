#!/bin/bash
#SBATCH --account=PDS0325
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --job-name=itermag-test
#SBATCH --output=itermag-test-%j.out

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

source /fs/ess/PDS0325/bioinformatic_tools/miniforge3/bin/activate iterMAG-test

# Preflight: fail fast before the 10h run
for tool in itermag snakemake megahit metabat2 coverm; do
    command -v "$tool" >/dev/null || { echo "Required tool not found: $tool" >&2; exit 1; }
done
[ -f /fs/ess/PDS0325/bioinformatic_tools/aviary-0.12.0/checkm2/uniref100.KO.1.dmnd ] || { echo "CheckM2 DB missing" >&2; exit 1; }
[ -d test_run/.snakemake/conda ] || { echo "Prebuilt conda prefix missing (CheckM2 env would rebuild)" >&2; exit 1; }

itermag run \
    -1 data/full-reads/2022_IncS_FSU.v40.S15_30d_r4.1.fq.gz \
    -2 data/full-reads/2022_IncS_FSU.v40.S15_30d_r4.2.fq.gz \
    -o full_test_2 \
    -t 40 \
    -it 5 \
    -g /fs/ess/PDS0325/bog-incubations/20251221_MAGs-and-Assemblies-renaming/mags-test/05_galah/galah-95-reps/ \
    --conda-prefix test_run/.snakemake/conda