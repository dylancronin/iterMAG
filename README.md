# iterMAG

**iterMAG** is a Snakemake-based bioinformatic pipeline for **Iterative Metagenome-Assembled Genome (MAG) Assembly**. 

**WARNING**: This pipeline is currently under development, use with caution.

It improves MAG recovery by iteratively assembling reads, binning contigs, and then re-assembling only the reads that did not map to the bins in the previous round. This allows for the discovery of lower-abundance organisms that might be masked by dominant species in a traditional single-pass assembly.

## Workflow Overview

The pipeline currently runs for **5 iterations** by default:

1.  **Assembly**: Assembly using **MEGAHIT** (`meta-large` preset) for short reads or **myloasm** (`--hifi`) for PacBio HiFi.
2.  **Mapping**: Read mapping and coverage calculation using **CoverM**.
3.  **Binning**: Metagenomic binning using **MetaBAT2**.
4.  **Extraction**: 
    *   Map original reads to the generated bins.
    *   Filter out mapped reads to extract "unmapped" reads.
    *   Convert unmapped reads back to FASTQ format.
5.  **Iteration**: Repeat steps 1-3 using the extracted unmapped reads as input for the next round.

## Installation

### Prerequisites
The following tools must be in your PATH:
- [MEGAHIT](https://github.com/voutcn/megahit) (short-read mode only)
- [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/)
- [CoverM](https://github.com/wwood/CoverM)
- [Samtools](http://www.htslib.org/)
- [Snakemake](https://snakemake.readthedocs.io/)
- Conda and a separate existing myloasm environment (HiFi mode only)
- minimap2, for CoverM mapping

### Setup
Clone the repository and create the Conda environment:
```bash
git clone https://github.com/dylancronin/iterMAG.git
cd iterMAG
conda env create -f environment.yml
conda activate itermag-0.0.0
```

## Usage

You can run the pipeline using the `itermag` command-line interface.

```bash
itermag run -1 <forward_reads.fq.gz> -2 <reverse_reads.fq.gz> -o <output_dir> -t <threads>
```

### Arguments:
- `-1`, `--forward`: Path to forward reads (FASTQ GZ).
- `-2`, `--reverse`: Path to reverse reads (FASTQ GZ).
- `-o`, `--output`: Output directory where results and iterations will be stored.
- `-t`, `--threads`: Number of CPU threads to use.
- `-g`, `--genomes`: (Optional) Path to a directory containing existing genomes/bins.
- `-it`, `--max_iterations`: Number of iterations to run (default: 5).

## Directory Structure
The output directory will contain subdirectories for each iteration:
```text
output_dir/
├── iter_1/
│   ├── assembly/          # MEGAHIT assembly results
│   ├── assembly_mapping/  # Coverage files
│   ├── binning/           # MetaBAT2 bins
│   └── bin_mapping/       # Unmapped reads for Round 2
└── iter_2/
    ├── assembly/
    ├── ...
```


## PacBio HiFi

Keep the existing short-read command above, or select HiFi explicitly:

```bash
itermag run --read-type hifi \
    --long-reads sample.hifi.fastq.gz \
    --output hifi_results --threads 40 --max_iterations 5 \
    --myloasm-env /fs/ess/PDS0325/bioinformatic_tools/miniforge3/envs/myloasm-0.7.0
```

HiFi mode requires `--long-reads` and rejects `-1/-2`. Short-read mode is the default
and requires both mates. FASTQ inputs can be plain or gzipped; plain files are
compressed into each iteration's reads directory, while gzipped inputs are linked.
Inputs must use four-line FASTQ records.

The existing myloasm 0.7.0 environment is used through
`conda run --no-capture-output --prefix ... myloasm ... --hifi`.
Its default prefix is the path shown above; override `--myloasm-env` on other systems.
Do not add the separate environment's bin directory to the pipeline's PATH.
HiFi execution does not require MEGAHIT. Verify Conda, Snakemake, CoverM,
Samtools, minimap2, and MetaBAT2 are on PATH; the existing CheckM2 environment
and database configuration also apply.

Each HiFi iteration retains myloasm outputs in `assembly/myloasm/` and copies
`assembly_primary.fa` to the shared `assembly/contigs.fa` path.
Both assembly coverage and MAG mapping use CoverM's
`--single ... --mapper minimap2-hifi`. Assembly coverage uses `--methods metabat`.
The thresholds remain **95% identity and 75% of the read aligned**, using per-read
options for HiFi. MAG extraction uses inverse filtering, so "residual" reads include
unmapped reads and reads whose primary alignment fails either threshold.
Secondary and supplementary alignments are excluded from FASTQ extraction.

Residual HiFi reads are written to `iter_N/bin_mapping/reads.fq.gz` and become
the next iteration's input. No mate file is created. Warm starts use `--genomes`
as in short-read mode. If no MAGs pass CheckM2 filtering, current reads are
preserved as residual reads and the workflow stops. If all reads are depleted,
the residual output is a valid empty gzip file and the workflow stops.
MetaBAT2 and CheckM2 quality criteria remain unchanged.

`run_metadata.json` records the sequencing mode, assembler, and myloasm prefix.
Use separate output directories for short and HiFi runs; the workflow rejects
resuming with a conflicting mode. Existing output directories without metadata
are treated as short-read runs unless they contain HiFi residual files.

A SLURM example accepting input and output paths is available:
```bash
sbatch tests/submit_hifi_test.sh /path/sample.hifi.fastq.gz /path/hifi_results
# Optional warm-start genome directory is the third argument.
```

## Validation

Run the tests in the pipeline environment with the updated package installed:
```bash
python -m pip install --no-deps -e .
python -m unittest discover -s tests -v
```

The suite covers CLI compatibility, input validation, gzip handling, mode conflicts,
two-iteration handoff, and early stopping. Snakemake dry-run tests cover short,
HiFi, and HiFi warm-start DAGs; they skip when Snakemake is unavailable.
The small alignment fixture exercises the real CoverM filter and Samtools FASTQ
conversion, including threshold boundaries, unmapped reads, reverse-strand reads,
and secondary/supplementary records; it skips when those tools are unavailable.
No full assembly or scheduler job is launched by these tests.
