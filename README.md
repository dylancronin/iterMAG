# iterMAG

**iterMAG** is a Snakemake-based bioinformatic pipeline for **Iterative Metagenome-Assembled Genome (MAG) Assembly**. 

It improves MAG recovery by iteratively assembling reads, binning contigs, and then re-assembling only the reads that did not map to the bins in the previous round. This allows for the discovery of lower-abundance organisms that might be masked by dominant species in a traditional single-pass assembly.

## Workflow Overview

The pipeline currently runs for **2 iterations** by default:

1.  **Assembly**: Initial assembly using **MEGAHIT** (`meta-large` preset).
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
- [MEGAHIT](https://github.com/voutcn/megahit)
- [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/)
- [CoverM](https://github.com/wwood/CoverM)
- [Samtools](http://www.htslib.org/)
- [Snakemake](https://snakemake.readthedocs.io/)

### Setup
Clone the repository and install the package:
```bash
git clone https://github.com/dylancronin/iterMAG.git
cd iterMAG
pip install .
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

