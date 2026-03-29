import subprocess
import sys
from pathlib import Path
import os

def run_workflow(forward, reverse, output, threads, genomes=None, max_iterations=5, checkm_db_path="/fs/ess/PDS0325/bioinformatic_tools/aviary-0.12.0/checkm2/uniref100.KO.1.dmnd"):
    # Resolve absolute paths for inputs
    forward = str(Path(forward).resolve())
    reverse = str(Path(reverse).resolve())
    output = str(Path(output).resolve())
    checkm_db_path = str(Path(checkm_db_path).resolve())
    if genomes:
        genomes = str(Path(genomes).resolve())

    print("Launching Snakemake workflow...")

    last_iter_mags = 0

    for current_iter in range(1, max_iterations + 1):
        print(f"Running iteration {current_iter}...")

        if current_iter != 1:
            prev_iter_dir = Path(output) / f"iter_{current_iter-1}"
            forward = str((prev_iter_dir / "bin_mapping/reads.1.fq.gz").resolve())
            reverse = str((prev_iter_dir / "bin_mapping/reads.2.fq.gz").resolve())
            
            # Check if directory exists before listing
            binning_filtered = prev_iter_dir / "binning_filtered"
            if binning_filtered.exists():
                last_iter_mags = len([f for f in os.listdir(binning_filtered) if f.endswith(".fa")])
            else:
                last_iter_mags = 0
        
        # Build the command
        cmd = [
            "snakemake",
            "--use-conda",
            "--conda-frontend", "mamba",
            "--snakefile", str(Path(__file__).parent / "workflow/Snakefile"),
            "--cores", str(threads),
            "--directory", output,
            "--config",
            f"forward={forward}", 
            f"reverse={reverse}",
            f"threads={threads}",
            f"iteration={current_iter}",
            f"checkm_db_path={checkm_db_path}"
        ]

        if genomes is not None:
            cmd.extend([
            f"genomes={genomes}"
        ])

        

        try:
            if (last_iter_mags == 0 and current_iter > 1):
                print("No MAGs recovered in the previous iteration.\nIterations complete.")
                break
            else:
                print(f"\nRunning command:\n {' '.join(cmd)}")
                subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Snakemake pipeline failed with error: {e}", file=sys.stderr)
