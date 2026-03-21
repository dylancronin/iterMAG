import subprocess
import sys
from pathlib import Path

def run_workflow(forward, reverse, output, threads, genomes=None):
    # Resolve absolute paths for inputs
    forward = Path(forward).resolve()
    reverse = Path(reverse).resolve()
    if genomes:
        genomes = str(Path(genomes).resolve())

    print("Launching Snakemake workflow...")

    max_iterations = 1
    
    for current_iter in range(1, max_iterations + 1):
        print(f"Running iteration {current_iter}...")

        
        
        # Build the command
        cmd = [
            "snakemake",
            "--snakefile", str(Path(__file__).parent / "workflow/Snakefile"),
            "--cores", str(threads),
            "--directory", output,
            "--config",
            f"forward={forward}", 
            f"reverse={reverse}",
            f"threads={threads}",
            f"iteration={current_iter}"
        ]

        if genomes is not None:
            cmd.extend([
            f"genomes={genomes}"
        ])

        print(f"\nRunning command:\n {' '.join(cmd)}")

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Snakemake pipeline failed with error: {e}", file=sys.stderr)
