import os
from pathlib import Path

rule warm_start_bins:
    input:
        config.get("genomes", [])
    output:
        bin_dir = directory("iter_{iteration}/binning/"),
        bininfo = "iter_{iteration}/binning/iter{iteration}_bins.BinInfo.txt"
    run:
        iteration = wildcards.iteration
        # Verify it's iterating 1 and genomes were provided
        if iteration == str(config['iteration']) and config.get("genomes"):
            os.makedirs(output.bin_dir, exist_ok=True)
            input_dir = config["genomes"]
            for f in os.listdir(input_dir):
                if f.endswith((".fa", ".fasta", ".fna")):
                    source = Path(input_dir) / f
                    dest = Path(output.bin_dir) / f
                    if dest.exists():
                        dest.unlink()
                    dest.symlink_to(source.resolve())
            
            # Create BinInfo identifying warm start
            with open(output.bininfo, "w") as f:
                f.write(f"Warm start MAGs imported from {input_dir}\n")
