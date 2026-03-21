from pathlib import Path

rule symlink_reads:
    input:
        # We grab the raw files handed to us by workflow.py
        r1 = config["forward"],
        r2 = config["reverse"]
    output:
        # We define the standardized output format for ANY iteration
        r1 = "iter_{iteration}/reads/reads.1.fq.gz",
        r2 = "iter_{iteration}/reads/reads.2.fq.gz"
    run:
        Path(output.r1).symlink_to(Path(input.r1).resolve())
        Path(output.r2).symlink_to(Path(input.r2).resolve())
        # Path().resolve() gets the absolute path
        # .symlink_to() creates the link