import subprocess
import sys
from pathlib import Path

from itermag.utils import paths
from itermag.utils.config import InputError, validate_run_args
from itermag.utils.logging import get_logger

DEFAULT_CHECKM_DB = "/fs/ess/PDS0325/bioinformatic_tools/aviary-0.12.0/checkm2/uniref100.KO.1.dmnd"


def run_workflow(
    forward,
    reverse,
    output,
    threads,
    genomes=None,
    max_iterations=5,
    checkm_db_path=DEFAULT_CHECKM_DB,
    conda_prefix=None,
):
    try:
        args = validate_run_args(
            forward=forward,
            reverse=reverse,
            output=output,
            threads=threads,
            genomes=genomes,
            checkm_db_path=checkm_db_path,
        )
    except InputError as e:
        print(f"[iterMAG] ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    forward = args["forward"]
    reverse = args["reverse"]
    output = args["output"]
    threads = args["threads"]
    genomes = args["genomes"]
    checkm_db_path = args["checkm_db_path"] or DEFAULT_CHECKM_DB

    out_root = Path(output)
    out_root.mkdir(parents=True, exist_ok=True)

    if conda_prefix:
        conda_prefix = str(Path(conda_prefix).resolve())
        Path(conda_prefix).mkdir(parents=True, exist_ok=True)

    if not Path(checkm_db_path).exists():
        print(
            f"[iterMAG] ERROR: CheckM2 database not found: {checkm_db_path}",
            file=sys.stderr,
        )
        sys.exit(1)

    logger = get_logger(log_file=out_root / "pipeline.log")
    logger.info(
        "Launching iterMAG: forward=%s reverse=%s output=%s threads=%s iterations=%s",
        forward,
        reverse,
        output,
        threads,
        max_iterations,
    )
    if genomes:
        logger.info("Warm start genomes directory: %s", genomes)

    snakefile = Path(__file__).parent / "workflow/Snakefile"
    if not snakefile.exists():
        logger.error("Snakefile not found: %s", snakefile)
        sys.exit(1)

    last_iter_mags = 0
    completed_iterations = 0

    for current_iter in range(1, int(max_iterations) + 1):
        logger.info("Running iteration %s.", current_iter)

        if current_iter != 1:
            prev_iter_dir = paths.iter_dir(output, current_iter - 1)
            forward = str(paths.reads(output, current_iter - 1, 1).resolve())
            reverse = str(paths.reads(output, current_iter - 1, 2).resolve())

            last_iter_mags = paths.magnitude_count(output, current_iter - 1)

            if not Path(forward).exists() or not Path(reverse).exists():
                logger.warning(
                    "Previous iteration produced no unmapped reads; iterations complete."
                )
                break

            if last_iter_mags == 0:
                logger.warning(
                    "No MAGs recovered in iteration %s; iterations complete.",
                    current_iter - 1,
                )
                break

        cmd = [
            "snakemake",
            "--use-conda",
            "--conda-frontend",
            "mamba",
            "--snakefile",
            str(snakefile),
            "--cores",
            str(threads),
            "--directory",
            output,
            "--rerun-incomplete",
            "--latency-wait",
            "30",
            "--restart-times",
            "1",
        ]
        if conda_prefix:
            cmd.extend(["--conda-prefix", conda_prefix])
        cmd.extend([
            "--config",
            f"forward={forward}",
            f"reverse={reverse}",
            f"threads={threads}",
            f"iteration={current_iter}",
            f"checkm_db_path={checkm_db_path}",
        ])
        if genomes is not None:
            cmd.extend([f"genomes={genomes}"])

        logger.info("Running command:\n  %s", " ".join(cmd))
        try:
            result = subprocess.run(cmd, check=False)
        except FileNotFoundError:
            logger.error("Could not launch snakemake. Is it installed and on PATH?")
            sys.exit(1)

        if result.returncode != 0:
            logger.error(
                "Snakemake failed during iteration %s (exit code %s). "
                "See logs under %s for details.",
                current_iter,
                result.returncode,
                out_root / ".snakemake/log",
            )
            sys.exit(result.returncode)

        reads_1 = paths.reads(output, current_iter, 1)
        reads_2 = paths.reads(output, current_iter, 2)
        if not reads_1.exists() or not reads_2.exists():
            logger.error(
                "Iteration %s finished but unmapped reads were not produced (%s).",
                current_iter,
                reads_1,
            )
            sys.exit(1)

        n_mags = paths.magnitude_count(output, current_iter)
        n_reads = _count_reads(reads_1)
        logger.info(
            "Iteration %s complete: %s high-quality MAG(s), %s unmapped read(s) remain.",
            current_iter,
            n_mags,
            n_reads,
        )
        completed_iterations = current_iter

        if n_mags == 0:
            logger.warning(
                "No new MAGs recovered in iteration %s; stopping early.",
                current_iter,
            )
            break

        if n_reads == 0:
            logger.warning("No unmapped reads remain; stopping early.")
            break

    _write_summary(output, completed_iterations)
    logger.info("iterMAG finished successfully.")


def _count_reads(fastq_gz):
    try:
        p1 = subprocess.Popen(
            ["gzip", "-cd", str(fastq_gz)],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        p2 = subprocess.Popen(
            ["wc", "-l"],
            stdin=p1.stdout,
            stdout=subprocess.PIPE,
            text=True,
        )
        p1.stdout.close()
        out, _ = p2.communicate()
        p1.wait()
        lines = int(out.split()[0])
    except Exception:
        return 0
    return lines // 4 if lines > 0 else 0


def _write_summary(output, completed_iterations):
    script = Path(__file__).parent / "workflow/scripts/summarize_bins.py"
    cmd = [
        sys.executable,
        str(script),
        str(output),
        "--completed",
        str(completed_iterations),
    ]
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        print(f"[iterMAG] WARNING: could not write summary: {e}", file=sys.stderr)