import json
import subprocess
import sys
from pathlib import Path

from itermag.utils import paths
from itermag.utils.config import InputError, validate_run_args, sanity_check_environment
from itermag.utils.fastq import count_reads
from itermag.utils.logging import get_logger

DEFAULT_CHECKM_DB = "/fs/ess/PDS0325/bioinformatic_tools/aviary-0.12.0/checkm2/uniref100.KO.1.dmnd"

DEFAULT_MYLOASM_ENV = "/fs/ess/PDS0325/bioinformatic_tools/miniforge3/envs/myloasm-0.7.0"


def run_workflow(
    forward=None,
    reverse=None,
    output=None,
    threads=1,
    genomes=None,
    max_iterations=5,
    checkm_db_path=DEFAULT_CHECKM_DB,
    conda_prefix=None,
    read_type="short",
    long_reads=None,
    myloasm_env=DEFAULT_MYLOASM_ENV,
):
    try:
        args = validate_run_args(
            forward=forward,
            reverse=reverse,
            output=output,
            threads=threads,
            genomes=genomes,
            checkm_db_path=checkm_db_path,
            read_type=read_type,
            long_reads=long_reads,
        )
        if int(max_iterations) < 1:
            raise InputError("Maximum iterations must be positive.")
        if read_type == "hifi":
            myloasm_env = str(Path(myloasm_env).expanduser().resolve())
        sanity_check_environment(read_type, myloasm_env)
        _check_run_metadata(Path(args["output"]), read_type)
        if read_type == "hifi" and count_reads(args["long_reads"]) == 0:
            raise InputError("HiFi input contains no reads.")
    except (InputError, OSError, ValueError, EOFError) as e:
        print(f"[iterMAG] ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    long_reads = args["long_reads"]
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

    metadata = {"read_type": read_type, "assembler": "myloasm" if read_type == "hifi" else "megahit",
                "myloasm_env": myloasm_env if read_type == "hifi" else None}
    (out_root / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    logger = get_logger(log_file=out_root / "pipeline.log")
    logger.info(
        "Launching iterMAG: mode=%s long_reads=%s myloasm_env=%s forward=%s reverse=%s output=%s threads=%s iterations=%s",
        read_type,
        long_reads,
        myloasm_env if read_type == "hifi" else None,
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
            if read_type == "hifi":
                long_reads = str(paths.reads(output, current_iter - 1).resolve())
            else:
                forward = str(paths.reads(output, current_iter - 1, 1).resolve())
                reverse = str(paths.reads(output, current_iter - 1, 2).resolve())

            last_iter_mags = paths.magnitude_count(output, current_iter - 1)

            required_reads = [long_reads] if read_type == "hifi" else [forward, reverse]
            if not all(Path(p).exists() for p in required_reads):
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
            f"read_type={read_type}",
            f"threads={threads}",
            f"iteration={current_iter}",
            f"checkm_db_path={checkm_db_path}",
        ])
        if read_type == "hifi":
            cmd.extend([f"long_reads={long_reads}", f"myloasm_env={myloasm_env}"])
        else:
            cmd.extend([f"forward={forward}", f"reverse={reverse}"])
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

        reads_1 = paths.reads(output, current_iter, None if read_type == "hifi" else 1)
        reads_2 = reads_1 if read_type == "hifi" else paths.reads(output, current_iter, 2)
        if not reads_1.exists() or not reads_2.exists():
            logger.error(
                "Iteration %s finished but unmapped reads were not produced (%s).",
                current_iter,
                reads_1,
            )
            sys.exit(1)

        n_mags = paths.magnitude_count(output, current_iter)
        try:
            n_reads = _count_reads(reads_1)
        except (OSError, ValueError, EOFError) as exc:
            logger.error("Cannot read residual FASTQ: %s", exc)
            sys.exit(1)
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


def _count_reads(fastq_path):
    return count_reads(fastq_path)


def _check_run_metadata(output, read_type):
    metadata_path = output / "run_metadata.json"
    if metadata_path.exists():
        try:
            previous_type = json.loads(metadata_path.read_text())["read_type"]
        except (OSError, ValueError, KeyError) as exc:
            raise InputError(f"Cannot read run metadata: {metadata_path}: {exc}")
    elif any(output.glob("iter_*/bin_mapping/reads.fq.gz")):
        previous_type = "hifi"
    elif any(output.glob("iter_*")):
        previous_type = "short"
    else:
        return
    if previous_type != read_type:
        raise InputError(f"Output directory contains a {previous_type} run; use a new directory for {read_type}.")


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
