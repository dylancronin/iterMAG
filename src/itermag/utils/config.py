import os
from pathlib import Path


class InputError(Exception):
    pass


def _require_file(path, what):
    p = Path(path)
    if not p.exists():
        raise InputError(f"{what} does not exist: {p}")
    if not p.is_file():
        raise InputError(f"{what} is not a regular file: {p}")
    return str(p.resolve())


def _require_dir(path, what, allow_empty=False):
    p = Path(path)
    if not p.exists():
        raise InputError(f"{what} does not exist: {p}")
    if not p.is_dir():
        raise InputError(f"{what} is not a directory: {p}")
    if not allow_empty and not any(p.iterdir()):
        raise InputError(f"{what} is empty: {p}")
    return str(p.resolve())


def validate_run_args(forward, reverse, output, threads, genomes=None, checkm_db_path=None):
    forward = _require_file(forward, "Forward reads")
    reverse = _require_file(reverse, "Reverse reads")

    out = Path(output)
    if out.exists() and not out.is_dir():
        raise InputError(f"Output path exists and is not a directory: {out}")

    try:
        threads = int(threads)
    except (TypeError, ValueError):
        raise InputError(f"Threads must be a positive integer, got: {threads!r}")
    if threads < 1:
        raise InputError(f"Threads must be a positive integer, got: {threads!r}")

    if genomes is not None:
        genomes = _require_dir(genomes, "Genome directory")

    if checkm_db_path is not None and Path(checkm_db_path).exists():
        checkm_db_path = str(Path(checkm_db_path).resolve())
    else:
        checkm_db_path = None

    return {
        "forward": forward,
        "reverse": reverse,
        "output": str(out.resolve()),
        "threads": threads,
        "genomes": genomes,
        "checkm_db_path": checkm_db_path,
    }


def sanity_check_environment():
    missing = []
    for tool in ["snakemake", "megahit", "metabat2", "coverm"]:
        if not _which(tool):
            missing.append(tool)
    if missing:
        raise InputError(
            "Required tools not found on PATH: " + ", ".join(missing)
        )


def _which(tool):
    for base in os.environ.get("PATH", "").split(os.pathsep):
        candidate = Path(base) / tool
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return None