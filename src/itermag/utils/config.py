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


def validate_read_inputs(forward, reverse, read_type="short", long_reads=None):
    if read_type not in ("short", "hifi"):
        raise InputError("Read type must be short or hifi.")
    if read_type == "hifi":
        if not long_reads or forward is not None or reverse is not None:
            raise InputError("HiFi mode requires --long-reads and does not accept -1/-2.")
    elif not forward or not reverse or long_reads is not None:
        raise InputError("Short mode requires -1 and -2 and does not accept --long-reads.")


def validate_run_args(forward, reverse, output, threads, genomes=None,
                      checkm_db_path=None, read_type="short", long_reads=None):
    validate_read_inputs(forward, reverse, read_type, long_reads)
    if read_type == "hifi":
        long_reads = _require_file(long_reads, "HiFi reads")
    else:
        forward = _require_file(forward, "Forward reads")
        reverse = _require_file(reverse, "Reverse reads")

    if output is None:
        raise InputError("An output directory is required.")
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
        "read_type": read_type,
        "long_reads": long_reads,
        "forward": forward,
        "reverse": reverse,
        "output": str(out.resolve()),
        "threads": threads,
        "genomes": genomes,
        "checkm_db_path": checkm_db_path,
    }


def sanity_check_environment(read_type="short", myloasm_env=None):
    missing = []
    tools = ["snakemake", "metabat2", "coverm", "samtools", "minimap2"]
    if read_type == "hifi":
        tools.append("conda")
        if not myloasm_env:
            raise InputError("A myloasm conda environment is required for HiFi mode.")
        prefix = Path(myloasm_env).expanduser().resolve()
        executable = prefix / "bin" / "myloasm"
        if not (prefix / "conda-meta").is_dir():
            raise InputError(f"Myloasm environment is not a conda prefix: {prefix}")
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise InputError(f"Myloasm executable not found or not executable: {executable}")
    else:
        tools.append("megahit")
    for tool in tools:
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
