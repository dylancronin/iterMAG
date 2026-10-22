"""Read plain or gzipped four-line FASTQ without relying on its filename."""

import gzip
import shutil
from pathlib import Path


def open_fastq(path):
    with open(path, "rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if compressed else open(path, "rt")


def count_reads(path):
    with open_fastq(path) as handle:
        count = 0
        while True:
            header = handle.readline()
            if not header:
                return count
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if not header.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
                raise ValueError(f"Invalid four-line FASTQ record {count + 1}: {path}")
            count += 1


def has_reads(path):
    with open_fastq(path) as handle:
        return bool(handle.read(1))


def stage_gzipped_reads(source, destination):
    """Link gzip inputs; compress plain FASTQ for tools expecting .fq.gz."""
    source = Path(source).resolve()
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink() or destination.exists():
        destination.unlink()
    with source.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    if compressed:
        destination.symlink_to(source)
    else:
        with source.open("rb") as reader, gzip.open(destination, "wb") as writer:
            shutil.copyfileobj(reader, writer)
