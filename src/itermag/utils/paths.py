from pathlib import Path


def iter_dir(output, iteration):
    return Path(output) / f"iter_{iteration}"


def reads_dir(output, iteration):
    return iter_dir(output, iteration) / "reads"


def assembly_dir(output, iteration):
    return iter_dir(output, iteration) / "assembly"


def assembly_mapping_dir(output, iteration):
    return iter_dir(output, iteration) / "assembly_mapping"


def binning_dir(output, iteration):
    return iter_dir(output, iteration) / "binning"


def checkm2_dir(output, iteration):
    return iter_dir(output, iteration) / "checkm2"


def binning_filtered_dir(output, iteration):
    return iter_dir(output, iteration) / "binning_filtered"


def bin_mapping_dir(output, iteration):
    return iter_dir(output, iteration) / "bin_mapping"


def reads(output, iteration, mate):
    return bin_mapping_dir(output, iteration) / f"reads.{mate}.fq.gz"


def bininfo(output, iteration):
    return binning_dir(output, iteration) / f"iter{iteration}_bins.BinInfo.txt"


def quality_report(output, iteration):
    return checkm2_dir(output, iteration) / "quality_report.tsv"


def filtered_report(output, iteration):
    return binning_filtered_dir(output, iteration) / "quality_report.filtered.tsv"


def magnitude_count(output, iteration):
    return sum(1 for f in binning_filtered_dir(output, iteration).glob("*.fa"))


def summary_file(output):
    return Path(output) / "summary.tsv"