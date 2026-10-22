"""Small real CoverM/Samtools fixture; skipped if these binaries are unavailable."""

import gzip
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import itermag


@unittest.skipUnless(shutil.which("coverm") and shutil.which("samtools"),
                     "requires CoverM and Samtools")
class MappingFixtureTests(unittest.TestCase):
    def test_inverse_filter_boundaries_and_fastq_preservation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sequence = "A" * 99 + "C"
            quality = "I" * 99 + "J"
            records = [
                ("passing", 0, "100M", 0),
                ("identity_boundary", 0, "100M", 5),
                ("low_identity", 0, "100M", 6),
                ("aligned_boundary", 0, "75M25S", 0),
                ("partial", 0, "74M26S", 0),
                ("reverse", 16, "100M", 6),
                ("low_identity", 256, "100M", 0),
                ("partial", 2048, "26S74M", 0),
            ]
            sam = root / "input.sam"
            with sam.open("w") as handle:
                handle.write("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:1000\n")
                for name, flag, cigar, nm in records:
                    seq = "G" + "T" * 99 if flag == 16 else sequence
                    qual = quality[::-1] if flag == 16 else quality
                    handle.write(f"{name}\t{flag}\tref\t1\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tNM:i:{nm}\n")
                handle.write(f"unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t{sequence}\t{quality}\n")
            bam = root / "sorted.bam"
            subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
            filtered = root / "filtered.bam"
            subprocess.run(["coverm", "filter", "-b", str(bam), "-o", str(filtered),
                            "--inverse", "--exclude-supplementary",
                            "--min-read-percent-identity", "95",
                            "--min-read-aligned-percent", "75"], check=True)
            result = subprocess.run(["samtools", "fastq", "-F", "2304", "-n", str(filtered)],
                                    check=True, text=True, capture_output=True)
            lines = result.stdout.splitlines()
            observed = {}
            for index in range(0, len(lines), 4):
                name = lines[index][1:]
                self.assertNotIn(name, observed)
                observed[name] = (lines[index + 1], lines[index + 3])
            self.assertEqual(set(observed), {"low_identity", "partial", "reverse", "unmapped"})
            for value in observed.values():
                self.assertEqual(value, (sequence, quality))


@unittest.skipUnless(shutil.which("snakemake") and shutil.which("samtools"),
                     "requires Snakemake and Samtools")
class NoBinsExtractionTests(unittest.TestCase):
    def test_no_bins_preserves_plain_input_as_valid_gzip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reads = root / "input.fastq"
            original = "@read1\nACGTT\n+\nIJKLM\n"
            reads.write_text(original)
            report = root / "iter_1/binning_filtered/quality_report.filtered.tsv"
            report.parent.mkdir(parents=True)
            report.write_text("Name\tCompleteness\tContamination\n")
            include = Path(itermag.__file__).parent / "workflow/rules/hifi.smk"
            snakefile = root / "Snakefile"
            snakefile.write_text('warm_start = True\ninclude: ' + repr(str(include)) + '\n')
            residual = root / "iter_1/bin_mapping/reads.fq.gz"
            subprocess.run(["snakemake", "--snakefile", str(snakefile),
                            "--directory", str(root), "--cores", "1",
                            "iter_1/bin_mapping/reads.fq.gz",
                            "--config", f"long_reads={reads}", "threads=1"],
                           check=True)
            with gzip.open(residual, "rt") as handle:
                self.assertEqual(handle.read(), original)
            self.assertFalse(residual.with_name("reads.2.fq.gz").exists())
            subprocess.run(["samtools", "quickcheck", "-u",
                            str(residual.with_name("coverm.filtered.bam"))], check=True)


if __name__ == "__main__":
    unittest.main()
