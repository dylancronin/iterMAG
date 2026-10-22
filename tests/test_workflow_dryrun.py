"""Build actual Snakemake DAGs without starting assemblers or CheckM2."""

from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import itermag


SNAKEFILE = Path(itermag.__file__).parent / "workflow/Snakefile"


@unittest.skipUnless(shutil.which("snakemake"), "requires Snakemake")
class DryRunTests(unittest.TestCase):
    def dry_run(self, root, mode, warm=False):
        reads = root / "input.fastq"
        reads.write_text("@r\nACGT\n+\nIIII\n")
        out = root / ("warm" if warm else mode)
        out.mkdir()
        config = ["iteration=1", "threads=2", f"read_type={mode}", "checkm_db_path=unused.dmnd"]
        if mode == "hifi":
            config.extend([f"long_reads={reads}", "myloasm_env=/existing/myloasm-env"])
        else:
            config.extend([f"forward={reads}", f"reverse={reads}"])
        if warm:
            genomes = root / "genomes"
            genomes.mkdir()
            (genomes / "bin.fa").write_text(">contig\nACGT\n")
            config.append(f"genomes={genomes}")
        result = subprocess.run(["snakemake", "--snakefile", str(SNAKEFILE),
                                 "--directory", str(out), "--cores", "2",
                                 "--dry-run", "--printshellcmds", "--config"] + config,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout + result.stderr

    def test_hifi_both_mappers_and_unpaired_filter(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = self.dry_run(Path(tmp), "hifi")
            self.assertIn("conda run --no-capture-output --prefix", output)
            self.assertIn("myloasm", output)
            self.assertIn("--hifi", output)
            self.assertEqual(output.count("--mapper minimap2-hifi"), 2)
            self.assertIn("--min-read-percent-identity 95", output)
            self.assertIn("--min-read-aligned-percent 75", output)
            self.assertNotIn("--proper-pairs-only", output)
            self.assertNotIn("--min-read-percent-identity-pair", output)
            self.assertNotIn("megahit -t", output)

    def test_hifi_warm_start_skips_assembly(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = self.dry_run(Path(tmp), "hifi", warm=True)
            self.assertNotIn("conda run --no-capture-output --prefix", output)
            self.assertNotIn("coverm contig", output)
            self.assertIn("coverm make", output)
            self.assertIn("--mapper minimap2-hifi", output)

    def test_short_read_commands_remain_paired(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = self.dry_run(Path(tmp), "short")
            self.assertIn("megahit -t", output)
            self.assertIn("--presets meta-large", output)
            self.assertIn("--proper-pairs-only", output)
            self.assertNotIn("--mapper minimap2-hifi", output)


if __name__ == "__main__":
    unittest.main()
