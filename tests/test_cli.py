import contextlib
import gzip
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

from itermag.cli import main
from itermag.utils.config import InputError, sanity_check_environment, validate_run_args
from itermag.utils.fastq import count_reads, stage_gzipped_reads
from itermag.workflow import _check_run_metadata, run_workflow


FASTQ = "@read1\nACGT\n+\nIJKL\n"


class HiFiInterfaceTests(unittest.TestCase):
    def test_cli_modes_and_forwarding(self):
        cases = [
            (["run", "-1", "r1", "-2", "r2", "-o", "out"], "short"),
            (["run", "--read-type", "hifi", "--long-reads", "reads", "-o", "out"], "hifi"),
        ]
        for argv, mode in cases:
            with self.subTest(mode=mode), patch("sys.argv", ["itermag"] + argv), \
                    patch("itermag.cli.run_workflow") as launch:
                main()
                self.assertEqual(launch.call_args.kwargs["read_type"], mode)

    def test_reject_missing_and_mixed_inputs(self):
        cases = [
            ["run", "-o", "out"],
            ["run", "-1", "r1", "-o", "out"],
            ["run", "--long-reads", "reads", "-o", "out"],
            ["run", "--read-type", "hifi", "-o", "out"],
            ["run", "--read-type", "hifi", "--long-reads", "reads", "-1", "r1", "-o", "out"],
        ]
        for argv in cases:
            with self.subTest(argv=argv), patch("sys.argv", ["itermag"] + argv), \
                    contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main()

    def test_plain_and_gzip_staging_counting_and_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for compressed in (False, True):
                source = root / ("reads.fastq.gz" if compressed else "reads.fastq")
                opener = gzip.open if compressed else open
                with opener(source, "wt") as out:
                    out.write(FASTQ)
                validated = validate_run_args(None, None, root / "out", 2,
                                              read_type="hifi", long_reads=source)
                self.assertEqual(validated["long_reads"], str(source.resolve()))
                dest = root / "staged.fq.gz"
                stage_gzipped_reads(source, dest)
                self.assertEqual(count_reads(dest), 1)
                with gzip.open(dest, "rt") as handle:
                    self.assertEqual(handle.read(), FASTQ)
                self.assertEqual(dest.is_symlink(), compressed)
                dest.unlink()
            empty = root / "empty.gz"
            with gzip.open(empty, "wt"):
                pass
            self.assertEqual(count_reads(empty), 0)
            malformed = root / "bad.fastq"
            malformed.write_text("@r\nACGT\n+\nIII\n")
            with self.assertRaises(ValueError):
                count_reads(malformed)

    def test_missing_environment_and_no_megahit_requirement(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp)
            with self.assertRaises(InputError):
                sanity_check_environment("hifi", prefix / "missing")
            (prefix / "conda-meta").mkdir()
            (prefix / "bin").mkdir()
            executable = prefix / "bin/myloasm"
            executable.write_text("#!/bin/sh\nexit 0\n")
            executable.chmod(0o755)
            checked = []
            def which(tool):
                checked.append(tool)
                return None if tool == "megahit" else "/bin/" + tool
            with patch("itermag.utils.config._which", side_effect=which):
                sanity_check_environment("hifi", prefix)
            self.assertNotIn("megahit", checked)

    def test_resume_mode_guard_and_legacy_short_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            _check_run_metadata(output, "hifi")
            (output / "run_metadata.json").write_text(json.dumps({"read_type": "hifi"}))
            _check_run_metadata(output, "hifi")
            with self.assertRaises(InputError):
                _check_run_metadata(output, "short")
            (output / "run_metadata.json").unlink()
            (output / "iter_1").mkdir()
            _check_run_metadata(output, "short")
            with self.assertRaises(InputError):
                _check_run_metadata(output, "hifi")

    def test_driver_two_rounds_and_stopping_conditions(self):
        for mode, first_mags, first_reads, expected in [
            ("hifi", 1, 1, 2), ("hifi", 0, 1, 1),
            ("hifi", 1, 0, 1), ("short", 1, 1, 2),
        ]:
            with self.subTest(mode=mode, mags=first_mags, reads=first_reads), \
                    tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                source = root / "input.fastq"
                source.write_text(FASTQ)
                db = root / "db.dmnd"
                db.touch()
                output = root / "output"
                commands = []

                def launch(cmd, **kwargs):
                    commands.append(cmd)
                    iteration = len(commands)
                    mapping = output / f"iter_{iteration}/bin_mapping"
                    bins = output / f"iter_{iteration}/binning_filtered"
                    mapping.mkdir(parents=True)
                    bins.mkdir()
                    if iteration == 1 and first_mags:
                        (bins / "bin.fa").write_text(">c\nACGT\n")
                    names = ["reads.fq.gz"] if mode == "hifi" else ["reads.1.fq.gz", "reads.2.fq.gz"]
                    for name in names:
                        with gzip.open(mapping / name, "wt") as handle:
                            handle.write(FASTQ if iteration == 1 and first_reads else "")
                    return SimpleNamespace(returncode=0)

                inputs = {"long_reads": source} if mode == "hifi" else {"forward": source, "reverse": source}
                with patch("itermag.workflow.sanity_check_environment"), \
                        patch("itermag.workflow.subprocess.run", side_effect=launch), \
                        patch("itermag.workflow._write_summary") as summary:
                    run_workflow(output=output, threads=2, max_iterations=3,
                                 read_type=mode, checkm_db_path=db, **inputs)
                self.assertEqual(len(commands), expected)
                summary.assert_called_once_with(str(output.resolve()), expected)
                self.assertEqual(json.loads((output / "run_metadata.json").read_text())["read_type"], mode)
                if expected == 2:
                    field = "long_reads" if mode == "hifi" else "forward"
                    name = "reads.fq.gz" if mode == "hifi" else "reads.1.fq.gz"
                    self.assertIn(f"{field}={output}/iter_1/bin_mapping/{name}", commands[1])


if __name__ == "__main__":
    unittest.main()
