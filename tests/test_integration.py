from __future__ import annotations

import gzip
import subprocess
import sys
import importlib.util
import shutil
import tempfile
import unittest
from pathlib import Path


class IntegrationTests(unittest.TestCase):
    @unittest.skipUnless(
        importlib.util.find_spec("cutadapt") is not None
        and importlib.util.find_spec("umi_tools") is not None
        and shutil.which("cutadapt") is not None
        and shutil.which("umi_tools") is not None,
        "requires the TagForge Conda environment",
    )
    def test_example_pipeline(self):
        root = Path(__file__).parents[1]
        subprocess.run([sys.executable, str(root / "examples/generate_example_data.py")], check=True)
        result = subprocess.run([sys.executable, "-m", "tagforge", "validate-config", "--config", str(root / "configs/config.example.yaml")], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        result = subprocess.run([sys.executable, "-m", "tagforge", "quick-test", "--config", str(root / "configs/config.example.yaml"), "--reads", "5", "--threads", "2"], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        quick_stats = root / "02_output/example/07_report/example.quick_test.tsv"
        self.assertIn("sample\treads_examined\t5", quick_stats.read_text(encoding="utf-8"))
        sampled_ids = root / "02_output/example/07_report/example.quick_test.sampled_read_ids.txt"
        first_ids = sampled_ids.read_text(encoding="utf-8")
        self.assertEqual(len(first_ids.strip().splitlines()), 6)
        self.assertIn("sample\treads_scanned\t5", quick_stats.read_text(encoding="utf-8"))
        quick_stats.unlink(missing_ok=True)
        result = subprocess.run([sys.executable, "-m", "tagforge", "run", "--config", str(root / "configs/config.example.yaml"), "--threads", "2", "--overwrite"], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("quick-test: sampled=9 reads, scanned=9", result.stdout)
        sample = root / "02_output/example"
        self.assertTrue((sample / "04_matrix/example.raw_count_matrix.tsv.gz").is_file())
        self.assertTrue((sample / "07_report/example.report.xlsx").is_file())
        self.assertTrue((sample / "02_extracted/example.extraction_stats.tsv").is_file())
        quick_stats = sample / "07_report/example.quick_test.tsv"
        self.assertTrue(quick_stats.is_file())
        self.assertIn("sample\treads_examined\t9", quick_stats.read_text(encoding="utf-8"))
        self.assertTrue((root / "00_report/TagForge_batch_report.html").is_file())
        log = (sample / "00_logs/example.pipeline.log").read_text(encoding="utf-8")
        self.assertIn("backend=cutadapt-python-api", log)
        self.assertIn("workers=2", log)
        self.assertIn("tagforge\tversion=0.1.3", log)
        with gzip.open(sample / "02_extracted/example.extracted.tsv.gz", "rt", encoding="utf-8") as handle:
            extracted = handle.read()
        self.assertIn('read009', extracted)
        self.assertIn('""CELL"":""fixed""', extracted)
        stats = (sample / "02_extracted/example.extraction_stats.tsv").read_text(encoding="utf-8")
        self.assertIn("CELL\tbarcode1\tR1\tlinker_fixed\t9\t9\t8\t1", stats)


if __name__ == "__main__":
    unittest.main()
