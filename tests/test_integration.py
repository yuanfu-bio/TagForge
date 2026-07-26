from __future__ import annotations

import gzip
import subprocess
import sys
import importlib.util
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tagforge.config import load_config
from tagforge.barcode_correct import correct_sample
from tagforge.extract import extract_sample
from tagforge.fastq import paired_fastq_batches
from tagforge.columnar import open_records, record_batches


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
        self.assertTrue((sample / "04_matrix/example.raw_molecule_pairs.parquet").is_file())
        self.assertTrue((sample / "04_matrix/example.raw_matrix_summary.json").is_file())
        self.assertTrue((sample / "07_report/example.report.xlsx").is_file())
        self.assertTrue((sample / "02_extracted/example.extraction_stats.tsv").is_file())
        quick_stats = sample / "07_report/example.quick_test.tsv"
        self.assertTrue(quick_stats.is_file())
        self.assertIn("sample\treads_examined\t9", quick_stats.read_text(encoding="utf-8"))
        self.assertTrue((root / "00_report/TagForge_batch_report.html").is_file())
        log = (sample / "00_logs/example.pipeline.log").read_text(encoding="utf-8")
        self.assertIn("backend=cutadapt-python-api", log)
        self.assertIn("workers=2", log)
        self.assertIn("tagforge\tversion=0.1.13", log)
        self.assertIn("correction_progress", log)
        self.assertIn("correction_summary", log)
        self.assertIn("correction_parallel_start", log)
        self.assertIn("dedup_summary", log)
        self.assertIn("dedup_progress", log)
        self.assertIn("workers=2", log)
        extracted_rows = list(open_records(sample / "02_extracted/example.extracted.parquet"))
        extracted = "\n".join("\t".join(row.values()) for row in extracted_rows)
        header = "\t".join(extracted_rows[0])
        self.assertEqual(
            header,
            "read_id\tCELL_segments\tFB_segments\tUMI_segments\tmethods\tstatus\tfailure_reason",
        )
        self.assertIn('read009', extracted)
        self.assertIn("FFF", extracted)
        self.assertNotIn("raw_barcode1", extracted)
        self.assertNotIn("barcode1_segments", header)
        self.assertNotIn("[", extracted)
        preview = sample / "02_extracted/example.extracted.preview.tsv"
        self.assertTrue(preview.is_file())
        self.assertEqual(len(preview.read_text(encoding="utf-8").splitlines()), 10)
        progress = sample / "00_logs/example.extraction_progress.tsv"
        self.assertIn("completed\t9\t100.00", progress.read_text(encoding="utf-8"))
        dedup_progress = sample / "00_logs/example.dedup_progress.tsv"
        self.assertTrue(dedup_progress.is_file())
        self.assertIn("completed\t9", dedup_progress.read_text(encoding="utf-8"))
        downsample_progress = sample / "00_logs/example.downsample_progress.tsv"
        self.assertTrue(downsample_progress.is_file())
        self.assertIn("completed", downsample_progress.read_text(encoding="utf-8"))
        self.assertTrue((sample / "06_downsample/example.downsample.html").is_file())
        stats = (sample / "02_extracted/example.extraction_stats.tsv").read_text(encoding="utf-8")
        self.assertIn("CELL\tCELL\tR1\tlinker_fixed\t9\t9\t8\t1", stats)
        valid_rows = list(open_records(sample / "05_detail/example.valid_reads.parquet"))
        self.assertEqual(list(valid_rows[0]), ["read_id", "CELL", "FB_sequence", "FB_name", "UMI", "correction_summary"])
        molecule_rows = list(open_records(sample / "05_detail/example.molecule_detail.parquet"))
        self.assertEqual(list(molecule_rows[0]), ["CELL", "FB_name", "corrected_umi", "reads_count", "raw_umi_count"])
        self.assertFalse((sample / "04_matrix/example.raw_count_matrix.tsv.gz").exists())

        # Commit one three-read gzip member, interrupt, then resume without
        # duplicating the header or completed records.
        with tempfile.TemporaryDirectory() as td:
            config = load_config(root / "configs/config.example.yaml")
            config.output_dir = Path(td) / "output"
            config.chunk_size = 3
            config.threads = 1
            real_batches = paired_fastq_batches

            def interrupted_batches(r1, r2, batch_size):
                yield next(iter(real_batches(r1, r2, batch_size)))
                raise RuntimeError("simulated interruption")

            with patch("tagforge.extract.paired_fastq_batches", side_effect=interrupted_batches):
                with self.assertRaisesRegex(RuntimeError, "simulated interruption"):
                    extract_sample(config, "example", resume=False)
            resume_file = config.output_dir / "example/01_checkpoint/example.extract.resume.json"
            self.assertTrue(resume_file.is_file())
            extract_sample(config, "example", resume=True)
            resumed_rows = list(open_records(config.output_dir / "example/02_extracted/example.extracted.parquet"))
            self.assertEqual(len(resumed_rows), 9)
            self.assertEqual(sum(row["read_id"] == "read001" for row in resumed_rows), 1)
            self.assertFalse(resume_file.exists())

            # Correction commits complete Parquet fragments for both outputs and
            # restores counters and output fragment boundaries after interruption.
            real_record_batches = record_batches

            def interrupted_correction_batches(path, batch_size):
                yield next(iter(real_record_batches(path, batch_size)))
                raise RuntimeError("simulated correction interruption")

            with patch(
                "tagforge.barcode_correct.record_batches",
                side_effect=interrupted_correction_batches,
            ):
                with self.assertRaisesRegex(RuntimeError, "simulated correction interruption"):
                    correct_sample(config, "example", resume=False)
            correction_resume = (
                config.output_dir / "example/01_checkpoint/example.correct.resume.json"
            )
            self.assertTrue(correction_resume.is_file())
            correct_sample(config, "example", resume=True)
            valid_rows = list(open_records(config.output_dir / "example/05_detail/example.valid_reads.parquet"))
            self.assertEqual(len(valid_rows), 9)
            self.assertEqual(sum(row["read_id"] == "read001" for row in valid_rows), 1)
            trace_rows = list(open_records(config.output_dir / "example/03_corrected/example.barcode_correction_trace.parquet"))
            self.assertEqual(len(trace_rows), 18)
            self.assertFalse(correction_resume.exists())
            correction_progress = (
                config.output_dir / "example/00_logs/example.correction_progress.tsv"
            ).read_text(encoding="utf-8")
            self.assertIn("completed\t9\t9\t9\t100.00", correction_progress)


if __name__ == "__main__":
    unittest.main()
