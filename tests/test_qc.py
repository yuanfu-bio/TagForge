from collections import Counter
from dataclasses import replace
import gzip
from itertools import product
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from tagforge.barcode_correct import (WhitelistCorrector, _build_correctors,
                                     evaluate_extracted_read, load_fb_annotation)
from tagforge.config import (ConfigError, CorrectionConfig, SampleConfig, SegmentConfig,
                            TagForgeConfig, load_config)
from tagforge.extract import _extract_record
from tagforge.fastq import paired_fastq
from tagforge.io_utils import open_tsv
from tagforge.quick_test import enforce_quick_test, quick_test_sample


class QCTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        disabled = CorrectionConfig(enabled=False)
        segments = [SegmentConfig(name, target, read, "fixed", 2, start=start,
                                  correction=disabled)
                    for name, target, read, start in (
                        ("PB1", "barcode1", "R1", 0), ("PB2", "barcode1", "R1", 2),
                        ("FB1", "barcode2", "R1", 0), ("FB2", "barcode2", "R2", 4),
                        ("UMI1", "umi", "R2", 0), ("UMI2", "umi", "R2", 2))]
        self.config = TagForgeConfig(self.root / "config.yaml", self.root,
                                    self.root / "output",
                                    [SampleConfig("sample", self.root / "R1.fq.gz",
                                                  self.root / "R2.fq.gz")],
                                    segments, None, barcode2_sequence_only=True, quick_test_read_details=True)

    def write_reads(self, records):
        for mate, path in ((1, self.config.samples[0].r1), (2, self.config.samples[0].r2)):
            text = "".join(f"@{name}/{mate}\n{row[mate]}\n+\n{'I' * len(row[mate])}\n"
                           for row in records for name in [row[0]])
            path.write_bytes(gzip.compress(text.encode(), mtime=0))

    def compare_formal_qc(self, records):
        self.write_reads(records)
        with patch("tagforge.quick_test.check_external_tools"):
            result = quick_test_sample(self.config, "sample", len(records))
        correctors = _build_correctors(self.config)
        annotation = load_fb_annotation(self.config)
        columns = {role: self.config.segment_column(role) for role in ("barcode1", "barcode2", "umi")}
        expected = []
        for record in paired_fastq(self.config.samples[0].r1, self.config.samples[0].r2):
            row, _ = _extract_record(record, self.config.segments, columns)
            summary = Counter()
            valid, _ = evaluate_extracted_read(row, self.config, correctors, annotation,
                {s.name: Counter() for s in self.config.segments if s.target != "umi"}, summary)
            reason = next((key[8:] for key, value in summary.items()
                           if key.startswith("failure_") and value), "")
            expected.append((record.read_id, valid is not None, reason))
        observed = [(row["read_id"], row["eligible"] == "true", row["rejection_reason"])
                    for row in open_tsv(result["read_qc"])]
        self.assertEqual(observed, expected)
        return result, observed

    def test_quick_test_read_details_can_be_disabled_without_changing_qc(self):
        result, _ = self.compare_formal_qc([("valid", "AACC", "TTAAGG"), ("invalid", "AA", "TTAAGG")])
        self.config.quick_test_read_details = False
        with patch("tagforge.quick_test.check_external_tools"):
            compact = quick_test_sample(self.config, "sample", 2)
        self.assertEqual(compact["valid_rate"], result["valid_rate"])
        self.assertIsNone(compact["read_qc"])
        self.assertIsNone(compact["read_ids"])
        self.assertFalse(result["read_qc"].exists())
        self.assertFalse(result["read_ids"].exists())
        self.assertTrue(compact["stats"].is_file())

    def test_partial_segments_and_invalid_umi_agree_per_read(self):
        records = [("valid", "AACC", "TTAAGG"),
                   ("partial_pb", "AA", "TTAAGG"),
                   ("partial_fb", "AACC", "TTAA"),
                   ("partial_umi", "AACC", "TT"),
                   ("umi_n", "AACC", "TNAGGG"),
                   ("umi_non_dna", "AACC", "TRAGGG")]
        result, observed = self.compare_formal_qc(records)
        self.assertEqual([row[1:] for row in observed], [(True, ""),
            (False, "extraction_failed"), (False, "extraction_failed"),
            (False, "extraction_failed"), (False, "umi_invalid"), (False, "umi_invalid")])
        self.assertEqual(result["valid_rate"], 1 / 6)

    def test_annotation_missing_and_barcode_ambiguity_are_distinct(self):
        whitelist = self.root / "PB1.txt"
        whitelist.write_text("AA\nAC\n")
        self.config.segments[0] = replace(self.config.segments[0], whitelist=whitelist,
            correction=CorrectionConfig(max_shift=0, allow_shift=False))
        self.config.barcode2_sequence_only = False
        self.config.fb_info = self.root / "annotation.tsv"
        self.config.fb_info.write_text("FB_ID\tsequence\tantibody_name\nF1\tAAGG\tfeature1\n")
        result, observed = self.compare_formal_qc([
            ("valid", "AACC", "TTAAGG"), ("ambiguous", "ATCC", "TTAAGG"),
            ("missing_annotation", "AACC", "TTAACC")])
        self.assertEqual([r[2] for r in observed], ["", "barcode_ambiguous", "annotation_missing"])
        self.assertEqual(result["valid_rate"], 1 / 3)

    def test_zero_valid_is_diagnostic_unless_threshold_configured(self):
        records = [("bad", "AA", "TTAAGG")]
        result, _ = self.compare_formal_qc(records)
        self.assertEqual(result["qc_status"], "diagnostic")
        enforce_quick_test(result)
        self.config.quick_test_min_valid_rate = 0.5
        result, _ = self.compare_formal_qc(records)
        self.assertEqual(result["qc_status"], "failed")
        self.assertTrue(result["stats"].is_file())
        with self.assertRaisesRegex(ConfigError, "min_valid_rate"):
            enforce_quick_test(result)

    def test_duplication_uses_the_same_percent_unit_as_downsample(self):
        result, _ = self.compare_formal_qc([
            ("one", "AACC", "TTAAGG"), ("two", "AACC", "TTAAGG")])
        values = {row["metric"]: row["value"] for row in open_tsv(result["stats"])
                  if row["scope"] == "sample"}
        self.assertEqual(float(values["duplication_ratio"]), 50)
        self.assertEqual(values["duplication_ratio_unit"], "percent")
        self.assertEqual(values["rate_unit"], "fraction")
        self.assertEqual(float(values["combined_valid_rate"]), 1)
        self.assertIn('<div class="value">50.00%</div>', result["html"].read_text())

    def test_n_base_threshold_is_explicit(self):
        self.config.quick_test_max_n_base_rate = 0
        result, _ = self.compare_formal_qc([("bad", "AACC", "TNAGGG")])
        with self.assertRaisesRegex(ConfigError, "max_n_base_rate"):
            enforce_quick_test(result)

    def test_disabled_correction_rejects_truncated_segment(self):
        corrector = WhitelistCorrector([], CorrectionConfig(enabled=False), 4)
        self.assertFalse(corrector.correct("AA").success)

    def test_original_exact_hit_skips_distance_search(self):
        values = ["".join(chars) for chars in product("ACGT", repeat=3)]
        corrector = WhitelistCorrector(values, CorrectionConfig(max_shift=2), 3)
        with patch("tagforge.barcode_correct.hamming", side_effect=AssertionError("distance search")):
            result = corrector.correct("TACGAA", anchor=1)
        self.assertEqual((result.corrected_sequence, result.shift_distance, result.mismatch_distance),
                         ("ACG", 0, 0))

    def test_threshold_validation(self):
        path = self.root / "threshold.yaml"
        original = Path(__file__).parents[1] / "configs/config.example.yaml"
        for key in ("min_valid_rate", "max_n_base_rate"):
            text = original.read_text().replace("workdir: ..", f"workdir: {original.parents[1]}")
            path.write_text(text.replace("quick_test:\n", f"quick_test:\n  {key}: 1.1\n"))
            with self.assertRaisesRegex(ConfigError, key):
                load_config(path, check_files=False)


if __name__ == "__main__":
    unittest.main()
