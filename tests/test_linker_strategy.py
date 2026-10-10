from dataclasses import replace
import json
from pathlib import Path
import unittest

from tagforge.config import SegmentConfig
from tagforge.extract import extract_segment


class LinkerStrategyTests(unittest.TestCase):
    def test_length_position_against_independent_tag_truth(self):
        cases = json.loads((Path(__file__).parent / "fixtures/linker_truth/cases.json").read_text())
        for case in cases:
            with self.subTest(case=case["name"]):
                segment = SegmentConfig("PB", "barcode1", "R1", "linker", case["length"],
                    start=case["start"], left_linker="GGG", right_linker="TTT",
                    linker_max_mismatch=case["max_mismatch"], linker_strategy="length_position")
                result = extract_segment(case["sequence"], segment)
                self.assertEqual(result.success, case["expected_success"])
                if result.success:
                    self.assertEqual(result.raw_sequence, case["expected_sequence"])
                else:
                    self.assertTrue(result.failure_reason.startswith(case["expected_reason"]), result.failure_reason)

    def test_legacy_shortest_reproduces_short_decoy_failure(self):
        segment = SegmentConfig("PB", "barcode1", "R1", "linker", 4, start=9,
            left_linker="GGG", right_linker="TTT", linker_strategy="legacy_shortest")
        legacy = extract_segment("GGGTTTGGGACACTTT", segment)
        modern = extract_segment("GGGTTTGGGACACTTT", replace(segment, linker_strategy="length_position"))
        self.assertFalse(legacy.success)
        self.assertTrue(modern.success)
        self.assertEqual(modern.raw_sequence, "ACAC")

    def test_ambiguous_linker_keeps_explicit_fixed_fallback(self):
        # Both equal-length candidates are equally distant from fixed start 8.
        segment = SegmentConfig("PB", "barcode1", "R1", "linker_fixed", 4, start=8,
            left_linker="GGG", right_linker="TTT", linker_strategy="length_position")
        result = extract_segment("GGGACACTTTGGGCGCGTTT", segment)
        self.assertTrue(result.success)
        self.assertEqual(result.extraction_method, "fixed")
        self.assertEqual(result.linker_failure_reason, "ambiguous_linker_pair")


if __name__ == "__main__":
    unittest.main()
