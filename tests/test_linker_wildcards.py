from dataclasses import replace
import unittest

from tagforge.config import SegmentConfig
from tagforge.extract import _find_all_linkers, _find_linker, extract_segment


class LinkerWildcardTests(unittest.TestCase):
    def test_n_matches_each_base_without_spending_mismatch_budget(self):
        for base in "ACGTN":
            for budget in (0, 1):
                with self.subTest(base=base, budget=budget):
                    self.assertEqual(_find_all_linkers(base + "ACG", "NACG", budget), [0])
                    self.assertEqual(_find_linker(base + "ACG", "NACG", budget), 0)

    def test_n_does_not_hide_errors_at_literal_positions(self):
        for motif in ("NACG", "NNACG"):
            prefix = "T" * motif.count("N")
            with self.subTest(motif=motif):
                self.assertEqual(_find_all_linkers(prefix + "ATG", motif, 0), [])
                self.assertEqual(_find_all_linkers(prefix + "ATG", motif, 1), [0])
                self.assertEqual(_find_linker(prefix + "ATG", motif, 1), 0)
                self.assertEqual(_find_all_linkers(prefix + "ATT", motif, 1), [])
                self.assertEqual(_find_linker(prefix + "ATT", motif, 1), None)
                self.assertEqual(_find_all_linkers(prefix + "ATT", motif, 2), [0])

    def test_other_iupac_codes_match_their_declared_bases(self):
        codes = {"R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT",
                 "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}
        for code, allowed in codes.items():
            for base in "ACGTN":
                with self.subTest(code=code, base=base):
                    expected = [0] if base in allowed else []
                    self.assertEqual(_find_all_linkers(base + "ACG", code + "ACG", 0), expected)
                    self.assertEqual(_find_all_linkers(base + "ACG", code + "ACG", 1), [0])

    def test_overlapping_poly_a_matches_are_retained(self):
        self.assertEqual(_find_all_linkers("T" + "A" * 12, "NAAAAAAAA", 0), [0, 1, 2, 3, 4])
        self.assertEqual(_find_all_linkers("AAAAA", "AAA", 0), [0, 1, 2])

    def test_matching_requires_full_length_and_disallows_indels(self):
        self.assertEqual(_find_all_linkers("TAC", "NACG", 1), [])
        self.assertEqual(_find_all_linkers("TATCG", "NACG", 0), [])
        self.assertEqual(_find_all_linkers("TANCG", "NACG", 0), [])
        # A read N is an error at an ordinary adapter base, not a read wildcard.
        self.assertEqual(_find_all_linkers("TANG", "NACG", 0), [])
        self.assertEqual(_find_all_linkers("TANG", "NACG", 1), [0])

    def test_pb5_spacer_is_part_of_right_linker(self):
        sequence = ("CAGCGAGACCGATTACGATGCTTCAGTCACAGGACACAACCATCATACCACAGCTGTCACACCTA"
                    "AGGTCGTCGGAATTCTGGAGT" + "A" * 40)
        segment = SegmentConfig("PB5", "barcode1", "R2", "linker_fixed", 8,
            start=77, left_linker="GGAA", right_linker="NAAAAAAAA",
            linker_strategy="length_position")
        result = extract_segment(sequence, segment)
        self.assertTrue(result.success)
        self.assertEqual(result.extraction_method, "linker")
        self.assertEqual(result.raw_sequence, "TTCTGGAG")
        self.assertEqual(result.selected_linker_gap, 8)
        # The original right linker includes the spacer in the extracted gap.
        original = extract_segment(sequence, replace(segment, right_linker="AAAAAAAA"))
        self.assertEqual(original.extraction_method, "fixed")
        self.assertEqual(original.linker_failure_reason, "unexpected_length:9")

    def test_pair_scoring_counts_real_errors_instead_of_wildcards(self):
        segment = SegmentConfig("PB", "barcode1", "R1", "linker", 4,
            start=14, left_linker="NCAGTC", right_linker="GTTACC",
            linker_max_mismatch=1, linker_strategy="length_position")
        # Starts 6 and 22 are equally distant from 14. Both left windows have
        # literal Hamming distance 1, but only one has a true substitution.
        for sequence, expected in (("CCAGTCACACGTTACCNCAGTACGCGGTTACC", "ACAC"),
                                   ("NCAGTAACACGTTACCCCAGTCCGCGGTTACC", "CGCG")):
            with self.subTest(sequence=sequence):
                result = extract_segment(sequence, segment)
                self.assertTrue(result.success)
                self.assertEqual(result.raw_sequence, expected)


if __name__ == "__main__":
    unittest.main()
