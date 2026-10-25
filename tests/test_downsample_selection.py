from collections import Counter
from itertools import product
import json
from pathlib import Path
import re
import tempfile
import unittest

from tagforge.downsample import choose_selection, curve_rows, expected_counts, _write_saturation_html


class DownsampleSelectionTests(unittest.TestCase):
    def rows(self):
        return [{"sample": "sample", "downsample_ratio": ratio,
                 "reads_sampled": 20, "umi_types": 10, "umi_detected_once": singles,
                 "sequencing_saturation": 50,
                 "non_singleton_molecule_fraction": fraction}
                for ratio, singles, fraction in ((0.1, 5, 50), (0.5, 3, 70))]

    def test_ratio_selected_by_observed_fraction(self):
        selected = choose_selection(self.rows())
        self.assertEqual(selected["downsample_ratio"], 0.5)
        self.assertEqual(selected["non_singleton_molecule_fraction"], 70)
        self.assertEqual(selected["selection_criterion"], "max_non_singleton_molecule_fraction")
        tied = self.rows()
        tied[0]["non_singleton_molecule_fraction"] = 70
        self.assertEqual(choose_selection(list(reversed(tied)))["downsample_ratio"], 0.1)

    def test_curve_keeps_one_observed_result_per_ratio(self):
        curves = curve_rows(list(reversed(self.rows())))
        self.assertEqual([row["downsample_ratio"] for row in curves], [0.1, 0.5])
        self.assertEqual([row["non_singleton_molecule_fraction"] for row in curves], [50, 70])
        with self.assertRaisesRegex(ValueError, "one downsampling result per ratio"):
            choose_selection(self.rows() + self.rows()[:1])

    def test_count_expectations_match_exhaustive_read_subsets(self):
        observed = []
        # Molecules have support [1,1,2]; each of four reads is kept with p=1/2.
        for bits in product((0, 1), repeat=4):
            counts = (bits[0], bits[1], bits[2] + bits[3])
            observed.append((sum(counts), sum(value > 0 for value in counts), sum(value == 1 for value in counts)))
        exact = tuple(sum(row[index] for row in observed) / len(observed) for index in range(3))
        self.assertEqual(expected_counts(Counter({1: 2, 2: 1}), 0.5), exact)
        self.assertEqual(exact, (2, 1.75, 1.5))

    def test_endpoints_and_large_support_remain_finite(self):
        histogram = Counter({1: 2, 1000000: 1})
        self.assertEqual(expected_counts(histogram, 0), (0, 0, 0))
        self.assertEqual(expected_counts(histogram, 1), (1000002, 3, 2))
        reads, molecules, singletons = expected_counts(histogram, 0.000001)
        self.assertAlmostEqual(reads, 1.000002)
        self.assertTrue(0 < singletons <= molecules <= 3)

    def test_report_embeds_single_sample_metrics(self):
        curves = curve_rows(self.rows())
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "saturation.html"
            _write_saturation_html(path, "sample", curves)
            html = path.read_text()
        payload = json.loads(re.search(r"const data = (.*);", html).group(1))
        self.assertEqual(payload["non_singleton"], [50, 70])
        self.assertNotIn("sd", payload)
        self.assertEqual(payload["sample"], "sample")
        self.assertIn("<b>sample</b>", html)
        self.assertIn("Saturation 50.00%", html)
        self.assertIn("Non-singleton Molecule Fraction 70.00%", html)
        self.assertEqual(payload["umi_types"], [10, 10])
        self.assertEqual(payload["singletons"], [5, 3])
        self.assertNotIn("<h1", html)
        self.assertNotIn("<p", html)
        self.assertIn("family:'Arial'", html)
        self.assertIn("name:'Saturation'", html)
        self.assertIn("name:'Non-singleton Molecule Fraction'", html)
        self.assertIn("yaxis:'y3'", html)
        self.assertIn("yaxis:'y2'", html)
        self.assertIn("domain:[0.04,0.48]", html)
        self.assertIn("domain:[0.56,1],range:[0,100]", html)
        self.assertNotIn("error_y", html)


if __name__ == "__main__":
    unittest.main()
