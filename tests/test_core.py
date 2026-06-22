from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from tagforge import __version__
from tagforge.barcode_correct import WhitelistCorrector
from tagforge.config import CorrectionConfig, SegmentConfig, default_downsample_ratios, load_config
from tagforge.downsample import calculate_metrics
from tagforge.extract import extract_segment
from tagforge.fastq import paired_fastq
from tagforge.slurm import make_slurm
from tagforge.quick_test import _take_leading_records
from tagforge.umi_correct import deduplicate_umis


class CoreTests(unittest.TestCase):
    def test_quick_test_takes_only_leading_records(self):
        selected, scanned = _take_leading_records(iter(range(1000)), 10)
        self.assertEqual(scanned, 10)
        self.assertEqual(selected, list(enumerate(range(10), 1)))

    def test_version_metadata_is_synchronized(self):
        root = Path(__file__).parents[1]
        pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
        setup = (root / "setup.py").read_text(encoding="utf-8")
        self.assertEqual(__version__, "0.1.3")
        self.assertIn(f'version = "{__version__}"', pyproject)
        self.assertIn(f'version="{__version__}"', setup)

    def test_fixed_and_linker_extraction(self):
        fixed = SegmentConfig("x", "umi", "R1", "fixed", 3, start=2)
        self.assertEqual(extract_segment("AACCCGG", fixed).raw_sequence, "CCC")
        linker = SegmentConfig("x", "barcode1", "R1", "linker", 4, right_linker="GG")
        combined = SegmentConfig(
            "x", "barcode1", "R1", "linker_fixed", 4,
            start=0, left_linker="GG", right_linker="CC",
            correction=CorrectionConfig(False, 0, 0, False, False),
        )
        def exact_linker(sequence, motif, max_mismatch, start=0):
            position = sequence.find(motif, start)
            return None if position < 0 else position
        def closest_pair(sequence, left, right, max_mismatch):
            lefts = [i for i in range(len(sequence)) if sequence.startswith(left, i)]
            rights = [i for i in range(len(sequence)) if sequence.startswith(right, i)]
            pairs = [(r - (l + len(left)), l, r) for l in lefts for r in rights if r >= l + len(left)]
            return (min(pairs) if pairs else None, len(lefts), len(rights), len(pairs))
        with patch("tagforge.extract._find_linker", side_effect=exact_linker), \
             patch("tagforge.extract._closest_linker_pair", side_effect=closest_pair):
            self.assertEqual(extract_segment("ACGTGG", linker).raw_sequence, "ACGT")
            # Linker succeeds: use ACGT between linkers, not fixed [0:4]=GGAC.
            linker_result = extract_segment("GGACGTCC", combined)
            self.assertEqual((linker_result.raw_sequence, linker_result.extraction_method), ("ACGT", "linker"))
            # Linker fails: fixed fallback [0:4] runs on the original read.
            fixed_result = extract_segment("ACGTZZZZ", combined)
            self.assertEqual((fixed_result.raw_sequence, fixed_result.extraction_method), ("ACGT", "fixed"))
            # Linkers are found but delimit five bases: wrong length also falls back.
            length_fallback = extract_segment("ACGTGGNNNNNCC", combined)
            self.assertEqual(length_fallback.extraction_method, "fixed")
            self.assertEqual(length_fallback.linker_failure_reason, "unexpected_length:5")
            short_linkers = SegmentConfig(
                "PB2", "barcode1", "R1", "linker_fixed", 8,
                start=16, left_linker="AACC", right_linker="ACAG",
                correction=CorrectionConfig(False, 0, 0, False, False),
            )
            multi = extract_segment("AACCxxxxAACC12345678ACAGzzACAG", short_linkers)
            self.assertEqual(multi.raw_sequence, "12345678")
            self.assertEqual(multi.selected_linker_gap, 8)
            self.assertEqual(multi.linker_candidate_pairs, 4)

    def test_corrections(self):
        corr = WhitelistCorrector(["ACGT"], CorrectionConfig(True, 1, 1, True, True), 4)
        self.assertEqual(corr.correct("ACGT").correction_type, "exact")
        self.assertEqual(corr.correct("ACGA").correction_type, "mismatch_only")
        self.assertEqual(corr.correct("NACGT").correction_type, "shift_only")
        self.assertEqual(corr.correct("NACGA").correction_type, "shift_and_mismatch")
        ambiguous = WhitelistCorrector(["AAAA", "AACC"], CorrectionConfig(True, 0, 1, False, True), 4)
        self.assertTrue(ambiguous.correct("AAAC").ambiguous)

    def test_umi_directional(self):
        class FakeClusterer:
            def __call__(self, counts, threshold):
                self.assertions = (counts, threshold)
                return [[b"AAAA", b"AAAT"], [b"CCCC"]]
        with patch("tagforge.umi_correct._umi_clusterer", return_value=FakeClusterer()):
            result = deduplicate_umis({"AAAA": 10, "AAAT": 2, "CCCC": 3}, "directional", 1)
        self.assertEqual(result["AAAT"], "AAAA")
        self.assertEqual(result["CCCC"], "CCCC")

    def test_downsample_metrics(self):
        self.assertEqual(calculate_metrics([3, 1, 0])[:3], (4, 2, 1))
        ratios = default_downsample_ratios()
        self.assertEqual(len(ratios), 36)
        self.assertEqual(ratios[:3], [0.0001, 0.0002, 0.0003])
        self.assertEqual(ratios[-3:], [0.7, 0.8, 0.9])

    def test_config_accepts_composed_extraction(self):
        with tempfile.TemporaryDirectory() as td:
            config_path = Path(td) / "config.yaml"
            config_path.write_text("""
project:
  workdir: .
samples:
  - sample: s
    r1: r1.fq.gz
    r2: r2.fq.gz
barcode2_annotation:
  fb_info: fb.tsv
segments:
  - name: CELL
    target: barcode1
    read: R1
    methods: [linker, fixed]
    left_linker: GG
    right_linker: TT
    start: 2
    length: 4
    correction:
      enabled: false
  - name: FB
    target: barcode2
    read: R2
    method: fixed
    start: 0
    length: 4
    correction:
      enabled: false
  - name: UMI
    target: umi
    read: R2
    method: fixed
    start: 4
    length: 4
""", encoding="utf-8")
            config = load_config(config_path, check_files=False)
            self.assertEqual(config.segments[0].method, "linker_fixed")

    def test_slurm_scheduler_options(self):
        with tempfile.TemporaryDirectory() as td:
            config = SimpleNamespace(
                samples=[SimpleNamespace(sample="s1")], path=Path(td) / "config.yaml"
            )
            files = make_slurm(
                config, Path(td) / "jobs", 12, "32G", "08:00:00",
                partition="compute", account="lab", qos="normal",
                constraint="avx2", nodes=1, ntasks=1, mail_user="a@example.org",
                mail_type="END,FAIL", conda_env="tagforge",
                skip_quick_test=True,
                extra_sbatch=["--exclusive"],
            )
            script = files[0].read_text(encoding="utf-8")
            self.assertIn("#SBATCH --partition=compute", script)
            self.assertIn("#SBATCH --account=lab", script)
            self.assertIn("#SBATCH --cpus-per-task=12", script)
            self.assertIn("#SBATCH --exclusive", script)
            self.assertIn("conda activate tagforge", script)
            self.assertIn("--threads 12", script)
            self.assertIn("--skip-quick-test", script)

    def test_paired_fastq(self):
        with tempfile.TemporaryDirectory() as td:
            paths = [Path(td) / f"r{x}.fq.gz" for x in (1, 2)]
            for i, path in enumerate(paths, 1):
                with gzip.open(path, "wt") as h: h.write(f"@r/ {i}\nAC\n+\nII\n")
            reads = list(paired_fastq(*paths))
            self.assertEqual(reads[0].read_id, "r/")


if __name__ == "__main__":
    unittest.main()
