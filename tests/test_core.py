from __future__ import annotations

import gzip
import tempfile
import unittest
import zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from tagforge import __version__
from tagforge.barcode_correct import WhitelistCorrector
from tagforge.config import ConfigError, CorrectionConfig, SegmentConfig, default_downsample_ratios, load_config
from tagforge.downsample import _analysis_ratios, calculate_metrics
from tagforge.extract import decode_method_payload, decode_segment_payload, extract_segment
from tagforge.fastq import _physical_position, open_text, paired_fastq, paired_fastq_batches
from tagforge.io_utils import touch_checkpoint
from tagforge.reports import write_summary
from tagforge.slurm import make_slurm
from tagforge.quick_test import _take_leading_records
from tagforge.umi_correct import _aggregated_tsv_cursor, _dedup_batch, _external_sort_aggregate, _group_batches, deduplicate_umis


class CoreTests(unittest.TestCase):
    def test_quick_test_takes_only_leading_records(self):
        selected, scanned = _take_leading_records(iter(range(1000)), 10)
        self.assertEqual(scanned, 10)
        self.assertEqual(selected, list(enumerate(range(10), 1)))

    def test_version_metadata_is_synchronized(self):
        root = Path(__file__).parents[1]
        pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
        setup = (root / "setup.py").read_text(encoding="utf-8")
        self.assertEqual(__version__, "0.1.13")
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

    def test_symmetric_fixed_shift_and_compact_payload(self):
        correction = CorrectionConfig(True, 1, 0, True, False)
        segment = SegmentConfig(
            "PB", "barcode1", "R1", "fixed", 8, start=16, correction=correction
        )
        raw = extract_segment("N" * 15 + "L" + "ABCDEFGH" + "R", segment).raw_sequence
        self.assertEqual(raw, "LABCDEFGHR")
        self.assertEqual(
            WhitelistCorrector(["LABCDEFG"], correction, 8).correct(raw, anchor=1).shift_distance,
            -1,
        )
        self.assertEqual(
            WhitelistCorrector(["ABCDEFGH"], correction, 8).correct(raw, anchor=1).shift_distance,
            0,
        )
        self.assertEqual(
            WhitelistCorrector(["BCDEFGHR"], correction, 8).correct(raw, anchor=1).shift_distance,
            1,
        )
        self.assertEqual(decode_segment_payload('AAAA,CCCC', [segment, SegmentConfig(
            "PB2", "barcode1", "R1", "fixed", 4, start=0
        )]), {"PB": "AAAA", "PB2": "CCCC"})
        with self.assertRaises(ValueError):
            decode_segment_payload('{"PB":"AAAA"}', [segment])
        self.assertEqual(
            decode_method_payload("FL", [segment, SegmentConfig(
                "PB2", "barcode1", "R1", "fixed", 4, start=0
            )]),
            {"PB": "fixed", "PB2": "linker"},
        )

    def test_umi_directional(self):
        class FakeClusterer:
            def __call__(self, counts, threshold):
                self.assertions = (counts, threshold)
                return [[b"AAAA", b"AAAT"], [b"CCCC"]]
        with patch("tagforge.umi_correct._umi_clusterer", return_value=FakeClusterer()):
            result = deduplicate_umis({"AAAA": 10, "AAAT": 2, "CCCC": 3}, "directional", 1)
        self.assertEqual(result["AAAT"], "AAAA")
        self.assertEqual(result["CCCC"], "CCCC")

    def test_external_sort_aggregation_counts_cross_chunk_duplicates(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            source = root / "valid.tsv"
            source.write_text(
                "PB\tFB_name\tUMI\n"
                "PB2\tCD3\tTTTT\n"
                "PB1\tCD3\tAAAA\n"
                "PB1\tCD3\tAAAA\n"
                "PB1\tCD3\tAAAT\n"
                "PB2\tCD3\tTTTT\n",
                encoding="utf-8",
            )
            scratch = root / "scratch"; scratch.mkdir()
            events = []
            aggregated, reads, raw_umis, groups, version = _external_sort_aggregate(
                source, scratch, "PB", "FB_name", "UMI", 1, 1, lambda *event: events.append(event),
            )
            self.assertIn("GNU", version)
            self.assertEqual(events[-1], ("aggregate_complete", 5, 3, 2))
            self.assertIn(("sort_complete", 5, 0, 0), events)
            self.assertEqual((reads, raw_umis, groups), (5, 3, 2))
            self.assertEqual(list(_aggregated_tsv_cursor(aggregated)), [
                ("PB1", "CD3", "AAAA", 2), ("PB1", "CD3", "AAAT", 1),
                ("PB2", "CD3", "TTTT", 2),
            ])

    def test_umi_group_batches_are_bounded_without_splitting_groups(self):
        cursor = [
            ("b1", "f1", "AAAA", 4),
            ("b1", "f1", "AAAT", 1),
            ("b1", "f2", "CCCC", 2),
            ("b2", "f1", "GGGG", 3),
            ("b2", "f1", "GGGA", 1),
        ]
        batches = list(_group_batches(cursor, batch_umis=2))
        self.assertEqual([[key for key, _ in batch] for batch in batches], [
            [("b1", "f1")], [("b1", "f2")], [("b2", "f1")],
        ])
        class FakeClusterer:
            def __call__(self, counts, threshold):
                return [[umi] for umi in counts]
        with patch("tagforge.umi_correct._umi_clusterer", return_value=FakeClusterer()):
            rows, groups, reads, raw_umis = _dedup_batch(batches[0], "unique", 0)
        self.assertEqual((groups, reads, raw_umis), (1, 5, 2))
        self.assertEqual(len(rows), 2)

    def test_downsample_metrics(self):
        self.assertEqual(calculate_metrics([3, 1, 0])[:3], (4, 2, 1))
        ratios = default_downsample_ratios()
        self.assertEqual(len(ratios), 36)
        self.assertEqual(ratios[:3], [0.0001, 0.0002, 0.0003])
        self.assertEqual(ratios[-3:], [0.7, 0.8, 0.9])
        config = SimpleNamespace(downsample_ratios=ratios)
        self.assertEqual(len(_analysis_ratios(config)), 38)
        self.assertEqual(_analysis_ratios(config)[0], 0.0)
        self.assertEqual(_analysis_ratios(config)[-1], 1.0)

    def test_summary_rebuilds_completed_samples_in_config_order(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td) / "02_output"
            config = SimpleNamespace(output_dir=root, samples=[SimpleNamespace(sample="alpha"), SimpleNamespace(sample="zeta")], segments=[SimpleNamespace(name="CELL", target="barcode1"), SimpleNamespace(name="FB", target="barcode2")], target_name=lambda target: {"barcode1": "CELL", "barcode2": "FB"}[target])

            def complete(sample, feature):
                base = root / sample
                stats = base / "03_corrected" / f"{sample}.barcode_correction_stats.tsv"
                stats.parent.mkdir(parents=True)
                stats.write_text("scope\ttotal_reads\tvalid_rate\nCELL\t100\t0.9\nFB\t100\t0.8\ncombined\t100\t0.7\n", encoding="utf-8")
                point = base / "06_downsample" / f"{sample}.optimal_saturation_point.tsv"
                point.parent.mkdir(parents=True)
                point.write_text("sample\toptimal_downsample_ratio\tmax_sequencing_saturation\treads_sampled\tumi_types\tumi_detected_once\tduplication_ratio\n" + f"{sample}\t0.5\t60\t50\t20\t8\t30\n", encoding="utf-8")
                detail = base / "05_detail" / f"{sample}.optimal_saturation_molecule_detail.tsv.gz"
                detail.parent.mkdir(parents=True)
                with gzip.open(detail, "wt", encoding="utf-8") as handle:
                    handle.write(f"CELL\tFB_name\tcorrected_umi\treads_count_at_optimal_downsample\nA\t{feature}\tAAAA\t1\n")
                checkpoint = base / "01_checkpoint" / "downsample.done"
                checkpoint.parent.mkdir(parents=True)
                touch_checkpoint(checkpoint, __version__)

            complete("zeta", "Z_feature")
            _, samples = write_summary(config)
            self.assertEqual(samples, ["zeta"])
            complete("alpha", "A_feature")
            output, samples = write_summary(config)
            self.assertEqual(samples, ["alpha", "zeta"])
            with zipfile.ZipFile(output) as workbook:
                meta = workbook.read("xl/worksheets/sheet1.xml").decode("utf-8")
                count = workbook.read("xl/worksheets/sheet2.xml").decode("utf-8")
            self.assertLess(meta.index("alpha"), meta.index("zeta"))
            self.assertIn("CELL Valid rate", meta)
            self.assertIn("A_feature", count)
            self.assertIn("Z_feature", count)

    def test_config_accepts_grouped_barcode_and_umi_sections(self):
        with tempfile.TemporaryDirectory() as td:
            raw = Path(td) / "raw"
            for sample in ("A-44", "A"):
                (raw / sample).mkdir(parents=True)
            config_path = Path(td) / "config.yaml"
            config_path.write_text("""
project:
  workdir: .
samples:
  auto:
    raw_dir: raw
    r1: "{sample}_raw_1.fq.gz"
    r2: "{sample}_raw_2.fq.gz"
barcode2_annotation:
  fb_info: fb.tsv
correction_barcode:
  enabled: true
  max_shift: 1
  max_mismatch: 1
linker:
  max_mismatch: 1
correction_umi:
  method: adjacency
  max_distance: 2
performance:
  umi_aggregation_workers: 3
barcode1:
  name: PB
  segments:
    - segment: PB1
      read: R1
      methods: [linker, fixed]
      left_linker: GG
      right_linker: TT
      start: 2
      length: 4
    - segment: PB2
      read: R1
      method: fixed
      start: 8
      length: 4
      linker_max_mismatch: 0
      correction:
        enabled: false
barcode2:
  name: FB
  segments:
    - segment: FB
      read: R2
      method: fixed
      start: 0
      length: 4
      correction:
        max_shift: 0
umi:
  name: UMI
  segments:
    - segment: UMI
      read: R2
      method: fixed
      start: 4
      length: 4
      correction:
        enabled: false
""", encoding="utf-8")
            config = load_config(config_path, check_files=False)
            self.assertEqual([sample.sample for sample in config.samples], ["A", "A-44"])
            self.assertEqual(config.samples[0].r1, (raw / "A/A_raw_1.fq.gz").resolve())
            self.assertEqual(config.samples[1].r2, (raw / "A-44/A-44_raw_2.fq.gz").resolve())
            self.assertEqual([s.name for s in config.segments], ["PB1", "PB2", "FB", "UMI"])
            self.assertEqual(config.segments[0].method, "linker_fixed")
            self.assertEqual(config.segments[0].target, "barcode1")
            self.assertEqual(config.segments[0].target_name, "PB")
            self.assertEqual(config.segments[0].linker_max_mismatch, 1)
            self.assertTrue(config.segments[0].correction.enabled)
            self.assertFalse(config.segments[1].correction.enabled)
            self.assertEqual(config.segments[1].linker_max_mismatch, 0)
            self.assertEqual(config.segments[2].target_name, "FB")
            self.assertEqual(config.segments[2].correction.max_shift, 0)
            self.assertEqual(config.umi_method, "adjacency")
            self.assertEqual(config.umi_max_distance, 2)
            self.assertEqual(config.umi_aggregation_workers, 3)
            self.assertEqual(config.umi_aggregation_backend, "external_sort")
            self.assertEqual(config.umi_sort_memory_mb, 512)

    def test_config_accepts_legacy_segments(self):
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

    def test_config_auto_samples_reports_missing_fastq(self):
        with tempfile.TemporaryDirectory() as td:
            raw = Path(td) / "raw/A"
            raw.mkdir(parents=True)
            config_path = Path(td) / "config.yaml"
            config_path.write_text("""
project:
  workdir: .
samples:
  auto:
    raw_dir: raw
barcode2_annotation:
  fb_info: fb.tsv
segments:
  - name: CELL
    target: barcode1
    read: R1
    method: fixed
    start: 0
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
            with self.assertRaisesRegex(ConfigError, "A_raw_1.fq.gz"):
                load_config(config_path, check_files=True)

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
            self.assertEqual([path.name for path in files], ["samples.tsv", "tagforge_array.slurm"])
            self.assertIn("1\ts1", files[0].read_text(encoding="utf-8"))
            script = files[1].read_text(encoding="utf-8")
            self.assertIn("#SBATCH --array=1-1", script)
            self.assertIn("#SBATCH --partition=compute", script)
            self.assertIn("#SBATCH --account=lab", script)
            self.assertIn("#SBATCH --cpus-per-task=12", script)
            self.assertIn("#SBATCH --exclusive", script)
            self.assertIn("conda activate tagforge", script)
            self.assertIn('sample=$(awk -v i="$SLURM_ARRAY_TASK_ID"', script)
            self.assertIn('--sample "$sample"', script)
            self.assertIn("--threads 12", script)
            self.assertIn("--skip-quick-test", script)

    def test_slurm_per_sample_mode(self):
        with tempfile.TemporaryDirectory() as td:
            config = SimpleNamespace(
                samples=[SimpleNamespace(sample="s1"), SimpleNamespace(sample="s2")],
                path=Path(td) / "config.yaml",
            )
            files = make_slurm(
                config, Path(td) / "jobs", 4, "8G", "01:00:00",
                mode="per-sample", array_limit=None,
            )
            self.assertEqual([path.name for path in files], ["s1.slurm", "s2.slurm", "submit_all.sh"])
            self.assertIn("--sample s1", files[0].read_text(encoding="utf-8"))
            self.assertIn("sbatch", files[-1].read_text(encoding="utf-8"))

    def test_paired_fastq(self):
        with tempfile.TemporaryDirectory() as td:
            paths = [Path(td) / f"r{x}.fq.gz" for x in (1, 2)]
            for i, path in enumerate(paths, 1):
                with gzip.open(path, "wt") as h:
                    for n in range(3):
                        h.write(f"@r{n}/ {i}\nAC\n+\nII\n")
            reads = list(paired_fastq(*paths))
            self.assertEqual(reads[0].read_id, "r0/")
            batches = list(paired_fastq_batches(*paths, batch_size=1))
            self.assertEqual([len(batch) for batch, _ in batches], [1, 1, 1])
            self.assertLess(batches[0][1], 1.0)
            self.assertEqual(batches[-1][1], 1.0)
            with open_text(paths[0]) as handle:
                handle.readline()
                self.assertLessEqual(_physical_position(handle), paths[0].stat().st_size)


if __name__ == "__main__":
    unittest.main()
