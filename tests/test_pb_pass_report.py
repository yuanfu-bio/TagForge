import gzip
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from tagforge.columnar import write_parquet
from tagforge.config import ConfigError
from tagforge.io_utils import sample_dirs
from tagforge.matrix import pair_mapping_sample
from tagforge.reports import (_pair_mapping_pb_stats, _pair_mapping_batch_html,
                              _pair_mapping_dominance_qc, _assignment_panel, _assignment_plot_js)
from tagforge.dominance_qc import DominanceQC, validate_dominance_qc


class PbPassReportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.config = SimpleNamespace(
            output_dir=Path(temporary.name), compression_level=3,
            pb_cb_enabled=True, pb_cb_dominance_threshold=0.8,
            raw={"pb_cb": {"observed_correction": False}},
            target_name=lambda role: {"barcode1": "PB", "barcode2": "CB"}[role])
        self.dirs = sample_dirs(self.config.output_dir, "sample")

    def write_legacy_mapping(self, pbs):
        path = self.dirs["matrix"] / "sample.pb_cb_mapping.tsv.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("PB\tCB\tmolecule_count\treads_count\n")
            for pb, cb in pbs:
                handle.write(f"{pb}\t{cb}\t1\t1\n")
        return path

    def test_pass_rate_counts_each_pb_once_and_includes_threshold_equality(self):
        rows = [{"PB": pb, "CB_name": cb, "corrected_umi": "AAAA", "reads_count": reads}
                for pb, first, second in (("PB1", 9, 1), ("PB2", 8, 2), ("PB3", 5, 5))
                for cb, reads in (("CB1", first), ("CB2", second))]
        write_parquet(self.dirs["detail"] / "sample.molecule_detail.parquet", rows, list(rows[0]))
        pair_mapping_sample(self.config, "sample")
        summary = json.loads((self.dirs["matrix"] / "sample.pb_cb_summary.json").read_text())
        self.assertEqual((summary["total_pbs"], summary["retained_pbs"]), (3, 2))
        self.assertEqual(_pair_mapping_pb_stats(self.config, "sample", summary), (3, 2, 2 / 3))
        del summary["total_pbs"]
        self.assertEqual(_pair_mapping_pb_stats(self.config, "sample", summary), (3, 2, 2 / 3))

    def test_legacy_count_handles_repeated_pb_and_changed_file(self):
        self.write_legacy_mapping((("PB1", "CB1"), ("PB1", "CB2"), ("PB2", "CB1")))
        summary = {"retained_pbs": 1}
        self.assertEqual(_pair_mapping_pb_stats(self.config, "sample", summary), (2, 1, 0.5))
        path = self.write_legacy_mapping((("PB1", "CB1"), ("PB2", "CB1"), ("PB3", "CB1")))
        # Shared filesystems may give rapid same-size writes the same timestamp.
        stat = path.stat()
        os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
        self.assertEqual(_pair_mapping_pb_stats(self.config, "sample", summary), (3, 1, 1 / 3))

    def test_empty_summary_and_legacy_mapping_have_zero_pass_rate(self):
        self.assertEqual(_pair_mapping_pb_stats(self.config, "sample", {"total_pbs": 0, "retained_pbs": 0}), (0, 0, 0.0))
        self.write_legacy_mapping(())
        self.assertEqual(_pair_mapping_pb_stats(self.config, "sample", {"retained_pbs": 0}), (0, 0, 0.0))

    def test_invalid_counts_and_unsorted_legacy_mapping_fail_clearly(self):
        with self.assertRaisesRegex(ConfigError, "Invalid PB counts"):
            _pair_mapping_pb_stats(self.config, "sample", {"total_pbs": 1, "retained_pbs": 2})
        self.write_legacy_mapping((("PB2", "CB1"), ("PB1", "CB1")))
        with self.assertRaisesRegex(ConfigError, "not sorted by PB"):
            _pair_mapping_pb_stats(self.config, "sample", {"retained_pbs": 1})

    def test_batch_chart_separates_already_unique_resolved_and_unassigned(self):
        header = ["Sample", "Total PBs", "Retained PBs", "PB pass rate", "Single-CB PBs", "Multi-CB PBs",
                  "PBs assigned to dominant CB", "Unassigned PBs", "Multi-CB PB pass rate",
                  "Total reads", "Valid ratio", "Saturation", "Non-singleton molecule fraction", "Retained CBs"]
        overview = [header, ["A&B", 4, 2, 50.0, 1, 3, 1, 2, 100 / 3, 1000, 90.0, 55.5, 70.0, 2],
                    ["empty", 0, 0, 0.0, 0, 0, 0, 0, None, 0, 0.0, None, None, 0]]
        page = _pair_mapping_batch_html(self.config, ["A&B", "empty"], overview)
        self.assertIn('"x": [0.25, 0]', page)
        self.assertIn('"x": [0.5, 0]', page)
        for name in ("Single-CB assigned", "Assigned to dominant CB", "Unassigned"):
            self.assertIn(name, page)
        self.assertIn('<td>A&amp;B</td>', page)
        self.assertIn('<td>33.3%</td>', page)
        self.assertIn('<td>N/A</td>', page)
        self.assertNotIn('Molecule Retention', page)

    def test_outcomes_and_relation_removals_match_mapping_and_legacy(self):
        rows = [{"PB": pb, "CB_name": cb, "corrected_umi": "AAAA", "reads_count": reads}
                for pb, first, second in (("PB1", 3, 0), ("PB2", 9, 1), ("PB3", 8, 2), ("PB4", 5, 5))
                for cb, reads in (("CB1", first), ("CB2", second)) if reads]
        write_parquet(self.dirs["detail"] / "sample.molecule_detail.parquet", rows, list(rows[0]))
        pair_mapping_sample(self.config, "sample")
        summary = json.loads((self.dirs["matrix"] / "sample.pb_cb_summary.json").read_text())
        qc = _pair_mapping_dominance_qc(self.config, "sample", summary)
        self.assertEqual([qc[k] for k in ('total_pbs', 'single_cb_pbs', 'multi_cb_pbs',
            'resolved_multi_cb_pbs', 'unassigned_pbs', 'total_relations',
            'removed_secondary_relations', 'removed_unresolved_relations')], [4, 1, 3, 2, 1, 7, 2, 2])
        self.assertEqual(qc['dominant_ratio_histogram'][16]['passed_pbs'], 1)
        self.assertEqual(qc['dominant_ratio_histogram'][18]['passed_pbs'], 1)
        del summary['dominance_qc']
        self.assertEqual(_pair_mapping_dominance_qc(self.config, 'sample', summary), qc)

    def test_legacy_diagnostics_merge_cb_aliases_before_classifying(self):
        mapping = self.dirs['matrix'] / 'sample.pb_cb_mapping.tsv.gz'
        with gzip.open(mapping, 'wt') as handle:
            handle.write('PB\tCB\tmolecule_count\treads_count\n')
            handle.write('PB1\tAAA\t1\t4\nPB1\tAAT\t1\t5\nPB1\tCCC\t1\t1\n')
            handle.write('PB2\tAAA\t1\t1\nPB2\tAAT\t1\t2\n')
        with gzip.open(self.dirs['matrix'] / 'sample.cb_observed_correction.tsv.gz', 'wt') as handle:
            handle.write('raw_cb\tcorrected_cb\nAAA\tAAT\n')
        summary = {'retained_pbs': 2, 'retained_supporting_reads': 12, 'supporting_reads': 13,
                   'dominance_threshold': 0.8, 'cb_aliases': 1}
        qc = _pair_mapping_dominance_qc(self.config, 'sample', summary)
        self.assertEqual((qc['single_cb_pbs'], qc['resolved_multi_cb_pbs'], qc['unassigned_pbs']), (1, 1, 0))
        self.assertEqual((qc['total_relations'], qc['removed_secondary_relations']), (3, 1))
        # Changing alias support invalidates the cache and matches the new mapping summary.
        with gzip.open(self.dirs['matrix'] / 'sample.cb_observed_correction.tsv.gz', 'wt') as handle:
            handle.write('raw_cb\tcorrected_cb\n')
        summary.update(retained_pbs=0, retained_supporting_reads=0, cb_aliases=0)
        qc = _pair_mapping_dominance_qc(self.config, 'sample', summary)
        self.assertEqual((qc['total_relations'], qc['unassigned_pbs']), (5, 2))

    def test_no_multicb_pbs_show_resolution_as_not_applicable(self):
        qc = DominanceQC(0.8)
        qc.add(1, 2, 2)
        plots = _assignment_plot_js(qc.payload())
        self.assertIn('Multi-CB PB Pass Rate: <b>N/A</b>', plots)
        self.assertIn('No multi-CB PBs', plots)

    def test_threshold_one_retains_pure_singlecb_pbs(self):
        qc = DominanceQC(1)
        qc.add(1, 2, 2)
        payload = validate_dominance_qc(qc.payload(), {'retained_pbs': 1})
        self.assertEqual((payload['single_cb_pbs'], payload['single_cb_retained_pbs'], payload['unassigned_pbs']), (1, 1, 0))
