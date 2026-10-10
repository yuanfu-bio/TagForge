import json
from pathlib import Path
import tempfile
import unittest

from tagforge.config import ConfigError, load_config
from tagforge.columnar import open_records, write_parquet
from tagforge.io_utils import open_tsv, sample_dirs
from tagforge.umi_correct import _dedup_batch, _group_batches, dedup_sample


class UmiBudgetTests(unittest.TestCase):
    def test_connected_chain_is_not_split_by_batch_target(self):
        counts = {"AAAA": 20, "AAAT": 4, "AATT": 1, "CCCC": 3}
        stats = {}
        batches = list(_group_batches((("PB", "FB", u, n) for u, n in counts.items()), 1, 1, stats=stats))
        self.assertEqual(stats, {"peak_group_umis": 4, "peak_component_umis": 3})
        self.assertEqual([len(batch[0][1]) for batch in batches], [3, 1])
        expected = _dedup_batch([(("PB", "FB"), counts)], "directional", 1)[0]
        actual = [row for batch in batches for row in _dedup_batch(batch, "directional", 1)[0]]
        self.assertEqual(sorted(actual), sorted(expected))

    def test_group_budget_stops_loading_before_unbounded_growth(self):
        consumed = []
        def records():
            for umi in ("AAAA", "AAAT", "AATT", "TTTT"):
                consumed.append(umi)
                yield "PB-risk", "FB-risk", umi, 1
        stats = {}
        with self.assertRaisesRegex(ConfigError, "PB-risk.*FB-risk.*umi_max_group_umis=2"):
            list(_group_batches(records(), 1, 1, max_group_umis=2, stats=stats))
        self.assertEqual(len(consumed), 3)
        self.assertEqual(stats["peak_group_umis"], 2)

    def test_component_budget_rejects_connected_chain(self):
        rows = [("PB", "FB", u, 1) for u in ("AAAA", "AAAT", "AATT")]
        with self.assertRaisesRegex(ConfigError, "component.*PB.*FB.*3 UMIs.*umi_max_component_umis=2"):
            list(_group_batches(rows, 1, 1, max_component_umis=2))
        # A wider distance keeps the entire correction scope together.
        with self.assertRaisesRegex(ConfigError, "umi_max_component_umis=2"):
            list(_group_batches(rows, 1, 2, max_component_umis=2))

    def test_components_match_full_scope_for_supported_methods(self):
        counts = {"AAAA": 20, "AAAT": 2, "AATT": 1, "CCCC": 10, "CCCT": 2}
        for method in ("directional", "adjacency", "cluster", "unique"):
            with self.subTest(method=method):
                expected = _dedup_batch([(("PB", "FB"), counts)], method, 1)[0]
                batches = _group_batches((("PB", "FB", u, n) for u, n in counts.items()), 2, 1)
                actual = [row for batch in batches for row in _dedup_batch(batch, method, 1)[0]]
                self.assertEqual(sorted(actual), sorted(expected))

    def test_budget_configuration_rejects_fractional_and_negative_values(self):
        repo = Path(__file__).parents[1]
        text = (repo / "configs/config.example.yaml").read_text().replace("workdir: ..", f"workdir: {repo}")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.yaml"
            for value in ("0", "-1", "1.5", "true"):
                with self.subTest(value=value):
                    path.write_text(text.replace("  umi_batch_size: 5000", "  umi_batch_size: 5000\n  umi_max_group_umis: " + value))
                    with self.assertRaises(ConfigError):
                        load_config(path)
            path.write_text(text.replace("  umi_batch_size: 5000", "  umi_batch_size: 5000\n  umi_max_group_umis: 100\n  umi_max_component_umis: 50"))
            config = load_config(path)
            self.assertEqual((config.umi_max_group_umis, config.umi_max_component_umis), (100, 50))

    def test_both_backends_preserve_counts_and_report_scope_peaks(self):
        repo = Path(__file__).parents[1]
        with tempfile.TemporaryDirectory() as directory:
            for backend in ("external_sort", "sqlite"):
                with self.subTest(backend=backend):
                    config = load_config(repo / "configs/config.example.yaml")
                    config.output_dir = Path(directory) / backend
                    config.umi_aggregation_backend = backend
                    config.threads = config.umi_workers = config.umi_aggregation_workers = 1
                    config.umi_batch_size = 1
                    config.chunk_size = 3
                    dirs = sample_dirs(config.output_dir, "example")
                    rows = [{"CELL": "PB", "FB_name": "FB", "UMI": umi}
                            for umi, count in (("AAAA", 20), ("AAAT", 4), ("AATT", 1), ("CCCC", 3))
                            for _ in range(count)]
                    write_parquet(dirs["corrected"] / "example.valid_reads.parquet", rows, list(rows[0]))
                    output, stats = dedup_sample(config, "example")
                    self.assertEqual((stats["valid_reads"], stats["molecules"], stats["peak_group_umis"], stats["peak_component_umis"], stats["peak_batch_umis"]), (28, 2, 4, 3, 3))
                    self.assertEqual(sum(int(r["reads_count"]) for r in open_records(output)), 28)
                    self.assertEqual((stats["groups"], stats["total_raw_umis"], stats["components"]), (1, 4, 2))
                    self.assertFalse((config.output_dir / "example/08_tmp").exists())
                    self.assertFalse(list(dirs["detail"].glob(".umi-*")))
                    progress = json.loads((dirs["logs"] / "example.progress.json").read_text())
                    self.assertEqual((progress["peak_group_umis"], progress["peak_component_umis"]), (4, 3))
                    config.umi_max_group_umis = 2
                    with self.assertRaisesRegex(ConfigError, "umi_max_group_umis=2"):
                        dedup_sample(config, "example")
                    self.assertEqual(json.loads((dirs["logs"] / "example.progress.json").read_text())["status"], "failed")
                    self.assertFalse(list(dirs["detail"].glob(".umi-*")))


if __name__ == "__main__":
    unittest.main()
