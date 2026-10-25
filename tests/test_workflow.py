from dataclasses import replace
import gzip
import json
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from tagforge.checkpoint import STAGE_NUMBERS, STAGE_VERSIONS, stage_current, stage_plan
from tagforge.columnar import open_records, write_parquet
from tagforge.config import ConfigError, CorrectionConfig, load_config
from tagforge.external_tools import ExternalToolVersions
from tagforge.io_utils import open_tsv, sample_dirs
from tagforge.matrix import pair_mapping_sample
from tagforge.pipeline import STEP_FUNCS, run_pipeline, run_step
from tagforge.reports import _completed_summary_samples
from tagforge.selected import selection_point


def workbook_tables(path):
    from xml.etree import ElementTree
    from zipfile import ZipFile
    namespace = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    tables = {}
    with ZipFile(path) as workbook:
        sheets = ElementTree.fromstring(workbook.read("xl/workbook.xml")).find("s:sheets", namespace)
        for index, sheet in enumerate(sheets, 1):
            xml = ElementTree.fromstring(workbook.read(f"xl/worksheets/sheet{index}.xml"))
            tables[sheet.attrib["name"]] = [
                [cell.findtext("s:v", namespaces=namespace) if cell.find("s:v", namespace) is not None
                 else cell.findtext("s:is/s:t", default="", namespaces=namespace) for cell in row]
                for row in xml.findall("s:sheetData/s:row", namespace)]
    return tables


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        repository = Path(__file__).parents[1]
        fixture = self.root / "fixture"
        shutil.copytree(repository / "tests/fixtures/small_fastq", fixture)
        self.config = load_config(repository / "configs/config.example.yaml")
        self.config.workdir = self.root
        self.config.output_dir = self.root / "output"
        self.config.samples = [replace(self.config.samples[0],
            r1=fixture / "raw/example/example_raw_1.fq.gz",
            r2=fixture / "raw/example/example_raw_2.fq.gz")]
        self.config.segments = [replace(s, whitelist=fixture / "config" / s.whitelist.name)
                                if s.whitelist else s for s in self.config.segments]
        self.config.fb_info = fixture / "config/FB_info.tsv"
        self.config.threads = 1
        self.config.chunk_size = 3
        self.config.downsample_ratios = [0.5, 1]
        mocked = patch("tagforge.pipeline.check_external_tools", return_value=ExternalToolVersions("test", "test"))
        mocked.start(); self.addCleanup(mocked.stop)

    def test_compact_outputs_keep_curve_information_and_preview_is_optional(self):
        from tagforge.downsample import curve_rows, expected_counts
        from collections import Counter
        old = sample_dirs(self.config.output_dir, "example")
        for suffix in ("quick_test.read_qc.tsv", "quick_test.sampled_read_ids.txt"):
            (old["report"] / f"example.{suffix}").write_text("old debug output")
        run_pipeline(self.config, ["example"])
        self.assertFalse((self.config.output_dir / "example/08_tmp").exists())
        self.assertFalse((old["checkpoint"] / ".work").exists())
        self.assertFalse((old["report"] / "example.quick_test.read_qc.tsv").exists())
        dirs = sample_dirs(self.config.output_dir, "example")
        self.assertFalse((dirs["extracted"] / "example.extracted.preview.tsv").exists())

        self.assertTrue((dirs["corrected"] / "example.correction_qc.json").exists())
        self.assertFalse((dirs["corrected"] / "example.barcode2_annotation_qc.tsv").exists())
        progress = json.loads((dirs["logs"] / "example.progress.json").read_text())
        self.assertEqual((progress["stage"], progress["status"]), ("report", "completed"))
        self.assertFalse(list(dirs["logs"].glob("*_progress.tsv")))
        metrics = list(open_tsv(dirs["downsample"] / "example.downsample_metrics.tsv"))
        curves = list(open_tsv(dirs["downsample"] / "example.downsample_curve.tsv"))
        summaries = curve_rows(metrics)
        histogram = Counter(int(row["reads_count"]) for row in open_records(dirs["detail"] / "example.molecule_detail.parquet"))
        for curve, summary in zip(curves, summaries):
            for key, value in summary.items():
                if key not in {"sample", "downsample_ratio"}:
                    self.assertEqual(float(curve["observed_" + key]), float(value))
            observed = tuple(float(curve[key]) for key in ("expected_reads", "expected_molecules", "expected_singletons"))
            self.assertEqual(observed, expected_counts(histogram, float(curve["downsample_ratio"])))
        for suffix in ("selection_point.tsv", "optimal_saturation_point.tsv", "downsample_summary.tsv", "downsample_expectations.tsv"):
            self.assertFalse((dirs["downsample"] / ("example." + suffix)).exists())
        self.config.extraction_preview_reads = 2
        self.assertEqual(run_step(self.config, "example", "extract")[1], "skipped")
        self.assertEqual(len(list(open_tsv(dirs["extracted"] / "example.extracted.preview.tsv"))), 2)
        self.config.extraction_preview_reads = 0
        run_step(self.config, "example", "extract")
        self.assertFalse((dirs["extracted"] / "example.extracted.preview.tsv").exists())

    def test_pipeline_layout_and_main_result_links_for_each_library(self):
        from urllib.parse import unquote
        import re
        expected_names = {"standard": "example.raw_molecule_pairs.parquet",
                          "pi-seq": "example.molecule_detail.rmMP.parquet",
                          "pb-cb": "example.pb_cb_map.tsv.gz"}
        for library, filename in expected_names.items():
            with self.subTest(library=library):
                config = replace(self.config, output_dir=self.root / library,
                                 pi_seq_enabled=library == "pi-seq", pb_cb_enabled=library == "pb-cb",
                                 downsample_enabled=False)
                if library == "pb-cb":
                    config.barcode2_sequence_only = True
                    config.fb_info = None
                    config.segments = [replace(segment, whitelist=None, correction=CorrectionConfig(enabled=False))
                                       if segment.target == "barcode2" else segment for segment in config.segments]
                run_pipeline(config, ["example"])
                root = config.output_dir / "example"
                self.assertEqual({p.name for p in root.iterdir() if p.is_dir()},
                                 {"00_logs", "01_extracted", "02_corrected", "03_molecules",
                                  "04_matrix", "05_downsample", "06_report"})
                self.assertTrue((root / "02_corrected/example.valid_reads.parquet").is_dir())
                self.assertFalse((root / "03_molecules/example.valid_reads.parquet").exists())
                manifests = {p.name for p in (root / "00_logs/checkpoints").glob("*.manifest.json")}
                self.assertEqual(manifests, {"01_extract.manifest.json", "02_correct.manifest.json",
                                            "03_dedup.manifest.json", "06_report.manifest.json",
                                            "04_pair-map.manifest.json" if library == "pb-cb" else "04_matrix.manifest.json"})
                for report in (root / "06_report/example.report.html",
                               config.output_dir / "00_report/TagForge_batch_report.html"):
                    page = report.read_text()
                    links = re.findall(r'href="([^"]+)"', page)
                    primary = [report.parent / unquote(link) for link in links
                               if Path(unquote(link)).name == filename]
                    if report.name == "TagForge_batch_report.html":
                        self.assertEqual(primary, [])
                        self.assertNotIn("Main results", page)
                        self.assertNotIn("Downsample not run.", page)
                        self.assertNotIn("Not run. Counts use full-depth molecules.", page)
                    else:
                        self.assertEqual(primary, [])
                        self.assertNotIn("Main results", page)
                rows = workbook_tables(root / "06_report/example.report.xlsx")["Results"]
                self.assertEqual(rows[1][1], library)
                self.assertEqual(Path(rows[1][2]).name, filename)
                self.assertTrue((config.output_dir / rows[1][2]).exists())

    def test_disabled_downsample_uses_full_counts_without_selection(self):
        self.config.downsample_enabled = False
        run_pipeline(self.config, ["example"])
        dirs = sample_dirs(self.config.output_dir, "example")
        self.assertNotIn("downsample", stage_plan(self.config))
        self.assertFalse((dirs["checkpoint"] / "05_downsample.manifest.json").exists())
        self.assertFalse((dirs["downsample"] / "example.selection.json").exists())
        self.assertEqual(sum(int(row["molecule_count"]) for row in open_records(dirs["matrix"] / "example.raw_molecule_pairs.parquet")), 6)
        html = (dirs["report"] / "example.report.html").read_text()
        self.assertNotIn("Not run. Counts use full-depth molecules", html)
        self.assertNotIn('class="sub">Selected downsample ratio', html)
        self.assertEqual(_completed_summary_samples(self.config), ["example"])
        with self.assertRaisesRegex(ConfigError, "downsample.enabled"):
            run_step(self.config, "example", "downsample")

    def test_toggling_downsample_excludes_stale_selection(self):
        run_pipeline(self.config, ["example"])
        self.config.downsample_enabled = False
        self.assertFalse(stage_current(self.config, "example", "report"))
        self.assertFalse(stage_current(self.config, "example", "downsample"))
        with patch("tagforge.pipeline.downsample_sample", side_effect=AssertionError("disabled stage")):
            run_pipeline(self.config, ["example"])
        self.assertTrue(stage_current(self.config, "example", "report"))

    def test_resume_disabled_recomputes_every_planned_stage(self):
        self.config.downsample_enabled = False
        self.config.resume_enabled = False
        run_pipeline(self.config, ["example"])
        checkpoint = self.config.output_dir / "example/00_logs/checkpoints"
        before = {step: json.loads((checkpoint / f"{STAGE_NUMBERS[step]}_{step}.manifest.json").read_text())["generation"]
                  for step in stage_plan(self.config)}
        run_pipeline(self.config, ["example"])
        after = {step: json.loads((checkpoint / f"{STAGE_NUMBERS[step]}_{step}.manifest.json").read_text())["generation"] for step in before}
        self.assertTrue(all(before[step] != after[step] for step in before))
        self.assertEqual(len(list(open_records(self.config.output_dir / "example/02_corrected/example.valid_reads.parquet"))), 9)

    def test_threads_chunks_and_trace_preserve_counts_and_seeded_selection(self):
        reference = None
        for workers, chunk, trace in ((1, 3, True), (2, 1, False), (4, 5, True)):
            with self.subTest(workers=workers, chunk=chunk, trace=trace):
                self.config.threads = workers
                self.config.chunk_size = chunk
                self.config.trace_enabled = trace
                self.config.workdir = self.root / f"run-{workers}"
                self.config.output_dir = self.config.workdir / "output"
                run_pipeline(self.config, ["example"])
                dirs = sample_dirs(self.config.output_dir, "example")
                rows = lambda path: sorted(tuple(sorted(row.items())) for row in open_records(path))
                actual = (rows(dirs["corrected"] / "example.valid_reads.parquet"),
                          rows(dirs["detail"] / "example.molecule_detail.parquet"),
                          rows(dirs["matrix"] / "example.raw_molecule_pairs.parquet"),
                          list(open_tsv(dirs["downsample"] / "example.downsample_metrics.tsv")),
                          selection_point(self.config, "example"))
                if reference is None:
                    reference = actual
                else:
                    self.assertEqual(actual, reference)
                self.assertEqual(len(actual[0]), 9)
                self.assertEqual(len(actual[1]), 6)

    def test_report_version_change_reuses_upstream_counts(self):
        self.config.downsample_enabled = False
        run_pipeline(self.config, ["example"])
        directory = self.config.output_dir / "example/00_logs/checkpoints"
        generations = lambda: {step: json.loads((directory / f"{STAGE_NUMBERS[step]}_{step}.manifest.json").read_text())["generation"]
                               for step in stage_plan(self.config)}
        before = generations()
        def upstream_ran(*_args, **_kwargs):
            self.fail("report version change recomputed an upstream stage")
        with patch.dict(STAGE_VERSIONS, {"report": STAGE_VERSIONS["report"] + 1}):
            with patch.dict(STEP_FUNCS, {step: upstream_ran for step in ("extract", "correct", "dedup", "matrix")}):
                run_pipeline(self.config, ["example"])
            after = generations()
            self.assertTrue(all(after[step] == before[step] for step in before if step != "report"))
            self.assertNotEqual(after["report"], before["report"])

    def test_batch_reports_stay_in_output_and_are_withdrawn_when_inputs_change(self):
        from tagforge.checkpoint import invalidate

        self.config.output_dir = self.root / "custom/nested/results"
        run_pipeline(self.config, ["example"])
        directory = self.config.output_dir / "00_report"
        names = ("TagForge_batch_report.xlsx", "TagForge_batch_report.html")
        self.assertTrue(all((directory / name).is_file() for name in names))
        self.assertFalse((self.root / "00_report").exists())
        self.assertFalse((self.config.output_dir / "00_summary.xlsx").exists())
        self.assertFalse((directory / "00_summary.xlsx").exists())
        self.assertFalse((directory / "01_FB_not_in_anno.xlsx").exists())
        counts = workbook_tables(directory / "TagForge_batch_report.xlsx")["Counts_raw"]
        self.assertEqual(counts[1][0], "example")
        self.assertEqual(sum(int(value) for value in counts[1][1:]), 6)
        invalidate(self.config, "example", "correct")
        self.assertTrue(all(not (directory / name).exists() for name in names))
        run_pipeline(self.config, ["example"])
        self.assertTrue(all((directory / name).is_file() for name in names))

    def test_merged_workbook_distinguishes_mixed_full_depth_and_selected_counts(self):
        from tagforge.downsample import choose_selection
        from tagforge.reports import write_summary
        from tagforge.selected import load_selection
        self.config.samples = [replace(self.config.samples[0], sample=name) for name in ("full", "subset")]
        def force_ratio(rows):
            ratio = 1.0 if rows[0]["sample"] == "full" else 0.5
            return choose_selection([row for row in rows if float(row["downsample_ratio"]) == ratio])
        with patch("tagforge.downsample.choose_selection", side_effect=force_ratio):
            run_pipeline(self.config, ["full", "subset"])
        output, samples = write_summary(self.config)
        self.assertEqual(samples, ["full", "subset"])
        tables = workbook_tables(output)
        self.assertEqual(list(tables), ["Summary", "Counts_raw", "Counts_final", "Downsample", "Warnings", "Results"])
        full = tables["Counts_raw"]
        selected = tables["Counts_final"]
        self.assertEqual(full[0], selected[0])
        self.assertEqual(full[1], selected[1])
        for row in full[1:]:
            self.assertEqual(sum(int(value) for value in row[1:]), 6)
        selection = load_selection(self.config, "subset")
        self.assertEqual(selection["ratio"], 0.5)
        self.assertEqual(sum(int(value) for value in selected[2][1:]), selection["metrics"]["umi_types"])
        self.assertNotEqual(full[2][1:], selected[2][1:])
        detail = self.config.output_dir / "subset" / selection["molecules"]["path"]
        next(detail.glob("*.parquet")).unlink()
        _, samples = write_summary(self.config)
        self.assertEqual(samples, ["full"])
        tables = workbook_tables(output)
        self.assertEqual(tables["Counts_final"], tables["Counts_raw"])

    def test_summary_command_updates_same_workbook_and_retires_duplicate_files(self):
        from tagforge.cli import main
        self.config.downsample_enabled = False
        run_pipeline(self.config, ["example"])
        directory = self.config.output_dir / "00_report"
        report = directory / "TagForge_batch_report.html"
        html = report.read_bytes()
        for name in ("00_summary.xlsx", "01_FB_not_in_anno.xlsx"):
            (directory / name).write_text("legacy")
        with patch("tagforge.cli.load_config", return_value=self.config):
            self.assertEqual(main(["summary", "--config", str(self.config.path)]), 0)
        self.assertEqual(report.read_bytes(), html)
        self.assertFalse((directory / "00_summary.xlsx").exists())
        self.assertFalse((directory / "01_FB_not_in_anno.xlsx").exists())
        tables = workbook_tables(directory / "TagForge_batch_report.xlsx")
        self.assertEqual(tables["Counts_final"], tables["Counts_raw"])
        self.assertEqual(tables["Downsample"][1][1:3], ["not run", ""])
        self.assertEqual(sum(int(value) for value in tables["Counts_raw"][1][1:]), 6)

    def test_merged_warnings_keep_annotation_details_above_warning_threshold(self):
        from tagforge.reports import write_summary
        self.config.downsample_enabled = False
        run_pipeline(self.config, ["example"])
        details = [{"barcode2_sequence": f"missing-{index}", "reads_count": 100 - index} for index in range(52)]
        with patch("tagforge.reports._annotation_rate_and_details", return_value=(0.95, details)):
            output, _ = write_summary(self.config)
        warnings = workbook_tables(output)["Warnings"]
        self.assertEqual(len(warnings), 51)
        self.assertEqual(warnings[1], ["example", "0.95", "missing-0", "100", "no"])

    def test_pb_cb_pipeline_reports_without_matrix(self):
        self.config.pb_cb_enabled = True
        self.config.barcode2_sequence_only = True
        self.config.fb_info = None
        self.config.segments = [replace(s, target_name="CB", whitelist=None,
            correction=CorrectionConfig(enabled=False)) if s.target == "barcode2" else s
            for s in self.config.segments]
        run_pipeline(self.config, ["example"])
        dirs = sample_dirs(self.config.output_dir, "example")
        self.assertFalse((dirs["matrix"] / "example.raw_matrix_summary.json").exists())
        self.assertFalse((dirs["matrix"] / "example.optimal_matrix_summary.json").exists())
        payload = json.loads((dirs["matrix"] / "example.pb_cb_summary.json").read_text())
        self.assertEqual(payload["supporting_reads"], 9)
        page = (dirs["report"] / "example.report.html").read_text()
        for expected in ("Barcode Validation", "PB Assignment Outcomes", "Multi-CB PB Pass Rate"):
            self.assertIn(expected, page)
        self.assertNotIn("Molecule Retention", page)
        self.assertNotIn("PBs per CB", page)
        self.assertNotIn("Mapping Summary", page)
        total_pbs = len({row[self.config.target_name("barcode1")] for row in
                         open_tsv(dirs["matrix"] / "example.pb_cb_mapping.tsv.gz")})
        self.assertEqual(payload["total_pbs"], total_pbs)
        rate = payload["retained_pbs"] / total_pbs
        self.assertIn(f'<div class="value">{total_pbs:,}</div>', page)
        self.assertIn(f'{rate:.1%}', page)
        overview = workbook_tables(dirs["report"] / "example.report.xlsx")["Overview"]
        self.assertEqual(int(overview[1][overview[0].index("Total PBs")]), total_pbs)
        self.assertAlmostEqual(float(overview[1][overview[0].index("PB pass rate")]), 100 * rate)
        batch_page = (self.config.output_dir / "00_report/TagForge_batch_report.html").read_text()
        self.assertIn("PB Assignment Outcomes", batch_page)
        self.assertNotIn("Molecule Retention", batch_page)
        self.assertTrue((dirs['report'] / 'example.dominance_qc.json').exists())
        self.assertIn('Dominance Distribution', workbook_tables(dirs['report'] / 'example.report.xlsx'))
        self.assertEqual(payload['dominance_qc']['retained_pbs'], payload['retained_pbs'])
        self.assertNotIn("Main results", page)
        self.assertEqual(_completed_summary_samples(self.config), ["example"])
        self.assertTrue((self.config.output_dir / "00_report/TagForge_batch_report.xlsx").is_file())
        self.assertTrue((self.config.output_dir / "00_report/TagForge_batch_report.html").is_file())
        generations = {p.name: p.read_text() for p in dirs["checkpoint"].glob("*.manifest.json")}
        run_pipeline(self.config, ["example"])
        self.assertEqual(generations, {p.name: p.read_text() for p in dirs["checkpoint"].glob("*.manifest.json")})

    def test_aliases_merge_shared_umi_and_preserve_reads(self):
        config = SimpleNamespace(output_dir=self.root / "alias", compression_level=3,
            pb_cb_enabled=True, pb_cb_dominance_threshold=0.8, raw={},
            target_name=lambda role: {"barcode1": "PB", "barcode2": "CB"}[role])
        dirs = sample_dirs(config.output_dir, "sample")
        rows = [{"PB": "PB1", "CB_name": cb, "corrected_umi": umi, "reads_count": reads, "raw_umi_count": 1}
                for cb, umi, reads in (("AACC", "AAAA", 20), ("AAAC", "AAAA", 1), ("AAAC", "CCCC", 1))]
        write_parquet(dirs["detail"] / "sample.molecule_detail.parquet", rows, list(rows[0]))
        result = pair_mapping_sample(config, "sample")
        mapping = list(open_tsv(result))
        self.assertEqual(len(mapping), 1)
        self.assertEqual((mapping[0]["CB"], mapping[0]["molecule_count"], mapping[0]["reads_count"]), ("AACC", "2", "22"))
        summary = json.loads((dirs["matrix"] / "sample.pb_cb_summary.json").read_text())
        self.assertEqual((summary["raw_molecules"], summary["resolved_molecules"], summary["supporting_reads"]), (3, 2, 22))

    def load_modified_example(self, original, replacement):
        repository = Path(__file__).parents[1]
        text = (repository / "configs/config.example.yaml").read_text()
        text = text.replace("workdir: ..", f"workdir: {repository}").replace(original, replacement)
        path = self.root / "modified.yaml"
        path.write_text(text)
        return load_config(path)

    def test_false_strings_are_parsed_as_booleans(self):
        config = self.load_modified_example("enabled: true", 'enabled: "false"')
        self.assertFalse(config.resume_enabled)
        self.assertFalse(config.downsample_enabled)
        self.assertFalse(config.quick_test_enabled)

    def test_parquet_compression_defaults_normalizes_and_rejects_unknown_codecs(self):
        config = self.load_modified_example("  parquet_compression: snappy  # snappy (default) or zstd\n", "")
        self.assertEqual(config.parquet_compression, "snappy")
        config = self.load_modified_example("parquet_compression: snappy", "parquet_compression: ZSTD")
        self.assertEqual(config.parquet_compression, "zstd")
        for codec in ("gzip", "null", "17"):
            with self.subTest(codec=codec):
                with self.assertRaisesRegex(ConfigError, "parquet_compression must be"):
                    self.load_modified_example("parquet_compression: snappy", "parquet_compression: " + codec)

    def test_fractional_coordinates_resources_and_nonfinite_alias_ratio_are_rejected(self):
        changes = (("length: 8", "length: 7.5"), ("start: 0", "start: 0.5"),
                   ("threads: 2", "threads: 2.5"), ("chunk_size: 10000", "chunk_size: 3.5"),
                   ("downsample:", "downsample:\n  repeats: 1"), ("random_seed: 12345", "random_seed: .inf"),
                   ("reads: 10000", "reads: true"),
                   ("performance:", "pb_cb:\n  min_parent_reads: 2.5\nperformance:"),
                   ("performance:", "pb_cb:\n  parent_child_ratio: .nan\nperformance:"))
        for original, replacement in changes:
            with self.subTest(value=replacement):
                with self.assertRaises(ConfigError):
                    self.load_modified_example(original, replacement)
    def test_inapplicable_umi_settings_and_conflicting_segment_roles_are_rejected(self):
        cases = [("segment: UMI", "segment: UMI\n      whitelist: ignored.txt", "UMI segments"),
                 ("segment: UMI", "segment: UMI\n      target: barcode1", "parent role"),
                 ("enabled: false", "enabled: true", "UMI segments")]
        for old, new, message in cases:
            with self.subTest(value=new), self.assertRaisesRegex(ConfigError, message):
                self.load_modified_example(old, new)

    def test_matrix_export_rejects_stale_config_and_preserves_existing_output(self):
        import io
        from tagforge.cli import main
        self.config.downsample_enabled = False
        run_pipeline(self.config, ["example"])
        output = self.root / "export.tsv"
        command = ["export-matrix", "--config", "mocked.yaml", "--sample", "example",
                   "--format", "tsv", "--out", str(output)]
        with patch("tagforge.cli.load_config", return_value=self.config), patch("sys.stdout", io.StringIO()):
            self.assertEqual(main(command), 0)
            original = output.read_bytes()
            self.config.umi_max_distance = 0
            with patch("sys.stderr", io.StringIO()) as error:
                self.assertEqual(main(command), 2)
                self.assertIn("current matrix stage", error.getvalue())
            self.assertEqual(output.read_bytes(), original)
            run_pipeline(self.config, ["example"])
            self.assertEqual(main(command), 0)
        self.assertEqual(sum(int(value) for line in output.read_text().splitlines()[1:]
                             for value in line.split("\t")[1:]), 7)

    def test_cli_parquet_export_uses_configured_codec(self):
        import io
        import pyarrow.parquet as pq
        from tagforge.cli import main
        self.config.downsample_enabled = False
        for codec in ("snappy", "zstd"):
            with self.subTest(codec=codec):
                self.config.parquet_compression = codec
                run_pipeline(self.config, ["example"])
                output = self.root / f"export-{codec}.parquet"
                command = ["export-matrix", "--config", "mocked.yaml", "--sample", "example",
                           "--format", "parquet", "--out", str(output)]
                with patch("tagforge.cli.load_config", return_value=self.config), patch("sys.stdout", io.StringIO()):
                    self.assertEqual(main(command), 0)
                table = pq.read_table(output)
                self.assertEqual(sum(table["molecule_count"].to_pylist()), 6)
                metadata = pq.ParquetFile(output).metadata
                for i in range(metadata.num_row_groups):
                    for j in range(metadata.num_columns):
                        self.assertEqual(metadata.row_group(i).column(j).compression, codec.upper())

    def test_matrix_migration_applies_and_validates_cli_thread_override(self):
        import io
        from tagforge.cli import main
        path = Path(__file__).parents[1] / "configs/config.example.yaml"
        with patch("tagforge.cli.run_step") as stage, patch("sys.stdout", io.StringIO()):
            self.assertEqual(main(["migrate-parquet", "--config", str(path), "--threads", "3"]), 0)
            self.assertEqual(stage.call_args.args[0].threads, 3)
        with patch("tagforge.cli.run_step") as stage, patch("sys.stderr", io.StringIO()):
            self.assertEqual(main(["migrate-parquet", "--config", str(path), "--threads", "0"]), 2)
            stage.assert_not_called()

    def test_unknown_negative_and_conflicting_configuration_is_rejected(self):
        cases = [("chunk_size: 10000", "chunk_szie: 10000", "Unknown"),
                 ("max_shift: 1", "max_shift: -1", "nonnegative"),
                 ("max_distance: 1", "max_distance: -1", "nonnegative"),
                 ("max_shift: 1", "max_shift: -0.5", "integer"),
                 ("name: UMI", "name: CELL", "conflict"),
                 ("name: CELL", "name: reads_count", "conflict")]
        for old, new, error in cases:
            with self.subTest(new=new), self.assertRaisesRegex(ConfigError, error):
                self.load_modified_example(old, new)

    def test_wrong_whitelist_length_is_rejected_before_run(self):
        repository = Path(__file__).parents[1]
        wrong = self.root / "wrong.txt"
        wrong.write_text("AA\n")
        with self.assertRaisesRegex(ConfigError, "whitelist entries.*length 8"):
            self.load_modified_example("tests/fixtures/small_fastq/config/CELL_WL.txt", str(wrong))

    def test_zero_valid_pipeline_has_readable_empty_matrix_and_report(self):
        sample = self.config.samples[0]
        for path, sequence in ((sample.r1, "AA"), (sample.r2, "TT")):
            path.write_bytes(gzip.compress(f"@invalid/1\n{sequence}\n+\nII\n".encode(), mtime=0))
        # Pair identifiers strip /1 or /2; both mates intentionally have one ID.
        run_pipeline(self.config, ["example"])
        dirs = sample_dirs(self.config.output_dir, "example")
        self.assertEqual(list(open_records(dirs["detail"] / "example.molecule_detail.parquet")), [])
        self.assertEqual(list(open_records(dirs["matrix"] / "example.raw_molecule_pairs.parquet")), [])
        self.assertTrue((dirs["report"] / "example.report.xlsx").is_file())
        self.assertTrue(stage_current(self.config, "example", "report"))


if __name__ == "__main__":
    unittest.main()
