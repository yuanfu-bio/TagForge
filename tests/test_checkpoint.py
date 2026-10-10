from __future__ import annotations

from dataclasses import replace
from contextlib import contextmanager
import gzip
from itertools import product
from pathlib import Path
import json
import os
import shutil
import tempfile
import unittest
from unittest.mock import patch

from tagforge.checkpoint import STAGE_NUMBERS
from tagforge.config import load_config
from tagforge.columnar import open_records
from tagforge.external_tools import ExternalToolVersions
from tagforge.io_utils import open_tsv
from tagforge.pipeline import run_step


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        repository = Path(__file__).parents[1]
        fixture = self.root / "fixture"
        shutil.copytree(repository / "tests/fixtures/small_fastq", fixture)
        self.config = load_config(repository / "configs/config.example.yaml")
        self.config.workdir = self.root
        self.config.output_dir = self.root / "output"
        self.config.samples = [replace(self.config.samples[0],
            r1=fixture / "raw/example/example_raw_1.fq.gz",
            r2=fixture / "raw/example/example_raw_2.fq.gz")]
        self.config.segments = [replace(segment, whitelist=fixture / "config" / segment.whitelist.name)
                                if segment.whitelist else segment for segment in self.config.segments]
        self.config.fb_info = fixture / "config/FB_info.tsv"
        self.config.threads = 1
        self.config.chunk_size = 3
        self.addCleanup(patch.stopall)
        patch("tagforge.pipeline.check_external_tools", return_value=ExternalToolVersions("test", "test")).start()

    def run_through(self, last):
        for step in ("extract", "correct", "dedup", "matrix", "downsample", "report"):
            run_step(self.config, "example", step)
            if step == last:
                break

    def test_codec_change_invalidates_parquet_stages_and_recomputes_outputs(self):
        import pyarrow.parquet as pq
        from tagforge.checkpoint import stage_current
        self.run_through("matrix")
        self.config.parquet_compression = "zstd"
        for step in ("extract", "correct", "dedup", "matrix"):
            self.assertFalse(stage_current(self.config, "example", step))
        self.run_through("matrix")
        for path in (self.config.output_dir / "example").rglob("*.parquet"):
            if not path.is_file():
                continue
            metadata = pq.ParquetFile(path).metadata
            for i in range(metadata.num_row_groups):
                for j in range(metadata.num_columns):
                    self.assertEqual(metadata.row_group(i).column(j).compression, "ZSTD")
        self.config.compression_level = 1
        self.assertFalse(stage_current(self.config, "example", "extract"))

    def test_whitelist_change_invalidates_only_affected_stages(self):
        from tagforge.checkpoint import stage_current
        self.run_through("matrix")
        self.config.threads = 2
        self.assertEqual(run_step(self.config, "example", "extract")[1], "skipped")
        self.assertEqual(run_step(self.config, "example", "correct")[1], "skipped")
        whitelist = next(s.whitelist for s in self.config.segments if s.target == "barcode1")
        whitelist.write_text(whitelist.read_text() + "\nAAAAAAAA\n")
        self.assertTrue(stage_current(self.config, "example", "extract"))
        self.assertFalse(stage_current(self.config, "example", "correct"))
        with self.assertRaisesRegex(ValueError, "current.*correct"):
            run_step(self.config, "example", "dedup")
        run_step(self.config, "example", "correct")
        for step in ("dedup", "matrix", "downsample", "report"):
            self.assertFalse((self.config.output_dir / "example/00_logs/checkpoints" / f"{STAGE_NUMBERS[step]}_{step}.manifest.json").exists())

    def test_content_checksum_detects_annotation_change_with_same_metadata(self):
        from tagforge.checkpoint import stage_current
        self.run_through("correct")
        annotation = self.config.fb_info
        stat = annotation.stat()
        annotation.write_text(annotation.read_text().replace("CD4", "CD3"))
        os.utime(annotation, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        self.assertTrue(stage_current(self.config, "example", "extract"))
        self.assertFalse(stage_current(self.config, "example", "correct"))

    def test_coordinate_change_invalidates_extraction_and_downstream(self):
        from tagforge.checkpoint import stage_current
        self.run_through("matrix")
        self.config.segments = [replace(s, start=s.start + 1) if s.target == "umi" else s
                                for s in self.config.segments]
        for step in ("extract", "correct", "dedup", "matrix"):
            self.assertFalse(stage_current(self.config, "example", step))
        with self.assertRaisesRegex(ValueError, "current.*extract"):
            run_step(self.config, "example", "correct")
        self.assertEqual(run_step(self.config, "example", "extract")[1], "completed")
        self.assertTrue(stage_current(self.config, "example", "extract"))
        self.assertFalse(stage_current(self.config, "example", "correct"))

    def test_umi_settings_change_does_not_repeat_correction(self):
        from tagforge.checkpoint import stage_current
        self.run_through("matrix")
        self.config.umi_max_distance = 0
        self.assertTrue(stage_current(self.config, "example", "correct"))
        self.assertFalse(stage_current(self.config, "example", "dedup"))
        self.assertFalse(stage_current(self.config, "example", "matrix"))

    def test_final_correction_stats_failure_can_resume(self):
        from tagforge.barcode_correct import write_json as real_write
        self.run_through("extract")
        def fail_stats(path, *args, **kwargs):
            if path.name.endswith(".correction_qc.json"):
                raise RuntimeError("correction stats interrupted")
            return real_write(path, *args, **kwargs)
        with patch("tagforge.barcode_correct.write_json", side_effect=fail_stats):
            with self.assertRaisesRegex(RuntimeError, "correction stats interrupted"):
                run_step(self.config, "example", "correct")
        # An interrupted old-layout stage must keep its committed batches.
        sample_root = self.config.output_dir / "example"
        work = sample_root / "00_logs/checkpoints/.work/stage-correct"
        legacy = sample_root / "08_tmp/stage-correct"
        legacy.parent.mkdir()
        work.rename(legacy)
        work.parent.rmdir()
        with patch("tagforge.barcode_correct._correct_batch", side_effect=AssertionError("committed read repeated")):
            run_step(self.config, "example", "correct")
        self.assertFalse(legacy.parent.exists())
        self.assertFalse(work.parent.exists())
        valid = self.config.output_dir / "example/02_corrected/example.valid_reads.parquet"
        self.assertEqual(len(list(open_records(valid))), 9)

    def test_missing_or_corrupt_fragment_invalidates_downstream(self):
        from tagforge.checkpoint import stage_current
        self.run_through("correct")
        dataset = self.config.output_dir / "example/01_extracted/example.extracted.parquet"
        next(dataset.glob("*.parquet")).write_bytes(b"broken parquet")
        self.assertFalse(stage_current(self.config, "example", "extract"))
        self.assertFalse(stage_current(self.config, "example", "correct"))
        run_step(self.config, "example", "extract")
        self.assertEqual(len(list(open_records(dataset))), 9)

    def test_final_extraction_stats_failure_can_resume(self):
        from tagforge.extract import write_tsv as real_write
        def fail_stats(path, *args, **kwargs):
            if path.name.endswith(".extraction_stats.tsv"):
                raise RuntimeError("stats interrupted")
            return real_write(path, *args, **kwargs)
        with patch("tagforge.extract.write_tsv", side_effect=fail_stats):
            with self.assertRaisesRegex(RuntimeError, "stats interrupted"):
                run_step(self.config, "example", "extract")
        self.assertFalse((self.config.output_dir / "example/00_logs/checkpoints/01_extract.manifest.json").exists())
        self.assertTrue((self.config.output_dir / "example/00_logs/checkpoints/example.extract.resume.json").exists())
        with patch("tagforge.extract._extract_record", side_effect=AssertionError("committed read repeated")):
            run_step(self.config, "example", "extract")
        self.assertEqual(len(list(open_records(self.config.output_dir / "example/01_extracted/example.extracted.parquet"))), 9)

    def test_uncommitted_fragment_and_preview_tail_are_discarded(self):
        self.config.extraction_preview_reads = 1000
        from tagforge.extract import atomic_text as real_atomic
        calls = 0
        @contextmanager
        def interrupt_second_commit(path, *args, **kwargs):
            nonlocal calls
            with real_atomic(path, *args, **kwargs) as handle:
                yield handle
                if path.name.endswith(".extract.resume.json"):
                    calls += 1
                    if calls == 3:
                        raise RuntimeError("commit interrupted")
        with patch("tagforge.extract.atomic_text", side_effect=interrupt_second_commit):
            with self.assertRaisesRegex(RuntimeError, "commit interrupted"):
                run_step(self.config, "example", "extract")
        run_step(self.config, "example", "extract")
        base = self.config.output_dir / "example/01_extracted"
        self.assertEqual(len(list(open_records(base / "example.extracted.parquet"))), 9)
        preview = list(open_tsv(base / "example.extracted.preview.tsv"))
        self.assertEqual(len(preview), 9)
        self.assertEqual(len({row["read_id"] for row in preview}), 9)

    def test_summary_excludes_stale_sample_before_upstream_recompute(self):
        from tagforge.reports import write_summary
        self.run_through("downsample")
        self.assertEqual(write_summary(self.config)[1], ["example"])
        self.config.umi_max_distance = 0
        self.assertEqual(write_summary(self.config)[1], [])

    def test_strong_fastq_checksum_detects_same_metadata_change(self):
        from tagforge.checkpoint import stage_current
        self.config.raw.setdefault("resume", {})["input_fingerprint"] = "sha256"
        self.run_through("extract")
        fastq = self.config.samples[0].r1
        stat = fastq.stat()
        data = bytearray(fastq.read_bytes())
        data[4] ^= 1  # Change only the valid gzip timestamp header.
        fastq.write_bytes(data)
        os.utime(fastq, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        self.assertFalse(stage_current(self.config, "example", "extract"))

    def test_interrupted_publication_recovers_without_recomputing(self):
        from tagforge import checkpoint
        real_publish = checkpoint._publish_artifact
        count = 0
        def fail_second(*args, **kwargs):
            nonlocal count
            count += 1
            if count == 2:
                raise RuntimeError("publication interrupted")
            return real_publish(*args, **kwargs)
        with patch("tagforge.checkpoint._publish_artifact", side_effect=fail_second):
            with self.assertRaisesRegex(RuntimeError, "publication interrupted"):
                run_step(self.config, "example", "extract")
        sample_root = self.config.output_dir / "example"
        work = sample_root / "00_logs/checkpoints/.work/stage-extract"
        self.assertTrue((work / "publication.json").is_file())
        legacy = sample_root / "08_tmp/stage-extract"
        legacy.parent.mkdir()
        work.rename(legacy)
        work.parent.rmdir()
        self.assertFalse(checkpoint.stage_current(self.config, "example", "extract"))
        with patch("tagforge.pipeline.extract_sample", side_effect=AssertionError("recomputed")):
            run_step(self.config, "example", "extract")
        self.assertTrue(checkpoint.stage_current(self.config, "example", "extract"))
        self.assertFalse(legacy.parent.exists())
        self.assertFalse(work.parent.exists())
        self.assertFalse((self.config.output_dir / "example/00_logs/checkpoints/example.extract.resume.json").exists())

    def test_legacy_checkpoint_is_recomputed(self):
        self.run_through("extract")
        checkpoint = self.config.output_dir / "example/00_logs/checkpoints/01_extract.manifest.json"
        checkpoint.unlink()
        legacy = checkpoint.with_name("extract.done")
        legacy.write_text("tagforge_version=0.1.13\n")
        self.assertEqual(run_step(self.config, "example", "extract")[1], "completed")
        self.assertEqual(json.loads(checkpoint.read_text())["step"], "extract")
        self.assertFalse(legacy.exists())

    def test_json_done_migrates_without_recompute_and_bad_current_manifest_blocks_fallback(self):
        from tagforge.checkpoint import stage_current
        self.run_through("matrix")
        directory = self.config.output_dir / "example/00_logs/checkpoints"
        generations = {}
        for step in ("extract", "correct", "dedup", "matrix"):
            current = directory / f"{STAGE_NUMBERS[step]}_{step}.manifest.json"
            generations[step] = json.loads(current.read_text())["generation"]
            current.rename(directory / f"{step}.done")
        self.assertTrue(stage_current(self.config, "example", "matrix"))
        for step in generations:
            with patch.dict("tagforge.pipeline.STEP_FUNCS", {step: lambda *args, **kwargs:
                            self.fail("valid old manifest caused recomputation")}):
                self.assertEqual(run_step(self.config, "example", step)[1], "skipped")
            current = directory / f"{STAGE_NUMBERS[step]}_{step}.manifest.json"
            self.assertEqual(json.loads(current.read_text())["generation"], generations[step])
            self.assertFalse((directory / f"{step}.done").exists())
        current = directory / "04_matrix.manifest.json"
        legacy = directory / "matrix.done"
        shutil.copy2(current, legacy)
        current.write_text("invalid JSON")
        self.assertFalse(stage_current(self.config, "example", "matrix"))
        self.assertEqual(run_step(self.config, "example", "matrix")[1], "completed")
        self.assertTrue(stage_current(self.config, "example", "matrix"))
        self.assertFalse(legacy.exists())

    def test_matrix_publication_removes_legacy_table_that_can_shadow_counts(self):
        self.run_through("dedup")
        legacy = self.config.output_dir / "example/04_matrix/example.raw_count_matrix.tsv.gz"
        with gzip.open(legacy, "wt") as handle:
            handle.write("CELL\tCD4\nold\t999\n")
        run_step(self.config, "example", "matrix")
        self.assertFalse(legacy.exists())
        from tagforge.reports import _matrix_summary
        self.assertEqual(sum(_matrix_summary(legacy)[1].values()), 6)

    def test_public_api_returns_published_paths(self):
        self.run_through("correct")
        (molecules, _), _ = run_step(self.config, "example", "dedup")
        self.assertIsInstance(molecules, Path)
        self.assertTrue(molecules.is_dir())
        (pairs, summary), _ = run_step(self.config, "example", "matrix")
        self.assertTrue(pairs.is_file())
        self.assertTrue(summary.is_file())

    def test_same_sample_lock_rejects_second_writer(self):
        from tagforge.checkpoint import sample_lock
        with sample_lock(self.config, "example"):
            with self.assertRaisesRegex(ValueError, "already running"):
                run_step(self.config, "example", "extract")

    def test_full_depth_selection_replaces_old_optimal_reference(self):
        from tagforge.selected import load_selection
        self.run_through("downsample")
        first = load_selection(self.config, "example")
        self.assertLess(first["ratio"], 1)
        self.config.downsample_ratios = [1.0]
        run_step(self.config, "example", "downsample")
        selected = load_selection(self.config, "example")
        self.assertEqual(selected["ratio"], 1)
        self.assertEqual(selected["molecules"]["path"], "03_molecules/example.molecule_detail.parquet")
        self.assertFalse((self.config.output_dir / "example/04_matrix/example.optimal_matrix_summary.json").exists())
        self.assertFalse((self.config.output_dir / "example/03_molecules/example.optimal_saturation_molecule_detail.parquet").exists())
        from tagforge.selected import selection_point
        point = selection_point(self.config, "example")
        self.assertEqual(selected["metrics"]["umi_types"], int(point["umi_types"]))

    def test_selection_json_counts_must_match_the_raw_result(self):
        from tagforge.selected import load_selection, selection_path
        self.run_through("downsample")
        path = selection_path(self.config, "example")
        selected = json.loads(path.read_text())
        selected["metrics"]["reads_sampled"] += 1
        path.write_text(json.dumps(selected))
        with self.assertRaisesRegex(ValueError, "selected raw result"):
            load_selection(self.config, "example")

    def test_selected_fragment_loss_requires_regeneration(self):
        from tagforge.checkpoint import stage_current
        from tagforge.selected import load_selection
        self.run_through("downsample")
        selection = load_selection(self.config, "example")
        dataset = self.config.output_dir / "example" / selection["molecules"]["path"]
        next(dataset.glob("*.parquet")).unlink()
        self.assertFalse(stage_current(self.config, "example", "downsample"))
        run_step(self.config, "example", "downsample")
        self.assertTrue(stage_current(self.config, "example", "downsample"))

    def test_zero_selection_is_a_valid_empty_dataset(self):
        from tagforge.selected import load_selection
        fastq = self.config.samples[0].r2
        with gzip.open(fastq, "rt") as handle:
            lines = handle.readlines()
        sequences = ("".join(bases) for bases in product("ACGT", repeat=6))
        for index in range(1, len(lines), 4):
            lines[index] = next(sequences) + lines[index][6:]
        with gzip.open(fastq, "wt") as handle:
            handle.writelines(lines)
        self.config.umi_method = "unique"
        self.config.downsample_ratios = [0]
        # All molecules are genuine singletons, so the legacy criterion ties
        # at zero and full depth and deterministically selects zero.
        self.run_through("downsample")
        selection = load_selection(self.config, "example")
        self.assertEqual(selection["ratio"], 0)
        self.assertEqual(selection["metrics"]["umi_types"], 0)
        dataset = self.config.output_dir / "example" / selection["molecules"]["path"]
        self.assertEqual(list(open_records(dataset)), [])

    def test_full_depth_then_subset_and_new_seed_have_consistent_counts(self):
        from tagforge.selected import load_selection
        self.config.downsample_ratios = [1]
        self.run_through("downsample")
        self.assertEqual(load_selection(self.config, "example")["ratio"], 1)
        self.config.downsample_ratios = [0.7]
        self.config.downsample_seed = 12345
        run_step(self.config, "example", "downsample")
        selected = load_selection(self.config, "example")
        self.assertLess(selected["ratio"], 1)
        self.config.downsample_seed = 42
        run_step(self.config, "example", "downsample")
        selected = load_selection(self.config, "example")
        rows = list(open_records(self.config.output_dir / "example" / selected["molecules"]["path"]))
        count_field = "reads_count" if selected["ratio"] == 1 else "reads_count_at_optimal_downsample"
        self.assertEqual(sum(int(row[count_field]) for row in rows), selected["metrics"]["reads_sampled"])


if __name__ == "__main__":
    unittest.main()
