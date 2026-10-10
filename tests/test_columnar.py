from pathlib import Path
import tempfile
import unittest

import pyarrow.parquet as pq

from tagforge.barcode_correct import correct_sample
from tagforge.columnar import (ParquetFragmentWriter, open_records, resolve_table,
                              table_schema, tsv_to_parquet, write_parquet)
from tagforge.config import (ConfigError, CorrectionConfig, SampleConfig,
                            SegmentConfig, TagForgeConfig)
from tagforge.extract import _extract_record
from tagforge.fastq import PairedRead
from tagforge.io_utils import open_tsv, sample_dirs
from tagforge.reports import _barcode2_annotation_gap


class ColumnarTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_file_fragments_and_empty_table_have_same_typed_schema(self):
        fields = ["PB", "reads_count", "raw_umi_count"]
        rows = [{"PB": f"PB{i}", "reads_count": str(i + 1), "raw_umi_count": 1} for i in range(5)]
        single = self.root / "single.parquet"
        write_parquet(single, rows, fields, batch_size=2)
        fragments = self.root / "parts.parquet"
        writer = ParquetFragmentWriter(fragments, fields, batch_size=2)
        writer.writerows(rows[:2]); writer.writerows(rows[2:4]); writer.writerows(rows[4:]); writer.close()
        self.assertEqual(list(open_records(single)), list(open_records(fragments)))
        empty = self.root / "empty.parquet"
        write_parquet(empty, [], fields)
        for path in (single, empty, next(fragments.glob("*.parquet"))):
            self.assertEqual(pq.ParquetFile(path).schema_arrow, table_schema(fields))
            self.assertEqual(str(pq.ParquetFile(path).schema_arrow.field("reads_count").type), "int64")
            metadata = pq.ParquetFile(path).metadata
            for i in range(metadata.num_row_groups):
                for j in range(metadata.num_columns):
                    self.assertEqual(metadata.row_group(i).column(j).compression, "SNAPPY")
        self.assertEqual(list(open_records(empty)), [])

    def test_zstd_single_fragments_and_empty_tables_preserve_codec_and_counts(self):
        fields = ["PB", "reads_count"]
        rows = [{"PB": "0001", "reads_count": 3}, {"PB": "0002", "reads_count": 2}]
        single = self.root / "zstd.parquet"
        write_parquet(single, rows, fields, compression="zstd")
        fragments = self.root / "zstd-parts.parquet"
        writer = ParquetFragmentWriter(fragments, fields, batch_size=1, compression="zstd")
        for row in rows:
            writer.writerow(row)
        writer.close()
        self.assertEqual(list(open_records(single)), list(open_records(fragments)))
        empty = self.root / "zstd-empty.parquet"
        writer = ParquetFragmentWriter(empty, fields, compression="zstd")
        writer.close()
        for path in [single, *fragments.glob("*.parquet"), *empty.glob("*.parquet")]:
            metadata = pq.ParquetFile(path).metadata
            for i in range(metadata.num_row_groups):
                for j in range(metadata.num_columns):
                    self.assertEqual(metadata.row_group(i).column(j).compression, "ZSTD")

    def test_legacy_migration_keeps_numeric_counts_and_empty_header(self):
        legacy = self.root / "legacy.tsv"
        legacy.write_text("PB\treads_count\nPB1\t3\n")
        migrated = tsv_to_parquet(legacy, batch_size=1)
        self.assertEqual(list(open_records(migrated)), [{"PB": "PB1", "reads_count": "3"}])
        legacy.write_text("PB\treads_count\n")
        migrated = tsv_to_parquet(legacy)
        self.assertEqual(pq.ParquetFile(migrated).metadata.num_rows, 0)
        self.assertEqual(pq.ParquetFile(migrated).schema_arrow, table_schema(["PB", "reads_count"]))

    def test_fractional_nonfinite_and_boolean_counts_are_rejected(self):
        path = self.root / "invalid.parquet"
        for value in (2.5, True, float("inf"), float("nan")):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ConfigError, "Invalid integer.*reads_count"):
                    write_parquet(path, [{"reads_count": value}], ["reads_count"])
                self.assertFalse(path.exists())
        write_parquet(path, [{"reads_count": 2.0}], ["reads_count"])
        self.assertEqual(pq.read_table(path)["reads_count"].to_pylist(), [2])

    def test_incomplete_current_table_does_not_fall_back_to_legacy(self):
        current = self.root / "current.parquet"
        current.mkdir()
        legacy = self.root / "legacy.tsv"
        legacy.write_text("PB\nPB1\n")
        with self.assertRaisesRegex(ConfigError, "Incomplete"):
            resolve_table(current, legacy)
        current.rmdir()
        self.assertEqual(resolve_table(current, legacy), legacy)

    def test_schema_mismatch_is_rejected_on_read_and_resume(self):
        dataset = self.root / "data.parquet"
        dataset.mkdir()
        write_parquet(dataset / "part-00000000.parquet", [{"PB": "AA"}], ["PB"])
        write_parquet(dataset / "part-00000001.parquet", [{"CB": "AA"}], ["CB"])
        with self.assertRaisesRegex(ConfigError, "inconsistent schemas"):
            list(open_records(dataset))
        with self.assertRaisesRegex(ConfigError, "schema does not match"):
            ParquetFragmentWriter(dataset, ["CB"], resume=True, parts=1)

    def annotation_run(self, trace_enabled):
        output = self.root / ("trace" if trace_enabled else "no_trace")
        annotation = self.root / "annotation.tsv"
        annotation.write_text("FB_ID\tsequence\tantibody_name\nF1\tAAGG\tfeature\n")
        disabled = CorrectionConfig(enabled=False)
        segments = [SegmentConfig(name, role, read, "fixed", length, start=start, correction=disabled)
                    for name, role, read, length, start in (
                        ("PB", "barcode1", "R1", 4, 0), ("FB1", "barcode2", "R1", 2, 0),
                        ("FB2", "barcode2", "R2", 2, 2), ("UMI", "umi", "R2", 2, 0))]
        config = TagForgeConfig(self.root / "config.yaml", self.root, output,
            [SampleConfig("sample", self.root / "R1", self.root / "R2")], segments, annotation,
            chunk_size=2, trace_enabled=trace_enabled)
        dirs = sample_dirs(output, "sample")
        columns = {role: config.segment_column(role) for role in ("barcode1", "barcode2", "umi")}
        rows = [_extract_record(PairedRead(f"read{i}", "AACC", "IIII", sequence, "IIII"),
                                config.segments, columns)[0]
                for i, sequence in enumerate(("TTGG", "TTCC", "AACC"))]
        write_parquet(dirs["extracted"] / "sample.extracted.parquet", rows, list(rows[0]))
        correct_sample(config, "sample", resume=False)
        from tagforge.barcode_correct import load_correction_stats
        stats = load_correction_stats(dirs["corrected"] / "sample.correction_qc.json")
        result = _barcode2_annotation_gap(config, "sample", stats, segments[1:3])
        return config, dirs, stats, result

    def test_annotation_counts_full_feature_with_and_without_trace(self):
        results = [self.annotation_run(enabled)[3] for enabled in (True, False)]
        self.assertEqual(results[0], results[1])
        self.assertEqual(results[0], (3, 2, [{"barcode2_sequence": "AACC", "reads_count": "2"}]))

    def test_legacy_directory_trace_recovers_multisegment_annotation_counts(self):
        config, dirs, stats, _ = self.annotation_run(True)
        (dirs["corrected"] / "sample.correction_qc.json").unlink()
        (dirs["corrected"] / "sample.barcode2_not_in_annotation.tsv.gz").unlink()
        result = _barcode2_annotation_gap(config, "sample", stats, config.segments[1:3])
        self.assertEqual(result, (3, 2, [{"barcode2_sequence": "AACC", "reads_count": 2}]))


if __name__ == "__main__":
    unittest.main()
