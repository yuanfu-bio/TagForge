import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest

import pyarrow.parquet as pq

from tagforge.columnar import open_records, write_parquet
from tagforge.sparse_matrix import build_sparse_matrix, export_matrix


class SparseMatrixTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def build(self, empty=False):
        source = self.root / "molecules.parquet"
        rows = [] if empty else [{"PB": pb, "FB_name": feature} for pb, feature in (
            ("0001", "A"), ("0001", "A"), ("0001", "B"), ("0002", "B"))]
        write_parquet(source, rows, ["PB", "FB_name"])
        return build_sparse_matrix(source, self.root, "sample", "PB", "FB_name", threads=2, memory_mb=64)

    def test_aggregation_and_all_export_formats_preserve_counts(self):
        pairs, summary = self.build()
        payload = json.loads(summary.read_text())
        self.assertEqual(payload["feature_counts"], {"A": 2, "B": 2})
        self.assertEqual(payload["barcode_count"], 2)
        self.assertEqual(payload["pi_seq"]["molecule_bins"][:3], [1, 0, 1])
        expected = [("0001", "A", "2"), ("0001", "B", "1"), ("0002", "B", "1")]
        self.assertEqual(sorted((r["barcode1"], r["feature"], r["molecule_count"]) for r in open_records(pairs)), expected)
        for name in ("matrix.tsv.gz", "matrix.tsv"):
            output = self.root / name
            export_matrix(pairs, output, "tsv", "PB", memory_mb=64)
            opener = gzip.open if name.endswith(".gz") else open
            with opener(output, "rt") as handle:
                self.assertEqual(list(csv.reader(handle, delimiter="\t")), [["PB", "A", "B"], ["0001", "2", "1"], ["0002", "0", "1"]])
        output = self.root / "matrix.mtx.gz"
        _, barcodes, features = export_matrix(pairs, output, "mtx", memory_mb=64)
        with gzip.open(output, "rt") as handle:
            self.assertEqual(handle.read().splitlines()[2:], ["2 2 3", "1 1 2", "1 2 1", "2 2 1"])
        with gzip.open(barcodes, "rt") as handle:
            self.assertEqual(handle.read().splitlines(), ["0001", "0002"])
        with gzip.open(features, "rt") as handle:
            self.assertEqual(handle.read().splitlines(), ["A", "B"])
        parquet = self.root / "export.parquet"
        export_matrix(pairs, parquet, "parquet", memory_mb=64)
        self.assertEqual(sorted((r["barcode1"], r["feature"], r["molecule_count"]) for r in open_records(parquet)), expected)
        for path in (pairs, parquet):
            metadata = pq.ParquetFile(path).metadata
            for i in range(metadata.num_row_groups):
                for j in range(metadata.num_columns):
                    self.assertEqual(metadata.row_group(i).column(j).compression, "SNAPPY")

    def test_zstd_matrix_build_and_export_preserve_codec_and_counts(self):
        self.build()
        pairs, _ = build_sparse_matrix(self.root / "molecules.parquet", self.root / "zstd", "sample",
                                       "PB", "FB_name", memory_mb=64, compression="zstd")
        output = self.root / "zstd-export.parquet"
        export_matrix(pairs, output, "parquet", memory_mb=64, compression="zstd")
        self.assertEqual(list(open_records(pairs)), list(open_records(output)))
        self.assertEqual(sum(int(row["molecule_count"]) for row in open_records(output)), 4)
        for path in (pairs, output):
            metadata = pq.ParquetFile(path).metadata
            for i in range(metadata.num_row_groups):
                for j in range(metadata.num_columns):
                    self.assertEqual(metadata.row_group(i).column(j).compression, "ZSTD")


    def test_empty_matrix_and_exports_are_valid(self):
        pairs, summary = self.build(empty=True)
        self.assertEqual(json.loads(summary.read_text())["barcode_count"], 0)
        for fmt, name in (("tsv", "empty.tsv"), ("mtx", "empty.mtx")):
            output = self.root / name
            export_matrix(pairs, output, fmt, memory_mb=64)
            self.assertIn("0 0 0" if fmt == "mtx" else "Barcode1", output.read_text())

    def test_legacy_numeric_barcodes_keep_leading_zeroes(self):
        source = self.root / "legacy.tsv"
        source.write_text("PB\tFB_name\n0001\tA\n0002\tA\n")
        pairs, _ = build_sparse_matrix(source, self.root, "legacy", "PB", "FB_name", memory_mb=64)
        self.assertEqual({r["barcode1"] for r in open_records(pairs)}, {"0001", "0002"})


if __name__ == "__main__":
    unittest.main()
