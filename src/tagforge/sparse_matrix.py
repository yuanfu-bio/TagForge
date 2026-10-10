"""Sparse matrix materialisation and bounded-memory exports."""
from __future__ import annotations

import csv
import gzip
import heapq
import json
from collections import Counter
from pathlib import Path
import tempfile

from .columnar import sql_reader
from .config import ConfigError, validate_parquet_compression
from .io_utils import atomic_text


def pair_path(matrix_dir: Path, sample: str, prefix: str = "raw") -> Path:
    return matrix_dir / f"{sample}.{prefix}_molecule_pairs.parquet"


def summary_path(matrix_dir: Path, sample: str, prefix: str = "raw") -> Path:
    return matrix_dir / f"{sample}.{prefix}_matrix_summary.json"


def _duckdb():
    try:
        import duckdb
    except ImportError as exc:
        raise ConfigError("Sparse matrices require duckdb") from exc
    return duckdb


def _connect(temporary, *, threads=1, memory_mb=512):
    con = _duckdb().connect()
    try:
        con.execute(f"SET threads={int(threads)}")
        con.execute(f"SET memory_limit='{int(memory_mb)}MB'")
        con.execute("SET temp_directory='" + str(temporary).replace("'", "''") + "'")
        con.execute("SET preserve_insertion_order=false")
    except Exception:
        con.close()
        raise
    return con


def _query_rows(con, query, batch_size=10000):
    cursor = con.execute(query)
    while True:
        rows = cursor.fetchmany(batch_size)
        if not rows:
            return
        yield from rows


def _parquet_options(compression, compression_level):
    compression = validate_parquet_compression(compression)
    options = f"FORMAT PARQUET, COMPRESSION {compression.upper()}"
    if compression == "zstd":
        options += f", COMPRESSION_LEVEL {int(compression_level)}"
    return options


def build_sparse_matrix(source: Path, matrix_dir: Path, sample: str, barcode_col: str,
                        feature_col: str, *, prefix: str = "raw", threads=1,
                        memory_mb=512, compression_level=3, compression="snappy"):
    """Aggregate and summarize using a configured disk-spilling SQL engine."""
    matrix_dir.mkdir(parents=True, exist_ok=True)
    pairs, summary = pair_path(matrix_dir, sample, prefix), summary_path(matrix_dir, sample, prefix)
    b1, b2 = barcode_col.replace('"', '""'), feature_col.replace('"', '""')
    output_sql = str(pairs).replace("'", "''")
    molecule_bins, protein_bins, molecule_heap, protein_heap = [0] * 8, [0] * 8, [], []
    options = _parquet_options(compression, compression_level)
    with tempfile.TemporaryDirectory(prefix="matrix-spill-", dir=matrix_dir) as temporary:
        con = _connect(temporary, threads=threads, memory_mb=memory_mb)
        try:
            con.execute(f'CREATE TEMP TABLE pairs AS SELECT "{b1}" AS barcode1, "{b2}" AS feature, COUNT(*)::BIGINT AS molecule_count FROM {sql_reader(source)} GROUP BY 1, 2')
            con.execute(f"COPY pairs TO '{output_sql}' ({options})")
            feature_counts = dict(_query_rows(con, "SELECT feature, SUM(molecule_count) FROM pairs GROUP BY feature"))
            for barcode, molecules, proteins in _query_rows(con, "SELECT barcode1, SUM(molecule_count), COUNT(*) FROM pairs GROUP BY barcode1"):
                molecules, proteins = int(molecules), int(proteins)
                molecule_bins[min(molecules, 8) - 1] += 1
                protein_bins[min(proteins, 8) - 1] += 1
                for heap, score in ((molecule_heap, molecules), (protein_heap, proteins)):
                    item = (score, barcode, molecules, proteins)
                    if len(heap) < 20:
                        heapq.heappush(heap, item)
                    elif item > heap[0]:
                        heapq.heapreplace(heap, item)
        finally:
            con.close()
    top = lambda heap: [{"PB": row[1], "Molecules": row[2], "Proteins": row[3]} for row in sorted(heap, reverse=True)]
    payload = {"schema_version": 1, "pair_table": pairs.name, "barcode_count": sum(molecule_bins), "feature_counts": feature_counts,
               "pi_seq": {"total": sum(molecule_bins), "molecule_bins": molecule_bins, "protein_bins": protein_bins,
                          "top_molecules": top(molecule_heap), "top_proteins": top(protein_heap)}}
    with atomic_text(summary) as handle:
        json.dump(payload, handle, separators=(",", ":"), sort_keys=True)
    return pairs, summary


def load_summary(path: Path):
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    return int(payload["barcode_count"]), Counter({k: int(v) for k, v in payload["feature_counts"].items()}), payload.get("pi_seq")


def _text_output(path):
    return (gzip.open(path, "wt", encoding="utf-8", newline="", compresslevel=3)
            if str(path).endswith(".gz") else path.open("w", encoding="utf-8", newline=""))


def export_matrix(pairs: Path, output: Path, fmt: str, barcode_header: str = "Barcode1", *, threads=1, memory_mb=512, compression="snappy", compression_level=3):
    if fmt not in {"parquet", "mtx", "tsv"}:
        raise ConfigError("matrix export format must be parquet, mtx, or tsv")
    options = _parquet_options(compression, compression_level)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="export-spill-", dir=output.parent) as temporary:
        con = _connect(temporary, threads=threads, memory_mb=memory_mb)
        try:
            reader = sql_reader(pairs)
            if fmt == "parquet":
                target = str(output).replace("'", "''")
                con.execute(f"COPY (SELECT * FROM {reader}) TO '{target}' ({options})")
                return [output]
            if fmt == "mtx":
                con.execute(f"CREATE TEMP TABLE barcodes AS SELECT barcode1, ROW_NUMBER() OVER (ORDER BY barcode1) AS idx FROM (SELECT DISTINCT barcode1 FROM {reader})")
                con.execute(f"CREATE TEMP TABLE features AS SELECT feature, ROW_NUMBER() OVER (ORDER BY feature) AS idx FROM (SELECT DISTINCT feature FROM {reader})")
                dimensions = (con.execute("SELECT COUNT(*) FROM barcodes").fetchone()[0],
                              con.execute("SELECT COUNT(*) FROM features").fetchone()[0],
                              con.execute(f"SELECT COUNT(*) FROM {reader}").fetchone()[0])
                prefix = output.name.removesuffix(".mtx.gz").removesuffix(".mtx")
                barcode_path = output.with_name(f"{prefix}.barcodes.tsv.gz")
                feature_path = output.with_name(f"{prefix}.features.tsv.gz")
                for path, query in ((barcode_path, "SELECT barcode1 FROM barcodes ORDER BY idx"),
                                    (feature_path, "SELECT feature FROM features ORDER BY idx")):
                    with _text_output(path) as handle:
                        for (value,) in _query_rows(con, query):
                            handle.write(value + "\n")
                with _text_output(output) as handle:
                    handle.write("%%MatrixMarket matrix coordinate integer general\n% Generated by TagForge\n")
                    handle.write(" ".join(map(str, dimensions)) + "\n")
                    query = f"SELECT b.idx, f.idx, p.molecule_count FROM {reader} p JOIN barcodes b USING (barcode1) JOIN features f USING (feature) ORDER BY 1, 2"
                    for row in _query_rows(con, query):
                        handle.write(" ".join(map(str, row)) + "\n")
                return [output, barcode_path, feature_path]
            features = [row[0] for row in _query_rows(con, f"SELECT DISTINCT feature FROM {reader} ORDER BY feature")]
            current_barcode, counts = None, {}
            with _text_output(output) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([barcode_header] + features)
                for barcode, feature, count in _query_rows(con, f"SELECT barcode1, feature, molecule_count FROM {reader} ORDER BY 1, 2"):
                    if current_barcode is not None and barcode != current_barcode:
                        writer.writerow([current_barcode] + [counts.get(name, 0) for name in features])
                        counts = {}
                    current_barcode = barcode
                    counts[feature] = count
                if current_barcode is not None:
                    writer.writerow([current_barcode] + [counts.get(name, 0) for name in features])
            return [output]
        finally:
            con.close()
