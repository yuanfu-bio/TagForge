from __future__ import annotations

import csv
import gzip
import json
import shutil
import sqlite3
import tempfile
from pathlib import Path

from .columnar import open_records, resolve_table
from .config import TagForgeConfig
from .dominance_qc import DominanceQC
from .io_utils import atomic_text, open_tsv, sample_dirs


def _observed_cb_aliases(cb_counts, min_parent_reads: int, parent_child_ratio: float):
    """Find conservative CB aliases without searching a whitelist."""
    parents = [(cb, reads) for cb, reads in cb_counts if "N" not in cb and reads >= min_parent_reads]
    index = {}
    for cb, reads in parents:
        for position in range(len(cb)):
            key = cb[:position] + "*" + cb[position + 1:]
            index.setdefault(key, []).append((cb, reads))
    aliases = []
    for cb, child_reads in cb_counts:
        n_count = cb.count("N")
        if n_count > 1:
            continue
        keys = [cb.replace("N", "*", 1)] if n_count else [
            cb[:position] + "*" + cb[position + 1:] for position in range(len(cb))]
        candidates = {}
        for key in keys:
            for parent, parent_reads in index.get(key, ()):
                if parent != cb and parent_reads > child_reads and parent_reads >= child_reads * parent_child_ratio:
                    candidates[parent] = parent_reads
        if len(candidates) == 1:
            parent, parent_reads = next(iter(candidates.items()))
            aliases.append((cb, parent, child_reads, parent_reads, "n_wildcard" if n_count else "hamming_1"))
    return aliases


def matrix_from_molecules(
    source: Path,
    output: Path,
    compression_level: int = 3,
    *,
    barcode_col: str = "barcode1",
    feature_col: str = "barcode2_name",
    row_header: str = "Barcode1",
):
    db = Path(str(output) + ".sqlite3.tmp")
    db.unlink(missing_ok=True)
    con = sqlite3.connect(db)
    try:
        con.execute("CREATE TABLE m (b1 TEXT, b2 TEXT, n INTEGER, PRIMARY KEY(b1,b2)) WITHOUT ROWID")
        batch = []
        for row in open_records(source):
            batch.append((row[barcode_col], row[feature_col], 1))
            if len(batch) >= 100000:
                con.executemany("INSERT INTO m VALUES(?,?,?) ON CONFLICT DO UPDATE SET n=n+1", batch); con.commit(); batch.clear()
        if batch:
            con.executemany("INSERT INTO m VALUES(?,?,?) ON CONFLICT DO UPDATE SET n=n+1", batch); con.commit()
        features = [x[0] for x in con.execute("SELECT DISTINCT b2 FROM m ORDER BY b2")]
        with atomic_text(output, compression_level) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([row_header] + features)
            cursor = con.execute("SELECT b1,b2,n FROM m ORDER BY b1,b2")
            current = None; values = {}
            for b1, b2, n in cursor:
                if current is not None and b1 != current:
                    writer.writerow([current] + [values.get(f, 0) for f in features]); values = {}
                current = b1; values[b2] = n
            if current is not None:
                writer.writerow([current] + [values.get(f, 0) for f in features])
        return {"barcode1_count": con.execute("SELECT COUNT(DISTINCT b1) FROM m").fetchone()[0],
                "feature_count": len(features), "nonzero_count": con.execute("SELECT COUNT(*) FROM m").fetchone()[0]}
    finally:
        con.close(); db.unlink(missing_ok=True)


def matrix_sample(config: TagForgeConfig, sample_name: str):
    dirs = sample_dirs(config.output_dir, sample_name)
    preferred = dirs["detail"] / (f"{sample_name}.molecule_detail.rmMP.parquet" if config.pi_seq_enabled else f"{sample_name}.molecule_detail.parquet")
    source = resolve_table(preferred, preferred.with_name(preferred.name.replace(".parquet", ".tsv.gz")))
    from .sparse_matrix import build_sparse_matrix
    (dirs["matrix"] / f"{sample_name}.raw_count_matrix.tsv.gz").unlink(missing_ok=True)
    return build_sparse_matrix(
        source, dirs["matrix"], sample_name,
        barcode_col=config.target_name("barcode1"),
        feature_col=f"{config.target_name('barcode2')}_name",
        threads=config.threads, memory_mb=config.matrix_memory_mb, compression_level=config.compression_level,
        compression=config.parquet_compression,
    )


def pair_mapping_sample(config: TagForgeConfig, sample_name: str):
    """Aggregate UMI-deduplicated molecules into a compact PB--CB mapping.

    For ``library.type: pb-cb``, also select the sole CB for each PB using
    read support after UMI deduplication. A PB is retained only if its leading
    CB has support at least the configured dominance threshold.
    """
    dirs = sample_dirs(config.output_dir, sample_name)
    source = resolve_table(dirs["detail"] / f"{sample_name}.molecule_detail.parquet",
                           dirs["detail"] / f"{sample_name}.molecule_detail.tsv.gz")
    output = dirs["matrix"] / f"{sample_name}.pb_cb_mapping.tsv.gz"
    map_output = dirs["matrix"] / f"{sample_name}.pb_cb_map.tsv.gz"
    cb_counts_output = dirs["matrix"] / f"{sample_name}.cb_pb_counts.tsv.gz"
    distribution_output = dirs["matrix"] / f"{sample_name}.cb_pb_count_distribution.tsv"
    correction_output = dirs["matrix"] / f"{sample_name}.cb_observed_correction.tsv.gz"
    pb_column = config.target_name("barcode1")
    cb_column = f"{config.target_name('barcode2')}_name"
    import duckdb
    import pyarrow as pa
    from .columnar import sql_reader
    quote = lambda name: '"' + name.replace('"', '""') + '"'
    scratch = Path(tempfile.mkdtemp(prefix=".pair-map.", dir=dirs["matrix"]))
    try:
        con = duckdb.connect()
        try:
            con.execute(f"PRAGMA threads={max(1, getattr(config, 'threads', 1))}")
            con.execute(f"PRAGMA memory_limit='{int(getattr(config, 'matrix_memory_mb', 16384))}MB'")
            con.execute("PRAGMA temp_directory=?", [str(scratch / "spill")])
            con.execute("PRAGMA preserve_insertion_order=false")
            con.execute(
                f"CREATE TEMP TABLE molecules AS SELECT {quote(pb_column)} AS pb, {quote(cb_column)} AS cb, "
                f"corrected_umi AS umi, SUM(reads_count::BIGINT)::BIGINT AS reads_count "
                f"FROM {sql_reader(source)} GROUP BY 1, 2, 3")
            con.execute("CREATE TEMP TABLE mapping AS SELECT pb, cb, COUNT(*)::BIGINT AS molecule_count, "
                        "SUM(reads_count)::BIGINT AS reads_count FROM molecules GROUP BY pb, cb")

            def rows(query, *params):
                cursor = con.execute(query, list(params))
                while True:
                    chunk = cursor.fetchmany(100000)
                    if not chunk:
                        return
                    yield from chunk

            with atomic_text(output, config.compression_level) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([pb_column, config.target_name("barcode2"), "molecule_count", "reads_count"])
                writer.writerows(rows("SELECT pb, cb, molecule_count, reads_count FROM mapping ORDER BY pb, cb"))
            if not getattr(config, "pb_cb_enabled", False):
                return output

            pb_cb_config = getattr(config, "raw", {}).get("pb_cb", {})
            aliases = []
            if bool(pb_cb_config.get("observed_correction", True)):
                aliases = _observed_cb_aliases(
                    [(cb, int(reads)) for cb, reads in con.execute(
                        "SELECT cb, SUM(reads_count) FROM mapping GROUP BY cb ORDER BY cb").fetchall()],
                    int(pb_cb_config.get("min_parent_reads", 5)),
                    float(pb_cb_config.get("parent_child_ratio", 10)),
                )
            with atomic_text(correction_output, config.compression_level) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["raw_cb", "corrected_cb", "raw_reads_count", "parent_reads_count", "correction_type"])
                writer.writerows(aliases)
            alias_table = pa.table({
                "raw_cb": pa.array([row[0] for row in aliases], type=pa.string()),
                "corrected_cb": pa.array([row[1] for row in aliases], type=pa.string())})
            con.register("cb_alias_arrow", alias_table)
            con.execute("CREATE TEMP TABLE cb_alias AS SELECT * FROM cb_alias_arrow")
            con.execute("""
                CREATE TEMP TABLE resolved_mapping AS
                SELECT pb, cb, COUNT(*)::BIGINT AS molecule_count, SUM(reads_count)::BIGINT AS reads_count FROM (
                    SELECT molecules.pb, COALESCE(cb_alias.corrected_cb, molecules.cb) AS cb,
                           molecules.umi, SUM(molecules.reads_count) AS reads_count
                    FROM molecules LEFT JOIN cb_alias ON molecules.cb = cb_alias.raw_cb
                    GROUP BY molecules.pb, COALESCE(cb_alias.corrected_cb, molecules.cb), molecules.umi
                ) GROUP BY pb, cb
            """)

            threshold = getattr(config, "pb_cb_dominance_threshold", 0.8)
            ranked = """
                WITH ranked AS (
                    SELECT pb, cb, molecule_count, reads_count,
                           SUM(molecule_count) OVER (PARTITION BY pb) AS total_molecule_count,
                           SUM(reads_count) OVER (PARTITION BY pb) AS total_reads_count,
                           COUNT(*) OVER (PARTITION BY pb) AS cb_count,
                           ROW_NUMBER() OVER (
                               PARTITION BY pb ORDER BY reads_count DESC, cb ASC
                           ) AS rn
                    FROM resolved_mapping
                )
                SELECT pb, cb, molecule_count, reads_count, total_molecule_count,
                       total_reads_count, CAST(reads_count AS DOUBLE) / total_reads_count AS dominant_ratio, cb_count
                FROM ranked WHERE rn = 1 ORDER BY pb
            """
            selected = [tuple(int(v) if i in (2, 3, 4, 5) else v for i, v in enumerate(row))
                        for row in con.execute(ranked).fetchall()]
            qc = DominanceQC(threshold)
            for row in selected:
                qc.add(row[7], row[5], row[3])
            retained = [row for row in selected if row[6] >= threshold]
            with atomic_text(map_output, config.compression_level) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([
                    pb_column, config.target_name("barcode2"), "molecule_count", "reads_count",
                    "total_molecule_count", "total_reads_count", "dominant_ratio",
                ])
                writer.writerows(row[:7] for row in retained)
            con.register("retained_arrow", pa.table({
                "pb": pa.array([row[0] for row in retained], type=pa.string()),
                "cb": pa.array([row[1] for row in retained], type=pa.string())}))
            con.execute("CREATE TEMP TABLE retained_mapping AS SELECT * FROM retained_arrow")
            cb_counts = con.execute("""
                SELECT cb, COUNT(*) AS pb_count
                FROM retained_mapping GROUP BY cb ORDER BY cb
            """).fetchall()
            with atomic_text(cb_counts_output, config.compression_level) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([config.target_name("barcode2"), "pb_count"])
                writer.writerows(cb_counts)
            distribution = con.execute("""
                SELECT pb_count, COUNT(*) AS cb_count
                FROM (
                    SELECT cb, COUNT(*) AS pb_count
                    FROM retained_mapping GROUP BY cb
                ) GROUP BY pb_count ORDER BY pb_count
            """).fetchall()
            with atomic_text(distribution_output, config.compression_level) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["pb_count", "cb_count"])
                writer.writerows(distribution)
            raw_molecules, raw_reads = con.execute("SELECT COALESCE(SUM(molecule_count),0), COALESCE(SUM(reads_count),0) FROM mapping").fetchone()
            resolved_molecules = con.execute("SELECT COALESCE(SUM(molecule_count),0) FROM resolved_mapping").fetchone()[0]
            payload = {"schema": 1, "sample": sample_name, "dominance_basis": "read_support",
                       "dominance_threshold": threshold, "raw_molecules": raw_molecules,
                       "resolved_molecules": resolved_molecules, "supporting_reads": raw_reads,
                       "total_pbs": len(selected),
                       "retained_pbs": len(retained), "retained_cbs": len(cb_counts),
                       "retained_molecules": sum(row[2] for row in retained),
                       "retained_supporting_reads": sum(row[3] for row in retained),
                       "cb_aliases": len(aliases), "dominance_qc": qc.payload()}
            with atomic_text(dirs["matrix"] / f"{sample_name}.pb_cb_summary.json") as handle:
                json.dump(payload, handle, sort_keys=True)
            return map_output
        finally:
            con.close()
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
