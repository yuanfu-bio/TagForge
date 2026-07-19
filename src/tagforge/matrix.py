from __future__ import annotations

import csv
import gzip
import sqlite3
from pathlib import Path

from .config import TagForgeConfig
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
                if parent != cb and parent_reads >= child_reads * parent_child_ratio:
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
        for row in open_tsv(source):
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
    source = dirs["detail"] / (f"{sample_name}.molecule_detail.rmMP.tsv.gz" if config.pi_seq_enabled else f"{sample_name}.molecule_detail.tsv.gz")
    output = dirs["matrix"] / f"{sample_name}.raw_count_matrix.tsv.gz"
    return output, matrix_from_molecules(
        source, output, config.compression_level,
        barcode_col=config.target_name("barcode1"),
        feature_col=f"{config.target_name('barcode2')}_name",
        row_header=config.target_name("barcode1"),
    )


def pair_mapping_sample(config: TagForgeConfig, sample_name: str):
    """Aggregate UMI-deduplicated molecules into a compact PB--CB mapping.

    For ``library.type: pb-cb``, also select the sole CB for each PB using
    read support after UMI deduplication. A PB is retained only if its leading
    CB has support strictly greater than the configured dominance threshold.
    """
    dirs = sample_dirs(config.output_dir, sample_name)
    source = dirs["detail"] / f"{sample_name}.molecule_detail.tsv.gz"
    output = dirs["matrix"] / f"{sample_name}.pb_cb_mapping.tsv.gz"
    map_output = dirs["matrix"] / f"{sample_name}.pb_cb_map.tsv.gz"
    cb_counts_output = dirs["matrix"] / f"{sample_name}.cb_pb_counts.tsv.gz"
    distribution_output = dirs["matrix"] / f"{sample_name}.cb_pb_count_distribution.tsv"
    correction_output = dirs["matrix"] / f"{sample_name}.cb_observed_correction.tsv.gz"
    db = Path(str(output) + ".sqlite3.tmp")
    db.unlink(missing_ok=True)
    pb_column = config.target_name("barcode1")
    cb_column = f"{config.target_name('barcode2')}_name"
    try:
        con = sqlite3.connect(db)
        try:
            con.execute(
                "CREATE TABLE mapping (pb TEXT, cb TEXT, molecule_count INTEGER, reads_count INTEGER, "
                "PRIMARY KEY(pb, cb)) WITHOUT ROWID"
            )
            batch = []
            for row in open_tsv(source):
                batch.append((row[pb_column], row[cb_column], int(row["reads_count"])))
                if len(batch) >= 100000:
                    con.executemany(
                        "INSERT INTO mapping VALUES (?, ?, 1, ?) "
                        "ON CONFLICT(pb, cb) DO UPDATE SET molecule_count=molecule_count+1, "
                        "reads_count=reads_count+excluded.reads_count", batch
                    )
                    con.commit(); batch.clear()
            if batch:
                con.executemany(
                    "INSERT INTO mapping VALUES (?, ?, 1, ?) "
                    "ON CONFLICT(pb, cb) DO UPDATE SET molecule_count=molecule_count+1, "
                    "reads_count=reads_count+excluded.reads_count", batch
                )
                con.commit()
            with atomic_text(output, config.compression_level) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([pb_column, config.target_name("barcode2"), "molecule_count", "reads_count"])
                writer.writerows(con.execute(
                    "SELECT pb, cb, molecule_count, reads_count FROM mapping ORDER BY pb, cb"
                ))
            if not getattr(config, "pb_cb_enabled", False):
                return output

            pb_cb_config = getattr(config, "raw", {}).get("pb_cb", {})
            aliases = []
            if bool(pb_cb_config.get("observed_correction", True)):
                aliases = _observed_cb_aliases(
                    con.execute("SELECT cb, SUM(reads_count) FROM mapping GROUP BY cb").fetchall(),
                    int(pb_cb_config.get("min_parent_reads", 5)),
                    float(pb_cb_config.get("parent_child_ratio", 10)),
                )
            with atomic_text(correction_output, config.compression_level) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["raw_cb", "corrected_cb", "raw_reads_count", "parent_reads_count", "correction_type"])
                writer.writerows(aliases)
            con.execute("CREATE TEMP TABLE cb_alias (raw_cb TEXT PRIMARY KEY, corrected_cb TEXT NOT NULL) WITHOUT ROWID")
            con.executemany("INSERT INTO cb_alias VALUES (?, ?)", [(row[0], row[1]) for row in aliases])
            con.execute(
                "CREATE TEMP TABLE resolved_mapping (pb TEXT, cb TEXT, molecule_count INTEGER, reads_count INTEGER, "
                "PRIMARY KEY(pb, cb)) WITHOUT ROWID"
            )
            con.execute("""
                INSERT INTO resolved_mapping
                SELECT mapping.pb, COALESCE(cb_alias.corrected_cb, mapping.cb),
                       SUM(mapping.molecule_count), SUM(mapping.reads_count)
                FROM mapping LEFT JOIN cb_alias ON mapping.cb = cb_alias.raw_cb
                GROUP BY mapping.pb, COALESCE(cb_alias.corrected_cb, mapping.cb)
            """)

            threshold = getattr(config, "pb_cb_dominance_threshold", 0.8)
            ranked = """
                WITH ranked AS (
                    SELECT pb, cb, molecule_count, reads_count,
                           SUM(molecule_count) OVER (PARTITION BY pb) AS total_molecule_count,
                           SUM(reads_count) OVER (PARTITION BY pb) AS total_reads_count,
                           ROW_NUMBER() OVER (
                               PARTITION BY pb ORDER BY reads_count DESC, cb ASC
                           ) AS rank
                    FROM resolved_mapping
                )
                SELECT pb, cb, molecule_count, reads_count, total_molecule_count,
                       total_reads_count, CAST(reads_count AS REAL) / total_reads_count AS dominant_ratio
                FROM ranked WHERE rank = 1
            """
            selected = list(con.execute(ranked))
            retained = [row for row in selected if row[6] > threshold]
            with atomic_text(map_output, config.compression_level) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([
                    pb_column, config.target_name("barcode2"), "molecule_count", "reads_count",
                    "total_molecule_count", "total_reads_count", "dominant_ratio",
                ])
                writer.writerows(retained)
            con.execute(
                "CREATE TEMP TABLE retained_mapping (pb TEXT PRIMARY KEY, cb TEXT NOT NULL) WITHOUT ROWID"
            )
            con.executemany("INSERT INTO retained_mapping VALUES (?, ?)", [(row[0], row[1]) for row in retained])
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
            return map_output
        finally:
            con.close()
    finally:
        db.unlink(missing_ok=True)
