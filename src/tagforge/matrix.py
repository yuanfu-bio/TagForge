from __future__ import annotations

import csv
import gzip
import sqlite3
from pathlib import Path

from .config import TagForgeConfig
from .io_utils import atomic_text, open_tsv, sample_dirs


def matrix_from_molecules(source: Path, output: Path, compression_level: int = 3):
    db = Path(str(output) + ".sqlite3.tmp")
    db.unlink(missing_ok=True)
    con = sqlite3.connect(db)
    try:
        con.execute("CREATE TABLE m (b1 TEXT, b2 TEXT, n INTEGER, PRIMARY KEY(b1,b2)) WITHOUT ROWID")
        batch = []
        for row in open_tsv(source):
            batch.append((row["barcode1"], row["barcode2_name"], 1))
            if len(batch) >= 100000:
                con.executemany("INSERT INTO m VALUES(?,?,?) ON CONFLICT DO UPDATE SET n=n+1", batch); con.commit(); batch.clear()
        if batch:
            con.executemany("INSERT INTO m VALUES(?,?,?) ON CONFLICT DO UPDATE SET n=n+1", batch); con.commit()
        features = [x[0] for x in con.execute("SELECT DISTINCT b2 FROM m ORDER BY b2")]
        with atomic_text(output, compression_level) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["Barcode1"] + features)
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
    source = dirs["detail"] / f"{sample_name}.molecule_detail.tsv.gz"
    output = dirs["matrix"] / f"{sample_name}.raw_count_matrix.tsv.gz"
    return output, matrix_from_molecules(source, output, config.compression_level)

