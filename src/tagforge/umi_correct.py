from __future__ import annotations

import csv
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Dict

from .config import TagForgeConfig
from .io_utils import atomic_text, open_tsv, sample_dirs


def _umi_clusterer(method: str):
    try:
        from umi_tools import UMIClusterer
    except ImportError as exc:
        raise RuntimeError(
            "umi_tools is required for UMI deduplication. Activate the TagForge Conda environment."
        ) from exc
    return UMIClusterer(cluster_method=method)


def deduplicate_umis(counts: Dict[str, int], method: str = "directional", max_distance: int = 1):
    """Return raw-to-corrected assignments using UMI-tools' UMIClusterer."""
    encoded_counts = {umi.encode("ascii"): count for umi, count in counts.items()}
    clusterer = _umi_clusterer(method)
    groups = clusterer(encoded_counts, threshold=max_distance)
    assignments = {}
    for group in groups:
        representative = group[0].decode("ascii")
        for raw_umi in group:
            assignments[raw_umi.decode("ascii")] = representative
    missing = set(counts) - set(assignments)
    if missing:
        raise RuntimeError(f"UMI-tools did not return {len(missing)} input UMI(s)")
    return assignments


MOLECULE_FIELDS = ["barcode1", "barcode2_name", "corrected_umi", "reads_count", "raw_umi_count"]


def dedup_sample(config: TagForgeConfig, sample_name: str):
    dirs = sample_dirs(config.output_dir, sample_name)
    source = dirs["detail"] / f"{sample_name}.valid_reads.tsv.gz"
    output = dirs["detail"] / f"{sample_name}.molecule_detail.tsv.gz"
    db_path = dirs["tmp"] / f"{sample_name}.umi_counts.sqlite3"
    db_path.unlink(missing_ok=True)
    connection = sqlite3.connect(db_path)
    try:
        connection.execute("PRAGMA journal_mode=WAL")
        connection.execute("PRAGMA synchronous=NORMAL")
        connection.execute("CREATE TABLE counts (b1 TEXT, b2 TEXT, umi TEXT, n INTEGER, PRIMARY KEY (b1,b2,umi)) WITHOUT ROWID")
        batch = []
        reads = 0
        for row in open_tsv(source):
            batch.append((row["barcode1"], row["barcode2_name"], row["umi"], 1)); reads += 1
            if len(batch) >= config.chunk_size:
                connection.executemany("INSERT INTO counts VALUES (?,?,?,?) ON CONFLICT DO UPDATE SET n=n+1", batch)
                connection.commit(); batch.clear()
        if batch:
            connection.executemany("INSERT INTO counts VALUES (?,?,?,?) ON CONFLICT DO UPDATE SET n=n+1", batch)
            connection.commit()
        cursor = connection.execute("SELECT b1,b2,umi,n FROM counts ORDER BY b1,b2")
        molecules = 0
        with atomic_text(output, config.compression_level) as handle:
            writer = csv.DictWriter(handle, fieldnames=MOLECULE_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            key = None; group = {}
            def flush(group_key, umi_counts):
                nonlocal molecules
                if group_key is None: return
                assignments = deduplicate_umis(umi_counts, config.umi_method, config.umi_max_distance)
                merged = defaultdict(lambda: [0, 0])
                for raw_umi, n in umi_counts.items():
                    merged[assignments[raw_umi]][0] += n; merged[assignments[raw_umi]][1] += 1
                for corrected in sorted(merged):
                    n, raw_count = merged[corrected]; molecules += 1
                    writer.writerow({"barcode1": group_key[0], "barcode2_name": group_key[1],
                                     "corrected_umi": corrected, "reads_count": n, "raw_umi_count": raw_count})
            for b1, b2, umi, n in cursor:
                current = (b1, b2)
                if key is not None and current != key:
                    flush(key, group); group = {}
                key = current; group[umi] = n
            flush(key, group)
    finally:
        connection.close(); db_path.unlink(missing_ok=True)
        Path(str(db_path) + "-wal").unlink(missing_ok=True); Path(str(db_path) + "-shm").unlink(missing_ok=True)
    return output, {"valid_reads": reads, "molecules": molecules, "duplicates": reads - molecules}
