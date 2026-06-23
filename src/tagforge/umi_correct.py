from __future__ import annotations

import csv
import sqlite3
import time
from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, Iterator

from .config import TagForgeConfig
from .io_utils import atomic_text, open_tsv, sample_dirs
from .logging_utils import sample_logger


def _umi_clusterer(method: str):
    try:
        from umi_tools import UMIClusterer
    except ImportError as exc:
        raise RuntimeError(
            "umi_tools is required for UMI deduplication. Activate the TagForge Conda environment."
        ) from exc
    return UMIClusterer(cluster_method=method)


def _assign_umis(counts: Dict[str, int], clusterer, max_distance: int):
    encoded_counts = {umi.encode("ascii"): count for umi, count in counts.items()}
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


def deduplicate_umis(counts: Dict[str, int], method: str = "directional", max_distance: int = 1):
    """Return raw-to-corrected assignments using UMI-tools' UMIClusterer."""
    return _assign_umis(counts, _umi_clusterer(method), max_distance)


def molecule_fields(config: TagForgeConfig):
    return [
        config.target_name("barcode1"),
        f"{config.target_name('barcode2')}_name",
        "corrected_umi", "reads_count", "raw_umi_count",
    ]


def _group_batches(cursor: Iterable[tuple], batch_umis: int) -> Iterator[list]:
    """Batch complete barcode groups without ever splitting a correction scope."""
    batch = []
    batch_size = 0
    key = None
    group = {}
    for b1, b2, umi, n in cursor:
        current = (b1, b2)
        if key is not None and current != key:
            if batch and batch_size + len(group) > batch_umis:
                yield batch
                batch = []
                batch_size = 0
            batch.append((key, group))
            batch_size += len(group)
            if batch_size >= batch_umis:
                yield batch
                batch = []
                batch_size = 0
            group = {}
        key = current
        group[umi] = n
    if key is not None:
        if batch and batch_size + len(group) > batch_umis:
            yield batch
            batch = []
        batch.append((key, group))
    if batch:
        yield batch


def _dedup_batch(batch: list, method: str, max_distance: int):
    """Process one bounded batch inside a worker process."""
    clusterer = _umi_clusterer(method)
    rows = []
    input_reads = 0
    raw_umis = 0
    for group_key, umi_counts in batch:
        assignments = _assign_umis(umi_counts, clusterer, max_distance)
        merged = defaultdict(lambda: [0, 0])
        for raw_umi, n in umi_counts.items():
            input_reads += n
            raw_umis += 1
            merged[assignments[raw_umi]][0] += n
            merged[assignments[raw_umi]][1] += 1
        for corrected in sorted(merged):
            n, raw_count = merged[corrected]
            rows.append((group_key[0], group_key[1], corrected, n, raw_count))
    return rows, len(batch), input_reads, raw_umis


def _worker_ready():
    """Small startup probe used before consuming the SQLite cursor."""
    return True


def dedup_sample(config: TagForgeConfig, sample_name: str):
    dirs = sample_dirs(config.output_dir, sample_name)
    logger = sample_logger(sample_name, dirs["logs"] / f"{sample_name}.pipeline.log")
    source = dirs["detail"] / f"{sample_name}.valid_reads.tsv.gz"
    output = dirs["detail"] / f"{sample_name}.molecule_detail.tsv.gz"
    db_path = dirs["tmp"] / f"{sample_name}.umi_counts.sqlite3"
    db_path.unlink(missing_ok=True)
    started = time.monotonic()
    aggregation_started = started
    connection = sqlite3.connect(db_path)
    executor = None
    requested_workers = config.umi_workers or config.threads
    workers = requested_workers
    barcode1_col = config.target_name("barcode1")
    barcode2_name_col = f"{config.target_name('barcode2')}_name"
    umi_col = config.target_name("umi")
    fields = molecule_fields(config)
    try:
        # This database is a disposable aggregation scratch file. Keeping its
        # cache bounded and temporary data on disk avoids competing with UMI
        # worker processes for RAM.
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("PRAGMA temp_store=FILE")
        connection.execute(f"PRAGMA cache_size=-{config.umi_sqlite_cache_mb * 1024}")
        connection.execute(
            "CREATE TABLE counts (b1 TEXT, b2 TEXT, umi TEXT, n INTEGER, "
            "PRIMARY KEY (b1,b2,umi)) WITHOUT ROWID"
        )
        batch = []
        reads = 0
        for row in open_tsv(source):
            batch.append((row[barcode1_col], row[barcode2_name_col], row[umi_col], 1))
            reads += 1
            if len(batch) >= config.chunk_size:
                connection.executemany(
                    "INSERT INTO counts VALUES (?,?,?,?) ON CONFLICT DO UPDATE SET n=n+1", batch
                )
                connection.commit()
                batch.clear()
        if batch:
            connection.executemany(
                "INSERT INTO counts VALUES (?,?,?,?) ON CONFLICT DO UPDATE SET n=n+1", batch
            )
            connection.commit()
        aggregation_seconds = time.monotonic() - aggregation_started
        cursor = connection.execute("SELECT b1,b2,umi,n FROM counts ORDER BY b1,b2")
        group_batches = _group_batches(cursor, config.umi_batch_size)
        molecules = 0
        groups_processed = 0
        raw_umis_processed = 0
        peak_batch_umis = 0
        clustering_started = time.monotonic()

        if workers > 1:
            try:
                executor = ProcessPoolExecutor(max_workers=workers)
                executor.submit(_worker_ready).result()
            except (PermissionError, NotImplementedError, OSError) as exc:
                if executor is not None:
                    executor.shutdown(wait=True, cancel_futures=True)
                executor = None
                workers = 1
                logger.warning(
                    "dedup_parallel_fallback\trequested_workers=%s\tworkers=1\treason=%s",
                    requested_workers, type(exc).__name__,
                )

        logger.info(
            "dedup_parallel_start\tbackend=umi_tools-UMIClusterer\trequested_workers=%s\tworkers=%s\t"
            "batch_umis=%s\tmax_pending_batches=%s\tsqlite_cache_mb=%s\treads=%s\t"
            "aggregation_seconds=%.3f",
            requested_workers, workers, config.umi_batch_size, workers, config.umi_sqlite_cache_mb,
            reads, aggregation_seconds,
        )

        with atomic_text(output, config.compression_level) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)

            def consume(result):
                nonlocal molecules, groups_processed, raw_umis_processed
                rows, group_count, _batch_reads, raw_umi_count = result
                writer.writerows(rows)
                molecules += len(rows)
                groups_processed += group_count
                raw_umis_processed += raw_umi_count

            if workers == 1:
                for group_batch in group_batches:
                    peak_batch_umis = max(peak_batch_umis, sum(len(group) for _, group in group_batch))
                    consume(_dedup_batch(group_batch, config.umi_method, config.umi_max_distance))
            else:
                # At most one submitted batch per worker: enough to keep all
                # CPUs busy without an unbounded multiprocessing task queue.
                pending = deque()
                with executor:
                    for group_batch in group_batches:
                        batch_umi_count = sum(len(group) for _, group in group_batch)
                        peak_batch_umis = max(peak_batch_umis, batch_umi_count)
                        pending.append(executor.submit(
                            _dedup_batch, group_batch, config.umi_method, config.umi_max_distance
                        ))
                        if len(pending) >= workers:
                            consume(pending.popleft().result())
                    while pending:
                        consume(pending.popleft().result())

        clustering_seconds = time.monotonic() - clustering_started
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
        connection.close()
        db_path.unlink(missing_ok=True)
        Path(str(db_path) + "-journal").unlink(missing_ok=True)
    wall_seconds = time.monotonic() - started
    return output, {
        "valid_reads": reads,
        "molecules": molecules,
        "duplicates": reads - molecules,
        "groups": groups_processed,
        "raw_umis": raw_umis_processed,
        "requested_workers": requested_workers,
        "workers": workers,
        "umi_batch_size": config.umi_batch_size,
        "peak_batch_umis": peak_batch_umis,
        "sqlite_cache_mb": config.umi_sqlite_cache_mb,
        "aggregation_seconds": aggregation_seconds,
        "clustering_seconds": clustering_seconds,
        "wall_seconds": wall_seconds,
    }
