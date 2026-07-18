from __future__ import annotations

import csv
import shutil
import sqlite3
import subprocess
import tempfile
import time
from collections import defaultdict, deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Iterable, Iterator

from .config import ConfigError, TagForgeConfig
from .io_utils import atomic_text, open_tsv, sample_dirs, write_tsv
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


def _umi_components(counts: Dict[str, int], max_distance: int) -> Iterator[Dict[str, int]]:
    """Yield independent Hamming-1 UMI components in deterministic order.

    Directional correction can only merge UMI sequences connected by an edge at
    the configured edit distance. For the standard one-base setting, separate
    connected components are therefore safe to correct independently. This
    avoids million-UMI calls to UMI-tools when nearly all UMIs are isolated.
    """
    bases = "ACGTN"
    if max_distance != 1 or len(counts) < 2 or any(set(umi) - set(bases) for umi in counts):
        yield counts
        return
    seen = set()
    for seed in counts:
        if seed in seen:
            continue
        component = {seed}
        seen.add(seed)
        pending = [seed]
        while pending:
            umi = pending.pop()
            for index, base in enumerate(umi):
                prefix, suffix = umi[:index], umi[index + 1:]
                for replacement in bases:
                    if replacement == base:
                        continue
                    neighbor = prefix + replacement + suffix
                    if neighbor in counts and neighbor not in seen:
                        seen.add(neighbor)
                        component.add(neighbor)
                        pending.append(neighbor)
        yield {umi: counts[umi] for umi in sorted(component)}


def molecule_fields(config: TagForgeConfig):
    return [
        config.target_name("barcode1"),
        f"{config.target_name('barcode2')}_name",
        "corrected_umi", "reads_count", "raw_umi_count",
    ]


def _group_batches(cursor: Iterable[tuple], batch_umis: int, component_max_distance: int | None = None) -> Iterator[list]:
    """Batch independent UMI components without splitting correction edges."""
    batch = []
    batch_size = 0
    key = None
    group = {}

    def add_components(group_key, umi_counts):
        if component_max_distance is None:
            return [(group_key, umi_counts)]
        return [(group_key, component) for component in _umi_components(umi_counts, component_max_distance)]

    def add_to_batches(items):
        nonlocal batch, batch_size
        for item in items:
            component_size = len(item[1])
            if batch and batch_size + component_size > batch_umis:
                yield batch
                batch = []
                batch_size = 0
            batch.append(item)
            batch_size += component_size
            if batch_size >= batch_umis:
                yield batch
                batch = []
                batch_size = 0

    for b1, b2, umi, n in cursor:
        current = (b1, b2)
        if key is not None and current != key:
            yield from add_to_batches(add_components(key, group))
            group = {}
        key = current
        group[umi] = n
    if key is not None:
        yield from add_to_batches(add_components(key, group))
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


def _aggregate_batch(rows: list[tuple[str, str, str]]):
    counts = defaultdict(int)
    for row in rows:
        counts[row] += 1
    return [(b1, b2, umi, count) for (b1, b2, umi), count in counts.items()], len(rows)


def _worker_ready():
    """Small startup probe used before consuming the SQLite cursor."""
    return True


def _external_sort_aggregate(
    source: Path, temp_dir: Path, barcode1_col: str, barcode2_name_col: str, umi_col: str,
    workers: int, memory_mb: int, progress,
):
    """Sort raw UMI tuples externally, then stream exact counts to a TSV file."""
    sort = shutil.which("sort")
    if sort is None:
        raise ConfigError(
            "UMI external_sort aggregation requires GNU sort in PATH. "
            "Install coreutils in the TagForge environment or set "
            "performance.umi_aggregation_backend: sqlite."
        )
    try:
        version = subprocess.run(
            [sort, "--version"], check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).stdout.splitlines()[0]
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ConfigError(f"Could not run required external tool 'sort': {exc}") from exc
    if "GNU" not in version:
        raise ConfigError("UMI external_sort aggregation requires GNU sort (coreutils)")

    sorted_path = temp_dir / "raw_umi_tuples.sorted.tsv"
    aggregated_path = temp_dir / "raw_umi_counts.tsv"
    command = [
        sort, "--parallel", str(workers), "--buffer-size", f"{memory_mb}M", "--temporary-directory", str(temp_dir),
        "--field-separator", "\t", "--key", "1,1", "--key", "2,2", "--key", "3,3",
    ]
    reads = 0
    process = None
    try:
        with sorted_path.open("w", encoding="utf-8", newline="") as sorted_handle:
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=sorted_handle, stderr=subprocess.PIPE, text=True)
            assert process.stdin is not None
            try:
                for row in open_tsv(source):
                    values = (row[barcode1_col], row[barcode2_name_col], row[umi_col])
                    if any("\t" in value or "\n" in value or "\r" in value for value in values):
                        raise ConfigError("PB, FB name, and UMI values must not contain tabs or newlines for external_sort aggregation")
                    process.stdin.write("\t".join(values) + "\n")
                    reads += 1
                    if reads % 1_000_000 == 0:
                        progress("sort_input", reads, 0, 0)
            finally:
                process.stdin.close()
            while True:
                try:
                    process.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    progress("sort_running", reads, 0, 0)
            stderr = process.stderr.read() if process.stderr is not None else ""
            if process.stderr is not None:
                process.stderr.close()
            if process.returncode != 0:
                raise RuntimeError(f"GNU sort failed during UMI aggregation: {stderr.strip() or 'unknown error'}")
        progress("sort_complete", reads, 0, 0)

        raw_umis = groups = 0
        previous = None
        count = 0
        with sorted_path.open("r", encoding="utf-8", newline="") as source_handle, \
             aggregated_path.open("w", encoding="utf-8", newline="") as output_handle:
            for line in source_handle:
                key = tuple(line.rstrip("\n").split("\t"))
                if len(key) != 3:
                    raise RuntimeError("GNU sort produced an invalid UMI tuple record")
                if previous is not None and key != previous:
                    output_handle.write("\t".join((*previous, str(count))) + "\n")
                    raw_umis += 1
                    if previous[:2] != key[:2]:
                        groups += 1
                    if raw_umis % 1_000_000 == 0:
                        progress("count_output", reads, raw_umis, groups)
                    previous, count = key, 1
                elif previous is None:
                    previous, count = key, 1
                else:
                    count += 1
            if previous is not None:
                output_handle.write("\t".join((*previous, str(count))) + "\n")
                raw_umis += 1
                groups += 1
        sorted_path.unlink(missing_ok=True)
        progress("aggregate_complete", reads, raw_umis, groups)
        return aggregated_path, reads, raw_umis, groups, version
    except Exception:
        if process is not None and process.poll() is None:
            process.terminate()
            process.wait()
        if process is not None and process.stderr is not None and not process.stderr.closed:
            process.stderr.close()
        sorted_path.unlink(missing_ok=True)
        aggregated_path.unlink(missing_ok=True)
        raise


def _aggregated_tsv_cursor(path: Path) -> Iterator[tuple[str, str, str, int]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            b1, b2, umi, count = line.rstrip("\n").split("\t")
            yield b1, b2, umi, int(count)



def _dedup_sample_external_sort(config: TagForgeConfig, sample_name: str):
    """Deduplicate after GNU sort has made each PB-FB UMI scope contiguous."""
    dirs = sample_dirs(config.output_dir, sample_name)
    logger = sample_logger(sample_name, dirs["logs"] / f"{sample_name}.pipeline.log")
    source = dirs["detail"] / f"{sample_name}.valid_reads.tsv.gz"
    output = dirs["detail"] / f"{sample_name}.molecule_detail.tsv.gz"
    progress_path = dirs["logs"] / f"{sample_name}.dedup_progress.tsv"
    started = time.monotonic()
    aggregation_workers = config.umi_aggregation_workers or min(config.threads, 4)
    requested_workers = config.umi_workers or config.threads
    workers = requested_workers
    b1_col, b2_col, umi_col = config.target_name("barcode1"), f"{config.target_name('barcode2')}_name", config.target_name("umi")
    reads = total_raw_umis = total_groups = molecules = groups_processed = raw_umis_processed = 0
    batches_submitted = batches_completed = peak_batch_umis = 0
    executor = None
    sort_dir = Path(tempfile.mkdtemp(prefix=f"{sample_name}.umi-sort.", dir=dirs["tmp"]))
    fields = molecule_fields(config)
    progress_fields = ["status", "valid_reads", "groups", "total_groups", "raw_umis", "total_raw_umis", "molecules", "batches_submitted", "batches_completed", "requested_workers", "workers", "aggregation_workers", "aggregation_records", "elapsed_seconds"]

    def progress(status: str):
        elapsed = time.monotonic() - started
        row = {"status": status, "valid_reads": reads, "groups": total_groups, "total_groups": total_groups, "components": groups_processed, "raw_umis": raw_umis_processed, "total_raw_umis": total_raw_umis, "molecules": molecules, "batches_submitted": batches_submitted, "batches_completed": batches_completed, "requested_workers": requested_workers, "workers": workers, "aggregation_workers": aggregation_workers, "aggregation_records": total_raw_umis, "elapsed_seconds": f"{elapsed:.2f}"}
        write_tsv(progress_path, progress_fields, [row])
        logger.info("dedup_progress\tstatus=%s\tvalid_reads=%s\tgroups=%s\ttotal_groups=%s\traw_umis=%s\ttotal_raw_umis=%s\tmolecules=%s\tbatches_completed=%s\taggregation_records=%s\telapsed_seconds=%s", status, reads, groups_processed, total_groups, raw_umis_processed, total_raw_umis, molecules, batches_completed, total_raw_umis, row["elapsed_seconds"])

    def aggregation_progress(status: str, input_reads: int, raw_umis: int, groups: int):
        nonlocal reads, total_raw_umis, total_groups
        reads, total_raw_umis, total_groups = input_reads, raw_umis, groups
        progress(status)

    try:
        logger.info("dedup_aggregation_start\tbackend=gnu-sort\trequested_workers=%s\tworkers=%s\tsort_memory_mb=%s\ttmp_dir=%s", aggregation_workers, aggregation_workers, config.umi_sort_memory_mb, sort_dir)
        aggregated_path, reads, total_raw_umis, total_groups, sort_version = _external_sort_aggregate(source, sort_dir, b1_col, b2_col, umi_col, aggregation_workers, config.umi_sort_memory_mb, aggregation_progress)
        aggregation_seconds = time.monotonic() - started
        logger.info("dedup_aggregation_complete\tbackend=gnu-sort\tsort_version=%s\treads=%s\tunique_raw_umis=%s\tgroups=%s\tseconds=%.3f", sort_version, reads, total_raw_umis, total_groups, aggregation_seconds)
        progress("aggregated")
        if workers > 1:
            try:
                executor = ProcessPoolExecutor(max_workers=workers)
                executor.submit(_worker_ready).result()
            except (PermissionError, NotImplementedError, OSError) as exc:
                if executor is not None:
                    executor.shutdown(wait=True, cancel_futures=True)
                executor = None
                workers = 1
                logger.warning("dedup_parallel_fallback\trequested_workers=%s\tworkers=1\treason=%s", requested_workers, type(exc).__name__)
        logger.info("dedup_parallel_start\tbackend=umi_tools-UMIClusterer\trequested_workers=%s\tworkers=%s\tbatch_umis=%s\treads=%s\ttotal_raw_umis=%s\taggregation_workers=%s\taggregation_seconds=%.3f", requested_workers, workers, config.umi_batch_size, reads, total_raw_umis, aggregation_workers, aggregation_seconds)
        clustering_started = time.monotonic()
        group_batches = _group_batches(
            _aggregated_tsv_cursor(aggregated_path), config.umi_batch_size, config.umi_max_distance,
        )
        with atomic_text(output, config.compression_level) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)
            def consume(result):
                nonlocal molecules, groups_processed, raw_umis_processed, batches_completed
                rows, group_count, _batch_reads, raw_umi_count = result
                writer.writerows(rows); molecules += len(rows); groups_processed += group_count; raw_umis_processed += raw_umi_count; batches_completed += 1; progress("running")
            if workers == 1:
                for group_batch in group_batches:
                    peak_batch_umis = max(peak_batch_umis, sum(len(group) for _, group in group_batch)); batches_submitted += 1
                    consume(_dedup_batch(group_batch, config.umi_method, config.umi_max_distance))
            else:
                pending = deque()
                with executor:
                    for group_batch in group_batches:
                        peak_batch_umis = max(peak_batch_umis, sum(len(group) for _, group in group_batch)); batches_submitted += 1
                        pending.append(executor.submit(_dedup_batch, group_batch, config.umi_method, config.umi_max_distance))
                        if len(pending) >= workers: consume(pending.popleft().result())
                    while pending: consume(pending.popleft().result())
        clustering_seconds = time.monotonic() - clustering_started
        progress("completed")
        return output, {"valid_reads": reads, "molecules": molecules, "duplicates": reads - molecules, "groups": total_groups, "total_groups": total_groups, "components": groups_processed, "raw_umis": raw_umis_processed, "total_raw_umis": total_raw_umis, "requested_workers": requested_workers, "workers": workers, "requested_aggregation_workers": aggregation_workers, "aggregation_workers": aggregation_workers, "batches_submitted": batches_submitted, "batches_completed": batches_completed, "umi_batch_size": config.umi_batch_size, "peak_batch_umis": peak_batch_umis, "aggregation_backend": "external_sort", "sort_memory_mb": config.umi_sort_memory_mb, "aggregation_seconds": aggregation_seconds, "clustering_seconds": clustering_seconds, "wall_seconds": time.monotonic() - started}
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
        shutil.rmtree(sort_dir, ignore_errors=True)

def dedup_sample(config: TagForgeConfig, sample_name: str):
    if config.umi_aggregation_backend == "external_sort":
        return _dedup_sample_external_sort(config, sample_name)
    dirs = sample_dirs(config.output_dir, sample_name)
    logger = sample_logger(sample_name, dirs["logs"] / f"{sample_name}.pipeline.log")
    source = dirs["detail"] / f"{sample_name}.valid_reads.tsv.gz"
    output = dirs["detail"] / f"{sample_name}.molecule_detail.tsv.gz"
    progress_path = dirs["logs"] / f"{sample_name}.dedup_progress.tsv"
    db_path = dirs["tmp"] / f"{sample_name}.umi_counts.sqlite3"
    db_path.unlink(missing_ok=True)
    started = time.monotonic()
    aggregation_started = started
    connection = sqlite3.connect(db_path)
    aggregation_executor = None
    executor = None
    requested_aggregation_workers = config.umi_aggregation_workers or min(config.threads, 4)
    aggregation_workers = requested_aggregation_workers
    requested_workers = config.umi_workers or config.threads
    workers = requested_workers
    barcode1_col = config.target_name("barcode1")
    barcode2_name_col = f"{config.target_name('barcode2')}_name"
    umi_col = config.target_name("umi")
    fields = molecule_fields(config)
    reads = 0
    total_groups = 0
    total_raw_umis = 0
    molecules = 0
    groups_processed = 0
    raw_umis_processed = 0
    batches_submitted = 0
    batches_completed = 0
    clustering_started = started
    progress_fields = [
        "status", "valid_reads", "groups", "total_groups", "group_percent",
        "raw_umis", "total_raw_umis", "molecules", "batches_submitted",
        "batches_completed", "pending_batches", "requested_workers", "workers",
        "aggregation_workers", "sqlite_rows_written",
        "elapsed_seconds", "valid_reads_per_second", "groups_per_second",
        "eta_seconds", "estimated_finish",
    ]

    def record_progress(status: str, pending_batches: int = 0):
        now = time.monotonic()
        phase_start = aggregation_started if status == "aggregating" else clustering_started
        elapsed = now - phase_start
        valid_rate = reads / (now - started) if now > started else 0.0
        rate = groups_processed / elapsed if elapsed else 0.0
        fraction = (
            groups_processed / total_groups
            if total_groups else 0.0 if status in {"aggregating", "aggregated"} else 1.0
        )
        eta = (
            elapsed * (1.0 - fraction) / fraction
            if 0 < fraction < 1 else 0.0
        )
        finish = (
            (datetime.now().astimezone() + timedelta(seconds=eta)).isoformat(timespec="seconds")
            if eta else ""
        )
        row = {
            "status": status, "valid_reads": reads, "groups": groups_processed,
            "total_groups": total_groups, "group_percent": f"{fraction * 100:.2f}",
            "raw_umis": raw_umis_processed, "total_raw_umis": total_raw_umis,
            "molecules": molecules, "batches_submitted": batches_submitted,
            "batches_completed": batches_completed, "pending_batches": pending_batches,
            "requested_workers": requested_workers, "workers": workers,
            "aggregation_workers": aggregation_workers,
            "sqlite_rows_written": total_raw_umis,
            "elapsed_seconds": f"{elapsed:.2f}",
            "valid_reads_per_second": f"{valid_rate:.2f}",
            "groups_per_second": f"{rate:.3f}",
            "eta_seconds": f"{eta:.2f}" if eta else "",
            "estimated_finish": finish,
        }
        write_tsv(progress_path, progress_fields, [row])
        logger.info(
            "dedup_progress\tstatus=%s\tvalid_reads=%s\tgroups=%s\t"
            "total_groups=%s\tgroup_percent=%s\traw_umis=%s\ttotal_raw_umis=%s\t"
            "molecules=%s\tbatches_submitted=%s\tbatches_completed=%s\t"
            "pending_batches=%s\trequested_workers=%s\tworkers=%s\t"
            "aggregation_workers=%s\tsqlite_rows_written=%s\t"
            "elapsed_seconds=%s\tvalid_reads_per_second=%s\tgroups_per_second=%s\t"
            "eta_seconds=%s\testimated_finish=%s",
            row["status"], row["valid_reads"], row["groups"], row["total_groups"],
            row["group_percent"], row["raw_umis"], row["total_raw_umis"],
            row["molecules"], row["batches_submitted"], row["batches_completed"],
            row["pending_batches"], row["requested_workers"], row["workers"],
            row["aggregation_workers"], row["sqlite_rows_written"],
            row["elapsed_seconds"], row["valid_reads_per_second"],
            row["groups_per_second"], row["eta_seconds"] or "NA",
            row["estimated_finish"] or "NA",
        )
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
        if aggregation_workers > 1:
            try:
                aggregation_executor = ProcessPoolExecutor(max_workers=aggregation_workers)
                aggregation_executor.submit(_worker_ready).result()
            except (PermissionError, NotImplementedError, OSError) as exc:
                if aggregation_executor is not None:
                    aggregation_executor.shutdown(wait=True, cancel_futures=True)
                aggregation_executor = None
                aggregation_workers = 1
                logger.warning(
                    "dedup_aggregation_fallback\trequested_workers=%s\tworkers=1\treason=%s",
                    requested_aggregation_workers, type(exc).__name__,
                )
        logger.info(
            "dedup_aggregation_start\tbackend=%s\trequested_workers=%s\tworkers=%s\t"
            "chunk_size=%s\tsqlite_cache_mb=%s",
            "process-pool-preaggregate" if aggregation_workers > 1 else "serial",
            requested_aggregation_workers, aggregation_workers, config.chunk_size,
            config.umi_sqlite_cache_mb,
        )

        def write_aggregated(result):
            nonlocal total_raw_umis
            rows_to_write, _input_reads = result
            total_raw_umis += len(rows_to_write)
            if rows_to_write:
                connection.executemany(
                    "INSERT INTO counts VALUES (?,?,?,?) "
                    "ON CONFLICT DO UPDATE SET n=n+excluded.n",
                    rows_to_write,
                )
                connection.commit()
            record_progress("aggregating", 0)

        batch = []
        pending_aggregation = set()
        for row in open_tsv(source):
            batch.append((row[barcode1_col], row[barcode2_name_col], row[umi_col]))
            reads += 1
            if len(batch) >= config.chunk_size:
                if aggregation_executor is None:
                    write_aggregated(_aggregate_batch(batch))
                else:
                    while len(pending_aggregation) >= aggregation_workers:
                        done, pending_aggregation = wait(
                            pending_aggregation, return_when=FIRST_COMPLETED
                        )
                        for future in done:
                            write_aggregated(future.result())
                    pending_aggregation.add(
                        aggregation_executor.submit(
                            _aggregate_batch,
                            batch,
                        )
                    )
                batch = []
        if batch:
            if aggregation_executor is None:
                write_aggregated(_aggregate_batch(batch))
            else:
                pending_aggregation.add(
                    aggregation_executor.submit(
                        _aggregate_batch,
                        batch,
                    )
                )
        while pending_aggregation:
            done, pending_aggregation = wait(
                pending_aggregation, return_when=FIRST_COMPLETED
            )
            for future in done:
                write_aggregated(future.result())
        aggregation_seconds = time.monotonic() - aggregation_started
        total_groups = connection.execute("SELECT COUNT(*) FROM (SELECT 1 FROM counts GROUP BY b1,b2)").fetchone()[0]
        cursor = connection.execute("SELECT b1,b2,umi,n FROM counts ORDER BY b1,b2")
        group_batches = _group_batches(cursor, config.umi_batch_size)
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
            "total_raw_umis=%s\taggregation_workers=%s\taggregation_seconds=%.3f",
            requested_workers, workers, config.umi_batch_size, workers, config.umi_sqlite_cache_mb,
            reads, total_raw_umis, aggregation_workers, aggregation_seconds,
        )
        record_progress("aggregated", 0)

        with atomic_text(output, config.compression_level) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)

            def consume(result, pending_batches: int = 0):
                nonlocal molecules, groups_processed, raw_umis_processed, batches_completed
                rows, group_count, _batch_reads, raw_umi_count = result
                writer.writerows(rows)
                molecules += len(rows)
                groups_processed += group_count
                raw_umis_processed += raw_umi_count
                batches_completed += 1
                record_progress("running", pending_batches)

            if workers == 1:
                for group_batch in group_batches:
                    peak_batch_umis = max(peak_batch_umis, sum(len(group) for _, group in group_batch))
                    batches_submitted += 1
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
                        batches_submitted += 1
                        if len(pending) >= workers:
                            consume(pending.popleft().result(), len(pending))
                    while pending:
                        consume(pending.popleft().result(), len(pending))
        record_progress("completed", 0)

        clustering_seconds = time.monotonic() - clustering_started
    except Exception:
        record_progress("failed", 0)
        raise
    finally:
        if aggregation_executor is not None:
            aggregation_executor.shutdown(wait=True, cancel_futures=True)
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
        "total_groups": total_groups,
        "raw_umis": raw_umis_processed,
        "total_raw_umis": total_raw_umis,
        "requested_workers": requested_workers,
        "workers": workers,
        "requested_aggregation_workers": requested_aggregation_workers,
        "aggregation_workers": aggregation_workers,
        "batches_submitted": batches_submitted,
        "batches_completed": batches_completed,
        "umi_batch_size": config.umi_batch_size,
        "peak_batch_umis": peak_batch_umis,
        "sqlite_cache_mb": config.umi_sqlite_cache_mb,
        "aggregation_seconds": aggregation_seconds,
        "clustering_seconds": clustering_seconds,
        "wall_seconds": wall_seconds,
    }
