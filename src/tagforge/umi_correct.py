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

from .columnar import ParquetFragmentWriter, open_records, resolve_table
from .config import ConfigError, TagForgeConfig
from .io_utils import atomic_text, open_tsv, sample_dirs, write_progress, write_tsv
from .logging_utils import ProgressThrottle, sample_logger
from .pi_seq import filter_pi_seq_molecules


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


def _group_batches(
    cursor: Iterable[tuple], batch_umis: int, component_max_distance: int | None = None,
    *, max_group_umis: int | None = None, max_component_umis: int | None = None,
    stats: dict | None = None,
) -> Iterator[list]:
    """Batch independent components; a connected component may exceed batch_umis.

    Optional budgets reject oversized scopes rather than splitting correction
    edges. The group budget is checked while consuming the aggregation cursor.
    """
    batch = []
    batch_size = 0
    key = None
    group = {}
    stats = stats if stats is not None else {}
    stats.setdefault("peak_group_umis", 0)
    stats.setdefault("peak_component_umis", 0)

    def add_components(group_key, umi_counts):
        components = (umi_counts,) if component_max_distance is None else _umi_components(umi_counts, component_max_distance)
        for component in components:
            size = len(component)
            stats["peak_component_umis"] = max(stats["peak_component_umis"], size)
            if max_component_umis is not None and size > max_component_umis:
                raise ConfigError(
                    f"UMI component for barcode1={group_key[0]!r}, feature={group_key[1]!r} "
                    f"has {size} UMIs, exceeding performance.umi_max_component_umis={max_component_umis}. "
                    "Increase the budget with sufficient memory; connected UMIs cannot be split safely."
                )
            yield group_key, component

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
        if umi not in group and max_group_umis is not None and len(group) >= max_group_umis:
            raise ConfigError(
                f"UMI group for barcode1={b1!r}, feature={b2!r} exceeds "
                f"performance.umi_max_group_umis={max_group_umis} while loading UMI {umi!r}. "
                "Increase the budget with sufficient memory."
            )
        group[umi] = n
        stats["peak_group_umis"] = max(stats["peak_group_umis"], len(group))
    if key is not None:
        yield from add_to_batches(add_components(key, group))
    if batch:
        yield batch


def _dedup_batch(batch: list, method: str, max_distance: int):
    """Process a batch; connected components are never split to fit its target."""
    clusterer = _umi_clusterer(method)
    rows = []
    input_reads = 0
    raw_umis = 0
    for group_key, umi_counts in batch:
        if len(umi_counts) == 1:
            # A lone UMI is its own cluster; skip the UMI-tools call.
            assignments = {next(iter(umi_counts)): next(iter(umi_counts))}
        else:
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
                for row in open_records(source):
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
    source = resolve_table(dirs["corrected"] / f"{sample_name}.valid_reads.parquet",
                           dirs["corrected"] / f"{sample_name}.valid_reads.tsv.gz")
    parquet_output = dirs["detail"] / f"{sample_name}.molecule_detail.parquet"
    progress_path = dirs["logs"] / f"{sample_name}.progress.json"
    started = time.monotonic()
    aggregation_workers = config.umi_aggregation_workers or min(config.threads, 4)
    requested_workers = config.umi_workers or config.threads
    workers = requested_workers
    b1_col, b2_col, umi_col = config.target_name("barcode1"), f"{config.target_name('barcode2')}_name", config.target_name("umi")
    reads = total_raw_umis = total_groups = molecules = groups_processed = raw_umis_processed = 0
    batches_submitted = batches_completed = peak_batch_umis = 0
    group_stats = {"peak_group_umis": 0, "peak_component_umis": 0}
    executor = None
    sort_dir = Path(tempfile.mkdtemp(prefix=".umi-sort.", dir=dirs["detail"]))
    fields = molecule_fields(config)
    progress_throttle = ProgressThrottle()

    def progress(status: str):
        elapsed = time.monotonic() - started
        row = {"status": status, "valid_reads": reads, "groups": total_groups, "total_groups": total_groups, "components": groups_processed, "raw_umis": raw_umis_processed, "total_raw_umis": total_raw_umis, "molecules": molecules, "batches_submitted": batches_submitted, "batches_completed": batches_completed, "requested_workers": requested_workers, "workers": workers, "aggregation_workers": aggregation_workers, "aggregation_records": total_raw_umis, "elapsed_seconds": f"{elapsed:.2f}"}
        if progress_throttle.allow(status):
            row.update(group_stats)
            write_progress(progress_path, sample_name, "dedup", row)
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
            max_group_umis=config.umi_max_group_umis, max_component_umis=config.umi_max_component_umis,
            stats=group_stats,
        )
        output_writer = ParquetFragmentWriter(parquet_output, fields, config.chunk_size, compression_level=getattr(config, "compression_level", 3), compression=getattr(config, "parquet_compression", "snappy"))
        def consume(result):
            nonlocal molecules, groups_processed, raw_umis_processed, batches_completed
            rows, group_count, _batch_reads, raw_umi_count = result
            output_writer.writerows(dict(zip(fields, row)) for row in rows)
            molecules += len(rows); groups_processed += group_count; raw_umis_processed += raw_umi_count; batches_completed += 1; progress("running")
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
        output_writer.close()
        clustering_seconds = time.monotonic() - clustering_started
        progress("completed")
        logger.info("dedup_umi_scope\tpeak_group_umis=%s\tpeak_component_umis=%s\tpeak_batch_umis=%s", group_stats["peak_group_umis"], group_stats["peak_component_umis"], peak_batch_umis)
        if config.pi_seq_enabled:
            pi_stats = filter_pi_seq_molecules(config, sample_name, parquet_output)
            logger.info("pi_seq_summary\ttotal_fb_umis=%s\tmulti_pi_ratio=%s\tmp_ratio=%s", pi_stats["total_fb_umis"], pi_stats["multi_pi_ratio"], pi_stats["mp_ratio"])
        return parquet_output, {"valid_reads": reads, "molecules": molecules, "duplicates": reads - molecules, "groups": total_groups, "total_groups": total_groups, "components": groups_processed, "raw_umis": raw_umis_processed, "total_raw_umis": total_raw_umis, "requested_workers": requested_workers, "workers": workers, "requested_aggregation_workers": aggregation_workers, "aggregation_workers": aggregation_workers, "batches_submitted": batches_submitted, "batches_completed": batches_completed, "umi_batch_size": config.umi_batch_size, "peak_batch_umis": peak_batch_umis, **group_stats, "aggregation_backend": "external_sort", "sort_memory_mb": config.umi_sort_memory_mb, "aggregation_seconds": aggregation_seconds, "clustering_seconds": clustering_seconds, "wall_seconds": time.monotonic() - started}
    except Exception:
        progress("failed")
        raise
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
        shutil.rmtree(sort_dir, ignore_errors=True)


def _quote_ident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def _bucket_dedup_worker(task):
    """Cluster every UMI scope of one hash bucket and commit its Parquet fragment."""
    (bucket_glob, output, fields, method, max_distance, batch_umis, max_group_umis,
     max_component_umis, compression, compression_level) = task
    import duckdb
    from .columnar import write_parquet_columns
    connection = duckdb.connect()
    try:
        connection.execute("PRAGMA threads=1")
        cursor = connection.execute(
            "SELECT b1, b2, umi, n FROM read_parquet(?) ORDER BY b1, b2, umi", [bucket_glob])
        raw_umis = components = groups = 0
        group_stats = {"peak_group_umis": 0, "peak_component_umis": 0}
        last_key = [None]

        def counted():
            for row in cursor.fetchall():
                key = (row[0], row[1])
                if key != last_key[0]:
                    nonlocal_groups[0] += 1
                    last_key[0] = key
                yield row

        nonlocal_groups = [0]
        rows = []
        peak_batch = 0
        for batch in _group_batches(
            counted(), batch_umis, max_distance, max_group_umis=max_group_umis,
            max_component_umis=max_component_umis, stats=group_stats,
        ):
            peak_batch = max(peak_batch, sum(len(group) for _, group in batch))
            batch_rows, batch_components, _reads, batch_raw = _dedup_batch(batch, method, max_distance)
            rows.extend(batch_rows)
            components += batch_components
            raw_umis += batch_raw
        groups = nonlocal_groups[0]
        columns = {field: [row[index] for row in rows] for index, field in enumerate(fields)}
        write_parquet_columns(Path(output), columns, fields, compression_level=compression_level,
                              compression=compression)
        return len(rows), components, raw_umis, groups, peak_batch, group_stats
    finally:
        connection.close()


def _dedup_sample_duckdb(config: TagForgeConfig, sample_name: str):
    """DuckDB aggregation into hash buckets, then bucket-parallel UMI clustering.

    Every (barcode1, feature) scope is hashed into exactly one bucket, so workers
    cluster independent scopes and write their own fragments; the coordinating
    process only schedules buckets and sums counters.
    """
    import duckdb
    from .columnar import parquet_files, write_parquet_columns, _pa
    dirs = sample_dirs(config.output_dir, sample_name)
    logger = sample_logger(sample_name, dirs["logs"] / f"{sample_name}.pipeline.log")
    source = resolve_table(dirs["corrected"] / f"{sample_name}.valid_reads.parquet",
                           dirs["corrected"] / f"{sample_name}.valid_reads.tsv.gz")
    parquet_output = dirs["detail"] / f"{sample_name}.molecule_detail.parquet"
    progress_path = dirs["logs"] / f"{sample_name}.progress.json"
    started = time.monotonic()
    requested_workers = config.umi_workers or config.threads
    workers = max(1, requested_workers)
    b1_col, b2_col, umi_col = (config.target_name("barcode1"),
        f"{config.target_name('barcode2')}_name", config.target_name("umi"))
    fields = molecule_fields(config)
    scratch = Path(tempfile.mkdtemp(prefix=".umi-duckdb.", dir=dirs["detail"]))
    progress_throttle = ProgressThrottle()
    reads = total_raw_umis = total_groups = molecules = components = raw_umis_processed = 0
    buckets_done = 0
    peak_batch_umis = 0
    group_stats = {"peak_group_umis": 0, "peak_component_umis": 0}
    bucket_count = max(64, workers * 4)
    compression = getattr(config, "parquet_compression", "snappy")
    compression_level = getattr(config, "compression_level", 3)

    def progress(status: str):
        row = {"status": status, "valid_reads": reads, "groups": total_groups, "total_groups": total_groups,
               "components": components, "raw_umis": raw_umis_processed, "total_raw_umis": total_raw_umis,
               "molecules": molecules, "batches_submitted": bucket_count, "batches_completed": buckets_done,
               "requested_workers": requested_workers, "workers": workers,
               "aggregation_records": total_raw_umis, "elapsed_seconds": f"{time.monotonic() - started:.2f}"}
        if progress_throttle.allow(status):
            row.update(group_stats)
            write_progress(progress_path, sample_name, "dedup", row)
            logger.info("dedup_progress\tstatus=%s\tvalid_reads=%s\tgroups=%s\ttotal_groups=%s\traw_umis=%s\t"
                        "total_raw_umis=%s\tmolecules=%s\tbatches_completed=%s\taggregation_records=%s\telapsed_seconds=%s",
                        status, reads, components, total_groups, raw_umis_processed, total_raw_umis,
                        molecules, buckets_done, total_raw_umis, row["elapsed_seconds"])

    executor = None
    try:
        progress("aggregating")
        if parquet_output.exists():
            shutil.rmtree(parquet_output) if parquet_output.is_dir() else parquet_output.unlink()
        parquet_output.mkdir(parents=True)
        if not (source.is_dir() or source.suffix == ".parquet"):
            raise ConfigError("duckdb UMI aggregation requires a Parquet valid-reads table; "
                              "use performance.umi_aggregation_backend: external_sort for legacy TSV input")
        files = [str(path) for path in parquet_files(source)]
        _, pq = _pa()
        reads = sum(pq.ParquetFile(path).metadata.num_rows for path in files)
        aggregated = scratch / "aggregated"
        connection = duckdb.connect()
        try:
            connection.execute(f"PRAGMA threads={max(1, config.threads)}")
            connection.execute(f"PRAGMA memory_limit='{config.umi_duckdb_memory_mb}MB'")
            connection.execute("PRAGMA temp_directory=?", [str(scratch / "spill")])
            connection.execute("PRAGMA preserve_insertion_order=false")
            quoted = {name: _quote_ident(name) for name in (b1_col, b2_col, umi_col)}
            logger.info("dedup_aggregation_start\tbackend=duckdb\tthreads=%s\tmemory_mb=%s\tbuckets=%s\ttmp_dir=%s",
                        config.threads, config.umi_duckdb_memory_mb, bucket_count, scratch)
            connection.execute(
                f"COPY (SELECT (hash({quoted[b1_col]}, {quoted[b2_col]}) % {bucket_count})::INTEGER AS bucket, "
                f"{quoted[b1_col]} AS b1, {quoted[b2_col]} AS b2, {quoted[umi_col]} AS umi, count(*)::BIGINT AS n "
                f"FROM read_parquet(?) GROUP BY ALL) TO '{str(aggregated).replace(chr(39), chr(39) * 2)}' "
                f"(FORMAT parquet, PARTITION_BY (bucket), COMPRESSION zstd)",
                [files])
        finally:
            connection.close()
        bucket_dirs = sorted(aggregated.glob("bucket=*"), key=lambda path: int(path.name.split("=")[1]))
        for bucket_dir in bucket_dirs:
            total_raw_umis += sum(pq.ParquetFile(path).metadata.num_rows for path in bucket_dir.glob("*.parquet"))
        aggregation_seconds = time.monotonic() - started
        logger.info("dedup_aggregation_complete\tbackend=duckdb\treads=%s\tunique_raw_umis=%s\tseconds=%.3f",
                    reads, total_raw_umis, aggregation_seconds)
        progress("aggregated")
        clustering_started = time.monotonic()
        tasks = [
            (str(bucket_dir / "*.parquet"),
             str(parquet_output / f"part-{int(bucket_dir.name.split('=')[1]):08d}.parquet"),
             fields, config.umi_method, config.umi_max_distance, config.umi_batch_size,
             config.umi_max_group_umis, config.umi_max_component_umis, compression, compression_level)
            for bucket_dir in bucket_dirs
        ]
        logger.info("dedup_parallel_start\tbackend=umi_tools-UMIClusterer\trequested_workers=%s\tworkers=%s\t"
                    "buckets=%s\treads=%s\ttotal_raw_umis=%s\taggregation_seconds=%.3f",
                    requested_workers, workers, len(tasks), reads, total_raw_umis, aggregation_seconds)

        def consume(result):
            nonlocal molecules, components, raw_umis_processed, total_groups, buckets_done, peak_batch_umis
            rows, bucket_components, bucket_raw, bucket_groups, bucket_peak, bucket_stats = result
            molecules += rows
            components += bucket_components
            raw_umis_processed += bucket_raw
            total_groups += bucket_groups
            peak_batch_umis = max(peak_batch_umis, bucket_peak)
            for key, value in bucket_stats.items():
                group_stats[key] = max(group_stats[key], value)
            buckets_done += 1
            progress("running")

        if workers == 1 or len(tasks) <= 1:
            for task in tasks:
                consume(_bucket_dedup_worker(task))
        else:
            executor = ProcessPoolExecutor(max_workers=workers)
            pending = deque()
            for task in tasks:
                pending.append(executor.submit(_bucket_dedup_worker, task))
                if len(pending) >= workers * 2:
                    consume(pending.popleft().result())
            while pending:
                consume(pending.popleft().result())
        if not any(parquet_output.glob("*.parquet")):
            write_parquet_columns(parquet_output / "part-00000000.parquet", {field: [] for field in fields},
                                  fields, compression_level=compression_level, compression=compression)
        clustering_seconds = time.monotonic() - clustering_started
        progress("completed")
        logger.info("dedup_umi_scope\tpeak_group_umis=%s\tpeak_component_umis=%s\tpeak_batch_umis=%s",
                    group_stats["peak_group_umis"], group_stats["peak_component_umis"], peak_batch_umis)
        if config.pi_seq_enabled:
            pi_stats = filter_pi_seq_molecules(config, sample_name, parquet_output)
            logger.info("pi_seq_summary\ttotal_fb_umis=%s\tmulti_pi_ratio=%s\tmp_ratio=%s",
                        pi_stats["total_fb_umis"], pi_stats["multi_pi_ratio"], pi_stats["mp_ratio"])
        return parquet_output, {
            "valid_reads": reads, "molecules": molecules, "duplicates": reads - molecules,
            "groups": total_groups, "total_groups": total_groups, "components": components,
            "raw_umis": raw_umis_processed, "total_raw_umis": total_raw_umis,
            "requested_workers": requested_workers, "workers": workers,
            "requested_aggregation_workers": config.threads, "aggregation_workers": config.threads,
            "batches_submitted": len(tasks), "batches_completed": buckets_done,
            "umi_batch_size": config.umi_batch_size, "peak_batch_umis": peak_batch_umis, **group_stats,
            "aggregation_backend": "duckdb", "duckdb_memory_mb": config.umi_duckdb_memory_mb,
            "aggregation_seconds": aggregation_seconds, "clustering_seconds": clustering_seconds,
            "wall_seconds": time.monotonic() - started,
        }
    except Exception:
        progress("failed")
        raise
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
        shutil.rmtree(scratch, ignore_errors=True)


def dedup_sample(config: TagForgeConfig, sample_name: str):
    if config.umi_aggregation_backend == "duckdb":
        return _dedup_sample_duckdb(config, sample_name)
    if config.umi_aggregation_backend == "external_sort":
        return _dedup_sample_external_sort(config, sample_name)
    dirs = sample_dirs(config.output_dir, sample_name)
    logger = sample_logger(sample_name, dirs["logs"] / f"{sample_name}.pipeline.log")
    source = resolve_table(dirs["corrected"] / f"{sample_name}.valid_reads.parquet",
                           dirs["corrected"] / f"{sample_name}.valid_reads.tsv.gz")
    parquet_output = dirs["detail"] / f"{sample_name}.molecule_detail.parquet"
    progress_path = dirs["logs"] / f"{sample_name}.progress.json"
    scratch = Path(tempfile.mkdtemp(prefix=".umi-sqlite.", dir=dirs["detail"]))
    db_path = scratch / "umi_counts.sqlite3"
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
    group_stats = {"peak_group_umis": 0, "peak_component_umis": 0}
    clustering_started = started
    progress_throttle = ProgressThrottle()

    def record_progress(status: str, pending_batches: int = 0):
        now = time.monotonic()
        phase_start = aggregation_started if status == "aggregating" else clustering_started
        elapsed = now - phase_start
        valid_rate = reads / (now - started) if now > started else 0.0
        rate = groups_processed / elapsed if elapsed else 0.0
        fraction = (
            raw_umis_processed / total_raw_umis
            if total_raw_umis else 0.0 if status in {"aggregating", "aggregated"} else 1.0
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
            "status": status, "valid_reads": reads, "groups": total_groups,
            "total_groups": total_groups, "components": groups_processed, "raw_umi_percent": f"{fraction * 100:.2f}",
            "raw_umis": raw_umis_processed, "total_raw_umis": total_raw_umis,
            "molecules": molecules, "batches_submitted": batches_submitted,
            "batches_completed": batches_completed, "pending_batches": pending_batches,
            "requested_workers": requested_workers, "workers": workers,
            "aggregation_workers": aggregation_workers,
            "sqlite_rows_written": total_raw_umis,
            "elapsed_seconds": f"{elapsed:.2f}",
            "valid_reads_per_second": f"{valid_rate:.2f}",
            "components_per_second": f"{rate:.3f}",
            "eta_seconds": f"{eta:.2f}" if eta else "",
            "estimated_finish": finish,
        }
        if not progress_throttle.allow(status):
            return
        row.update(group_stats)
        write_progress(progress_path, sample_name, "dedup", row)
        logger.info(
            "dedup_progress\tstatus=%s\tvalid_reads=%s\tgroups=%s\t"
            "total_groups=%s\traw_umi_percent=%s\traw_umis=%s\ttotal_raw_umis=%s\t"
            "molecules=%s\tbatches_submitted=%s\tbatches_completed=%s\t"
            "pending_batches=%s\trequested_workers=%s\tworkers=%s\t"
            "aggregation_workers=%s\tsqlite_rows_written=%s\t"
            "elapsed_seconds=%s\tvalid_reads_per_second=%s\tcomponents_per_second=%s\t"
            "eta_seconds=%s\testimated_finish=%s",
            row["status"], row["valid_reads"], row["groups"], row["total_groups"],
            row["raw_umi_percent"], row["raw_umis"], row["total_raw_umis"],
            row["molecules"], row["batches_submitted"], row["batches_completed"],
            row["pending_batches"], row["requested_workers"], row["workers"],
            row["aggregation_workers"], row["sqlite_rows_written"],
            row["elapsed_seconds"], row["valid_reads_per_second"],
            row["components_per_second"], row["eta_seconds"] or "NA",
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
        for row in open_records(source):
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
        total_raw_umis = connection.execute("SELECT COUNT(*) FROM counts").fetchone()[0]
        total_groups = connection.execute("SELECT COUNT(*) FROM (SELECT 1 FROM counts GROUP BY b1,b2)").fetchone()[0]
        cursor = connection.execute("SELECT b1,b2,umi,n FROM counts ORDER BY b1,b2")
        group_batches = _group_batches(
            cursor, config.umi_batch_size, config.umi_max_distance,
            max_group_umis=config.umi_max_group_umis, max_component_umis=config.umi_max_component_umis,
            stats=group_stats,
        )
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

        output_writer = ParquetFragmentWriter(parquet_output, fields, config.chunk_size, compression_level=getattr(config, "compression_level", 3), compression=getattr(config, "parquet_compression", "snappy"))

        def consume(result, pending_batches: int = 0):
            nonlocal molecules, groups_processed, raw_umis_processed, batches_completed
            rows, group_count, _batch_reads, raw_umi_count = result
            output_writer.writerows(dict(zip(fields, row)) for row in rows)
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
        output_writer.close()
        record_progress("completed", 0)

        clustering_seconds = time.monotonic() - clustering_started
        logger.info("dedup_umi_scope\tpeak_group_umis=%s\tpeak_component_umis=%s\tpeak_batch_umis=%s", group_stats["peak_group_umis"], group_stats["peak_component_umis"], peak_batch_umis)
    except Exception:
        record_progress("failed", 0)
        raise
    finally:
        if aggregation_executor is not None:
            aggregation_executor.shutdown(wait=True, cancel_futures=True)
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
        connection.close()
        shutil.rmtree(scratch)
    wall_seconds = time.monotonic() - started
    if config.pi_seq_enabled:
        pi_stats = filter_pi_seq_molecules(config, sample_name, parquet_output)
        logger.info("pi_seq_summary\ttotal_fb_umis=%s\tmulti_pi_ratio=%s\tmp_ratio=%s", pi_stats["total_fb_umis"], pi_stats["multi_pi_ratio"], pi_stats["mp_ratio"])
    return parquet_output, {
        "valid_reads": reads,
        "molecules": molecules,
        "duplicates": reads - molecules,
        "groups": total_groups,
        "components": groups_processed,
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
        **group_stats,
        "sqlite_cache_mb": config.umi_sqlite_cache_mb,
        "aggregation_seconds": aggregation_seconds,
        "clustering_seconds": clustering_seconds,
        "wall_seconds": wall_seconds,
    }
