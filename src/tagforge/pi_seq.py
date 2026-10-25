from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

from .io_utils import sample_dirs, write_tsv


def filter_pi_seq_molecules(config, sample_name: str, source: Path):
    """Enforce one PB per FB/UMI molecule and emit PI-seq QC metrics.

    Multi-PB molecules retain their highest-read PB only when it exceeds the
    configured dominance threshold; ties on read count keep the smaller PB.
    The grouping runs in DuckDB, which parallelises and spills to disk.
    """
    import duckdb
    from .columnar import parquet_files, write_parquet_columns
    dirs = sample_dirs(config.output_dir, sample_name)
    qc_path = dirs["detail"] / f"{sample_name}.pi_seq_qc.tsv"
    scratch = Path(tempfile.mkdtemp(prefix=".pi-seq.", dir=dirs["detail"]))
    filtered_parquet = dirs["detail"] / f"{sample_name}.molecule_detail.rmMP.parquet"
    pb_col = config.target_name("barcode1")
    fb_col = f"{config.target_name('barcode2')}_name"
    fields = [pb_col, fb_col, "corrected_umi", "reads_count", "raw_umi_count"]
    quote = lambda name: '"' + name.replace('"', '""') + '"'
    compression = getattr(config, "parquet_compression", "snappy")
    compression_level = getattr(config, "compression_level", 3)
    connection = duckdb.connect()
    try:
        threads = max(1, getattr(config, "threads", 1))
        memory_mb = getattr(config, "umi_duckdb_memory_mb", 65536)
        connection.execute(f"PRAGMA threads={threads}")
        connection.execute(f"PRAGMA memory_limit='{memory_mb}MB'")
        connection.execute("PRAGMA temp_directory=?", [str(scratch / "spill")])
        connection.execute("PRAGMA preserve_insertion_order=false")
        if source.suffix == ".parquet" or source.is_dir():
            files = [str(path) for path in parquet_files(source)]
        else:  # explicit legacy TSV source
            from .columnar import tsv_to_parquet
            files = [str(tsv_to_parquet(source, scratch / "legacy.parquet"))]
        if filtered_parquet.exists():
            shutil.rmtree(filtered_parquet) if filtered_parquet.is_dir() else filtered_parquet.unlink()
        filtered_parquet.mkdir(parents=True)
        threshold = float(config.pi_seq_dominance_threshold)
        # Hash-partition by (FB, UMI) so each group lives in one bucket and every
        # bucket's aggregation fits comfortably in memory.
        bucket_count = 64
        partitioned = scratch / "by_molecule"
        partitioned_sql = str(partitioned).replace("'", "''")
        connection.execute(
            f"COPY (SELECT (hash({quote(fb_col)}, corrected_umi) % {bucket_count})::INTEGER AS bucket, "
            f"{quote(pb_col)} AS pb, {quote(fb_col)} AS fb, corrected_umi AS umi, "
            f"reads_count::BIGINT AS reads, raw_umi_count::BIGINT AS raw FROM read_parquet(?)) "
            f"TO '{partitioned_sql}' (FORMAT parquet, PARTITION_BY (bucket), COMPRESSION zstd)", [files])
        total_molecules = total_fb_umis = multi = dominant = retained = 0
        written = 0
        for bucket in sorted(int(path.name.split("=")[1]) for path in partitioned.glob("bucket=*")):
            glob = str(partitioned / f"bucket={bucket}" / "*.parquet")
            # Winner: highest read count, then the smallest PB (the filter's long-standing order).
            connection.execute(
                "CREATE OR REPLACE TEMP TABLE kept AS SELECT cnt, tot, best.pb AS pb, fb, umi, "
                "best.reads AS reads, best.raw AS raw FROM ("
                "SELECT fb, umi, count(*) AS cnt, sum(reads) AS tot, "
                "arg_min({'pb': pb, 'reads': reads, 'raw': raw}, {'r': -reads, 'p': pb}) AS best "
                "FROM read_parquet(?) GROUP BY fb, umi)", [glob])
            molecules, groups, bucket_multi, bucket_dominant, bucket_retained = connection.execute(
                "SELECT coalesce(sum(cnt), 0), count(*), "
                "count(*) FILTER (WHERE cnt > 1), "
                "count(*) FILTER (WHERE cnt > 1 AND reads::DOUBLE / tot::DOUBLE > ?), "
                "count(*) FILTER (WHERE cnt = 1 OR reads::DOUBLE / tot::DOUBLE > ?) FROM kept",
                [threshold, threshold]).fetchone()
            total_molecules += int(molecules)
            total_fb_umis += groups
            multi += bucket_multi
            dominant += bucket_dominant
            retained += bucket_retained
            selected = connection.execute(
                "SELECT pb, fb, umi, reads, raw FROM kept WHERE cnt = 1 OR reads::DOUBLE / tot::DOUBLE > ? "
                "ORDER BY fb, umi", [threshold])
            while True:
                chunk = selected.fetchmany(getattr(config, "chunk_size", 200000))
                if not chunk:
                    break
                write_parquet_columns(
                    filtered_parquet / f"part-{written:08d}.parquet",
                    {field: [row[index] for row in chunk] for index, field in enumerate(fields)},
                    fields, compression_level=compression_level, compression=compression)
                written += 1
        connection.execute("DROP TABLE IF EXISTS kept")
        if not written:
            write_parquet_columns(filtered_parquet / "part-00000000.parquet", {field: [] for field in fields},
                                  fields, compression_level=compression_level, compression=compression)
        removed_fb_umis = total_fb_umis - retained
        row = {
            "sample": sample_name, "total_molecules": total_molecules, "total_fb_umis": total_fb_umis,
            "multi_pi_fb_umis": multi,
            "multi_pi_ratio": f"{multi / total_fb_umis:.6f}" if total_fb_umis else "0.000000",
            "mp_ratio": f"{multi / total_fb_umis:.6f}" if total_fb_umis else "0.000000",
            "dominant_multi_pi_fb_umis": dominant,
            "dominant_ratio": f"{dominant / multi:.6f}" if multi else "0.000000",
            "retained_molecules": retained, "removed_fb_umis": removed_fb_umis,
            "removed_ratio": f"{removed_fb_umis / total_fb_umis:.6f}" if total_fb_umis else "0.000000",
            "dominance_threshold": config.pi_seq_dominance_threshold,
        }
        write_tsv(qc_path, list(row), [row])
        return row
    finally:
        connection.close()
        shutil.rmtree(scratch, ignore_errors=True)
