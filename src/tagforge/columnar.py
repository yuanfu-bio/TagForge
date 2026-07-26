"""Columnar storage helpers for TagForge intermediate tables."""
from __future__ import annotations

from pathlib import Path
import shutil
from typing import Iterator

from .config import ConfigError
from .io_utils import open_tsv


def _pa():
    try:
        import pyarrow as pa  # type: ignore
        import pyarrow.parquet as pq  # type: ignore
    except ImportError as exc:
        raise ConfigError("Parquet intermediates require pyarrow") from exc
    return pa, pq


def parquet_path(tsv_path: Path) -> Path:
    name = tsv_path.name
    if name.endswith(".tsv.gz"):
        name = name[:-7] + ".parquet"
    elif name.endswith(".tsv"):
        name = name[:-4] + ".parquet"
    return tsv_path.with_name(name)


def write_parquet(path: Path, rows, fields: list[str], batch_size: int = 200000):
    pa, pq = _pa()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.unlink(missing_ok=True)
    writer = None
    batch = []
    try:
        for row in rows:
            batch.append({field: row.get(field, "") for field in fields})
            if len(batch) >= batch_size:
                table = pa.Table.from_pylist(batch)
                writer = writer or pq.ParquetWriter(temporary, table.schema, compression="zstd")
                writer.write_table(table); batch.clear()
        if batch:
            table = pa.Table.from_pylist(batch)
            writer = writer or pq.ParquetWriter(temporary, table.schema, compression="zstd")
            writer.write_table(table)
        elif writer is None:
            writer = pq.ParquetWriter(temporary, pa.schema([(field, pa.string()) for field in fields]), compression="zstd")
    finally:
        if writer is not None: writer.close()
    temporary.replace(path)



def is_parquet_table(path: Path) -> bool:
    """Whether ``path`` is a Parquet file or a nonempty Parquet dataset."""
    return path.is_file() or (path.is_dir() and any(path.glob("*.parquet")))


def parquet_size(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size
    return sum(part.stat().st_size for part in path.glob("*.parquet")) if path.is_dir() else 0


class ParquetFragmentWriter:
    """Crash-safe, append-only Parquet dataset writer.

    Each committed batch is an atomic Zstandard Parquet file.  Unlike a
    single ParquetWriter this remains readable after an interrupted run, and
    resume only needs to discard fragments after the saved commit count.
    """
    def __init__(self, path: Path, fields: list[str], batch_size: int = 200000, *, resume: bool = False, parts: int = 0):
        self.path, self.fields, self.batch_size = path, fields, batch_size
        self.buffer: list[dict] = []
        if resume:
            if not path.is_dir():
                raise ConfigError(f"Parquet resume output is missing: {path}")
            existing = sorted(path.glob("*.parquet"))
            for stale in existing[parts:]:
                stale.unlink()
            self.part = parts
        else:
            if path.is_dir():
                shutil.rmtree(path)
            elif path.exists():
                path.unlink()
            path.mkdir(parents=True, exist_ok=True)
            self.part = 0

    def writerow(self, row: dict):
        self.buffer.append(row)
        if len(self.buffer) >= self.batch_size:
            self.flush()

    def writerows(self, rows):
        self.buffer.extend(rows)
        if len(self.buffer) >= self.batch_size:
            self.flush()

    def flush(self):
        if not self.buffer:
            return
        write_parquet(self.path / f"part-{self.part:08d}.parquet", self.buffer, self.fields, self.batch_size)
        self.buffer.clear()
        self.part += 1

    def close(self):
        self.flush()


def tsv_to_parquet(source: Path, destination: Path | None = None, batch_size: int = 200000) -> Path:
    """Stream a legacy TSV table directly to Zstandard Parquet with DuckDB."""
    destination = destination or parquet_path(source)
    temporary = destination.with_name(destination.name + ".tmp")
    temporary.unlink(missing_ok=True)
    try:
        import duckdb
        connection = duckdb.connect()
        try:
            source_sql = str(source).replace("'", "''")
            target_sql = str(temporary).replace("'", "''")
            connection.execute(
                f"COPY (SELECT * FROM read_csv('{source_sql}', delim=chr(9), header=true, compression='auto', all_varchar=true)) TO '{target_sql}' (FORMAT PARQUET, COMPRESSION ZSTD)"
            )
        finally:
            connection.close()
    except ImportError:
        iterator = open_tsv(source)
        try:
            first = next(iterator)
        except StopIteration:
            write_parquet(destination, [], [])
            return destination
        write_parquet(temporary, _prepend(first, iterator), list(first), batch_size)
    temporary.replace(destination)
    return destination


def _prepend(first, iterator):
    yield first
    yield from iterator


def parquet_batches(path: Path, batch_size: int = 200000) -> Iterator[tuple[list[dict], float]]:
    pa, pq = _pa()
    if path.is_dir():
        import pyarrow.dataset as ds  # type: ignore
        dataset = ds.dataset(path, format="parquet")
        total = dataset.count_rows() or 1
        batches = dataset.scanner(batch_size=batch_size).to_batches()
    else:
        parquet = pq.ParquetFile(path)
        total = parquet.metadata.num_rows or 1
        batches = parquet.iter_batches(batch_size=batch_size)
    completed = 0
    for table in batches:
        rows = [{field: "" if value is None else str(value) for field, value in row.items()} for row in table.to_pylist()]
        completed += len(rows)
        yield rows, min(1.0, completed / total)


def open_records(path: Path):
    if path.suffix == ".parquet" or path.is_dir():
        for rows, _ in parquet_batches(path):
            yield from rows
    else:
        yield from open_tsv(path)



def record_batches(path: Path, batch_size: int = 200000):
    """Yield record batches from either a Parquet intermediate or legacy TSV.gz."""
    if path.suffix == ".parquet" or path.is_dir():
        yield from parquet_batches(path, batch_size)
    else:
        from .io_utils import tsv_batches
        yield from tsv_batches(path, batch_size)
