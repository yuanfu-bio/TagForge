"""Columnar storage helpers for TagForge intermediate tables."""
from __future__ import annotations

import csv
from pathlib import Path
import shutil
from typing import Iterator

from .config import ConfigError, _integer, validate_parquet_compression
from .io_utils import open_tsv
from .fastq import open_text


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


INTEGER_FIELDS = frozenset({"reads_count", "raw_umi_count", "molecule_count",
                            "shift_distance", "mismatch_distance", "reads_count_at_optimal_downsample"})


def table_schema(fields):
    """Stable table schema, including empty intermediates and legacy migration."""
    pa, _ = _pa()
    if len(fields) != len(set(fields)):
        raise ConfigError("Intermediate table has duplicate field names")
    return pa.schema([(field, pa.int64() if field in INTEGER_FIELDS else pa.string()) for field in fields])


def _typed_row(row, fields):
    result = {}
    for field in fields:
        value = row.get(field, "")
        if field in INTEGER_FIELDS:
            try:
                value = _integer(value, f"Intermediate column {field}")
            except ConfigError as exc:
                raise ConfigError(f"Invalid integer in intermediate column {field}: {value!r}") from exc
        else:
            value = "" if value is None else str(value)
        result[field] = value
    return result


def parquet_files(path):
    """Identify and validate every fragment; incomplete directories fail clearly."""
    files = sorted(path.glob("*.parquet")) if path.is_dir() else [path]
    if not files:
        raise ConfigError(f"Incomplete Parquet dataset: {path}")
    _, pq = _pa()
    schemas = []
    for part in files:
        try:
            schemas.append(pq.ParquetFile(part).schema_arrow)
        except (OSError, ValueError) as exc:
            raise ConfigError(f"Invalid Parquet artifact: {part}") from exc
    if any(schema != schemas[0] for schema in schemas[1:]):
        raise ConfigError(f"Parquet fragments have inconsistent schemas: {path}")
    return files


def resolve_table(preferred, legacy=None):
    """Use the current table, or an explicit legacy source when it is absent."""
    if preferred.exists():
        if preferred.suffix == ".parquet" or preferred.is_dir():
            parquet_files(preferred)
        return preferred
    if legacy is not None and legacy.is_file():
        return legacy
    raise FileNotFoundError(f"Intermediate table is missing: {preferred}")


def sql_reader(source):
    """DuckDB reader sharing the same explicit fragment/legacy convention."""
    quote = lambda value: "'" + str(value).replace("'", "''") + "'"
    if source.suffix == ".parquet" or source.is_dir():
        return "read_parquet([" + ",".join(quote(part) for part in parquet_files(source)) + "])"
    return f"read_csv_auto({quote(source)}, delim=chr(9), header=true, compression='auto', all_varchar=true)"


def write_parquet(path: Path, rows, fields: list[str], batch_size: int = 200000, *, compression_level: int = 3, compression: str = "snappy"):
    pa, pq = _pa()
    compression = validate_parquet_compression(compression)
    options = {"compression": compression}
    if compression == "zstd":
        options["compression_level"] = compression_level
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.unlink(missing_ok=True)
    writer = None
    schema = table_schema(fields)
    batch = []
    try:
        for row in rows:
            batch.append(_typed_row(row, fields))
            if len(batch) >= batch_size:
                table = pa.Table.from_pylist(batch, schema=schema)
                writer = writer or pq.ParquetWriter(temporary, table.schema, **options)
                writer.write_table(table); batch.clear()
        if batch:
            table = pa.Table.from_pylist(batch, schema=schema)
            writer = writer or pq.ParquetWriter(temporary, table.schema, **options)
            writer.write_table(table)
        elif writer is None:
            writer = pq.ParquetWriter(temporary, schema, **options)
    finally:
        if writer is not None: writer.close()
    temporary.replace(path)



def is_parquet_table(path: Path) -> bool:
    """Whether ``path`` is a Parquet file or a nonempty Parquet dataset."""
    return (path.suffix == ".parquet" and path.is_file()) or (path.is_dir() and any(path.glob("*.parquet")))


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
    def __init__(self, path: Path, fields: list[str], batch_size: int = 200000, *, resume: bool = False, parts: int = 0, compression_level: int = 3, compression: str = "snappy"):
        self.path, self.fields, self.batch_size = path, fields, batch_size
        self.compression_level = compression_level
        self.compression = validate_parquet_compression(compression)
        self.buffer: list[dict] = []
        self.rows = 0
        if resume:
            if parts < 0:
                raise ConfigError("Committed Parquet part count cannot be negative")
            if not path.is_dir():
                raise ConfigError(f"Parquet resume output is missing: {path}")
            existing = sorted(path.glob("*.parquet"))
            _, pq = _pa()
            for index in range(parts):
                part = path / f"part-{index:08d}.parquet"
                try:
                    parquet = pq.ParquetFile(part)
                except (OSError, ValueError) as exc:
                    raise ConfigError(f"Committed Parquet fragment is missing or invalid: {part}") from exc
                if parquet.schema_arrow != table_schema(fields):
                    raise ConfigError(f"Committed Parquet schema does not match current table: {part}")
                self.rows += parquet.metadata.num_rows
            committed = {f"part-{index:08d}.parquet" for index in range(parts)}
            for stale in existing:
                if stale.name not in committed:
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
        write_parquet(self.path / f"part-{self.part:08d}.parquet", self.buffer, self.fields, self.batch_size, compression_level=self.compression_level, compression=self.compression)
        self.rows += len(self.buffer)
        self.buffer.clear()
        self.part += 1

    def close(self):
        self.flush()
        if self.part == 0:
            write_parquet(self.path / "part-00000000.parquet", [], self.fields, self.batch_size, compression_level=self.compression_level, compression=self.compression)
            self.part = 1


def tsv_to_parquet(source: Path, destination: Path | None = None, batch_size: int = 200000) -> Path:
    """Stream an explicit legacy TSV source using the current typed schema."""
    destination = destination or parquet_path(source)
    with open_text(source) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ConfigError(f"Legacy table has no header: {source}")
        write_parquet(destination, reader, reader.fieldnames, batch_size)
    return destination


def _prepend(first, iterator):
    yield first
    yield from iterator


def parquet_batches(path: Path, batch_size: int = 200000) -> Iterator[tuple[list[dict], float]]:
    pa, pq = _pa()
    if path.is_dir():
        import pyarrow.dataset as ds  # type: ignore
        dataset = ds.dataset([str(part) for part in parquet_files(path)], format="parquet")
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
