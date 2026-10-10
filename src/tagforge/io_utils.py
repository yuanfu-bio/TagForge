from __future__ import annotations

import csv
import gzip
import json
from datetime import datetime, timezone
import os
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Iterator, Mapping


SUBDIRS = {
    "logs": "00_logs", "checkpoint": "00_logs/checkpoints", "extracted": "01_extracted",
    "corrected": "02_corrected", "matrix": "04_matrix", "detail": "03_molecules",
    "downsample": "05_downsample", "report": "06_report",
}


def batch_report_dir(output_dir: Path) -> Path:
    """Keep all batch reports within the configured output root."""
    return output_dir / "00_report"


def sample_paths(output_dir: Path, sample: str):
    root = output_dir / sample
    return {key: root / value for key, value in SUBDIRS.items()}


def sample_dirs(output_dir: Path, sample: str):
    paths = sample_paths(output_dir, sample)
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


@contextmanager
def atomic_text(path: Path, gzip_level: int = 3):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(path) + ".tmp")
    if str(path).endswith(".gz"):
        handle = gzip.open(tmp, "wt", encoding="utf-8", newline="", compresslevel=gzip_level)
    else:
        handle = open(tmp, "w", encoding="utf-8", newline="")
    try:
        yield handle
        handle.flush()
        handle.close()
        os.replace(tmp, path)
    except BaseException:
        handle.close()
        tmp.unlink(missing_ok=True)
        raise


def open_tsv(path: Path) -> Iterator[dict]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def _physical_position(handle) -> int:
    """Return compressed bytes consumed for gzip, normal bytes otherwise."""
    binary = getattr(handle, "buffer", handle)
    physical = getattr(binary, "fileobj", None)
    return physical.tell() if physical is not None else binary.tell()


def tsv_batches(path: Path, batch_size: int):
    """Yield bounded TSV row batches and an approximate physical input fraction."""
    if batch_size < 1:
        raise ValueError("batch_size must be >= 1")
    opener = gzip.open if str(path).endswith(".gz") else open
    total_bytes = path.stat().st_size
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        pending = None
        while True:
            batch = []
            if pending is not None:
                batch.append(pending)
                pending = None
            while len(batch) < batch_size:
                row = next(rows, None)
                if row is None:
                    break
                batch.append(row)
            if not batch:
                return
            pending = next(rows, None)
            consumed = _physical_position(handle)
            fraction = 1.0 if pending is None else (
                min(0.9999, consumed / total_bytes) if total_bytes else 0.9999
            )
            yield batch, fraction


def write_tsv(path: Path, fields: list[str], rows: Iterable[Mapping], gzip_level: int = 3):
    with atomic_text(path, gzip_level) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def touch_checkpoint(path: Path, version: str):
    with atomic_text(path) as handle:
        handle.write(f"tagforge_version={version}\n")


def checkpoint_version(version: str, step: str) -> str:
    """Version output schemas independently so metric changes rerun only affected stages."""
    schema = {
        "downsample": "downsample-schema-2",
        "report": "report-schema-2",
    }.get(step)
    return f"{version}:{schema}" if schema else version


def output_ready(path: Path) -> bool:
    """Accept normal files and nonempty Parquet dataset directories."""
    return (path.is_file() and path.stat().st_size > 0) or (
        path.is_dir() and any(path.glob("*.parquet"))
    )


def step_complete(checkpoint: Path, outputs: list[Path], overwrite: bool, version: str) -> bool:
    return (
        not overwrite
        and checkpoint.is_file()
        and checkpoint.read_text(encoding="utf-8").strip() == f"tagforge_version={version}"
        and all(output_ready(path) for path in outputs)
    )


def write_json(path: Path, value):
    with atomic_text(path) as handle:
        json.dump(value, handle, ensure_ascii=False)
        handle.write("\n")


def write_progress(path: Path, sample: str, stage: str, values):
    """Publish the latest stage snapshot; pipeline.log retains the history."""
    write_json(path, {"schema": 1, "sample": sample, "stage": stage,
                      "updated_at": datetime.now(timezone.utc).isoformat(), **values})
    for name in ("extraction", "correction", "dedup", "downsample"):
        (path.parent / f"{sample}.{name}_progress.tsv").unlink(missing_ok=True)
