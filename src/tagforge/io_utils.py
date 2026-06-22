from __future__ import annotations

import csv
import gzip
import json
import os
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Iterator, Mapping


SUBDIRS = {
    "logs": "00_logs", "checkpoint": "01_checkpoint", "extracted": "02_extracted",
    "corrected": "03_corrected", "matrix": "04_matrix", "detail": "05_detail",
    "downsample": "06_downsample", "report": "07_report", "tmp": "08_tmp",
}


def sample_dirs(output_dir: Path, sample: str):
    root = output_dir / sample
    paths = {key: root / value for key, value in SUBDIRS.items()}
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


def write_tsv(path: Path, fields: list[str], rows: Iterable[Mapping], gzip_level: int = 3):
    with atomic_text(path, gzip_level) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def json_compact(value) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def touch_checkpoint(path: Path):
    with atomic_text(path) as handle:
        handle.write("done\n")


def step_complete(checkpoint: Path, outputs: list[Path], overwrite: bool) -> bool:
    return not overwrite and checkpoint.is_file() and all(path.is_file() and path.stat().st_size > 0 for path in outputs)

