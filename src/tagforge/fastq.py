from __future__ import annotations

import gzip
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


class FastqError(ValueError):
    pass


@dataclass(frozen=True)
class PairedRead:
    read_id: str
    r1_seq: str
    r1_qual: str
    r2_seq: str
    r2_qual: str


@contextmanager
def open_text(path: Path, mode: str = "rt", compresslevel: int = 3):
    if str(path).endswith(".gz"):
        handle = gzip.open(path, mode, encoding="utf-8" if "t" in mode else None, compresslevel=compresslevel)
    else:
        handle = open(path, mode, encoding="utf-8" if "t" in mode else None)
    try:
        yield handle
    finally:
        handle.close()


def _records(handle: TextIO, path: Path):
    line_number = 0
    while True:
        header = handle.readline()
        if not header:
            return
        seq = handle.readline(); plus = handle.readline(); qual = handle.readline()
        line_number += 4
        if not seq or not plus or not qual:
            raise FastqError(f"Truncated FASTQ record ending near line {line_number}: {path}")
        if not header.startswith("@") or not plus.startswith("+"):
            raise FastqError(f"Malformed FASTQ record ending near line {line_number}: {path}")
        seq, qual = seq.rstrip("\r\n").upper(), qual.rstrip("\r\n")
        if len(seq) != len(qual):
            raise FastqError(f"Sequence/quality length mismatch near line {line_number}: {path}")
        read_id = header[1:].split()[0]
        if read_id.endswith("/1") or read_id.endswith("/2"):
            read_id = read_id[:-2]
        yield read_id, seq, qual


def paired_fastq(r1: Path, r2: Path) -> Iterator[PairedRead]:
    with open_text(r1) as h1, open_text(r2) as h2:
        it1, it2 = _records(h1, r1), _records(h2, r2)
        while True:
            a = next(it1, None); b = next(it2, None)
            if a is None and b is None:
                return
            if a is None or b is None:
                raise FastqError(f"Paired FASTQ files have different record counts: {r1}, {r2}")
            if a[0] != b[0]:
                raise FastqError(f"Read ID mismatch: R1={a[0]!r}, R2={b[0]!r}")
            yield PairedRead(a[0], a[1], a[2], b[1], b[2])
