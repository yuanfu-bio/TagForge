"""Create a paired, reproducible prefix benchmark without reading the full run."""
import argparse
import gzip
import hashlib
from itertools import islice
import json
from pathlib import Path
import shutil

import yaml

from tagforge.fastq import paired_fastq


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reads", type=int, default=1000000)
    args = parser.parse_args()
    if args.reads < 1:
        parser.error("--reads must be positive")
    output = args.output.resolve()
    if output.exists():
        parser.error("output already exists; use a new benchmark directory")
    sample = "100-1-N1"
    raw = output / "01_raw" / sample
    raw.mkdir(parents=True)
    paths = [raw / f"{sample}_raw_{mate}.fq.gz" for mate in (1, 2)]
    counts = 0
    with paths[0].open("wb") as f1, paths[1].open("wb") as f2:
        with gzip.GzipFile(fileobj=f1, mode="wb", filename="", mtime=0, compresslevel=3) as h1, gzip.GzipFile(fileobj=f2, mode="wb", filename="", mtime=0, compresslevel=3) as h2:
            for read in islice(paired_fastq(args.r1, args.r2), args.reads):
                h1.write(f"@{read.read_id}\n{read.r1_seq}\n+\n{read.r1_qual}\n".encode())
                h2.write(f"@{read.read_id}\n{read.r2_seq}\n+\n{read.r2_qual}\n".encode())
                counts += 1
    if counts != args.reads:
        raise RuntimeError(f"Expected {args.reads} pairs; source contains only {counts}")
    shutil.copytree(args.template.parent, output / "00_config")
    config = yaml.safe_load(args.template.read_text())
    config["project"]["workdir"] = ".."
    config["samples"] = [{"sample": sample, "r1": str(paths[0].relative_to(output)), "r2": str(paths[1].relative_to(output))}]
    (output / "00_config/config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    source = [{"path": str(p.resolve()), "size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns} for p in (args.r1, args.r2)]
    hashes = {}
    for path in paths:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        hashes[str(path.relative_to(output))] = digest.hexdigest()
    manifest = {"paired_records": counts, "sampling": "prefix", "source": source, "sha256": hashes,
                "limitation": "One sample prefix tests scale; it is not a random population sample or accuracy truth."}
    (output / "fixture.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest), flush=True)


if __name__ == "__main__":
    main()
