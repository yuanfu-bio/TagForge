"""Measure a full run, including sampled process-tree RSS and scratch disk use.

GNU time provides its independent resource record. Sampled peaks can miss short
spikes and include output files staged for publication under 00_logs/checkpoints/.work.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from tagforge.config import load_config


def tree_rss(pid):
    seen, pending, total = set(), [pid], 0
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        try:
            status = Path(f"/proc/{current}/status").read_text()
            total += next((int(line.split()[1]) * 1024 for line in status.splitlines() if line.startswith("VmRSS:")), 0)
            children = Path(f"/proc/{current}/task/{current}/children").read_text()
            pending.extend(int(child) for child in children.split())
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return total


def scratch_bytes(roots):
    apparent = allocated = 0
    for root in roots:
        for directory, _, files in os.walk(root, followlinks=False):
            for name in files:
                path = Path(directory) / name
                try:
                    if path.is_symlink():
                        continue
                    stat = path.stat()
                    apparent += stat.st_size
                    allocated += stat.st_blocks * 512
                except FileNotFoundError:
                    continue
    return apparent, allocated


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--record-dir", type=Path, required=True)
    args = parser.parse_args()
    source, config_path, record = args.source.resolve(), args.config.resolve(), args.record_dir.resolve()
    record.mkdir(parents=True, exist_ok=True)
    config = load_config(config_path)
    scratch = [config.output_dir / s.sample / "00_logs/checkpoints/.work" for s in config.samples]
    environment = dict(os.environ, PYTHONPATH=str(source / "src"), PYTHONDONTWRITEBYTECODE="1")
    environment["PATH"] = str(Path(sys.executable).parent) + os.pathsep + environment.get("PATH", "")
    command = ["/usr/bin/time", "-v", "-o", str(record / "time.txt"), sys.executable,
               "-m", "tagforge", "run", "--config", str(config_path), "--overwrite"]
    start = time.monotonic()
    rss = apparent = allocated = observations = 0
    with (record / "run.log").open("w") as log:
        process = subprocess.Popen(command, cwd=source, env=environment, stdout=log, stderr=subprocess.STDOUT)
        while process.poll() is None:
            rss = max(rss, tree_rss(process.pid))
            current_apparent, current_allocated = scratch_bytes(scratch)
            apparent, allocated = max(apparent, current_apparent), max(allocated, current_allocated)
            observations += 1
            time.sleep(0.5)
        exit_code = process.wait()
    summary = {"exit_code": exit_code, "wall_seconds": time.monotonic() - start,
               "sampled_process_tree_rss_peak_bytes": rss,
               "sampled_scratch_apparent_peak_bytes": apparent, "sampled_scratch_allocated_peak_bytes": allocated,
               "sampling_interval_seconds": 0.5, "observations": observations,
               "source": str(source), "config": str(config_path), "command": command,
               "limitations": "Sampled peaks can miss short spikes; summed RSS counts shared pages more than once. Scratch excludes symlink targets."}
    (record / "resources.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
