"""Run with the installed environment's Python -I, outside its source checkout."""
import argparse
import copy
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys

import pyarrow.parquet as pq
import yaml
import tagforge
from tagforge.columnar import open_records
from tagforge.config import load_config
from tagforge.checkpoint import stage_current


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=args.resume)
    executable = Path(sys.executable).parent / "tagforge"
    package_path = Path(tagforge.__file__).resolve()
    assert "site-packages" in package_path.parts, package_path
    assert tagforge.__version__ == importlib.metadata.version("tagforge")
    environment = dict(os.environ, PYTHONNOUSERSITE="1")
    environment.pop("PYTHONPATH", None)
    environment["PATH"] = str(executable.parent) + os.pathsep + environment.get("PATH", "")
    commands = json.loads((root / "commands.json").read_text()) if args.resume else []
    def run(arguments, expected=0):
        command = [str(executable), *map(str, arguments)]
        result = subprocess.run(command, cwd=root, env=environment, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        commands.append({"command": command, "returncode": result.returncode, "output": result.stdout})
        (root / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        assert result.returncode == expected, result.stdout
        return result
    run(["--version"])
    demo = root / "demo"
    if not args.resume:
        run(["init-example", "--out", demo])
    first = hashlib.sha256((demo / "00_config/config.yaml").read_bytes()).hexdigest()
    run(["init-example", "--out", demo], expected=2)
    assert first == hashlib.sha256((demo / "00_config/config.yaml").read_bytes()).hexdigest()
    template = yaml.safe_load((demo / "00_config/config.yaml").read_text())
    cases = []
    for case, library, downsample, trace in (("ordinary", "standard", True, True), ("pi", "pi-seq", False, False), ("mapping", "pb-cb", True, False), ("multi", "standard", False, True)):
        config = copy.deepcopy(template)
        config["project"]["workdir"] = str(demo)
        config["project"]["output_dir"] = str(root / case / "02_output")
        config["library"] = {"type": library}
        config["output"] = {"correction_trace": trace}
        config["downsample"]["enabled"] = downsample
        config["downsample"]["ratios"] = [0.5, 1]
        config["downsample"]["repeats"] = 2
        if library == "pb-cb":
            config["barcode2"]["name"] = "CB"
            config["barcode2"]["sequence_only"] = True
            del config["barcode2_annotation"]
            for segment in config["barcode2"]["segments"]:
                segment.pop("whitelist", None)
                segment["correction"] = {"enabled": False}
            config["pb_cb"] = {"min_parent_reads": 2, "parent_child_ratio": 2}
        samples = ["example"]
        if case == "multi":
            samples = ["demo-A", "demo-B"]
            config["samples"] = [{"sample": sample, "r1": str(demo / "01_raw/example/example_raw_1.fq.gz"), "r2": str(demo / "01_raw/example/example_raw_2.fq.gz")} for sample in samples]
        directory = root / case
        directory.mkdir(exist_ok=args.resume)
        path = directory / "config.yaml"
        path.write_text(yaml.safe_dump(config, sort_keys=False))
        run(["validate-config", "--config", path])
        loaded = load_config(path)
        if not args.resume or not all(stage_current(loaded, sample, "report") for sample in samples):
            run(["run", "--config", path])
        for sample in samples:
            output = directory / "02_output" / sample
            molecules = list(open_records(output / "03_molecules" / f"{sample}.molecule_detail.parquet"))
            assert len(molecules) == 6
            assert sum(int(row["reads_count"]) for row in molecules) == 9
            if library == "pb-cb":
                summary = json.loads((output / "04_matrix" / f"{sample}.pb_cb_summary.json").read_text())
                assert (summary["resolved_molecules"], summary["retained_molecules"], summary["supporting_reads"]) == (6, 6, 9)
            else:
                for format in ("tsv", "mtx", "parquet"):
                    exported = directory / f"{sample}.{format}"
                    run(["export-matrix", "--config", path, "--sample", sample, "--format", format, "--out", exported])
                    if format == "parquet":
                        assert sum(pq.read_table(exported)["molecule_count"].to_pylist()) == 6
                    elif format == "tsv":
                        lines = exported.read_text().splitlines()
                        assert sum(int(value) for line in lines[1:] for value in line.split("\t")[1:]) == 6
                    else:
                        lines = [line for line in exported.read_text().splitlines() if not line.startswith('%')]
                        assert sum(int(line.split()[2]) for line in lines[1:]) == 6
        for mode in ("array", "per-sample"):
            scripts = directory / ("slurm-" + mode)
            run(["make-slurm", "--config", path, "--out", scripts, "--mode", mode])
            files = [p for p in scripts.rglob("*") if p.suffix in (".sh", ".slurm")]
            assert files, scripts
            for script in files:
                subprocess.run(["bash", "-n", str(script)], check=True)
        cases.append({"case": case, "samples": samples, "valid_reads_per_sample": 9, "molecules_per_sample": 6, "result": "OK"})
    versions = {name: importlib.metadata.version(name) for name in ("tagforge", "cutadapt", "umi-tools", "duckdb", "pyarrow", "numpy", "pandas", "PyYAML")}
    summary = {"python": sys.version, "package_path": str(package_path), "versions": versions, "cases": cases,
               "commands": len(commands), "result": "OK", "scope": "Installed wheel CLI, all library types, multi-sample, flags, exports and Slurm syntax"}
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
