"""Exercise release gates through an installed wheel (run with Python -I)."""
import argparse
import copy
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

import openpyxl
import yaml
import tagforge
from tagforge.checkpoint import STAGE_NUMBERS, stage_plan
from tagforge.columnar import open_records
from tagforge.config import load_config
from tagforge.io_utils import open_tsv
from tagforge.selected import load_selection


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.out.resolve()
    root.mkdir(parents=True, exist_ok=False)
    assert "site-packages" in Path(tagforge.__file__).parts
    assert tagforge.__version__ == importlib.metadata.version("tagforge")
    executable = Path(sys.executable).parent / "tagforge"
    environment = dict(os.environ, PYTHONNOUSERSITE="1")
    environment.pop("PYTHONPATH", None)
    environment["PATH"] = str(executable.parent) + os.pathsep + environment.get("PATH", "")
    commands, checks = [], []

    def run(arguments, expected=0, message=None):
        command = [str(executable), *map(str, arguments)]
        result = subprocess.run(command, cwd=root, env=environment, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        commands.append({"command": command, "returncode": result.returncode, "output": result.stdout})
        (root / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        assert result.returncode == expected, result.stdout
        if message:
            assert message in result.stdout, result.stdout
        return result

    def save_config(case, data):
        data = copy.deepcopy(data)
        data["project"]["output_dir"] = str(root / case / "02_output")
        directory = root / case
        directory.mkdir()
        path = directory / "config.yaml"
        path.write_text(yaml.safe_dump(data, sort_keys=False))
        return path

    def generations(config):
        directory = config.output_dir / "example/00_logs/checkpoints"
        return {step: json.loads((directory / f"{STAGE_NUMBERS[step]}_{step}.manifest.json").read_text())["generation"]
                for step in stage_plan(config)}

    def molecule_counts(config):
        directory = config.output_dir / "example"
        rows = list(open_records(directory / "03_molecules/example.molecule_detail.parquet"))
        return len(rows), sum(int(row["reads_count"]) for row in rows)

    run(["init-example", "--out", root / "demo"])
    template = yaml.safe_load((root / "demo/00_config/config.yaml").read_text())
    template["project"]["workdir"] = str(root / "demo")
    template["downsample"]["ratios"] = [0.5, 1]
    template["downsample"]["repeats"] = 3
    ordinary_config = None
    for case in ("ordinary", "pi", "mapping"):
        data = copy.deepcopy(template)
        if case == "pi":
            data["library"] = {"type": "pi-seq"}
        elif case == "mapping":
            data["library"] = {"type": "pb-cb"}
            data["barcode2"]["name"] = "CB"
            data["barcode2"]["sequence_only"] = True
            data.pop("barcode2_annotation")
            for segment in data["barcode2"]["segments"]:
                segment.pop("whitelist", None)
                segment["correction"] = {"enabled": False}
            data["pb_cb"] = {"min_parent_reads": 2, "parent_child_ratio": 2}
        path = save_config(case, data)
        config = load_config(path)
        run(["report", "--config", path], expected=2, message="requires a current")
        for step in stage_plan(config):
            run([step, "--config", path])
        assert molecule_counts(config) == (6, 9), case
        before = generations(config)
        for step in stage_plan(config):
            run([step, "--config", path])
        assert generations(config) == before
        run(["summary", "--config", path])
        if case == "mapping":
            directory = config.output_dir / "example/04_matrix"
            assert not list(directory.glob("*.parquet"))
            payload = json.loads((directory / "example.pb_cb_summary.json").read_text())
            assert (payload["resolved_molecules"], payload["retained_molecules"], payload["supporting_reads"]) == (6, 6, 9)
        else:
            selected = load_selection(config, "example")
            rows = list(open_records(config.output_dir / "example" / selected["molecules"]["path"]))
            field = "reads_count" if selected["ratio"] == 1 else "reads_count_at_optimal_downsample"
            assert sum(int(row[field]) for row in rows) == selected["metrics"]["reads_sampled"]
        checks.append({"gate": "independent_stages_and_repeat", "library": case,
                       "stages": list(stage_plan(config)), "molecules": 6, "reads": 9, "result": "OK"})
        if case == "ordinary":
            ordinary_config = config

    # Inspect numbers and units across actual TSV, JSON, workbook and HTML.
    config = ordinary_config
    directory = config.output_dir / "example"
    metric_rows = list(open_tsv(directory / "05_downsample/example.downsample_metrics.tsv"))
    for row in metric_rows:
        reads, molecules, singles = (int(row[key]) for key in ("reads_sampled", "umi_types", "umi_detected_once"))
        expected = 100 * (reads - molecules) / reads if reads else 0
        assert math.isclose(float(row["sequencing_saturation"]), expected, abs_tol=1e-6)
        assert "duplication_ratio" not in row
        expected = 100 * (molecules - singles) / molecules if molecules else 0
        assert math.isclose(float(row["non_singleton_molecule_fraction"]), expected, abs_tol=1e-6)
    workbook = openpyxl.load_workbook(directory / "06_report/example.report.xlsx", data_only=True)
    values = list(workbook["Downsample Metrics"].values)
    for observed, expected in zip(values[1:], metric_rows):
        for key, value in zip(values[0], observed):
            if key == "sample":
                assert value == expected[key]
            else:
                assert math.isclose(float(value), float(expected[key]), abs_tol=1e-6)
    assert len(values) == len(metric_rows) + 1
    workbook.close()
    summaries = list(open_tsv(directory / "05_downsample/example.downsample_curve.tsv"))
    html = (directory / "05_downsample/example.downsample.html").read_text()
    payload = json.loads(re.search(r"const data = (.*);", html).group(1))
    assert payload["saturation"] == [float(row["observed_sequencing_saturation"]) for row in summaries]
    assert payload["non_singleton"] == [float(row["observed_non_singleton_molecule_fraction"]) for row in summaries]
    run(["quick-test", "--config", config.path, "--reads", "9"])
    metrics = {row["metric"]: row["value"] for row in open_tsv(directory / "06_report/example.quick_test.tsv") if row["scope"] == "sample"}
    assert math.isclose(float(metrics["duplication_ratio"]), 100 / 3, abs_tol=1e-8)
    assert metrics["duplication_ratio_unit"] == "percent"
    assert float(metrics["combined_valid_rate"]) == 1
    checks.append({"gate": "metric_units_and_artifact_values", "result": "OK"})

    # Old completion state is rejected by the installed downstream CLI, then rebuilt.
    checkpoint = directory / "00_logs/checkpoints/01_extract.manifest.json"
    checkpoint.unlink()
    checkpoint.with_name("extract.done").write_text("tagforge_version=0.1.13\n")
    run(["report", "--config", config.path], expected=2, message="requires a current")
    run(["extract", "--config", config.path])
    for step in stage_plan(config)[1:]:
        assert not (directory / "00_logs/checkpoints" / f"{STAGE_NUMBERS[step]}_{step}.manifest.json").exists()
        run([step, "--config", config.path])
    assert molecule_counts(config) == (6, 9)
    checks.append({"gate": "old_completion_rejected_and_rebuilt", "result": "OK"})

    # Export must reject old counts after a parameter changes the molecule result.
    data = yaml.safe_load(config.path.read_text())
    data["correction_umi"]["max_distance"] = 0
    stale_path = root / "stale-export.yaml"
    stale_path.write_text(yaml.safe_dump(data, sort_keys=False))
    rejected = root / "stale-export.tsv"
    run(["export-matrix", "--config", stale_path, "--sample", "example", "--format", "tsv", "--out", rejected],
        expected=2, message="current matrix stage")
    assert not rejected.exists()
    data["downsample"]["enabled"] = False
    rebuilt = save_config("export-rebuilt", data)
    run(["run", "--config", rebuilt])
    assert molecule_counts(load_config(rebuilt)) == (7, 9)
    exported = root / "rebuilt-export.tsv"
    run(["export-matrix", "--config", rebuilt, "--sample", "example", "--format", "tsv", "--out", exported])
    assert sum(int(value) for line in exported.read_text().splitlines()[1:]
               for value in line.split("\t")[1:]) == 7
    checks.append({"gate": "stale_export_rejected_then_recomputed_counts_exported", "result": "OK"})

    # Legacy flat configuration must execute, with equivalent molecule counts.
    data = copy.deepcopy(template)
    data["segments"] = []
    for role in ("barcode1", "barcode2", "umi"):
        section = data.pop(role)
        data["segments"].extend(dict(segment, target=role) for segment in section["segments"])
    path = save_config("legacy", data)
    run(["run", "--config", path])
    assert molecule_counts(load_config(path)) == (6, 9)
    checks.append({"gate": "legacy_flat_configuration", "result": "OK"})

    # Configured quality failure stops both explicit quick-test and full run.
    data = copy.deepcopy(template)
    data["umi"]["segments"][0]["start"] = 100
    data["quick_test"]["min_valid_rate"] = 0.5
    data["downsample"]["enabled"] = False
    path = save_config("threshold", data)
    run(["quick-test", "--config", path], expected=2, message="min_valid_rate")
    run(["run", "--config", path], expected=2, message="min_valid_rate")
    assert not (root / "threshold/02_output/example/00_logs/checkpoints/01_extract.manifest.json").exists()
    checks.append({"gate": "configured_threshold_stops_cli_pipeline", "result": "OK"})
    data["quick_test"].pop("min_valid_rate")
    path = save_config("zero", data)
    run(["run", "--config", path])
    assert molecule_counts(load_config(path)) == (0, 0)
    checks.append({"gate": "zero_valid_prefix_remains_diagnostic_without_threshold", "result": "OK"})

    # Separate samples publish concurrently; the last writer must retain both.
    data = copy.deepcopy(template)
    data["downsample"]["enabled"] = False
    data["samples"] = [{"sample": name,
        "r1": str(root / "demo/01_raw/example/example_raw_1.fq.gz"),
        "r2": str(root / "demo/01_raw/example/example_raw_2.fq.gz")}
        for name in ("parallel-A", "parallel-B")]
    path = save_config("parallel", data)
    jobs = []
    for sample in ("parallel-A", "parallel-B"):
        command = [str(executable), "run", "--config", str(path), "--sample", sample]
        jobs.append((command, subprocess.Popen(command, cwd=root, env=environment,
                     text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)))
    for command, process in jobs:
        output, _ = process.communicate(timeout=120)
        commands.append({"command": command, "returncode": process.returncode, "output": output})
        (root / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        assert process.returncode == 0, output
    workbook = openpyxl.load_workbook(root / "parallel/02_output/00_report/TagForge_batch_report.xlsx", data_only=True)
    count_rows = list(workbook["Counts_raw"].values)[1:]
    assert [row[0] for row in count_rows] == ["parallel-A", "parallel-B"]
    assert sum(float(value) for row in count_rows for value in row[1:]) == 12
    workbook.close()
    checks.append({"gate": "different_samples_run_concurrently_without_losing_summary_rows", "result": "OK"})
    summary = {"version": tagforge.__version__, "python": sys.version, "package_path": str(tagforge.__file__),
               "result": "OK", "commands": len(commands), "checks": checks}
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
