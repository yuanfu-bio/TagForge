"""Explicit references to a single downsampling selection."""
from __future__ import annotations

import json
from pathlib import Path

from .config import ConfigError
from .io_utils import atomic_text, open_tsv


def selection_path(config, sample):
    return config.output_dir / sample / "05_downsample" / f"{sample}.selection.json"


def write_selection(config, sample, optimal, molecules, matrix_paths, owned):
    from .checkpoint import artifact_info
    root = config.output_dir / sample
    def reference(path):
        return {"path": str(path.relative_to(root)), "artifact": artifact_info(path)} if path else None
    payload = {"schema": 2, "sample": sample, "ratio": float(optimal["downsample_ratio"]),
        "criterion": optimal["selection_criterion"],
        "non_singleton_molecule_fraction": float(optimal["non_singleton_molecule_fraction"]),
        "seed": config.downsample_seed, "molecules": reference(molecules),
        "matrix_pairs": reference(matrix_paths[0]) if matrix_paths else None,
        "matrix_summary": reference(matrix_paths[1]) if matrix_paths else None,
        "owned_artifacts": [str(path.relative_to(root)) for path in owned],
        "metrics": {key: int(optimal[key]) for key in ("reads_sampled", "umi_types", "umi_detected_once")}}
    with atomic_text(selection_path(config, sample)) as handle:
        json.dump(payload, handle, sort_keys=True)
    return payload


def load_selection(config, sample):
    from .checkpoint import artifact_info
    path = selection_path(config, sample)
    try:
        selected = json.loads(path.read_text(encoding="utf-8"))
        if selected["schema"] != 2 or selected["sample"] != sample:
            raise ValueError("unsupported selection schema/sample")
        root = config.output_dir / sample
        for key in ("molecules", "matrix_pairs", "matrix_summary"):
            reference = selected[key]
            if reference is None:
                continue
            relative = Path(reference["path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("selection path must be relative to the sample")
            if artifact_info(root / relative) != reference["artifact"]:
                raise ValueError(f"selected {key} has changed")
        if selected["molecules"]["artifact"]["rows"] != selected["metrics"]["umi_types"]:
            raise ValueError("selected molecule count does not match its metrics")
        metrics = selected["metrics"]
        if not (0 <= selected["ratio"] <= 1 and
                0 <= metrics["umi_detected_once"] <= metrics["umi_types"] <= metrics["reads_sampled"]):
            raise ValueError("selection ratio/counts are inconsistent")
        candidates = [row for row in open_tsv(root / "05_downsample" / f"{sample}.downsample_metrics.tsv")
                      if float(row["downsample_ratio"]) == selected["ratio"]]
        if len(candidates) != 1:
            raise ValueError("expected one raw result for the selected ratio")
        actual = candidates[0]
        if any(int(actual[key]) != value for key, value in metrics.items()):
            raise ValueError("selection metrics disagree with the selected raw result")
        fraction = round(100 * (1 - metrics["umi_detected_once"] / metrics["umi_types"]), 6) if metrics["umi_types"] else 0.0
        if abs(float(selected["non_singleton_molecule_fraction"]) - fraction) > 1e-6:
            raise ValueError("selected non-singleton fraction disagrees with its counts")
        if selected["matrix_summary"]:
            summary = json.loads((root / selected["matrix_summary"]["path"]).read_text())
            if sum(summary["feature_counts"].values()) != selected["metrics"]["umi_types"]:
                raise ValueError("selected matrix and molecule counts disagree")
        return selected
    except (OSError, ValueError, KeyError, TypeError, StopIteration) as exc:
        raise ConfigError(f"Invalid or missing selected artifacts for {sample}: {path}: {exc}. Rerun downsample") from exc


def selected_summary_path(config, sample):
    selected = load_selection(config, sample)
    reference = selected["matrix_summary"]
    if reference is None:
        raise ConfigError(f"Selection for {sample} does not contain a feature matrix")
    return config.output_dir / sample / reference["path"]


def selection_point(config, sample):
    """A report view derived from the sole authoritative selection record."""
    selected = load_selection(config, sample)
    reads = selected["metrics"]["reads_sampled"]
    duplication = round(100 * (reads - selected["metrics"]["umi_types"]) / reads, 6) if reads else 0.0
    return {"sample": sample, "optimal_downsample_ratio": selected["ratio"],
            "selected_ratio": selected["ratio"],
            **selected["metrics"], "sequencing_saturation": duplication, "duplication_ratio": duplication,
            "max_non_singleton_molecule_fraction": selected["non_singleton_molecule_fraction"],
            "non_singleton_molecule_fraction": selected["non_singleton_molecule_fraction"],
            "selection_criterion": selected["criterion"]}
