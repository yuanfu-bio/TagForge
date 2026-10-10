"""Stage identities, dependency validation and recoverable publication."""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, replace
import fcntl
import hashlib
from importlib.metadata import version, PackageNotFoundError
import json
import os
from pathlib import Path
import shutil
import sys
import uuid

from . import __version__
from .config import ConfigError
from .io_utils import SUBDIRS, atomic_text, batch_report_dir


MANIFEST_SCHEMA = 1
# Bump only the affected stage when its counting logic or output contract changes.
STAGE_VERSIONS = {step: 1 for step in (
    "extract", "correct", "dedup", "matrix", "pair-map", "downsample", "report")}
STAGE_VERSIONS.update(extract=2, correct=4, dedup=3, downsample=5, matrix=2, report=4, **{"pair-map": 3})
STAGE_VERSIONS["report"] = 6
STAGE_NUMBERS = {"extract": "01", "correct": "02", "dedup": "03", "matrix": "04",
                 "pair-map": "04", "downsample": "05", "report": "06"}


def _digest(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_info(path: Path, *, content: bool = True):
    """Validate files/Parquet footers and identify every committed fragment.

    Large Parquet payloads use fragment metadata and footer/schema validation;
    FASTQ can explicitly opt into a full content checksum. Small artifacts are
    content hashed, including whitelist/annotation files.
    """
    if path.is_dir():
        parts = sorted(path.glob("*.parquet"))
        if not parts:
            raise ConfigError(f"Incomplete Parquet dataset: {path}")
        records = [{"name": p.name, **artifact_info(p, content=False)} for p in parts]
        schemas = {json.dumps(p["columns"]) for p in records}
        if len(schemas) != 1:
            raise ConfigError(f"Parquet fragments have inconsistent schemas: {path}")
        return {"kind": "dataset", "parts": records, "rows": sum(p["rows"] for p in records)}
    stat = path.stat()
    result = {"kind": "file", "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if path.suffix == ".parquet":
        import pyarrow.parquet as pq
        try:
            parquet = pq.ParquetFile(path)
            result.update(rows=parquet.metadata.num_rows,
                          columns=[(f.name, str(f.type)) for f in parquet.schema_arrow])
            # Lists survive JSON round trips, unlike tuples.
            result["columns"] = [list(column) for column in result["columns"]]
        except (OSError, ValueError) as exc:
            raise ConfigError(f"Invalid Parquet artifact: {path}: {exc}") from exc
    elif content:
        result["sha256"] = _digest(path)
    return result


def _package(name):
    try:
        return version(name)
    except PackageNotFoundError:
        return "unavailable"


def _root(config, sample):
    return config.output_dir / sample


def _dirs(config, sample):
    root = _root(config, sample)
    return {key: root / name for key, name in SUBDIRS.items()}


def molecule_source(config, sample):
    suffix = ".molecule_detail.rmMP.parquet" if config.pi_seq_enabled else ".molecule_detail.parquet"
    return _dirs(config, sample)["detail"] / (sample + suffix)


def dependencies(config, step):
    return {
        "extract": (), "correct": ("extract",), "dedup": ("correct",),
        "matrix": ("dedup",), "pair-map": ("dedup",),
        "downsample": ("dedup",) if config.pb_cb_enabled else ("dedup", "matrix"),
        "report": ("extract", "correct", "dedup", "pair-map" if config.pb_cb_enabled else "matrix")
                  + (("downsample",) if config.downsample_enabled else ()),
    }[step]


def stage_plan(config):
    steps = ["extract", "correct", "dedup"]
    if not config.pb_cb_enabled:
        steps.append("matrix")
    if config.downsample_enabled:
        steps.append("downsample")
    if config.pb_cb_enabled:
        steps.append("pair-map")
    steps.append("report")
    return tuple(steps)


def expected_outputs(config, sample, step):
    d = _dirs(config, sample)
    outputs = {
        "extract": [d["extracted"] / f"{sample}.extracted.parquet", d["extracted"] / f"{sample}.extraction_stats.tsv"],
        "correct": [d["corrected"] / f"{sample}.valid_reads.parquet", d["corrected"] / f"{sample}.correction_qc.json",
                    d["corrected"] / f"{sample}.barcode2_not_in_annotation.tsv.gz"] + (
            [d["corrected"] / f"{sample}.barcode_correction_trace.parquet"] if config.trace_enabled else []),
        "dedup": [d["detail"] / f"{sample}.molecule_detail.parquet"] + (
            [d["detail"] / f"{sample}.molecule_detail.rmMP.parquet", d["detail"] / f"{sample}.pi_seq_qc.tsv"] if config.pi_seq_enabled else []),
        "matrix": [d["matrix"] / f"{sample}.raw_molecule_pairs.parquet", d["matrix"] / f"{sample}.raw_matrix_summary.json"],
        "pair-map": [d["matrix"] / f"{sample}.pb_cb_mapping.tsv.gz"] + (
            [d["matrix"] / f"{sample}.{suffix}" for suffix in (
                "pb_cb_map.tsv.gz", "cb_pb_counts.tsv.gz", "cb_pb_count_distribution.tsv", "cb_observed_correction.tsv.gz")]
            + [d["matrix"] / f"{sample}.pb_cb_summary.json"]
            if config.pb_cb_enabled else []),
        "downsample": [d["downsample"] / f"{sample}.{suffix}" for suffix in (
            "downsample_metrics.tsv", "downsample_curve.tsv", "downsample.html", "selection.json")],
        "report": [d["report"] / f"{sample}.report.xlsx", d["report"] / f"{sample}.report.html"],
    }[step]
    if step == "downsample":
        selection = _json(d["downsample"] / f"{sample}.selection.json")
        if selection:
            outputs += [_root(config, sample) / relative for relative in selection["owned_artifacts"]]
    return outputs


def _settings(config, step):
    settings = _stage_settings(config, step)
    if step in {"extract", "correct", "dedup", "matrix", "downsample"}:
        compression = getattr(config, "parquet_compression", "snappy")
        settings["parquet_compression"] = compression
        if compression == "zstd":
            settings["parquet_compression_level"] = config.compression_level
    return settings


def _stage_settings(config, step):
    names = {role: config.target_name(role) for role in ("barcode1", "barcode2", "umi")}
    if step == "extract":
        segments = []
        for segment in config.segments:
            item = {key: getattr(segment, key) for key in (
                "name", "target", "target_name", "read", "method", "length", "start",
                "left_linker", "right_linker", "linker_max_mismatch", "linker_strategy")}
            item["shift_window"] = (segment.correction.max_shift if segment.target != "umi"
                and segment.correction.enabled and segment.correction.allow_shift else 0)
            segments.append(item)
        return {"segments": segments, "cutadapt": _package("cutadapt")}
    if step == "correct":
        return {"names": names, "segments": [{"name": s.name, "target": s.target,
            "length": s.length, "start": s.start, "correction": asdict(s.correction)}
            for s in config.segments if s.target != "umi"], "sequence_only": config.barcode2_sequence_only,
            "annotation_columns": [config.fb_id_column, config.fb_sequence_column, config.fb_name_column],
            "allow_duplicate_names": config.allow_duplicate_names, "trace": config.trace_enabled}
    if step == "dedup":
        return {"names": names, "method": config.umi_method, "distance": config.umi_max_distance,
            "max_group_umis": config.umi_max_group_umis, "max_component_umis": config.umi_max_component_umis,
            "pi_seq": config.pi_seq_enabled, "dominance": config.pi_seq_dominance_threshold,
            "umi_tools": _package("umi-tools")}
    if step == "matrix":
        return {"names": names, "pi_seq": config.pi_seq_enabled}
    if step == "downsample":
        return {"names": names, "pi_seq": config.pi_seq_enabled, "pb_cb": config.pb_cb_enabled,
            "enabled": config.downsample_enabled,
            "ratios": sorted(set(config.downsample_ratios)), "seed": config.downsample_seed,
            "python": list(sys.version_info[:3])}
    if step == "pair-map":
        raw = config.raw.get("pb_cb", {})
        return {"names": names, "pb_cb": config.pb_cb_enabled, "dominance": config.pb_cb_dominance_threshold,
            "observed_correction": bool(raw.get("observed_correction", True)),
            "min_parent_reads": int(raw.get("min_parent_reads", 5)),
            "parent_child_ratio": float(raw.get("parent_child_ratio", 10))}
    return {"names": names, "pi_seq": config.pi_seq_enabled,
            "pb_cb": config.pb_cb_enabled, "downsample_enabled": config.downsample_enabled,
            "dominance": config.pi_seq_dominance_threshold, "report_version": __version__}


def _inputs(config, sample, step):
    d = _dirs(config, sample)
    if step == "extract":
        reads = config.sample(sample)
        return [reads.r1, reads.r2]
    if step == "correct":
        return [d["extracted"] / f"{sample}.extracted.parquet"] + [s.whitelist for s in config.segments
            if s.target != "umi" and s.whitelist] + ([config.fb_info] if not config.barcode2_sequence_only and config.fb_info else [])
    if step == "dedup":
        return [d["corrected"] / f"{sample}.valid_reads.parquet"]
    if step in {"matrix", "downsample"}:
        paths = [molecule_source(config, sample)]
        if step == "downsample" and not config.pb_cb_enabled:
            paths += expected_outputs(config, sample, "matrix")
        return paths
    if step == "pair-map":
        return [d["detail"] / f"{sample}.molecule_detail.parquet"]
    return [path for parent in dependencies(config, step) for path in expected_outputs(config, sample, parent)] + (
        [config.fb_info] if config.fb_info else []) + [s.whitelist for s in config.segments if s.whitelist]


def _identity(config, sample, step):
    mode = config.raw.get("resume", {}).get("input_fingerprint", "metadata")
    if mode not in {"metadata", "sha256"}:
        raise ConfigError("resume.input_fingerprint must be metadata or sha256")
    inputs = [{"path": str(path.resolve()), "artifact": artifact_info(path,
        content=step != "extract" or mode == "sha256")} for path in _inputs(config, sample, step)]
    return {"stage_version": STAGE_VERSIONS[step], "settings": _settings(config, step), "inputs": inputs}


def _identity_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def resume_fingerprint(config, sample, step):
    return _identity_hash(_identity(config, sample, step))


def _json(path):
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else None
    except (OSError, ValueError):
        return None


def manifest_path(config, sample, step):
    return _dirs(config, sample)["checkpoint"] / f"{STAGE_NUMBERS[step]}_{step}.manifest.json"


def _legacy_manifest_path(config, sample, step):
    return _dirs(config, sample)["checkpoint"] / f"{step}.done"


def _read_manifest(config, sample, step):
    path = manifest_path(config, sample, step)
    # A damaged current manifest must never fall back to an older completion.
    if not path.exists() and not path.is_symlink():
        path = path.with_name(f"{step}.manifest.json")
        if not path.exists() and not path.is_symlink():
            path = _legacy_manifest_path(config, sample, step)
    return _json(path)


def migrate_manifest(config, sample, step):
    """Rename a validated legacy manifest while the caller holds the sample lock."""
    destination = manifest_path(config, sample, step)
    if not destination.exists() and not destination.is_symlink():
        for legacy in (destination.with_name(f"{step}.manifest.json"),
                       _legacy_manifest_path(config, sample, step)):
            if legacy.is_file():
                os.replace(legacy, destination)
                break
    return destination


def stage_request(config, sample, step):
    parents = {}
    for parent in dependencies(config, step):
        manifest = _read_manifest(config, sample, parent)
        if not manifest:
            raise ConfigError(f"Stage {step} requires a current {parent} stage for {sample}; run tagforge {parent} first")
        parents[parent] = manifest["generation"]
    identity = _identity(config, sample, step)
    inputs = []
    for item in identity["inputs"]:
        info = item["artifact"]
        # Upstream manifests retain all fragment details. Record a compact
        # direct-input identity here instead of duplicating that entire ledger.
        audit = {"kind": info["kind"], "identity_sha256": _identity_hash(info)}
        audit.update({key: info[key] for key in ("rows", "size", "mtime_ns", "sha256") if key in info})
        if "parts" in info:
            audit["fragment_count"] = len(info["parts"])
        inputs.append({"path": item["path"], "artifact": audit})
    return {"step": step, "sample": sample, "stage_version": STAGE_VERSIONS[step],
            "fingerprint": _identity_hash(identity), "settings": identity["settings"],
            "inputs": inputs, "parents": parents}


def stage_current(config, sample, step, _memo=None):
    memo = {} if _memo is None else _memo
    if step in memo:
        return memo[step]
    memo[step] = False
    manifest = _read_manifest(config, sample, step)
    if not manifest or manifest.get("schema") != MANIFEST_SCHEMA or not manifest.get("generation"):
        return False
    try:
        if not all(stage_current(config, sample, parent, memo) for parent in dependencies(config, step)):
            return False
        if manifest["request"] != stage_request(config, sample, step):
            return False
        paths = expected_outputs(config, sample, step)
        root = _root(config, sample)
        records = {record["path"]: record["artifact"] for record in manifest["outputs"]}
        if set(records) != {str(path.relative_to(root)) for path in paths}:
            return False
        if not all(artifact_info(path) == records[str(path.relative_to(root))] for path in paths):
            return False
        if step == "downsample":
            from .selected import load_selection
            load_selection(config, sample)
        memo[step] = True
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def require_dependencies(config, sample, step):
    for parent in dependencies(config, step):
        if not stage_current(config, sample, parent):
            raise ConfigError(f"Stage {step} requires a current {parent} stage for {sample}; run tagforge {parent} first")


@contextmanager
def project_lock(config):
    path = config.output_dir / ".summary.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


@contextmanager
def sample_lock(config, sample):
    if (_root(config, sample) / "01_checkpoint").exists():
        raise ConfigError(f"Sample {sample} uses the previous output directory layout. "
                          "Archive its sample directory and rerun to use the current layout; "
                          "do not copy old completion manifests into new directories.")
    path = _dirs(config, sample)["checkpoint"] / ".sample.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ConfigError(f"Sample {sample} is already running in {config.output_dir}") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def invalidate(config, sample, step):
    """Withdraw completion before any canonical output is changed."""
    affected = {step}
    for _ in STAGE_VERSIONS:
        affected.update(child for child in STAGE_VERSIONS if affected.intersection(dependencies(config, child)))
    with project_lock(config):
        for child in affected:
            manifest_path(config, sample, child).unlink(missing_ok=True)
            manifest_path(config, sample, child).with_name(f"{child}.manifest.json").unlink(missing_ok=True)
            _legacy_manifest_path(config, sample, child).unlink(missing_ok=True)
        # Previously exported batch snapshots must not continue to look current.
        if step != "report":
            directory = batch_report_dir(config.output_dir)
            for name in ("00_summary.xlsx", "01_FB_not_in_anno.xlsx",
                         "TagForge_batch_report.xlsx", "TagForge_batch_report.html"):
                path = directory / name
                path.unlink(missing_ok=True)


def prepare_stage(config, sample, step, request, *, overwrite=False):
    work = _dirs(config, sample)["checkpoint"] / ".work" / f"stage-{step}"
    # Preserve interrupted stages from the former output layout, including
    # publication journals and committed extraction/correction fragments.
    legacy = _root(config, sample) / "08_tmp" / f"stage-{step}"
    if legacy.exists() and not work.exists():
        work.parent.mkdir(parents=True, exist_ok=True)
        os.replace(legacy, work)
        if not any(legacy.parent.iterdir()):
            legacy.parent.rmdir()
    previous = _json(work / "request.json")
    if overwrite or previous != request:
        if work.exists():
            shutil.rmtree(work)
        work.mkdir(parents=True)
        with atomic_text(work / "request.json") as handle:
            json.dump(request, handle, sort_keys=True)
        (_dirs(config, sample)["checkpoint"] / f"{sample}.{step}.resume.json").unlink(missing_ok=True)
    staged = replace(config, output_dir=work / "output")
    root = _root(staged, sample)
    for key, directory in _dirs(staged, sample).items():
        directory.parent.mkdir(parents=True, exist_ok=True)
        if key in {"logs", "checkpoint"}:
            target = _dirs(config, sample)[key]
            target.mkdir(parents=True, exist_ok=True)
            if not directory.exists():
                directory.symlink_to(target.resolve(), target_is_directory=True)
        else:
            directory.mkdir(parents=True, exist_ok=True)
    canonical = _root(config, sample)
    for source in _inputs(config, sample, step):
        try:
            relative = source.relative_to(canonical)
        except ValueError:
            continue  # Raw FASTQ, whitelist and annotation paths are unchanged.
        link = root / relative
        link.parent.mkdir(parents=True, exist_ok=True)
        if not link.exists() and not link.is_symlink():
            link.symlink_to(source.resolve(), target_is_directory=source.is_dir())
    if step == "extract" and config.extraction_preview_reads > 0:
        # The explicitly labelled preview remains live, outside completed data.
        relative = Path("01_extracted") / f"{sample}.extracted.preview.tsv"
        preview = root / relative
        if not preview.exists() and not preview.is_symlink():
            preview.symlink_to((canonical / relative).resolve())
    return staged, work


def prepare_publication(config, staged, sample, step, request, work):
    root = _root(staged, sample)
    paths = expected_outputs(staged, sample, step)
    if step == "downsample":
        from .selected import load_selection
        load_selection(staged, sample)
    manifest = {"schema": MANIFEST_SCHEMA, "step": step, "sample": sample,
        "tagforge_version": __version__, "generation": uuid.uuid4().hex, "request": request,
        "dependencies": {name: _package(name) for name in ("cutadapt", "umi-tools", "pyarrow", "duckdb")},
        "outputs": [{"path": str(path.relative_to(root)), "artifact": artifact_info(path)} for path in paths]}
    with atomic_text(work / "publication.json") as handle:
        json.dump(manifest, handle, sort_keys=True)
    return manifest


def pending_publication(work):
    return _json(work / "publication.json")


def _publish_artifact(source, destination, backup, expected):
    if source.exists():
        if artifact_info(source) != expected:
            raise ConfigError(f"Staged artifact changed before publication: {source}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            backup.parent.mkdir(parents=True, exist_ok=True)
            if backup.exists():
                if backup.is_dir(): shutil.rmtree(backup)
                else: backup.unlink()
            os.replace(destination, backup)
        os.replace(source, destination)
    elif artifact_info(destination) != expected:
        raise ConfigError(f"Interrupted publication cannot recover artifact: {destination}")


def _obsolete_outputs(config, sample, step):
    """Known previous-generation files that must not shadow current results."""
    if step == "extract":
        return [f"01_extracted/{sample}.extracted.tsv.gz"] + (
            [f"01_extracted/{sample}.extracted.preview.tsv"] if not config.extraction_preview_reads else [])
    if step == "correct":
        return [f"02_corrected/{sample}.barcode_correction_stats.tsv",
                f"02_corrected/{sample}.barcode2_annotation_qc.tsv",
                f"02_corrected/{sample}.valid_reads.tsv.gz",
                f"02_corrected/{sample}.barcode_correction_trace.tsv.gz"] + (
                    [f"02_corrected/{sample}.barcode_correction_trace.parquet"] if not config.trace_enabled else [])
    if step == "matrix":
        return [f"04_matrix/{sample}.raw_count_matrix.tsv.gz"]
    if step == "downsample":
        return [f"05_downsample/{sample}.{suffix}" for suffix in (
            "selection_point.tsv", "optimal_saturation_point.tsv",
            "downsample_summary.tsv", "downsample_expectations.tsv")]
    return []


def publish_stage(config, staged, sample, step, request, work, manifest):
    require_dependencies(config, sample, step)
    if stage_request(config, sample, step) != request:
        raise ConfigError(f"Inputs changed while {step} was running for {sample}; rerun the stage")
    with project_lock(config):
        canonical, source_root = _root(config, sample), _root(staged, sample)
        for record in manifest["outputs"]:
            relative = Path(record["path"])
            _publish_artifact(source_root / relative, canonical / relative, work / "backup" / relative, record["artifact"])
        for relative in _obsolete_outputs(config, sample, step):
            obsolete = canonical / relative
            if obsolete.exists():
                backup = work / "backup" / relative
                backup.parent.mkdir(parents=True, exist_ok=True)
                os.replace(obsolete, backup)
        if step == "downsample":
            from .selected import load_selection
            selection = load_selection(config, sample)
            retained = set(selection["owned_artifacts"])
            for relative in (f"03_molecules/{sample}.optimal_saturation_molecule_detail.parquet",
                             f"04_matrix/{sample}.optimal_molecule_pairs.parquet",
                             f"04_matrix/{sample}.optimal_matrix_summary.json"):
                obsolete = canonical / relative
                if relative not in retained and obsolete.exists():
                    if obsolete.is_dir(): shutil.rmtree(obsolete)
                    else: obsolete.unlink()
        checkpoint = manifest_path(config, sample, step)
        with atomic_text(checkpoint) as handle:
            json.dump(manifest, handle, sort_keys=True)
        _legacy_manifest_path(config, sample, step).unlink(missing_ok=True)
        checkpoint.with_name(f"{step}.manifest.json").unlink(missing_ok=True)
    (_dirs(config, sample)["checkpoint"] / f"{sample}.{step}.resume.json").unlink(missing_ok=True)
    shutil.rmtree(work)
    if not any(work.parent.iterdir()):
        work.parent.rmdir()
