from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from . import simple_yaml


class ConfigError(ValueError):
    pass


def default_downsample_ratios() -> List[float]:
    """Return the logarithmic four-band downsampling grid."""
    return [round(base * multiplier, 4) for base in (0.0001, 0.001, 0.01, 0.1) for multiplier in range(1, 10)]


@dataclass(frozen=True)
class CorrectionConfig:
    enabled: bool = True
    max_shift: int = 1
    max_mismatch: int = 1
    allow_shift: bool = True
    allow_mismatch: bool = True


@dataclass(frozen=True)
class SegmentConfig:
    name: str
    target: str
    read: str
    method: str
    length: int
    start: Optional[int] = None
    left_linker: Optional[str] = None
    right_linker: Optional[str] = None
    linker_max_mismatch: int = 0
    whitelist: Optional[Path] = None
    correction: CorrectionConfig = field(default_factory=CorrectionConfig)


@dataclass(frozen=True)
class SampleConfig:
    sample: str
    r1: Path
    r2: Path


@dataclass
class TagForgeConfig:
    path: Path
    workdir: Path
    output_dir: Path
    samples: List[SampleConfig]
    segments: List[SegmentConfig]
    fb_info: Path
    fb_id_column: str = "FB_ID"
    fb_sequence_column: str = "sequence"
    fb_name_column: str = "antibody_name"
    allow_duplicate_names: bool = False
    umi_method: str = "directional"
    umi_max_distance: int = 1
    downsample_enabled: bool = True
    downsample_ratios: List[float] = field(default_factory=default_downsample_ratios)
    downsample_seed: int = 12345
    downsample_repeats: int = 1
    quick_test_enabled: bool = True
    quick_test_reads: int = 10000
    threads: int = 1
    barcode_workers: Optional[int] = None
    umi_workers: Optional[int] = None
    umi_batch_size: int = 5000
    umi_sqlite_cache_mb: int = 64
    chunk_size: int = 10000
    extraction_preview_reads: int = 1000
    compression_level: int = 3
    overwrite: bool = False
    trace_enabled: bool = True
    raw: Dict[str, Any] = field(default_factory=dict)

    def sample(self, name: str) -> SampleConfig:
        for sample in self.samples:
            if sample.sample == name:
                return sample
        raise ConfigError(f"Sample not found in configuration: {name}")


def _read_yaml(path: Path) -> Dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    try:
        import yaml  # type: ignore
    except ImportError:
        data = simple_yaml.loads(text)
    else:
        data = yaml.safe_load(text)
    if not isinstance(data, dict):
        raise ConfigError("Configuration root must be a YAML mapping")
    return data


def _resolve(base: Path, value: Any) -> Path:
    if value is None:
        raise ConfigError("A required path is missing")
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else (base / path).resolve()


def load_config(config_path: str | Path, check_files: bool = True) -> TagForgeConfig:
    path = Path(config_path).expanduser().resolve()
    if not path.is_file():
        raise ConfigError(f"Configuration file does not exist: {path}")
    raw = _read_yaml(path)
    project = raw.get("project") or {}
    workdir = _resolve(path.parent, project.get("workdir", path.parent))
    output_dir = _resolve(workdir, project.get("output_dir", "02_output"))

    samples_raw = raw.get("samples")
    if not isinstance(samples_raw, list) or not samples_raw:
        raise ConfigError("'samples' must be a non-empty list")
    samples: List[SampleConfig] = []
    seen_samples = set()
    for item in samples_raw:
        if not isinstance(item, dict) or not item.get("sample"):
            raise ConfigError("Each sample requires sample, r1, and r2 fields")
        name = str(item["sample"])
        if name in seen_samples:
            raise ConfigError(f"Duplicate sample name: {name}")
        seen_samples.add(name)
        samples.append(SampleConfig(name, _resolve(workdir, item.get("r1")), _resolve(workdir, item.get("r2"))))

    segments_raw = raw.get("segments")
    if not isinstance(segments_raw, list) or not segments_raw:
        raise ConfigError("'segments' must be a non-empty list")
    segments: List[SegmentConfig] = []
    seen_segments = set()
    for item in segments_raw:
        if not isinstance(item, dict):
            raise ConfigError("Each segment must be a mapping")
        name = str(item.get("name", ""))
        target = str(item.get("target", "")).lower()
        read = str(item.get("read", "")).upper()
        method_value = item.get("methods", item.get("method"))
        if not name or name in seen_segments:
            raise ConfigError(f"Segment names must be present and unique: {name!r}")
        if target not in {"barcode1", "barcode2", "umi"}:
            raise ConfigError(f"Segment {name}: target must be barcode1, barcode2, or umi")
        if read not in {"R1", "R2"}:
            raise ConfigError(f"Segment {name}: read must be R1 or R2")
        try:
            length = int(item["length"])
        except (KeyError, TypeError, ValueError):
            raise ConfigError(f"Segment {name}: length must be a positive integer") from None
        if length <= 0:
            raise ConfigError(f"Segment {name}: length must be positive")
        start = item.get("start")
        left = item.get("left_linker") or None
        right = item.get("right_linker") or None
        if isinstance(method_value, list):
            requested_methods = {str(value).strip().lower() for value in method_value}
        elif method_value is None:
            requested_methods = set()
        else:
            normalized = str(method_value).lower().replace("linker_fixed", "linker+fixed")
            for separator in (",", "|", "->"):
                normalized = normalized.replace(separator, "+")
            requested_methods = {value.strip() for value in normalized.split("+") if value.strip()}
        unknown_methods = requested_methods - {"fixed", "linker"}
        if unknown_methods:
            raise ConfigError(f"Segment {name}: unsupported extraction method(s): {', '.join(sorted(unknown_methods))}")
        # Presence of coordinates/linkers is itself an explicit method declaration.
        # This allows either ``methods: [linker, fixed]`` or the compact form of
        # specifying linker fields and ``start`` together.
        use_fixed = "fixed" in requested_methods or start is not None
        use_linker = "linker" in requested_methods or bool(left or right)
        if not use_fixed and not use_linker:
            raise ConfigError(f"Segment {name}: specify fixed coordinates, linker sequence(s), or both")
        if use_fixed and (start is None or int(start) < 0):
            raise ConfigError(f"Segment {name}: fixed extraction requires a 0-based start >= 0")
        if use_linker and not left and not right:
            raise ConfigError(f"Segment {name}: linker extraction requires left_linker and/or right_linker")
        method = "linker_fixed" if use_linker and use_fixed else "linker" if use_linker else "fixed"
        corr_raw = item.get("correction") or {}
        corr = CorrectionConfig(
            bool(corr_raw.get("enabled", True)), int(corr_raw.get("max_shift", 1)),
            int(corr_raw.get("max_mismatch", 1)), bool(corr_raw.get("allow_shift", True)),
            bool(corr_raw.get("allow_mismatch", True)),
        )
        whitelist = item.get("whitelist")
        segments.append(SegmentConfig(
            name=name, target=target, read=read, method=method, length=length,
            start=int(start) if start is not None else None,
            left_linker=str(left).upper() if left else None,
            right_linker=str(right).upper() if right else None,
            linker_max_mismatch=int(item.get("linker_max_mismatch", 0)),
            whitelist=_resolve(workdir, whitelist) if whitelist else None, correction=corr,
        ))
        seen_segments.add(name)

    targets = {s.target for s in segments}
    missing_targets = {"barcode1", "barcode2", "umi"} - targets
    if missing_targets:
        raise ConfigError(f"Missing segment target(s): {', '.join(sorted(missing_targets))}")
    ann = raw.get("barcode2_annotation") or {}
    fb_info = _resolve(workdir, ann.get("fb_info"))
    umi = raw.get("umi") or {}
    ds = raw.get("downsample") or {}
    perf = raw.get("performance") or {}
    resume = raw.get("resume") or {}
    output = raw.get("output") or {}
    quick = raw.get("quick_test") or {}
    ratios_value = ds.get("ratios", "auto")
    if ratios_value is None or (isinstance(ratios_value, str) and ratios_value.lower() in {"auto", "log_grid"}):
        ratios = default_downsample_ratios()
    elif isinstance(ratios_value, list):
        ratios = [float(x) for x in ratios_value]
    else:
        raise ConfigError("downsample.ratios must be 'auto' or a list of numbers")
    if not ratios or any(x <= 0 or x > 1 for x in ratios):
        raise ConfigError("downsample.ratios must contain values in (0, 1]")
    method = str(umi.get("correction_method", "directional"))
    if method not in {"unique", "cluster", "adjacency", "directional"}:
        raise ConfigError("umi.correction_method must be unique, cluster, adjacency, or directional")
    threads = int(perf.get("threads", 1))
    if threads < 1:
        raise ConfigError("performance.threads must be >= 1")
    barcode_workers_value = perf.get("barcode_workers")
    barcode_workers = int(barcode_workers_value) if barcode_workers_value is not None else None
    umi_workers_value = perf.get("umi_workers")
    umi_workers = int(umi_workers_value) if umi_workers_value is not None else None
    umi_batch_size = int(perf.get("umi_batch_size", 5000))
    umi_sqlite_cache_mb = int(perf.get("umi_sqlite_cache_mb", 64))
    if barcode_workers is not None and barcode_workers < 1:
        raise ConfigError("performance.barcode_workers must be >= 1")
    if umi_workers is not None and umi_workers < 1:
        raise ConfigError("performance.umi_workers must be >= 1")
    if umi_batch_size < 1:
        raise ConfigError("performance.umi_batch_size must be >= 1")
    if umi_sqlite_cache_mb < 1:
        raise ConfigError("performance.umi_sqlite_cache_mb must be >= 1")
    chunk_size = int(perf.get("chunk_size", 10000))
    preview_reads = int(perf.get("extraction_preview_reads", 1000))
    if chunk_size < 1:
        raise ConfigError("performance.chunk_size must be >= 1")
    if preview_reads < 0:
        raise ConfigError("performance.extraction_preview_reads must be >= 0")
    quick_test_reads = int(quick.get("reads", 10000))
    if quick_test_reads < 1:
        raise ConfigError("quick_test.reads must be >= 1")
    config = TagForgeConfig(
        path=path, workdir=workdir, output_dir=output_dir, samples=samples, segments=segments,
        fb_info=fb_info, fb_id_column=str(ann.get("id_column", "FB_ID")),
        fb_sequence_column=str(ann.get("sequence_column", "sequence")),
        fb_name_column=str(ann.get("name_column", "antibody_name")),
        allow_duplicate_names=bool(ann.get("allow_duplicate_names", False)),
        umi_method=method, umi_max_distance=int(umi.get("max_distance", 1)),
        downsample_enabled=bool(ds.get("enabled", True)), downsample_ratios=sorted(set(ratios)),
        downsample_seed=int(ds.get("random_seed", 12345)), downsample_repeats=int(ds.get("repeats", 1)),
        quick_test_enabled=bool(quick.get("enabled", True)), quick_test_reads=quick_test_reads,
        threads=threads, barcode_workers=barcode_workers,
        umi_workers=umi_workers, umi_batch_size=umi_batch_size,
        umi_sqlite_cache_mb=umi_sqlite_cache_mb,
        chunk_size=chunk_size, extraction_preview_reads=preview_reads,
        compression_level=int(perf.get("compression_level", 3)),
        overwrite=bool(resume.get("overwrite", False)), trace_enabled=bool(output.get("correction_trace", True)), raw=raw,
    )
    if check_files:
        missing = [p for s in samples for p in (s.r1, s.r2) if not p.is_file()]
        missing += [s.whitelist for s in segments if s.whitelist and not s.whitelist.is_file()]
        if not fb_info.is_file():
            missing.append(fb_info)
        if missing:
            raise ConfigError("Missing input file(s):\n  " + "\n  ".join(str(p) for p in missing))
    return config
