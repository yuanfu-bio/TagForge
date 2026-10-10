from __future__ import annotations

from dataclasses import dataclass, field
import math
from pathlib import Path
from typing import Any, Dict, List, Optional

from . import simple_yaml


class ConfigError(ValueError):
    pass


def validate_parquet_compression(value):
    if not isinstance(value, str) or value.strip().lower() not in {"snappy", "zstd"}:
        raise ConfigError("performance.parquet_compression must be snappy or zstd")
    return value.strip().lower()


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

    def __post_init__(self):
        if self.max_shift < 0 or self.max_mismatch < 0:
            raise ConfigError("Correction max_shift and max_mismatch must be nonnegative")


@dataclass(frozen=True)
class SegmentConfig:
    name: str
    target: str
    read: str
    method: str
    length: int
    start: Optional[int] = None
    target_name: Optional[str] = None
    left_linker: Optional[str] = None
    right_linker: Optional[str] = None
    linker_max_mismatch: int = 0
    whitelist: Optional[Path] = None
    correction: CorrectionConfig = field(default_factory=CorrectionConfig)
    linker_strategy: str = "legacy_shortest"


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
    fb_info: Optional[Path]
    barcode2_sequence_only: bool = False
    fb_id_column: str = "FB_ID"
    fb_sequence_column: str = "sequence"
    fb_name_column: str = "antibody_name"
    allow_duplicate_names: bool = False
    umi_method: str = "directional"
    umi_max_distance: int = 1
    downsample_enabled: bool = True
    downsample_ratios: List[float] = field(default_factory=default_downsample_ratios)
    downsample_seed: int = 12345
    quick_test_enabled: bool = True
    quick_test_reads: int = 10000
    quick_test_min_valid_rate: Optional[float] = None
    quick_test_max_n_base_rate: Optional[float] = None
    threads: int = 1
    barcode_workers: Optional[int] = None
    downsample_workers: Optional[int] = None
    umi_aggregation_workers: Optional[int] = None
    umi_aggregation_backend: str = "external_sort"
    umi_sort_memory_mb: int = 512
    umi_workers: Optional[int] = None
    umi_batch_size: int = 5000
    umi_max_group_umis: Optional[int] = None
    umi_max_component_umis: Optional[int] = None
    umi_sqlite_cache_mb: int = 64
    matrix_memory_mb: int = 512
    pi_seq_enabled: bool = False
    pi_seq_dominance_threshold: float = 0.8
    pb_cb_enabled: bool = False
    pb_cb_dominance_threshold: float = 0.8
    chunk_size: int = 10000
    extraction_preview_reads: int = 0
    compression_level: int = 3
    parquet_compression: str = "snappy"
    overwrite: bool = False
    resume_enabled: bool = True
    quick_test_read_details: bool = False
    trace_enabled: bool = True
    raw: Dict[str, Any] = field(default_factory=dict)

    def sample(self, name: str) -> SampleConfig:
        for sample in self.samples:
            if sample.sample == name:
                return sample
        raise ConfigError(f"Sample not found in configuration: {name}")

    def target_name(self, role: str) -> str:
        for segment in self.segments:
            if segment.target == role:
                return segment.target_name or role
        return role

    def segment_column(self, role: str) -> str:
        return f"{self.target_name(role)}_segments"


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


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "yes", "on", "1"}:
            return True
        if normalized in {"false", "no", "off", "0"}:
            return False
    if isinstance(value, int) and value in (0, 1):
        return bool(value)
    raise ConfigError(f"Expected a boolean, got {value!r}")


def _known_keys(value, allowed, context):
    if not isinstance(value, dict):
        raise ConfigError(f"{context} must be a mapping")
    unknown = set(value) - set(allowed)
    if unknown:
        raise ConfigError(f"Unknown configuration key(s) in {context}: {', '.join(sorted(map(str, unknown)))}")


_CORRECTION_KEYS = {"enabled", "max_shift", "max_mismatch", "allow_shift", "allow_mismatch"}
_SEGMENT_KEYS = {"segment", "name", "target", "read", "method", "methods", "length", "start",
                 "left_linker", "right_linker", "linker_max_mismatch", "whitelist", "correction", "linker_strategy"}
_ANNOTATION_KEYS = {"fb_info", "id_column", "sequence_column", "name_column", "allow_duplicate_names"}


def _validate_config_keys(raw):
    sections = {
        "project": {"name", "workdir", "output_dir"},
        "correction_barcode": _CORRECTION_KEYS,
        "correction_umi": {"method", "correction_method", "max_distance"},
        "linker": {"max_mismatch", "strategy"},
        "downsample": {"enabled", "ratios", "random_seed"},
        "quick_test": {"enabled", "reads", "min_valid_rate", "max_n_base_rate"},
        "performance": {"threads", "barcode_workers", "downsample_workers", "umi_aggregation_workers",
                        "umi_aggregation_backend", "umi_sort_memory_mb", "umi_workers", "umi_batch_size",
                        "umi_max_group_umis", "umi_max_component_umis",
                        "umi_sqlite_cache_mb", "chunk_size", "extraction_preview_reads", "compression_level", "parquet_compression", "matrix_memory_mb"},
        "resume": {"enabled", "overwrite", "input_fingerprint"},
        "output": {"correction_trace", "quick_test_read_details"}, "library": {"type"},
        "pi_seq": {"dominance_threshold"},
        "pb_cb": {"dominance_threshold", "observed_correction", "min_parent_reads", "parent_child_ratio"},
        "barcode2_annotation": _ANNOTATION_KEYS,
    }
    _known_keys(raw, set(sections) | {"samples", "barcode", "barcode1", "barcode2", "umi",
                                   "segments", "library_type", "linker_max_mismatch"}, "root")
    for name, allowed in sections.items():
        if name in raw and raw[name] is not None:
            if name == "library" and isinstance(raw[name], str):
                continue
            _known_keys(raw[name], allowed, name)
    samples = raw.get("samples")
    if isinstance(samples, dict):
        _known_keys(samples, {"auto"}, "samples")
        _known_keys(samples.get("auto"), {"raw_dir", "r1", "r2"}, "samples.auto")
    elif isinstance(samples, list):
        for sample in samples:
            _known_keys(sample, {"sample", "r1", "r2"}, "samples[]")
    for role in ("barcode1", "barcode2", "umi"):
        section = raw.get(role)
        if isinstance(section, dict):
            allowed = {"name", "target", "segments"}
            if role == "barcode2":
                allowed |= {"sequence_only", "annotation", "barcode2_annotation"}
            if role == "umi":
                allowed |= {"method", "correction_method", "max_distance"}
            _known_keys(section, allowed, role)
            for key in ("annotation", "barcode2_annotation"):
                if key in section:
                    _known_keys(section[key], _ANNOTATION_KEYS, f"{role}.{key}")
    for group in raw.get("barcode", []) or []:
        _known_keys(group, {"barcode1", "barcode2", "role", "name", "target", "segments", "annotation", "barcode2_annotation"}, "barcode[]")
        for key in ("annotation", "barcode2_annotation"):
            if key in group:
                _known_keys(group[key], _ANNOTATION_KEYS, f"barcode[].{key}")


def _integer(value, key):
    try:
        number = int(value)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ConfigError(f"{key} must be an integer") from exc
    if isinstance(value, bool) or not isinstance(value, str) and value != number:
        raise ConfigError(f"{key} must be an integer")
    return number


def _correction_from_mapping(raw: Dict[str, Any]) -> CorrectionConfig:
    return CorrectionConfig(
        enabled=_as_bool(raw.get("enabled", True)),
        max_shift=_integer(raw.get("max_shift", 1), "correction.max_shift"),
        max_mismatch=_integer(raw.get("max_mismatch", 1), "correction.max_mismatch"),
        allow_shift=_as_bool(raw.get("allow_shift", True)),
        allow_mismatch=_as_bool(raw.get("allow_mismatch", True)),
    )


def _merge_correction(default: CorrectionConfig, override: Dict[str, Any] | None) -> CorrectionConfig:
    if not override:
        return default
    return CorrectionConfig(
        enabled=_as_bool(override.get("enabled", default.enabled)),
        max_shift=_integer(override.get("max_shift", default.max_shift), "correction.max_shift"),
        max_mismatch=_integer(override.get("max_mismatch", default.max_mismatch), "correction.max_mismatch"),
        allow_shift=_as_bool(override.get("allow_shift", default.allow_shift)),
        allow_mismatch=_as_bool(override.get("allow_mismatch", default.allow_mismatch)),
    )


def _requested_methods(method_value: Any, segment_name: str) -> set[str]:
    if isinstance(method_value, list):
        requested = {str(value).strip().lower() for value in method_value}
    elif method_value is None:
        requested = set()
    else:
        normalized = str(method_value).lower().replace("linker_fixed", "linker+fixed")
        for separator in (",", "|", "->"):
            normalized = normalized.replace(separator, "+")
        requested = {value.strip() for value in normalized.split("+") if value.strip()}
    unknown = requested - {"fixed", "linker"}
    if unknown:
        raise ConfigError(
            f"Segment {segment_name}: unsupported extraction method(s): "
            + ", ".join(sorted(unknown))
        )
    return requested


def _parse_segment(
    item: Dict[str, Any],
    *,
    target: str,
    target_name: str,
    workdir: Path,
    default_correction: CorrectionConfig,
    default_linker_max_mismatch: int,
    seen_segments: set[str],
    default_linker_strategy: str = "legacy_shortest",
) -> SegmentConfig:
    _known_keys(item, _SEGMENT_KEYS, "segment")
    if "correction" in item:
        _known_keys(item["correction"], _CORRECTION_KEYS, "segment.correction")
    name = str(item.get("segment", item.get("name", "")))
    if "target" in item and str(item["target"]).lower() != target:
        raise ConfigError(f"Segment {name}: target conflicts with its parent role {target}")
    read = str(item.get("read", "")).upper()
    method_value = item.get("methods", item.get("method"))
    if not name or name in seen_segments:
        raise ConfigError(f"Segment names must be present and unique: {name!r}")
    if target not in {"barcode1", "barcode2", "umi"}:
        raise ConfigError(
            f"Segment {name}: internal target role must be barcode1, barcode2, or umi"
        )
    if read not in {"R1", "R2"}:
        raise ConfigError(f"Segment {name}: read must be R1 or R2")
    try:
        length = _integer(item["length"], f"Segment {name}: length")
    except (KeyError, TypeError, ValueError):
        raise ConfigError(f"Segment {name}: length must be a positive integer") from None
    if length <= 0:
        raise ConfigError(f"Segment {name}: length must be positive")
    start = item.get("start")
    if start is not None:
        start = _integer(start, f"Segment {name}: start")
    left = item.get("left_linker") or None
    right = item.get("right_linker") or None
    requested = _requested_methods(method_value, name)
    # Presence of coordinates/linkers is itself an explicit method declaration.
    # This allows either ``methods: [linker, fixed]`` or the compact form of
    # specifying linker fields and ``start`` together.
    use_fixed = "fixed" in requested or start is not None
    use_linker = "linker" in requested or bool(left or right)
    if not use_fixed and not use_linker:
        raise ConfigError(f"Segment {name}: specify fixed coordinates, linker sequence(s), or both")
    if use_fixed and (start is None or start < 0):
        raise ConfigError(f"Segment {name}: fixed extraction requires a 0-based start >= 0")
    if use_linker and not left and not right:
        raise ConfigError(f"Segment {name}: linker extraction requires left_linker and/or right_linker")
    method = "linker_fixed" if use_linker and use_fixed else "linker" if use_linker else "fixed"
    correction_default = (
        CorrectionConfig(False, 0, 0, False, False)
        if target == "umi" else default_correction
    )
    correction = _merge_correction(correction_default, item.get("correction") or {})
    if target == "umi" and (correction.enabled or item.get("whitelist")):
        raise ConfigError(
            f"Segment {name}: UMI segments do not use barcode correction or a whitelist; "
            "configure UMI clustering under correction_umi instead"
        )
    linker_strategy = str(item.get("linker_strategy", default_linker_strategy))
    if linker_strategy not in {"legacy_shortest", "length_position"}:
        raise ConfigError("linker strategy must be legacy_shortest or length_position")
    linker_max_mismatch = _integer(item.get("linker_max_mismatch", default_linker_max_mismatch), f"Segment {name}: linker_max_mismatch")
    if linker_max_mismatch < 0:
        raise ConfigError(f"Segment {name}: linker_max_mismatch must be >= 0")
    whitelist = item.get("whitelist")
    seen_segments.add(name)
    return SegmentConfig(
        name=name, target=target, read=read, method=method, length=length,
        start=start,
        target_name=target_name,
        left_linker=str(left).upper() if left else None,
        right_linker=str(right).upper() if right else None,
        linker_max_mismatch=linker_max_mismatch,
        whitelist=_resolve(workdir, whitelist) if whitelist else None,
        correction=correction,
        linker_strategy=linker_strategy,
    )


def _barcode_group_role_and_name(group: Dict[str, Any]) -> tuple[str, str]:
    role_keys = [key for key in ("barcode1", "barcode2") if key in group]
    if len(role_keys) == 1:
        role = role_keys[0]
        name = str(group[role] or "").strip()
    elif "role" in group:
        role = str(group["role"]).strip().lower()
        name = str(group.get("name", group.get("target", role))).strip()
    else:
        raise ConfigError(
            "Each barcode group must specify 'barcode1: NAME', 'barcode2: NAME', "
            "or both 'role' and 'name'."
        )
    if role not in {"barcode1", "barcode2"}:
        raise ConfigError("Barcode group role must be barcode1 or barcode2")
    if not name:
        raise ConfigError(f"Barcode group {role} requires a non-empty name")
    return role, name


def _section_name_and_segments(role: str, section: Any) -> tuple[str, list[Any]]:
    if isinstance(section, dict):
        name = str(section.get("name", section.get("target", role))).strip()
        segments = section.get("segments")
    elif isinstance(section, list):
        name = role
        segments = section
    else:
        raise ConfigError(
            f"'{role}' must be either a list of segments or a mapping with "
            "'name' and 'segments'."
        )
    if not name:
        raise ConfigError(f"'{role}.name' must not be empty")
    if not isinstance(segments, list) or not segments:
        raise ConfigError(f"'{role}' requires a non-empty 'segments' list")
    return name, segments


def _format_sample_template(template: str, sample: str, key: str) -> str:
    try:
        return template.format(sample=sample)
    except KeyError as exc:
        raise ConfigError(
            f"samples.auto.{key} contains unsupported template field {{{exc.args[0]}}}; "
            "only {sample} is supported"
        ) from None


def _parse_samples(samples_raw: Any, workdir: Path) -> List[SampleConfig]:
    samples: List[SampleConfig] = []
    seen_samples = set()
    if isinstance(samples_raw, list):
        if not samples_raw:
            raise ConfigError("'samples' must be a non-empty list or an auto-discovery mapping")
        for item in samples_raw:
            if not isinstance(item, dict) or not item.get("sample"):
                raise ConfigError("Each sample requires sample, r1, and r2 fields")
            name = str(item["sample"])
            if name in seen_samples:
                raise ConfigError(f"Duplicate sample name: {name}")
            seen_samples.add(name)
            samples.append(SampleConfig(
                name, _resolve(workdir, item.get("r1")), _resolve(workdir, item.get("r2"))
            ))
        return samples
    if isinstance(samples_raw, dict) and "auto" in samples_raw:
        auto = samples_raw.get("auto")
        if not isinstance(auto, dict):
            raise ConfigError("'samples.auto' must be a mapping")
        raw_dir = _resolve(workdir, auto.get("raw_dir", "01_raw"))
        if not raw_dir.is_dir():
            raise ConfigError(f"samples.auto.raw_dir does not exist or is not a directory: {raw_dir}")
        r1_template = str(auto.get("r1", "{sample}_raw_1.fq.gz"))
        r2_template = str(auto.get("r2", "{sample}_raw_2.fq.gz"))
        for sample_dir in sorted((p for p in raw_dir.iterdir() if p.is_dir()), key=lambda p: p.name):
            name = sample_dir.name
            if name in seen_samples:
                raise ConfigError(f"Duplicate sample name: {name}")
            seen_samples.add(name)
            samples.append(SampleConfig(
                name,
                sample_dir / _format_sample_template(r1_template, name, "r1"),
                sample_dir / _format_sample_template(r2_template, name, "r2"),
            ))
        if not samples:
            raise ConfigError(f"samples.auto.raw_dir contains no sample directories: {raw_dir}")
        return samples
    raise ConfigError("'samples' must be a non-empty list or contain an 'auto' mapping")


def load_config(config_path: str | Path, check_files: bool = True) -> TagForgeConfig:
    path = Path(config_path).expanduser().resolve()
    if not path.is_file():
        raise ConfigError(f"Configuration file does not exist: {path}")
    raw = _read_yaml(path)
    _validate_config_keys(raw)
    project = raw.get("project") or {}
    workdir = _resolve(path.parent, project.get("workdir", path.parent))
    output_dir = _resolve(workdir, project.get("output_dir", "02_output"))

    samples = _parse_samples(raw.get("samples"), workdir)
    for sample in samples:
        if sample.sample in {".", ".."} or any(c in sample.sample for c in "/\\\t\r\n"):
            raise ConfigError("Sample names must be directory names without path separators or control characters")

    barcode_correction_default = _correction_from_mapping(raw.get("correction_barcode") or {})
    linker_raw = raw.get("linker") or {}
    default_linker_max_mismatch = _integer(linker_raw.get("max_mismatch", raw.get("linker_max_mismatch", 0)), "linker.max_mismatch")
    default_linker_strategy = str(linker_raw.get("strategy", "legacy_shortest"))
    if default_linker_strategy not in {"legacy_shortest", "length_position"}:
        raise ConfigError("linker.strategy must be legacy_shortest or length_position")
    segments: List[SegmentConfig] = []
    seen_segments = set()

    barcode_groups_raw = raw.get("barcode")
    direct_barcode_roles = [role for role in ("barcode1", "barcode2") if role in raw]
    umi_section = raw.get("umi")
    has_direct_sections = bool(direct_barcode_roles) or (
        isinstance(umi_section, dict) and "segments" in umi_section
    )
    umi_segments_raw = umi_section if isinstance(umi_section, list) else None
    if has_direct_sections:
        for role in ("barcode1", "barcode2"):
            if role not in raw:
                continue
            target_name, group_segments = _section_name_and_segments(role, raw[role])
            for item in group_segments:
                if not isinstance(item, dict):
                    raise ConfigError(f"Each segment in {role} group {target_name} must be a mapping")
                segments.append(_parse_segment(
                    item, target=role, target_name=target_name, workdir=workdir,
                    default_correction=barcode_correction_default,
                    default_linker_max_mismatch=default_linker_max_mismatch,
                    default_linker_strategy=default_linker_strategy,
                    seen_segments=seen_segments,
                ))
        if "umi" not in raw:
            raise ConfigError("'umi' must be provided when using direct barcode1/barcode2 sections")
        umi_name, umi_group_segments = _section_name_and_segments("umi", raw["umi"])
        for item in umi_group_segments:
            if not isinstance(item, dict):
                raise ConfigError("Each UMI segment must be a mapping")
            segments.append(_parse_segment(
                item, target="umi", target_name=umi_name, workdir=workdir,
                default_correction=barcode_correction_default,
                default_linker_max_mismatch=default_linker_max_mismatch,
                    default_linker_strategy=default_linker_strategy,
                seen_segments=seen_segments,
            ))
    elif barcode_groups_raw is not None or umi_segments_raw is not None:
        if not isinstance(barcode_groups_raw, list) or not barcode_groups_raw:
            raise ConfigError("'barcode' must be a non-empty list of barcode groups")
        for group in barcode_groups_raw:
            if not isinstance(group, dict):
                raise ConfigError("Each barcode group must be a mapping")
            role, target_name = _barcode_group_role_and_name(group)
            group_segments = group.get("segments")
            if not isinstance(group_segments, list) or not group_segments:
                raise ConfigError(f"Barcode group {target_name} requires a non-empty 'segments' list")
            for item in group_segments:
                if not isinstance(item, dict):
                    raise ConfigError(f"Each segment in barcode group {target_name} must be a mapping")
                segments.append(_parse_segment(
                    item, target=role, target_name=target_name, workdir=workdir,
                    default_correction=barcode_correction_default,
                    default_linker_max_mismatch=default_linker_max_mismatch,
                    default_linker_strategy=default_linker_strategy,
                    seen_segments=seen_segments,
                ))
        if not isinstance(umi_segments_raw, list) or not umi_segments_raw:
            raise ConfigError("'umi' must be a non-empty list of UMI segments")
        for item in umi_segments_raw:
            if not isinstance(item, dict):
                raise ConfigError("Each UMI segment must be a mapping")
            segments.append(_parse_segment(
                item, target="umi", target_name="umi", workdir=workdir,
                default_correction=barcode_correction_default,
                default_linker_max_mismatch=default_linker_max_mismatch,
                    default_linker_strategy=default_linker_strategy,
                seen_segments=seen_segments,
            ))
    else:
        segments_raw = raw.get("segments")
        if not isinstance(segments_raw, list) or not segments_raw:
            raise ConfigError(
                "Configuration requires either new 'barcode'/'umi' sections "
                "or a legacy non-empty 'segments' list"
            )
        for item in segments_raw:
            if not isinstance(item, dict):
                raise ConfigError("Each segment must be a mapping")
            target = str(item.get("target", "")).lower()
            if target not in {"barcode1", "barcode2", "umi"}:
                raise ConfigError(
                    f"Segment {item.get('name', item.get('segment', ''))}: target must be "
                    "barcode1, barcode2, or umi in legacy 'segments'. For custom names, "
                    "use the new 'barcode' and 'umi' sections."
                )
            segments.append(_parse_segment(
                item, target=target, target_name=target, workdir=workdir,
                default_correction=barcode_correction_default,
                default_linker_max_mismatch=default_linker_max_mismatch,
                    default_linker_strategy=default_linker_strategy,
                seen_segments=seen_segments,
            ))

    targets = {s.target for s in segments}
    missing_targets = {"barcode1", "barcode2", "umi"} - targets
    if missing_targets:
        raise ConfigError(f"Missing segment target(s): {', '.join(sorted(missing_targets))}")
    role_names = {role: {segment.target_name for segment in segments if segment.target == role}
                  for role in targets}
    if any(len(names) != 1 for names in role_names.values()):
        raise ConfigError("Each target role must have a single output name")
    names = {role: next(iter(values)) for role, values in role_names.items()}
    reserved = {"read_id", "corrected_umi", "reads_count", "raw_umi_count", "molecule_count",
                "methods", "status", "failure_reason", "correction_summary", "shift_distance", "mismatch_distance",
                "reads_count_at_optimal_downsample"}
    generated = [names["barcode1"], names["umi"], names["barcode2"] + "_sequence", names["barcode2"] + "_name"]
    if len(set(names.values())) != 3 or len(set(generated)) != len(generated) or set(names.values()) & reserved or set(generated) & reserved:
        raise ConfigError("Target names conflict with another role or reserved output column")
    for name in [*names.values(), *(segment.name for segment in segments)]:
        if any(character in name for character in ",\t\r\n"):
            raise ConfigError("Target and segment names must not contain commas, tabs or newlines")
    if any(segment.linker_max_mismatch < 0 for segment in segments):
        raise ConfigError("linker_max_mismatch must be nonnegative")
    group_annotation = None
    if has_direct_sections:
        group_annotation = (
            raw.get("barcode2", {}).get("annotation", raw.get("barcode2", {}).get("barcode2_annotation"))
            if isinstance(raw.get("barcode2"), dict) else None
        )
    elif isinstance(barcode_groups_raw, list):
        for group in barcode_groups_raw:
            if isinstance(group, dict):
                try:
                    role, _ = _barcode_group_role_and_name(group)
                except ConfigError:
                    continue
                if role == "barcode2":
                    group_annotation = group.get("annotation", group.get("barcode2_annotation"))
                    if group_annotation:
                        break
    ann = group_annotation or raw.get("barcode2_annotation") or {}
    barcode2_raw = raw.get("barcode2") if isinstance(raw.get("barcode2"), dict) else {}
    barcode2_sequence_only = _as_bool(barcode2_raw.get("sequence_only", False))
    fb_info = None if barcode2_sequence_only else _resolve(workdir, ann.get("fb_info"))
    umi = raw.get("correction_umi") or (raw.get("umi") if isinstance(raw.get("umi"), dict) else {}) or {}
    ds = raw.get("downsample") or {}
    perf = raw.get("performance") or {}
    resume = raw.get("resume") or {}
    output = raw.get("output") or {}
    library = raw.get("library") or {}
    if not isinstance(library, dict):
        library = {"type": library}
    library_type = str(library.get("type", raw.get("library_type", ""))).strip().lower().replace("_", "-")
    if library_type not in {"", "standard", "barcode-feature", "pi-seq", "pb-cb"}:
        raise ConfigError("library.type must be standard, barcode-feature, pi-seq or pb-cb")
    pi_seq = raw.get("pi_seq") or {}
    if not isinstance(pi_seq, dict):
        raise ConfigError("pi_seq must be a mapping")
    pi_seq_enabled = library_type == "pi-seq"
    pi_seq_dominance_threshold = float(pi_seq.get("dominance_threshold", 0.8))
    if not 0 <= pi_seq_dominance_threshold <= 1:
        raise ConfigError("pi_seq.dominance_threshold must be in [0, 1]")
    pb_cb = raw.get("pb_cb") or {}
    if not isinstance(pb_cb, dict):
        raise ConfigError("pb_cb must be a mapping")
    pb_cb_enabled = library_type == "pb-cb"
    pb_cb_dominance_threshold = float(pb_cb.get("dominance_threshold", 0.8))
    if not 0 <= pb_cb_dominance_threshold <= 1:
        raise ConfigError("pb_cb.dominance_threshold must be in [0, 1]")
    if _integer(pb_cb.get("min_parent_reads", 5), "pb_cb.min_parent_reads") < 1:
        raise ConfigError("pb_cb.min_parent_reads must be >= 1")
    parent_child_ratio = float(pb_cb.get("parent_child_ratio", 10))
    if not math.isfinite(parent_child_ratio) or parent_child_ratio < 1:
        raise ConfigError("pb_cb.parent_child_ratio must be >= 1")
    if pb_cb_enabled:
        if not barcode2_sequence_only:
            raise ConfigError("PB-CB libraries require barcode2.sequence_only: true")
        cb_segments = [segment for segment in segments if segment.target == "barcode2"]
        if any(segment.correction.enabled for segment in cb_segments):
            raise ConfigError(
                "PB-CB libraries require barcode2 correction.enabled: false; "
                "CB is retained as the observed sequence rather than searched against its whitelist"
            )
        if any(segment.whitelist is not None for segment in cb_segments):
            raise ConfigError(
                "PB-CB libraries must not set a barcode2 whitelist; remove it to avoid loading the large CB whitelist"
            )
    quick = raw.get("quick_test") or {}
    ratios_value = ds.get("ratios", "auto")
    if ratios_value is None or (isinstance(ratios_value, str) and ratios_value.lower() in {"auto", "log_grid"}):
        ratios = default_downsample_ratios()
    elif isinstance(ratios_value, list):
        ratios = [float(x) for x in ratios_value]
    else:
        raise ConfigError("downsample.ratios must be 'auto' or a list of numbers")
    if not ratios or any(not math.isfinite(x) or x < 0 or x > 1 for x in ratios):
        raise ConfigError("downsample.ratios must contain values in [0, 1]")
    method = str(umi.get("method", umi.get("correction_method", "directional")))
    umi_max_distance = _integer(umi.get("max_distance", 1), "correction_umi.max_distance")
    if umi_max_distance < 0:
        raise ConfigError("correction_umi.max_distance must be nonnegative")
    if method not in {"unique", "cluster", "adjacency", "directional"}:
        raise ConfigError("umi.correction_method must be unique, cluster, adjacency, or directional")
    threads = _integer(perf.get("threads", 1), "performance.threads")
    if threads < 1:
        raise ConfigError("performance.threads must be >= 1")
    barcode_workers_value = perf.get("barcode_workers")
    barcode_workers = _integer(barcode_workers_value, "performance.barcode_workers") if barcode_workers_value is not None else None
    downsample_workers_value = perf.get("downsample_workers")
    downsample_workers = _integer(downsample_workers_value, "performance.downsample_workers") if downsample_workers_value is not None else None
    umi_aggregation_workers_value = perf.get("umi_aggregation_workers")
    umi_aggregation_workers = (
        _integer(umi_aggregation_workers_value, "performance.umi_aggregation_workers")
        if umi_aggregation_workers_value is not None else None
    )
    umi_aggregation_backend = str(perf.get("umi_aggregation_backend", "external_sort")).replace("-", "_")
    umi_sort_memory_mb = _integer(perf.get("umi_sort_memory_mb", 512), "performance.umi_sort_memory_mb")
    umi_workers_value = perf.get("umi_workers")
    umi_workers = _integer(umi_workers_value, "performance.umi_workers") if umi_workers_value is not None else None
    umi_batch_size = _integer(perf.get("umi_batch_size", 5000), "performance.umi_batch_size")
    umi_budgets = {}
    for name in ("umi_max_group_umis", "umi_max_component_umis"):
        value = perf.get(name)
        umi_budgets[name] = None if value is None else _integer(value, f"performance.{name}")
        if umi_budgets[name] is not None and umi_budgets[name] < 1:
            raise ConfigError(f"performance.{name} must be >= 1 or null")
    umi_sqlite_cache_mb = _integer(perf.get("umi_sqlite_cache_mb", 64), "performance.umi_sqlite_cache_mb")
    if barcode_workers is not None and barcode_workers < 1:
        raise ConfigError("performance.barcode_workers must be >= 1")
    if downsample_workers is not None and downsample_workers < 1:
        raise ConfigError("performance.downsample_workers must be >= 1")
    if umi_aggregation_workers is not None and umi_aggregation_workers < 1:
        raise ConfigError("performance.umi_aggregation_workers must be >= 1")
    if umi_aggregation_backend not in {"external_sort", "sqlite"}:
        raise ConfigError("performance.umi_aggregation_backend must be external_sort or sqlite")
    if umi_sort_memory_mb < 1:
        raise ConfigError("performance.umi_sort_memory_mb must be >= 1")
    if umi_workers is not None and umi_workers < 1:
        raise ConfigError("performance.umi_workers must be >= 1")
    if umi_batch_size < 1:
        raise ConfigError("performance.umi_batch_size must be >= 1")
    if umi_sqlite_cache_mb < 1:
        raise ConfigError("performance.umi_sqlite_cache_mb must be >= 1")
    chunk_size = _integer(perf.get("chunk_size", 10000), "performance.chunk_size")
    preview_reads = _integer(perf.get("extraction_preview_reads", 0), "performance.extraction_preview_reads")
    if chunk_size < 1:
        raise ConfigError("performance.chunk_size must be >= 1")
    if preview_reads < 0:
        raise ConfigError("performance.extraction_preview_reads must be >= 0")
    quick_test_reads = _integer(quick.get("reads", 10000), "quick_test.reads")
    if quick_test_reads < 1:
        raise ConfigError("quick_test.reads must be >= 1")
    qc_thresholds = {}
    for key in ("min_valid_rate", "max_n_base_rate"):
        value = quick.get(key)
        if value is not None:
            value = float(value)
            if not 0 <= value <= 1:
                raise ConfigError(f"quick_test.{key} must be in [0, 1]")
        qc_thresholds[key] = value
    if resume.get("input_fingerprint", "metadata") not in {"metadata", "sha256"}:
        raise ConfigError("resume.input_fingerprint must be metadata or sha256")
    compression_level = _integer(perf.get("compression_level", 3), "performance.compression_level")
    if not 0 <= compression_level <= 9:
        raise ConfigError("performance.compression_level must be in [0, 9]")
    parquet_compression = validate_parquet_compression(perf.get("parquet_compression", "snappy"))
    if "observed_correction" in pb_cb:
        pb_cb["observed_correction"] = _as_bool(pb_cb["observed_correction"])
    matrix_memory_mb = _integer(perf.get("matrix_memory_mb", 512), "performance.matrix_memory_mb")
    if matrix_memory_mb < 32:
        raise ConfigError("performance.matrix_memory_mb must be >= 32")
    config = TagForgeConfig(
        path=path, workdir=workdir, output_dir=output_dir, samples=samples, segments=segments,
        fb_info=fb_info, barcode2_sequence_only=barcode2_sequence_only,
        fb_id_column=str(ann.get("id_column", "FB_ID")),
        fb_sequence_column=str(ann.get("sequence_column", "sequence")),
        fb_name_column=str(ann.get("name_column", "antibody_name")),
        allow_duplicate_names=_as_bool(ann.get("allow_duplicate_names", False)),
        umi_method=method, umi_max_distance=umi_max_distance,
        downsample_enabled=_as_bool(ds.get("enabled", True)), downsample_ratios=sorted(set(ratios)),
        downsample_seed=_integer(ds.get("random_seed", 12345), "downsample.random_seed"),
        quick_test_enabled=_as_bool(quick.get("enabled", True)), quick_test_reads=quick_test_reads,
        quick_test_min_valid_rate=qc_thresholds["min_valid_rate"],
        quick_test_max_n_base_rate=qc_thresholds["max_n_base_rate"],
        threads=threads, barcode_workers=barcode_workers, downsample_workers=downsample_workers,
        umi_aggregation_workers=umi_aggregation_workers,
        umi_aggregation_backend=umi_aggregation_backend, umi_sort_memory_mb=umi_sort_memory_mb,
        umi_workers=umi_workers, umi_batch_size=umi_batch_size,
        **umi_budgets,
        umi_sqlite_cache_mb=umi_sqlite_cache_mb, matrix_memory_mb=matrix_memory_mb,
        pi_seq_enabled=pi_seq_enabled, pi_seq_dominance_threshold=pi_seq_dominance_threshold,
        pb_cb_enabled=pb_cb_enabled, pb_cb_dominance_threshold=pb_cb_dominance_threshold,
        chunk_size=chunk_size, extraction_preview_reads=preview_reads,
        compression_level=compression_level,
        parquet_compression=parquet_compression,
        overwrite=_as_bool(resume.get("overwrite", False)), resume_enabled=_as_bool(resume.get("enabled", True)),
        trace_enabled=_as_bool(output.get("correction_trace", True)),
        quick_test_read_details=_as_bool(output.get("quick_test_read_details", False)), raw=raw,
    )
    if check_files:
        missing = [p for s in samples for p in (s.r1, s.r2) if not p.is_file()]
        missing += [s.whitelist for s in segments if s.whitelist and not s.whitelist.is_file()]
        if fb_info and not fb_info.is_file():
            missing.append(fb_info)
        if missing:
            raise ConfigError("Missing input file(s):\n  " + "\n  ".join(str(p) for p in missing))
        for segment in segments:
            if segment.whitelist:
                sequences = [line.strip().split("\t")[-1] for line in segment.whitelist.read_text().splitlines()
                             if line.strip() and not line.startswith("#")]
                if not sequences or any(len(value) != segment.length for value in sequences):
                    raise ConfigError(f"Segment {segment.name}: whitelist entries must have length {segment.length}")
    return config
