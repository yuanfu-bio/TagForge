"""Small YAML subset loader used when PyYAML is unavailable.

It supports the mappings, lists, booleans, nulls, numbers, quoted strings and
inline lists used by TagForge configuration files. Production installations
may install the ``yaml`` extra for full YAML syntax.
"""

from __future__ import annotations

import ast
import re
from typing import Any


def _flow_parts(value: str):
    parts, current, quote, depth = [], [], None, 0
    for character in value:
        if quote:
            current.append(character)
            if character == quote:
                quote = None
            continue
        if character in {'"', "'"}:
            quote = character; current.append(character)
        elif character in "[{(":
            depth += 1; current.append(character)
        elif character in "]})":
            depth -= 1; current.append(character)
        elif character == "," and depth == 0:
            parts.append("".join(current).strip()); current = []
        else:
            current.append(character)
    if current:
        parts.append("".join(current).strip())
    return parts


def _scalar(value: str) -> Any:
    value = value.strip()
    if not value:
        return {}
    low = value.lower()
    if low in {"null", "none", "~"}:
        return None
    if low in {"true", "false"}:
        return low == "true"
    if value[0:1] in {'"', "'"}:
        return ast.literal_eval(value)
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1].strip()
        return [] if not inner else [_scalar(part) for part in _flow_parts(inner)]
    if value.startswith("{") and value.endswith("}"):
        result = {}
        for part in _flow_parts(value[1:-1]):
            if ":" not in part:
                raise ValueError(f"Malformed inline YAML mapping: {value}")
            key, item = part.split(":", 1)
            key = key.strip().strip("'\"")
            result[key] = _scalar(item)
        return result
    if re.fullmatch(r"[-+]?\d+", value):
        return int(value)
    if re.fullmatch(r"[-+]?(?:\d+\.\d*|\.\d+)(?:[eE][-+]?\d+)?", value):
        return float(value)
    return value


def loads(text: str) -> Any:
    raw_lines = []
    for lineno, raw in enumerate(text.splitlines(), 1):
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip(" "))
        if "\t" in raw[:indent]:
            raise ValueError(f"Tabs are not supported in YAML indentation (line {lineno})")
        content = raw.strip()
        raw_lines.append((indent, content, lineno))
    if not raw_lines:
        return {}

    def parse_block(index: int, indent: int):
        is_list = raw_lines[index][1].startswith("-")
        container: Any = [] if is_list else {}
        while index < len(raw_lines):
            level, content, lineno = raw_lines[index]
            if level < indent:
                break
            if level > indent:
                raise ValueError(f"Unexpected indentation on YAML line {lineno}")
            if is_list:
                if not content.startswith("-"):
                    break
                rest = content[1:].strip()
                if not rest:
                    if index + 1 >= len(raw_lines) or raw_lines[index + 1][0] <= indent:
                        container.append(None); index += 1; continue
                    item, index = parse_block(index + 1, raw_lines[index + 1][0])
                    container.append(item); continue
                if ":" in rest:
                    key, val = rest.split(":", 1)
                    item = {key.strip(): _scalar(val)}
                    index += 1
                    while index < len(raw_lines) and raw_lines[index][0] > indent:
                        child_indent = raw_lines[index][0]
                        child, next_index = parse_block(index, child_indent)
                        if not isinstance(child, dict):
                            raise ValueError(f"Expected mapping below list item on line {lineno}")
                        item.update(child); index = next_index
                    container.append(item); continue
                container.append(_scalar(rest)); index += 1
            else:
                if content.startswith("-") or ":" not in content:
                    break
                key, val = content.split(":", 1)
                key = key.strip()
                if not key:
                    raise ValueError(f"Empty YAML key on line {lineno}")
                index += 1
                if val.strip():
                    container[key] = _scalar(val)
                elif index < len(raw_lines) and raw_lines[index][0] > indent:
                    child, index = parse_block(index, raw_lines[index][0])
                    container[key] = child
                else:
                    container[key] = {}
        return container, index

    parsed, end = parse_block(0, raw_lines[0][0])
    if end != len(raw_lines):
        raise ValueError(f"Could not parse YAML near line {raw_lines[end][2]}")
    return parsed
