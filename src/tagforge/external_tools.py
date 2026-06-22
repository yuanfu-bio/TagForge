from __future__ import annotations

import importlib
import shutil
import subprocess
from dataclasses import dataclass

from .config import ConfigError


@dataclass(frozen=True)
class ExternalToolVersions:
    cutadapt: str
    umi_tools: str


def _command_version(command: str, arguments: list[str]) -> str:
    executable = shutil.which(command)
    if executable is None:
        raise ConfigError(
            f"Required external tool '{command}' was not found in PATH. "
            "Create and activate the Conda environment with: "
            "conda env create -f environment.yml && conda activate tagforge"
        )
    try:
        result = subprocess.run(
            [executable, *arguments], check=True, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ConfigError(f"Could not run required external tool '{command}': {exc}") from exc
    return result.stdout.strip().splitlines()[0]


def check_external_tools() -> ExternalToolVersions:
    missing_modules = []
    for module in ("cutadapt", "umi_tools"):
        try:
            importlib.import_module(module)
        except ImportError:
            missing_modules.append(module)
    if missing_modules:
        raise ConfigError(
            "Required Python package(s) missing from the active environment: "
            + ", ".join(missing_modules)
            + ". Create and activate the Conda environment with: "
              "conda env create -f environment.yml && conda activate tagforge"
        )
    return ExternalToolVersions(
        cutadapt=_command_version("cutadapt", ["--version"]),
        umi_tools=_command_version("umi_tools", ["--version"]),
    )
