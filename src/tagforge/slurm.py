from __future__ import annotations

import shlex
from pathlib import Path

from .config import TagForgeConfig
from .io_utils import atomic_text


def _safe(value: str, option: str) -> str:
    if "\n" in value or "\r" in value:
        raise ValueError(f"{option} must not contain newlines")
    return value


def make_slurm(
    config: TagForgeConfig,
    out: Path,
    threads: int,
    memory: str,
    walltime: str,
    *,
    partition: str | None = None,
    account: str | None = None,
    qos: str | None = None,
    constraint: str | None = None,
    gres: str | None = None,
    nodes: int = 1,
    ntasks: int = 1,
    mail_user: str | None = None,
    mail_type: str | None = None,
    conda_env: str | None = None,
    skip_quick_test: bool = False,
    extra_sbatch: list[str] | None = None,
):
    if threads < 1 or nodes < 1 or ntasks < 1:
        raise ValueError("threads, nodes, and ntasks must be >= 1")
    out.mkdir(parents=True, exist_ok=True); (out / "logs").mkdir(exist_ok=True)
    optional = [
        ("partition", partition), ("account", account), ("qos", qos),
        ("constraint", constraint), ("gres", gres),
        ("mail-user", mail_user), ("mail-type", mail_type),
    ]
    optional_directives = "".join(
        f"#SBATCH --{name}={_safe(value, name)}\n" for name, value in optional if value
    )
    extra_directives = "".join(
        f"#SBATCH {_safe(value, 'extra-sbatch')}\n" for value in (extra_sbatch or [])
    )
    activation = ""
    if conda_env:
        env = shlex.quote(_safe(conda_env, "conda-env"))
        activation = (
            'if ! command -v conda >/dev/null 2>&1; then\n'
            '  echo "conda is not available in the Slurm job PATH" >&2\n  exit 127\nfi\n'
            'source "$(conda info --base)/etc/profile.d/conda.sh"\n'
            f"conda activate {env}\n\n"
        )
    scripts = []
    for sample in config.samples:
        path = out / f"{sample.sample}.slurm"; scripts.append(path)
        command = f"tagforge run --config {shlex.quote(str(config.path))} --sample {shlex.quote(sample.sample)} --threads {threads}"
        if skip_quick_test:
            command += " --skip-quick-test"
        text = f"""#!/bin/bash
#SBATCH --job-name=tagforge_{sample.sample}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={memory}
#SBATCH --time={walltime}
#SBATCH --output={out / 'logs' / (sample.sample + '.%j.out')}
#SBATCH --error={out / 'logs' / (sample.sample + '.%j.err')}
{optional_directives}{extra_directives}

set -euo pipefail

{activation}{command}
"""
        with atomic_text(path) as handle: handle.write(text)
        path.chmod(0o755)
    submit = out / "submit_all.sh"
    with atomic_text(submit) as handle:
        handle.write("#!/bin/bash\nset -euo pipefail\n")
        for script in scripts: handle.write(f"sbatch {shlex.quote(str(script))}\n")
    submit.chmod(0o755)
    return scripts + [submit]
