from __future__ import annotations

import shlex
from pathlib import Path

from .config import TagForgeConfig
from .io_utils import atomic_text


def _safe(value: str, option: str) -> str:
    if "\n" in value or "\r" in value:
        raise ValueError(f"{option} must not contain newlines")
    return value


def _common_header(
    out: Path,
    job_name: str,
    threads: int,
    memory: str,
    walltime: str,
    *,
    nodes: int,
    ntasks: int,
    partition: str | None,
    account: str | None,
    qos: str | None,
    constraint: str | None,
    gres: str | None,
    mail_user: str | None,
    mail_type: str | None,
    extra_sbatch: list[str] | None,
    output_name: str,
    error_name: str,
    array: str | None = None,
) -> str:
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
    array_directive = f"#SBATCH --array={array}\n" if array else ""
    return f"""#!/bin/bash
#SBATCH --job-name={_safe(job_name, 'job-name')}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={memory}
#SBATCH --time={walltime}
{array_directive}#SBATCH --output={out / 'logs' / output_name}
#SBATCH --error={out / 'logs' / error_name}
{optional_directives}{extra_directives}
"""


def _activation(conda_env: str | None) -> str:
    if not conda_env:
        return ""
    env = shlex.quote(_safe(conda_env, "conda-env"))
    return (
        'if ! command -v conda >/dev/null 2>&1; then\n'
        '  echo "conda is not available in the Slurm job PATH" >&2\n  exit 127\nfi\n'
        'source "$(conda info --base)/etc/profile.d/conda.sh"\n'
        f"conda activate {env}\n\n"
    )


def _run_command(config_path: Path, sample: str, threads: int, skip_quick_test: bool) -> str:
    command = (
        f"tagforge run --config {shlex.quote(str(config_path))} "
        f"--sample {sample} --threads {threads}"
    )
    if skip_quick_test:
        command += " --skip-quick-test"
    return command


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
    mode: str = "array",
    array_limit: int | None = None,
):
    if threads < 1 or nodes < 1 or ntasks < 1:
        raise ValueError("threads, nodes, and ntasks must be >= 1")
    if mode not in {"array", "per-sample"}:
        raise ValueError("mode must be 'array' or 'per-sample'")
    if array_limit is not None and array_limit < 1:
        raise ValueError("array_limit must be >= 1")
    out.mkdir(parents=True, exist_ok=True); (out / "logs").mkdir(exist_ok=True)
    activation = _activation(conda_env)
    if mode == "array":
        samples_tsv = out / "samples.tsv"
        with atomic_text(samples_tsv) as handle:
            handle.write("index\tsample\n")
            for index, sample in enumerate(config.samples, 1):
                name = _safe(sample.sample, "sample")
                if "\t" in name:
                    raise ValueError("sample names must not contain tabs for samples.tsv")
                handle.write(f"{index}\t{name}\n")
        array_range = f"1-{len(config.samples)}"
        if array_limit is not None:
            array_range += f"%{array_limit}"
        array_script = out / "tagforge_array.slurm"
        command = _run_command(config.path, '"$sample"', threads, skip_quick_test)
        text = _common_header(
            out, "tagforge_array", threads, memory, walltime,
            nodes=nodes, ntasks=ntasks, partition=partition, account=account,
            qos=qos, constraint=constraint, gres=gres, mail_user=mail_user,
            mail_type=mail_type, extra_sbatch=extra_sbatch,
            output_name="array.%A_%a.out", error_name="array.%A_%a.err",
            array=array_range,
        ) + f"""
set -euo pipefail

{activation}samples_tsv={shlex.quote(str(samples_tsv))}
sample=$(awk -v i="$SLURM_ARRAY_TASK_ID" 'NR==i+1 {{print $2}}' "$samples_tsv")
if [[ -z "${{sample}}" ]]; then
  echo "No sample found for SLURM_ARRAY_TASK_ID=${{SLURM_ARRAY_TASK_ID}}" >&2
  exit 2
fi

{command}
"""
        with atomic_text(array_script) as handle:
            handle.write(text)
        array_script.chmod(0o755)
        return [samples_tsv, array_script]
    scripts = []
    for sample in config.samples:
        path = out / f"{sample.sample}.slurm"; scripts.append(path)
        command = _run_command(config.path, shlex.quote(sample.sample), threads, skip_quick_test)
        text = _common_header(
            out, f"tagforge_{sample.sample}", threads, memory, walltime,
            nodes=nodes, ntasks=ntasks, partition=partition, account=account,
            qos=qos, constraint=constraint, gres=gres, mail_user=mail_user,
            mail_type=mail_type, extra_sbatch=extra_sbatch,
            output_name=sample.sample + ".%j.out", error_name=sample.sample + ".%j.err",
        ) + f"""
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
