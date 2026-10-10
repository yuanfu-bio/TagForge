"""Public deterministic synthetic example, available in installed wheels."""
from __future__ import annotations

import gzip
from pathlib import Path

from .config import ConfigError


RECORDS = (
    ("read001", "ACGTACGTAGGTC", "AAAAAAAACCGG"),
    ("read002", "ACGTACGTAGGTC", "AAAAAAAACCGG"),
    ("read003", "ACGTACGTAGGTC", "AAAAATAACCGG"),
    ("read004", "TGCATGCAAGGTC", "CCCCCCTTGGCC"),
    ("read005", "TGCATGCAAGGTC", "CCCCCCTTGGCC"),
    ("read006", "ACGTACGAAGGTC", "GGGGGGAACCGG"),
    ("read007", "NACGTACGTAGGTC", "TTTTTTAACCGG"),
    ("read008", "TGCATGCAAGGTC", "ACACACTTGGCT"),
    ("read009", "ACGTACGTNNNNN", "TATATAAACCGG"),
)


def generate_example(directory: Path, config_text: str, *, force: bool = False):
    """Write a standalone project; reject existing targets before any writes."""
    directory = directory.resolve()
    texts = {
        "00_config/CELL_WL.txt": "ACGTACGT\nTGCATGCA\n",
        "00_config/FB_WL.txt": "AACCGG\nTTGGCC\n",
        "00_config/FB_info.tsv": "FB_ID\tsequence\tantibody_name\nFB0001\tAACCGG\tCD4\nFB0002\tTTGGCC\tCD8A\n",
        "00_config/config.yaml": config_text.replace("tests/fixtures/small_fastq/raw", "01_raw").replace("tests/fixtures/small_fastq/config/", "00_config/"),
    }
    fastqs = [directory / f"01_raw/example/example_raw_{mate}.fq.gz" for mate in (1, 2)]
    targets = [directory / name for name in texts] + fastqs
    existing = [str(path) for path in targets if path.exists()]
    if existing and not force:
        raise ConfigError("Refusing to overwrite existing example files (use --force): " + ", ".join(existing))
    for name, content in texts.items():
        path = directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    for mate, path in enumerate(fastqs, 1):
        payload = "".join(f"@{name}/{mate}\n{r1 if mate == 1 else r2}\n+\n{'I' * len(r1 if mate == 1 else r2)}\n" for name, r1, r2 in RECORDS)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(gzip.compress(payload.encode("ascii"), compresslevel=6, mtime=0))
    return directory / "00_config/config.yaml"
