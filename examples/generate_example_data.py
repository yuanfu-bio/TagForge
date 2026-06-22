"""Generate the tiny deterministic paired FASTQ example."""
from __future__ import annotations

import gzip
from pathlib import Path

ROOT = Path(__file__).parent / "small_fastq"
ROOT.mkdir(parents=True, exist_ok=True)

# CELL + linker in R1; six-base UMI + six-base feature in R2.
records = [
    ("read001", "ACGTACGTAGGTC", "AAAAAAAACCGG"),
    ("read002", "ACGTACGTAGGTC", "AAAAAAAACCGG"),
    ("read003", "ACGTACGTAGGTC", "AAAAATAACCGG"),
    ("read004", "TGCATGCAAGGTC", "CCCC CCTTGGCC".replace(" ", "")),
    ("read005", "TGCATGCAAGGTC", "CCCC CCTTGGCC".replace(" ", "")),
    ("read006", "ACGTACGAAGGTC", "GGGGGGAACCGG"),  # one barcode mismatch
    ("read007", "NACGTACGTAGGTC", "TTTTTTAACCGG"), # one-base forward shift
    ("read008", "TGCATGCAAGGTC", "ACACACTTGGCT"),  # one feature mismatch
    ("read009", "ACGTACGTNNNNN", "TATATAAACCGG"),  # linker fails; fixed fallback succeeds
]

for mate in (1, 2):
    path = ROOT / f"example_raw_{mate}.fq.gz"
    with gzip.open(path, "wt", encoding="ascii", compresslevel=6) as handle:
        for name, r1, r2 in records:
            sequence = r1 if mate == 1 else r2
            handle.write(f"@{name}/{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n")
