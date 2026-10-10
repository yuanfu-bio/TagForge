from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
from pathlib import Path

from .columnar import ParquetFragmentWriter, open_records
from .config import ConfigError
from .io_utils import atomic_text, open_tsv, sample_dirs, write_tsv


def filter_pi_seq_molecules(config, sample_name: str, source: Path):
    """Enforce one PB per FB/UMI molecule and emit PI-seq QC metrics.

    Multi-PB molecules retain their highest-read PB only when it exceeds the
    configured dominance threshold. Otherwise the molecule is removed. GNU
    sort keeps the operation disk-backed for multi-million-molecule samples.
    """
    dirs = sample_dirs(config.output_dir, sample_name)
    qc_path = dirs["detail"] / f"{sample_name}.pi_seq_qc.tsv"
    scratch = Path(tempfile.mkdtemp(prefix=".pi-seq.", dir=dirs["detail"]))
    raw_path = scratch / "molecules.tsv"
    sorted_path = scratch / "molecules.sorted.tsv"
    filtered_parquet = dirs["detail"] / f"{sample_name}.molecule_detail.rmMP.parquet"
    pb_col = config.target_name("barcode1")
    fb_col = f"{config.target_name('barcode2')}_name"
    fields = [pb_col, fb_col, "corrected_umi", "reads_count", "raw_umi_count"]
    total_molecules = 0
    try:
        with raw_path.open("w", encoding="utf-8", newline="") as handle:
            for row in open_records(source):
                values = (row[fb_col], row["corrected_umi"], row[pb_col], row["reads_count"], row["raw_umi_count"])
                if any("\t" in value or "\n" in value or "\r" in value for value in values):
                    raise ConfigError("PI-seq PB, FB name, and UMI values must not contain tabs or newlines")
                handle.write("\t".join(values) + "\n")
                total_molecules += 1
        sort = shutil.which("sort")
        if sort is None:
            raise ConfigError("PI-seq multi-PB filtering requires GNU sort in PATH")
        command = [sort, "--temporary-directory", str(scratch), "--field-separator", "\t",
                   "--key", "1,1", "--key", "2,2", "--key", "4,4nr", "--key", "3,3", str(raw_path)]
        with sorted_path.open("w", encoding="utf-8", newline="") as handle:
            subprocess.run(command, check=True, stdout=handle)

        total_fb_umis = multi_pi_fb_umis = dominant_multi_pi_fb_umis = retained_molecules = 0
        writer = ParquetFragmentWriter(filtered_parquet, fields, getattr(config, "chunk_size", 200000), compression_level=getattr(config, "compression_level", 3), compression=getattr(config, "parquet_compression", "snappy"))

        def consume(group):
            nonlocal total_fb_umis, multi_pi_fb_umis, dominant_multi_pi_fb_umis, retained_molecules
            if not group:
                return
            total_fb_umis += 1
            if len(group) == 1:
                kept = group[0]
            else:
                multi_pi_fb_umis += 1
                kept = group[0]
                if kept[3] / sum(row[3] for row in group) <= config.pi_seq_dominance_threshold:
                    return
                dominant_multi_pi_fb_umis += 1
            writer.writerow({pb_col: kept[2], fb_col: kept[0], "corrected_umi": kept[1],
                             "reads_count": kept[3], "raw_umi_count": kept[4]})
            retained_molecules += 1

        current_key = None
        group = []
        with sorted_path.open("r", encoding="utf-8", newline="") as handle:
            for line in handle:
                fb, umi, pb, reads, raw_umis = line.rstrip("\n").split("\t")
                row = (fb, umi, pb, int(reads), int(raw_umis))
                key = row[:2]
                if current_key is not None and key != current_key:
                    consume(group)
                    group = []
                current_key = key
                group.append(row)
        consume(group)
        writer.close()
        removed_fb_umis = total_fb_umis - retained_molecules
        row = {
            "sample": sample_name, "total_molecules": total_molecules, "total_fb_umis": total_fb_umis,
            "multi_pi_fb_umis": multi_pi_fb_umis,
            "multi_pi_ratio": f"{multi_pi_fb_umis / total_fb_umis:.6f}" if total_fb_umis else "0.000000",
            "mp_ratio": f"{multi_pi_fb_umis / total_fb_umis:.6f}" if total_fb_umis else "0.000000",
            "dominant_multi_pi_fb_umis": dominant_multi_pi_fb_umis,
            "dominant_ratio": f"{dominant_multi_pi_fb_umis / multi_pi_fb_umis:.6f}" if multi_pi_fb_umis else "0.000000",
            "retained_molecules": retained_molecules, "removed_fb_umis": removed_fb_umis,
            "removed_ratio": f"{removed_fb_umis / total_fb_umis:.6f}" if total_fb_umis else "0.000000",
            "dominance_threshold": config.pi_seq_dominance_threshold,
        }
        write_tsv(qc_path, list(row), [row])
        return row
    finally:
        shutil.rmtree(scratch)
