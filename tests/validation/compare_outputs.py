"""Compare logical Parquet rows across count-schema changes using EXCEPT ALL."""
import argparse
import csv
import json
from pathlib import Path

import duckdb
from tagforge.columnar import sql_reader


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--optimized", type=Path, required=True)
    parser.add_argument("--sample", action="append", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    connection = duckdb.connect()
    connection.execute("SET threads=1")
    connection.execute("SET memory_limit='512MB'")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    scratch = args.out.parent / "comparison-tmp"
    scratch.mkdir(exist_ok=True)
    connection.execute("SET temp_directory='" + str(scratch.resolve()).replace("'", "''") + "'")
    results = []
    def tsv(path):
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    for sample in args.sample:
        roots = [directory / "02_output" / sample for directory in (args.baseline, args.optimized)]
        result = {"sample": sample, "tables": {}}
        for folder, suffix in (("03_molecules", "valid_reads"), ("03_molecules", "molecule_detail"),
                               ("03_molecules", "molecule_detail.rmMP"), ("04_matrix", "raw_molecule_pairs")):
            paths = [root / folder / f"{sample}.{suffix}.parquet" for root in roots]
            if suffix.endswith(".rmMP") and not paths[0].exists():
                continue
            sources = [sql_reader(path) for path in paths]
            columns = [row[0] for row in connection.execute("DESCRIBE SELECT * FROM " + sources[0]).fetchall()]
            quoted = ['"' + col.replace('"', '""') + '"' for col in columns]
            projection = ", ".join(f"CAST({col} AS VARCHAR) AS {col}" for col in quoted)
            queries = ["SELECT " + projection + " FROM " + source for source in sources]
            difference = connection.execute("SELECT COUNT(*) FROM ((" + queries[0] + " EXCEPT ALL " + queries[1] + ") UNION ALL (" + queries[1] + " EXCEPT ALL " + queries[0] + "))").fetchone()[0]
            rows = connection.execute("SELECT COUNT(*) FROM " + sources[1]).fetchone()[0]
            result["tables"][suffix] = {"rows": rows, "symmetric_difference_rows": difference}
            assert difference == 0, (sample, suffix, difference)
        summaries = [root / "04_matrix" / f"{sample}.raw_matrix_summary.json" for root in roots]
        result["raw_matrix_summary_equal"] = json.loads(summaries[0].read_text()) == json.loads(summaries[1].read_text())
        assert result["raw_matrix_summary_equal"], sample
        for folder, suffix in (("03_molecules", "pi_seq_qc.tsv"), ("05_downsample", "downsample_metrics.tsv")):
            paths = [root / folder / f"{sample}.{suffix}" for root in roots]
            if not paths[0].exists():
                continue
            equal = tsv(paths[0]) == tsv(paths[1])
            result[suffix + "_equal"] = equal
            assert equal, (sample, suffix)
        results.append(result)
        args.out.write_text(json.dumps(results, indent=2) + "\n")
        print(json.dumps(result), flush=True)
    connection.close()
    scratch.rmdir()


if __name__ == "__main__":
    main()
