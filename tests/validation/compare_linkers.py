"""Compare linker policies against declared coordinates and known whitelists.

This is an evidence report, not a claim that whitelist or fixed-position
agreement is biological ground truth. Run from a checkout with PYTHONPATH=src.
"""
import argparse
from collections import Counter
import csv
from dataclasses import replace
from itertools import islice
import json
from pathlib import Path

from tagforge.barcode_correct import _build_correctors, evaluate_extracted_read, load_fb_annotation
from tagforge.config import load_config
from tagforge.extract import _extract_record, decode_method_payload, decode_segment_payload
from tagforge.fastq import paired_fastq
from tagforge.io_utils import atomic_text


def review_changes(config, configs, columns, correctors, annotation, output):
    """Audit every changed read against configured coordinates and trace evidence."""
    reviews = []
    totals = Counter()
    for sample in config.samples:
        with (output / f"{sample.sample}.changed_valid_reads.tsv").open() as handle:
            changed = {row["read_id"]: row["change"] for row in csv.DictReader(handle, delimiter="\t")}
        for record in paired_fastq(sample.r1, sample.r2):
            if record.read_id not in changed:
                continue
            outcomes = []
            for policy in configs:
                row, runtime = _extract_record(record, policy.segments, columns)
                counters = {s.name: Counter() for s in policy.segments if s.target != "umi"}
                summary = Counter()
                eligible, trace = evaluate_extracted_read(row, policy, correctors, annotation, counters, summary)
                outcomes.append({"eligible": eligible is not None, "failure": [key for key in summary if key.startswith("failure_")],
                                 "assignment": eligible, "trace": trace, "extraction": runtime})
            old = {r["segment_name"]: r for r in outcomes[0]["trace"]}
            new = {r["segment_name"]: r for r in outcomes[1]["trace"]}
            evidence = []
            for segment in config.segments:
                if segment.target == "umi":
                    continue
                before, after = old.get(segment.name, {}), new.get(segment.name, {})
                if before.get("corrected_sequence") == after.get("corrected_sequence"):
                    continue
                seq = record.r1_seq if segment.read == "R1" else record.r2_seq
                fixed = seq[segment.start:segment.start + segment.length] if segment.start is not None else ""
                exact = bool(fixed) and fixed in correctors[segment.name].value_set
                evidence.append({"segment": segment.name, "fixed_sequence": fixed, "fixed_in_whitelist": exact,
                                 "legacy_agrees_fixed": exact and before.get("corrected_sequence") == fixed,
                                 "new_agrees_fixed": exact and after.get("corrected_sequence") == fixed})
            review = {"sample": sample.sample, "read_id": record.read_id, "change": changed[record.read_id],
                      "r1": record.r1_seq, "r2": record.r2_seq, "legacy": outcomes[0], "new": outcomes[1], "coordinate_evidence": evidence}
            reviews.append(review)
            totals[review["change"]] += 1
            for item in evidence:
                if item["fixed_in_whitelist"]:
                    totals["changed_segments_with_fixed_exact_reference"] += 1
                    totals["legacy_agrees_fixed"] += item["legacy_agrees_fixed"]
                    totals["new_agrees_fixed"] += item["new_agrees_fixed"]
                else:
                    totals["changed_segments_without_fixed_exact_reference"] += 1
        expected = set(changed)
        actual = {row["read_id"] for row in reviews if row["sample"] == sample.sample}
        if actual != expected:
            raise RuntimeError(f"Missing changed read review for {sample.sample}: {expected - actual}")
    result = {"totals": dict(totals), "limitation": "Configured coordinate and whitelist agreement is independent supporting evidence, not biological truth.", "reads": reviews}
    (output / "changed-read-review.json").write_text(json.dumps(result, indent=2) + "\n")
    print(dict(totals), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--review-existing", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    output = Path(args.out)
    output.mkdir(parents=True, exist_ok=True)
    annotation, correctors = load_fb_annotation(config), _build_correctors(config)
    columns = {role: config.segment_column(role) for role in ("barcode1", "barcode2", "umi")}
    configs = [replace(config, segments=[replace(s, linker_strategy=strategy) for s in config.segments])
               for strategy in ("legacy_shortest", "length_position")]
    if args.review_existing:
        review_changes(config, configs, columns, correctors, annotation, output)
        return
    summary = []
    for sample in config.samples:
        totals, segments = Counter(), {s.name: Counter() for s in config.segments}
        correction_counters = [{s.name: Counter() for s in config.segments if s.target != "umi"} for _ in configs]
        correction_totals = [Counter(), Counter()]
        with atomic_text(output / f"{sample.sample}.changed_valid_reads.tsv") as handle:
            handle.write("read_id\tchange\tlegacy_barcode1\tnew_barcode1\tlegacy_feature\tnew_feature\n")
            records = paired_fastq(sample.r1, sample.r2)
            if args.limit is not None:
                records = islice(records, args.limit)
            for record in records:
                totals["reads"] += 1
                extracted, valid = [], []
                for index, policy in enumerate(configs):
                    row, runtime = _extract_record(record, policy.segments, columns)
                    eligible, trace = evaluate_extracted_read(row, policy, correctors, annotation,
                                                              correction_counters[index], correction_totals[index])
                    extracted.append((row, runtime))
                    valid.append(eligible)
                    label = "legacy" if index == 0 else "new"
                    totals[label + "_valid_reads"] += eligible is not None
                    for name, stat in runtime.items():
                        segments[name][label + "_linker_success"] += stat["success"] and stat["method"] == "linker"
                        segments[name][label + "_fixed_success"] += stat["success"] and stat["method"] == "fixed"
                if valid[0] != valid[1]:
                    change = "rescued" if valid[0] is None else "lost" if valid[1] is None else "changed_assignment"
                    totals[change] += 1
                    b1, fb = config.target_name("barcode1"), config.target_name("barcode2") + "_name"
                    handle.write("\t".join([record.read_id, change,
                        (valid[0] or {}).get(b1, ""), (valid[1] or {}).get(b1, ""),
                        (valid[0] or {}).get(fb, ""), (valid[1] or {}).get(fb, "")]) + "\n")
                for segment in config.segments:
                    if segment.target == "umi" or segment.start is None:
                        continue
                    sequence = record.r1_seq if segment.read == "R1" else record.r2_seq
                    fixed = sequence[segment.start:segment.start + segment.length]
                    if fixed not in correctors[segment.name].value_set:
                        continue
                    segments[segment.name]["fixed_exact_reference_reads"] += 1
                    for index, policy in enumerate(configs):
                        row, runtime = extracted[index]
                        role_segments = [s for s in policy.segments if s.target == segment.target]
                        raw = decode_segment_payload(row[columns[segment.target]], role_segments)[segment.name]
                        method = runtime[segment.name]["method"]
                        extra = segment.correction.max_shift if method == "fixed" and segment.correction.enabled and segment.correction.allow_shift else 0
                        anchor = min(extra, segment.start) if method == "fixed" else 0
                        corrected = correctors[segment.name].correct(raw, anchor)
                        label = "legacy" if index == 0 else "new"
                        segments[segment.name][label + "_fixed_reference_agreement"] += corrected.corrected_sequence == fixed
                        segments[segment.name][label + "_fixed_reference_conflict"] += corrected.success and corrected.corrected_sequence != fixed
        result = {"sample": sample.sample, "totals": dict(totals),
                  "segments": {name: dict(counts) for name, counts in segments.items()}}
        summary.append(result)
        (output / "comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(sample.sample, dict(totals), flush=True)


if __name__ == "__main__":
    main()
