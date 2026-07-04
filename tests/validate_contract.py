#!/usr/bin/env python3
"""Validate a poly-suite output contract (pgs_scores.tsv) — schema + gate
invariants. Catches regressions in the deliverable the whole pipeline produces
and that OmniGen consumes. No network/panel needed.

Importable (validate) + CLI: validate_contract.py [results_dir]  (exit 1 on problems).
"""
import sys, os, csv
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
from grade_pgs import CONTRACT_COLS

MATCH_RATE_MIN = 0.75


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def validate(rows, header):
    """rows: list of dicts; header: list of column names. -> list of problem strings."""
    problems = []
    missing = [c for c in CONTRACT_COLS if c not in header]
    if missing:
        problems.append(f"missing columns: {missing}")
    for i, r in enumerate(rows):
        tag = f"row {i} ({r.get('pgs_id','?')})"
        g = r.get("evidence_grade")
        if g not in ("A", "B", "C", "D"):
            problems.append(f"{tag}: bad evidence_grade {g!r}")
        if str(r.get("portability_flag")).lower() not in ("true", "false"):
            problems.append(f"{tag}: portability_flag not boolean ({r.get('portability_flag')!r})")
        if not (r.get("allowed_statement") or "").strip():
            problems.append(f"{tag}: empty allowed_statement")
        pct, ar, mr = _num(r.get("percentile")), _num(r.get("absolute_risk")), _num(r.get("match_rate"))
        if pct is not None and not (0 <= pct <= 100):
            problems.append(f"{tag}: percentile out of range ({pct})")
        if ar is not None and not (0 <= ar <= 1):
            problems.append(f"{tag}: absolute_risk out of range ({ar})")
        if ar is not None and pct is None:
            problems.append(f"{tag}: absolute_risk set but percentile is NA (can't have risk without calibration)")
        # QC gate: low coverage must be reflected as a low-coverage caveat
        if mr is not None and mr < MATCH_RATE_MIN and "coverage" not in (r.get("allowed_statement") or "").lower():
            problems.append(f"{tag}: match_rate {mr:.2f} < {MATCH_RATE_MIN} but no low-coverage caveat")
        # uncalibrated must not carry a percentile-derived absolute risk
        if pct is None and ar is not None:
            problems.append(f"{tag}: uncalibrated (no percentile) yet absolute_risk present")
    return problems


def validate_file(results_dir):
    p = os.path.join(results_dir, "pgs_scores.tsv")
    if not os.path.exists(p):
        return [f"no contract at {p}"], 0
    with open(p) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        rows = list(rdr)
        header = rdr.fieldnames or []
    return validate(rows, header), len(rows)


if __name__ == "__main__":
    rd = sys.argv[1] if len(sys.argv) > 1 else "results"
    problems, n = validate_file(rd)
    if problems:
        print(f"[validate_contract] {len(problems)} problem(s) in {rd}/pgs_scores.tsv ({n} rows):")
        for p in problems:
            print("  ✗", p)
        sys.exit(1)
    print(f"[validate_contract] OK — {n} rows, all {len(CONTRACT_COLS)} columns, gates consistent")
