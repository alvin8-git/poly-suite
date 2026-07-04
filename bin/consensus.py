#!/usr/bin/env python3
"""Multi-score consensus / robustness for poly-suite (review differentiator #3).

For a trait scored by >=2 published PGS, measure agreement across scores — the
pgx multi-caller-concordance idea, applied to PRS. Divergent scores => a low-
robustness flag on the finding, so a percentile that only one score supports is
not presented as if it were settled.

Concordance (calibrated inputs, percentiles present):
  concordance = max(0, 1 - (max_pct - min_pct)/100)     # 1 = identical, 0 = opposite ends
  tertile_concordant = all scores place the sample in the same tertile (<33 / 33-66 / >66)
Uncalibrated inputs (raw AVG only, not comparable across scores): concordance=None,
report n_scores only — real consensus needs calibrated percentiles.

Importable (consensus) and runnable (self-check).
"""
from collections import defaultdict


def _tertile(p):
    return 0 if p < (100 / 3) else (2 if p > (200 / 3) else 1)


def consensus(rows):
    """rows: dicts with trait, pgs_id, percentile (number or None).
    -> {trait: {n_scores, concordance, spread, tertile_concordant, members}}."""
    by_trait = defaultdict(list)
    for r in rows:
        by_trait[r["trait"]].append(r)
    out = {}
    for trait, rs in by_trait.items():
        pcts = [r["percentile"] for r in rs
                if isinstance(r.get("percentile"), (int, float))]
        n = len(rs)
        if len(pcts) >= 2:
            spread = max(pcts) - min(pcts)
            conc = round(max(0.0, 1.0 - spread / 100.0), 3)
            tert = len({_tertile(p) for p in pcts}) == 1
        else:
            spread = conc = tert = None
        out[trait] = dict(n_scores=n, concordance=conc, spread=spread,
                          tertile_concordant=tert,
                          members=[r["pgs_id"] for r in rs])
    return out


if __name__ == "__main__":
    demo = [
        {"trait": "CAD", "pgs_id": "A", "percentile": 82},
        {"trait": "CAD", "pgs_id": "B", "percentile": 88},      # concordant
        {"trait": "T2D", "pgs_id": "C", "percentile": 15},
        {"trait": "T2D", "pgs_id": "D", "percentile": 71},      # discordant
        {"trait": "LDL", "pgs_id": "E", "percentile": 60},      # single score
    ]
    c = consensus(demo)
    assert c["CAD"]["concordance"] > 0.9 and c["CAD"]["tertile_concordant"]
    assert c["T2D"]["concordance"] < 0.6 and not c["T2D"]["tertile_concordant"]
    assert c["LDL"]["n_scores"] == 1 and c["LDL"]["concordance"] is None
    print("consensus self-check OK: "
          f"CAD concordance={c['CAD']['concordance']} (robust), "
          f"T2D concordance={c['T2D']['concordance']} (divergent -> flag), "
          f"LDL n=1 (no consensus)")
