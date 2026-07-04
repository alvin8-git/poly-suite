#!/usr/bin/env python3
"""Ensemble / meta-PGS for poly-suite — combine >=2 validated published scores for
one trait into one combined predictor. The review's recommended novel-score route:
validated inputs, improved robustness, far cheaper to validate than a from-scratch
GWAS-derived score, and submittable back to the PGS Catalog.

Method: score-level linear combination of standardized PRSs (the metaGRS approach,
Inouye 2018 / Abraham 2016):
    meta_z = sum(w_i * z_i) / sum(w_i)
with weights w_i = evidence weight (GWAS N from the Catalog) or equal.
meta_percentile = Phi(meta_z)*100 — APPROXIMATE: the component scores are
correlated, so the exact percentile requires re-standardizing the meta_z
distribution in the reference panel (a calibration follow-on, not done here).

Importable (meta) + CLI: ensemble.py [results_dir]  (reads pgs_scores.tsv z_score).
"""
import sys, os, csv, json
from statistics import NormalDist
from collections import defaultdict

_N = NormalDist()


def meta(rows, weights=None):
    """rows: dicts with trait, pgs_id, z (number or None).
    -> {trait: {n_scores, meta_z, meta_percentile, members, weighting, note}}."""
    by_trait = defaultdict(list)
    for r in rows:
        by_trait[r["trait"]].append(r)
    out = {}
    for trait, rs in by_trait.items():
        zs = [(r["pgs_id"], r["z"], (weights or {}).get(r["pgs_id"], 1.0))
              for r in rs if isinstance(r.get("z"), (int, float))]
        if len(zs) < 2:
            continue
        wsum = sum(w for _, _, w in zs)
        mz = sum(z * w for _, z, w in zs) / wsum
        out[trait] = dict(
            n_scores=len(zs), meta_z=round(mz, 4),
            meta_percentile=round(_N.cdf(mz) * 100, 1),
            members=[p for p, _, _ in zs],
            weighting="evidence" if weights else "equal",
            note="approximate percentile — correlated components; re-standardize "
                 "meta_z in the reference panel for an exact value")
    return out


def meta_from_calibrated(pgs_rows, trait_of, weights=None):
    """Exact reference-based meta-PGS (ensemble v2). Removes the v1 'approximate'
    caveat: instead of assuming meta_z ~ N(0,1), rank the target's meta_z against
    the reference panel's meta_z distribution (both from the --run_ancestry output).

    pgs_rows: dicts with sampleset, IID, PGS, z. trait_of: pgs_id -> trait.
    -> {trait: {meta_z, exact_percentile, n_scores, n_reference, members}}."""
    by_sample, traits_pgs = defaultdict(dict), defaultdict(set)
    for r in pgs_rows:
        tr, z = trait_of.get(r["PGS"]), r.get("z")
        if tr is None or z is None:
            continue
        by_sample[(r["sampleset"], r["IID"])][r["PGS"]] = z
        traits_pgs[tr].add(r["PGS"])
    out = {}
    for trait, pgs_set in traits_pgs.items():
        pgs_list = sorted(pgs_set)
        if len(pgs_list) < 2:
            continue
        w = [(weights or {}).get(p, 1.0) for p in pgs_list]
        wsum = sum(w)

        def mz(zmap):
            vals = [zmap.get(p) for p in pgs_list]
            if any(v is None for v in vals):
                return None
            return sum(v * wi for v, wi in zip(vals, w)) / wsum

        ref = [m for (ss, _), z in by_sample.items()
               if ss.lower() == "reference" for m in [mz(z)] if m is not None]
        if not ref:
            continue
        for (ss, _), z in by_sample.items():
            if ss.lower() == "reference":
                continue
            tmz = mz(z)
            if tmz is None:
                continue
            pct = 100.0 * sum(1 for x in ref if x < tmz) / len(ref)
            out[trait] = dict(meta_z=round(tmz, 4), exact_percentile=round(pct, 1),
                              n_scores=len(pgs_list), n_reference=len(ref),
                              members=pgs_list)
    return out


def _load(results_dir):
    contract = os.path.join(results_dir, "pgs_scores.tsv")
    if not os.path.exists(contract):
        raise SystemExit(f"no contract at {contract}")
    rows = []
    with open(contract) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            z = r.get("z_score")
            try:
                z = float(z)
            except (TypeError, ValueError):
                z = None
            rows.append({"trait": r["trait"], "pgs_id": r["pgs_id"], "z": z})
    weights = None
    mp = os.path.join(results_dir, "pgs_catalog_meta.json")
    if not os.path.exists(mp):
        mp = "results/pgs_catalog_meta.json"
    if os.path.exists(mp):
        cat = json.load(open(mp))
        weights = {pid: (m.get("gwas_n") or 1.0) for pid, m in cat.items()}
    return rows, weights


if __name__ == "__main__":
    rd = sys.argv[1] if len(sys.argv) > 1 else "results"
    rows, weights = _load(rd)
    m = meta(rows, weights=weights)
    if not m:
        print("[ensemble] no trait has >=2 calibrated scores (z_score) — "
              "score >=2 PGS per trait, with ancestry calibration, to build a meta-PGS")
        sys.exit(0)
    out = os.path.join(rd, "meta_scores.tsv")
    with open(out, "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["trait", "n_scores", "meta_z", "meta_percentile", "weighting", "members"])
        for t, v in m.items():
            w.writerow([t, v["n_scores"], v["meta_z"], v["meta_percentile"],
                        v["weighting"], ",".join(v["members"])])
    print(f"[ensemble] meta-PGS -> {out}")
    for t, v in m.items():
        print(f"  {t}: meta_z={v['meta_z']} -> ~{v['meta_percentile']}th pct "
              f"({v['weighting']}-weighted over {v['n_scores']}: {','.join(v['members'])})")
