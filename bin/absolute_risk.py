#!/usr/bin/env python3
"""Absolute-risk conversion for poly-suite (review differentiator #4).

percentile -> Z -> relative risk (via published per-SD OR/HR) -> absolute risk
(x baseline incidence). The relative risk is calibrated so the POPULATION MEAN
RR = 1: for a log-normal PGS with per-SD log-effect beta, an individual at Z SDs
has RR = exp(beta*Z - beta^2/2), whose expectation over N(0,1) is 1. So a median
(50th-pct) person sits slightly BELOW baseline and the high-PGS tail carries the
excess — which is correct, not a bug.

Resources (illustrative launch values, sourced, verify per release):
  resources/pgs_effect.tsv        per-SD OR/HR per score
  resources/baseline_incidence.tsv  population absolute risk per trait/sex/ancestry

Importable (load_effects/load_baselines/lookup_baseline/risk) and runnable
(self-check).  Run: python3 bin/absolute_risk.py
"""
import os, csv, math
from statistics import NormalDist

_HERE = os.path.dirname(os.path.abspath(__file__))
_RES = os.path.join(os.path.dirname(_HERE), "resources")
EFFECT_TSV = os.path.join(_RES, "pgs_effect.tsv")
BASELINE_TSV = os.path.join(_RES, "baseline_incidence.tsv")

# sex-specific traits use the sex-specific baseline
TRAIT_SEX = {"breast cancer": "female", "prostate cancer": "male"}
_N = NormalDist()


def load_effects(path=EFFECT_TSV):
    """-> (by_pgs_id, by_trait) per-SD effect sizes."""
    by_pgs, by_trait = {}, {}
    if not os.path.exists(path):
        return by_pgs, by_trait
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                per = float(r["per_sd"])
            except (ValueError, KeyError):
                continue
            by_pgs[r["pgs_id"]] = per
            by_trait.setdefault(r["trait"].strip().lower(), per)
    return by_pgs, by_trait


def load_baselines(path=BASELINE_TSV):
    """-> {(trait_lower, sex, ancestry): (risk, source_pmid, note)}."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                out[(r["trait"].strip().lower(), r["sex"], r["ancestry"])] = (
                    float(r["baseline_risk"]), r.get("source_pmid", ""), r.get("note", ""))
            except (ValueError, KeyError):
                continue
    return out


def z_from_percentile(p):
    return _N.inv_cdf(min(max(p / 100.0, 1e-6), 1 - 1e-6))


def risk(percentile, per_sd, baseline):
    """(absolute_risk, relative_risk, z). RR calibrated so population mean = 1,
    then applied on the ODDS scale — the correct bounded model for common
    diseases (linear baseline*RR overshoots when baseline is high, e.g. CAD ~49%,
    giving impossible >90% risks). odds = baseline/(1-baseline) * RR; AR = odds/(1+odds)."""
    z = z_from_percentile(percentile)
    beta = math.log(per_sd)
    rr = math.exp(beta * z - beta * beta / 2.0)
    odds = (baseline / (1.0 - baseline)) * rr
    return odds / (1.0 + odds), rr, z


def lookup_baseline(baselines, trait, ancestry="overall", sex=None):
    """(risk, source, note), sex — resolves sex (explicit -> trait-implied ->
    'overall') and falls back on ancestry. Returns None for a sex-dimorphic
    trait when sex is unknown and no 'overall' baseline exists (honest)."""
    tl = trait.strip().lower()
    # trait-mandated sex wins (breast=female, prostate=male); otherwise the
    # sample's sex (for sex-dimorphic traits like CAD); else population-overall.
    sx = TRAIT_SEX.get(tl) or sex or "overall"
    for s in (sx, "overall"):
        for anc in (ancestry, "overall"):
            hit = baselines.get((tl, s, anc))
            if hit:
                return hit, s
    return None, sx


def estimate(trait, percentile, pgs_id=None, ancestry="overall", sex=None,
             eff=None, base=None):
    """High-level: -> dict(absolute_risk, relative_risk, baseline, sex, source)
    or None when effect or baseline is unavailable (report percentile only)."""
    if percentile is None:
        return None
    if eff is None:
        eff = load_effects()
    if base is None:
        base = load_baselines()
    by_pgs, by_trait = eff
    per = by_pgs.get(pgs_id) if pgs_id else None
    if per is None:
        per = by_trait.get(trait.strip().lower())
    if per is None:
        return None
    hit, sex = lookup_baseline(base, trait, ancestry, sex)
    if hit is None:
        return None
    baseline, src, note = hit
    ar, rr, z = risk(percentile, per, baseline)
    return dict(absolute_risk=ar, relative_risk=rr, baseline=baseline,
                per_sd=per, sex=sex, baseline_source=src, baseline_note=note)


if __name__ == "__main__":
    eff = load_effects()
    base = load_baselines()
    assert eff[1], "no effects loaded"
    assert base, "no baselines loaded"
    e = estimate("breast cancer", 90, pgs_id="PGS000004", eff=eff, base=base)
    b = e["baseline"]
    assert 0.15 < e["absolute_risk"] < 0.30, e            # ~21% at 90th pct
    assert estimate("breast cancer", 50, eff=eff, base=base)["absolute_risk"] < b  # median below mean
    assert estimate("breast cancer", 99, eff=eff, base=base)["absolute_risk"] > e["absolute_risk"]
    # CAD male, 95th pct
    c = estimate("coronary artery disease", 95, pgs_id="PGS000018", sex="male", eff=eff, base=base)
    assert c["sex"] == "male" and c["absolute_risk"] > c["baseline"]
    print(f"absolute_risk self-check OK: breast 90th-pct -> {e['absolute_risk']*100:.1f}% "
          f"(baseline {b*100:.1f}%, RR {e['relative_risk']:.2f}); "
          f"CAD 95th-pct male -> {c['absolute_risk']*100:.1f}% (baseline {c['baseline']*100:.0f}%)")
