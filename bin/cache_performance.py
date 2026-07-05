#!/usr/bin/env python3
"""Cache each score's PUBLICATION-REPORTED discrimination metric from the PGS
Catalog performance-metrics API, so the report can show predictive performance
(AUROC / C-index / R2) next to poly-suite's evidence grade.

A score has many performance records (different cohorts/ancestries). We pick ONE
representative family — AUROC, else C-index, else R2 — preferring European-ancestry
records (our calibration ancestry), and report the MEDIAN with the range across
reports (honest about the spread; these are reported, not computed here).

usage: cache_performance.py CONTRACT.tsv [OUT.tsv]
   OUT defaults to resources/pgs_performance.tsv (merged/updated, keyed by pgs_id).
Only prints a short summary — never floods context.
"""
import os, sys, csv, json, re, time, urllib.request
from statistics import median

UA = {"User-Agent": "poly-suite/0.1"}
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT = os.path.join(os.path.dirname(_HERE), "resources", "pgs_performance.tsv")
COLS = ["pgs_id", "metric", "value", "n_reports", "lo", "hi"]

# metric family -> (label, matcher on a metric record's name/name_short)
AUROC = re.compile(r"auroc|auc|area under", re.I)
CINDEX = re.compile(r"c-?index|harrell|concordance", re.I)
R2 = re.compile(r"\br2\b|r²|variance explained|r-squared", re.I)


def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
        return json.load(r)


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _collect(records, pat):
    """(values, values_eur) for metric records matching pat across performance records."""
    allv, eurv = [], []
    for rec in records:
        ss = rec.get("sampleset") or {}
        anc = " ".join((s.get("ancestry_broad") or "") for s in (ss.get("samples") or []))
        is_eur = "european" in anc.lower()
        pm = rec.get("performance_metrics") or {}
        for m in (pm.get("class_acc") or []) + (pm.get("othermetrics") or []):
            name = f"{m.get('name_short') or ''} {m.get('name') or ''}"
            if pat.search(name):
                est = _num(m.get("estimate"))
                if est is not None:
                    allv.append(est)
                    if is_eur:
                        eurv.append(est)
    return allv, eurv


# plausible ranges — the Catalog's othermetrics labels are noisy (R2 sometimes
# holds a percentage or liability-scale value >1, AUROC sometimes < .5), so drop
# anything physically out of range rather than rescale/guess.
BOUNDS = {"AUROC": (0.5, 1.0), "C-index": (0.5, 1.0), "R2": (0.0, 1.0)}


def representative(records):
    """-> dict(metric,value,n_reports,lo,hi) or None. Prefer AUROC>C-index>R2,
    EUR records if any report the chosen family, in-range values only."""
    for label, pat in (("AUROC", AUROC), ("C-index", CINDEX), ("R2", R2)):
        lo_b, hi_b = BOUNDS[label]
        allv, eurv = _collect(records, pat)
        allv = [v for v in allv if lo_b <= v <= hi_b]
        eurv = [v for v in eurv if lo_b <= v <= hi_b]
        vals = eurv or allv
        if vals:
            return dict(metric=label, value=round(median(vals), 3), n_reports=len(vals),
                        lo=round(min(vals), 3), hi=round(max(vals), 3))
    return None


def load(path):
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return {r["pgs_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    contract, out = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else DEFAULT_OUT)
    ids = sorted({re.match(r"(PGS\d+)", r["pgs_id"]).group(1)
                  for r in csv.DictReader(open(contract), delimiter="\t")
                  if re.match(r"(PGS\d+)", r.get("pgs_id", ""))})
    cache = load(out)
    got = miss = 0
    for pid in ids:
        if pid in cache and cache[pid].get("value") not in (None, "", "NA"):
            continue  # already cached
        try:
            j = get(f"https://www.pgscatalog.org/rest/performance/search?pgs_id={pid}")
        except Exception as e:
            print(f"  {pid}: API error {e}", file=sys.stderr)
            continue
        rep = representative(j.get("results") or [])
        if rep:
            cache[pid] = dict(pgs_id=pid, **rep)
            got += 1
        else:
            cache[pid] = dict(pgs_id=pid, metric="NA", value="NA",
                              n_reports=0, lo="NA", hi="NA")
            miss += 1
        time.sleep(0.2)  # be polite
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        for pid in sorted(cache):
            w.writerow({k: cache[pid].get(k, "") for k in COLS})
    print(f"performance cache -> {out}  ({len(cache)} scores; +{got} fetched, "
          f"{miss} with no reported metric)")


if __name__ == "__main__":
    main()
