#!/usr/bin/env python3
"""poly-suite grading + report stage (standalone).

Reads the pgsc_calc scoring-core output (`results/score/aggregated_scores.txt.gz`,
long format: sampleset,IID,PGS,SUM,DENOM,AVG[,percentile_MostSimilarPop,
Z_MostSimilarPop,MostSimilarPop]) and emits the STABLE contract OmniGen consumes:

  results/pgs_scores.tsv   (+ .json)   one row per (sample, trait)
  and a human-readable card to stdout.

This is poly-suite's OWN honesty layer — QC gate + ancestry/portability gate +
evidence grade + controlled-vocab caveat (pgs-pipeline-spec.md §2/§6). The scoring
math is pgsc_calc's; nothing is re-scored here. Absolute-risk / multi-score /
CI columns are contract placeholders (NA) until those [expand] stages exist.

Usage: grade_pgs.py [results_dir]   (default: ./results)
"""
import sys, os, glob, gzip, json, csv, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import absolute_risk, provenance, report_html, consensus

_EFF, _BASE = ({}, {}), {}   # per-SD effects, baseline incidence (loaded in main)
_SEX = None                  # sample sex (for sex-dimorphic absolute risk, e.g. CAD)


def sample_sex(results_dir):
    """Sample sex from env POLY_SUITE_SEX or {results_dir}/sample_sex.txt
    (written by infer_sex.py). None -> sex-dimorphic traits report percentile only."""
    s = os.environ.get("POLY_SUITE_SEX")
    if not s:
        p = os.path.join(results_dir, "sample_sex.txt")
        if os.path.exists(p):
            s = open(p).read().strip()
    s = (s or "").lower()
    return s if s in ("male", "female") else None

RESULTS = sys.argv[1] if len(sys.argv) > 1 else "results"
MATCH_RATE_MIN = 0.75          # pgsc_calc's own default overlap floor
GRADES = ["A", "B", "C", "D"]

# Starter-set PGS metadata. Ideally read from PGS Catalog API at run time (the
# evidence-selection [expand] stage); pinned here for the launch set.
PGS_META = {   # pgs_id: dict
    "PGS000018": dict(trait="coronary artery disease", efo="EFO_0001645",
                      n_variants=1745179, base_grade="A", training_ancestry="European",
                      pmid="30309464", source="Inouye 2018 metaGRS_CAD"),
    "PGS000004": dict(trait="breast cancer", efo="EFO_0000305",
                      n_variants=313, base_grade="A", training_ancestry="European",
                      pmid="30554720", source="Mavaddat 2019 BCAC 313-variant"),
    "PGS000662": dict(trait="prostate cancer", efo="EFO_0001663",
                      n_variants=269, base_grade="A", training_ancestry="Multi-ancestry",
                      pmid="33398198", source="Conti 2021 multi-ancestry"),
    "PGS000014": dict(trait="type 2 diabetes", efo="EFO_0001360",
                      n_variants=6917436, base_grade="B", training_ancestry="European",
                      pmid="30104762", source="Khera 2018 GPS_T2D"),
}

CONTRACT_COLS = [
    "sample", "trait", "efo_id", "pgs_id", "source_pmid", "source_doi", "n_variants",
    "training_ancestry", "n_matched", "match_rate", "inferred_ancestry",
    "most_similar_pop", "ancestry_distance", "percentile", "z_score",
    "ci_low", "ci_high", "absolute_risk", "baseline_incidence", "risk_ratio",
    "robustness_n_scores", "robustness_concordance", "evidence_grade",
    "portability_flag", "allowed_statement",
]


def downgrade(g, n=1):
    return GRADES[min(len(GRADES) - 1, GRADES.index(g) + n)] if g in GRADES else "D"


def find_scores(results_dir):
    for pat in (f"{results_dir}/score/aggregated_scores.txt.gz",
                f"{results_dir}/**/aggregated_scores*.txt.gz",
                f"{results_dir}/**/*_pgs.txt.gz",        # --run_ancestry calibrated output
                f"{results_dir}/**/*scores*.txt*"):
        hits = sorted(glob.glob(pat, recursive=True))
        if hits:
            return hits[0]
    return None


def load_pop(results_dir):
    """Most-similar population from pgsc_calc's pop_summary.csv (--run_ancestry).
    The sample's column is non-NA for exactly its assigned population."""
    for p in glob.glob(f"{results_dir}/**/pop_summary.csv", recursive=True):
        with open(p) as fh:
            r = csv.reader(fh)
            header = next(r, None)          # Most similar population, <sample>, reference
            for row in r:
                if len(row) >= 2 and row[1] and row[1].upper() != "NA":
                    return row[0]           # e.g. 'EUR'
    return None


def read_rows(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        # normalize header case for robust lookup
        return [{k.lower(): v for k, v in r.items()} for r in rdr]


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def norm_pgs(acc):
    """PGS000018_hmPOS_GRCh38 -> PGS000018 (pgsc_calc harmonization suffix)."""
    m = re.match(r"(PGS\d+)", str(acc))
    return m.group(1) if m else str(acc)


def parse_match_summary(results_dir):
    """Authoritative per-score coverage from pgsc_calc's match summary.
    The CSV is ragged (matched rows have 12 cols, unmatched/short rows ~6), so
    parse positionally: accession=col1, match_status=col3, count=col[-2].
    Returns {pgs_id: {"matched": n, "total": n, "match_rate": f}}."""
    hits = sorted(glob.glob(f"{results_dir}/**/match/*_summary.csv", recursive=True))
    if not hits:
        return {}
    agg = {}
    with open(hits[0]) as fh:
        rdr = csv.reader(fh)
        header = next(rdr, None)
        for row in rdr:
            if len(row) < 5:
                continue
            pgs = norm_pgs(row[1])
            status, count = row[3], num(row[-2])
            if count is None:
                continue
            a = agg.setdefault(pgs, {"matched": 0, "total": 0})
            a["total"] += int(count)
            if status == "matched":
                a["matched"] += int(count)
    for a in agg.values():
        a["match_rate"] = a["matched"] / a["total"] if a["total"] else None
    return agg


def grade_row(r, cov=None):
    pgs = norm_pgs(r.get("pgs") or r.get("pgs_id") or "?")
    m = PGS_META.get(pgs, dict(trait=pgs, efo="", n_variants=None,
                               base_grade="C", training_ancestry="unknown",
                               pmid="", doi="", source="unknown PGS"))
    c = (cov or {}).get(pgs)
    if c:  # authoritative coverage from the match summary
        n_matched = c["matched"]
        match_rate = c["match_rate"]
        nvar = m["n_variants"] or c["total"]
    else:  # fallback: derive from aggregated-scores DENOM
        denom = num(r.get("denom"))
        n_matched = int(denom // 2) if denom is not None else None
        nvar = m["n_variants"]
        match_rate = (n_matched / nvar) if (n_matched is not None and nvar) else None

    pct = num(r.get("percentile_mostsimilarpop") or r.get("percentile"))
    z = num(r.get("z_mostsimilarpop") or r.get("z_score"))
    pop = r.get("mostsimilarpop") or r.get("most_similar_pop") or None

    grade = m["base_grade"]
    notes, portability = [], False

    # QC gate — coverage / overlap
    if match_rate is not None and match_rate < MATCH_RATE_MIN:
        grade = downgrade(grade, 2)
        notes.append(f"low variant coverage ({match_rate:.0%} of {nvar:,} scoring "
                     "variants matched) — input is not full-genotype; raw score is "
                     "biased and NOT a valid PRS")

    # Calibration / ancestry gate
    if pct is None:
        grade = "D"
        notes.append("UNCALIBRATED — raw sum only; a percentile requires the "
                     "ancestry step (pgsc_calc --run_ancestry). No population "
                     "reference was applied")
    else:
        # portability: EUR-trained score read outside European ancestry
        if pop and pop.upper() != "EUR" and m["training_ancestry"] == "European":
            portability = True
            grade = downgrade(grade)
            notes.append(f"portability: score trained in {m['training_ancestry']} "
                         f"ancestry, sample most-similar to {pop} — reduced accuracy")

    # absolute-risk conversion (odds-scale), where a percentile + effect + baseline exist
    ar = (absolute_risk.estimate(m["trait"], pct, pgs_id=pgs,
                                 ancestry=(pop or "overall"), sex=_SEX, eff=_EFF, base=_BASE)
          if pct is not None else None)
    if ar:
        notes.append(f"absolute risk ~{ar['absolute_risk']*100:.0f}% for a {ar['sex']} "
                     f"vs ~{ar['baseline']*100:.0f}% population baseline "
                     f"(odds-scale, per-SD effect {ar['per_sd']}, illustrative — verify per release)")

    notes.append("explains only a fraction of trait variance; reduced accuracy "
                 "outside training ancestry; screening context only, not diagnostic")

    return {
        "sample": r.get("iid") or r.get("sample") or r.get("sampleset"),
        "trait": m["trait"], "efo_id": m["efo"], "pgs_id": pgs,
        # Provenance semantics (see docs/CHANGES.md):
        #   source_pmid non-empty            -> peer-reviewed pub, PubMed ID given
        #   source_pmid "" + source_doi set  -> provenance IS known but the source
        #                                       is a preprint/DOI-only record with
        #                                       no PMID upstream (NOT a fetch failure)
        #   both ""                          -> provenance could not be resolved
        "source_pmid": m["pmid"], "source_doi": m.get("doi", ""), "n_variants": nvar,
        "training_ancestry": m["training_ancestry"],
        "n_matched": n_matched,
        "match_rate": round(match_rate, 4) if match_rate is not None else "NA",
        "inferred_ancestry": pop or "NA", "most_similar_pop": pop or "NA",
        "ancestry_distance": "NA",
        "percentile": round(pct, 1) if pct is not None else "NA",
        "z_score": round(z, 4) if z is not None else "NA",
        "ci_low": "NA", "ci_high": "NA",
        "absolute_risk": round(ar["absolute_risk"], 4) if ar else "NA",
        "baseline_incidence": round(ar["baseline"], 4) if ar else "NA",
        "risk_ratio": round(ar["relative_risk"], 3) if ar else "NA",
        "robustness_n_scores": 1, "robustness_concordance": "NA",
        "evidence_grade": grade,
        "portability_flag": str(portability).lower(),
        "allowed_statement": " | ".join(notes),
        "_raw_sum": num(r.get("sum")), "_source": m["source"],
    }


def load_catalog_meta(results_dir):
    """Merge live evidence-graded metadata from select_pgs.py (if present) over
    the hardcoded launch-set table, so any scored PGS the selector covered gets
    its Catalog-derived grade/ancestry instead of a pinned guess."""
    for p in (os.path.join(results_dir, "pgs_catalog_meta.json"),
              "results/pgs_catalog_meta.json"):
        if os.path.exists(p):
            with open(p) as fh:
                raw = json.load(fh)
            for pid, m in raw.items():
                PGS_META[pid] = dict(
                    trait=m.get("trait", pid), efo=m.get("efo_id", ""),
                    n_variants=m.get("n_variants"),
                    base_grade=m.get("base_grade", "C"),
                    training_ancestry=m.get("training_ancestry", "unknown"),
                    pmid=str(m.get("pmid") or ""),
                    doi=str(m.get("doi") or ""),
                    source=f"{m.get('author','?')} {m.get('year','')}".strip())
            return len(raw)
    return 0


def main():
    global _EFF, _BASE, _SEX
    _EFF, _BASE = absolute_risk.load_effects(), absolute_risk.load_baselines()
    _SEX = sample_sex(RESULTS)
    load_catalog_meta(RESULTS)
    path = find_scores(RESULTS)
    if not path:
        print(f"[poly-suite] no scoring output under {RESULTS}/ — run pgsc_calc first")
        return 1
    cov_map = parse_match_summary(RESULTS)
    pop = load_pop(RESULTS)
    raw = [r for r in read_rows(path)
           if (r.get("sampleset") or "").strip().lower() != "reference"]  # drop panel rows
    if pop:
        for r in raw:
            r.setdefault("mostsimilarpop", pop)     # --run_ancestry: pop is in pop_summary.csv
    rows = [grade_row(r, cov_map) for r in raw]
    if not rows:
        print(f"[poly-suite] {path} has no score rows")
        return 1

    # multi-score consensus: robustness across >=2 scores per trait
    cons = consensus.consensus([
        {"trait": r["trait"], "pgs_id": r["pgs_id"],
         "percentile": (r["percentile"] if r["percentile"] != "NA" else None)}
        for r in rows])
    for r in rows:
        c = cons.get(r["trait"])
        if not c:
            continue
        r["robustness_n_scores"] = c["n_scores"]
        r["robustness_concordance"] = c["concordance"] if c["concordance"] is not None else "NA"
        if (c["n_scores"] >= 2 and c["concordance"] is not None
                and (not c["tertile_concordant"] or c["concordance"] < 0.67)):
            r["allowed_statement"] += (
                f" | multi-score robustness LOW: {c['n_scores']} scores disagree "
                f"(percentile spread {c['spread']:.0f}, concordance {c['concordance']}"
                f"{', different risk tertiles' if not c['tertile_concordant'] else ''}) "
                f"— treat as uncertain")

    out_tsv = os.path.join(RESULTS, "pgs_scores.tsv")
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CONTRACT_COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(RESULTS, "pgs_scores.json"), "w") as fh:
        json.dump([{k: v for k, v in r.items() if not k.startswith("_")}
                   for r in rows], fh, indent=2)

    # standalone deliverable: provenance record + self-contained HTML report
    sample = rows[0].get("sample") if rows else None
    provenance.write(RESULTS, sample=sample, sample_sex=_SEX,
                     ref_fasta="/data/alvin/ref/GRCh38/hg38.canonical.fa")
    report_html.render(RESULTS)

    # human card
    print(f"\n=== poly-suite report card ===  (source: {os.path.relpath(path)})")
    print(f"contract → {out_tsv}\n")
    for r in rows:
        pctstr = (f"percentile {r['percentile']}" if r["percentile"] != "NA"
                  else f"raw sum {r['_raw_sum']:.4g} (no percentile)")
        covstr = (f"{r['match_rate']:.0%}" if r["match_rate"] != "NA" else "NA")
        nvarstr = f"{r['n_variants']:,}" if r["n_variants"] else "?"
        print(f"  [{r['pgs_id']}  grade {r['evidence_grade']}]  {r['trait']}")
        print(f"      {pctstr}   coverage {covstr} of {nvarstr} variants"
              f"   pop {r['most_similar_pop']}")
        print(f"      {r['_source']}  (PMID {r['source_pmid']}, "
              f"trained {r['training_ancestry']})")
        print(f"      caveat: {r['allowed_statement']}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
