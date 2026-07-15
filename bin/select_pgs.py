#!/usr/bin/env python3
"""Evidence-graded PGS selection from the PGS Catalog API.

For each launch-set trait: fetch candidate scores + metadata (GWAS sample size,
ancestry composition, evaluation ancestry, method, n_variants, publication),
compute an evidence grade, and pick the best-evidence score — the defensible
alternative to picking a PGS ID blindly (review §3 differentiator #2).

Writes results/pgs_catalog_meta.json (pgs_id -> metadata + grade + rationale),
which grade_pgs.py loads instead of its hardcoded table. curl/wget are blocked
on this box, so the API is hit via urllib inside this script.

usage: select_pgs.py [OUT.json]
"""
import sys, json, time, urllib.request, urllib.error

OUT = sys.argv[1] if len(sys.argv) > 1 else "results/pgs_catalog_meta.json"
TOP_N = int(sys.argv[2]) if len(sys.argv) > 2 else 2   # scores per trait (>=2 -> consensus)
BASE = "https://www.pgscatalog.org/rest"

# launch-set traits (label -> EFO/MONDO id used by the Catalog)
TIER = sys.argv[3] if len(sys.argv) > 3 else "extended"   # core | extended | all | gated

# Tiered launch set: (tier, trait_label, PGS Catalog trait_id). IDs + score counts
# resolved via bin/resolve_traits.py; evidence grade is the final data-sufficiency
# filter, so thin traits self-downgrade rather than being hard-excluded.
#   core     ~25  grade-A-capable, actionable (headline)
#   extended ~39  grade>=B, informative (default suite = core+extended)
#   gated     ~6  sensitive/low-actionability (opt-in only)
LAUNCH_SET = [
    # --- core ---
    ("core", "coronary artery disease", "MONDO_0005010"),
    ("core", "type 2 diabetes", "MONDO_0005148"),
    ("core", "atrial fibrillation", "MONDO_0004981"),
    ("core", "ischemic stroke", "HP_0002140"),
    ("core", "LDL cholesterol", "EFO_0004611"),
    ("core", "HDL cholesterol", "EFO_0004612"),
    ("core", "triglycerides", "EFO_0004530"),
    ("core", "lipoprotein(a)", "EFO_0006925"),
    ("core", "body mass index", "EFO_0004340"),
    ("core", "chronic kidney disease", "MONDO_0005300"),
    ("core", "heart failure", "MONDO_0005252"),
    ("core", "venous thromboembolism", "MONDO_0005399"),
    ("core", "hypertension", "MONDO_0001134"),
    ("core", "breast cancer", "MONDO_0004989"),
    ("core", "prostate cancer", "MONDO_0005159"),
    ("core", "colorectal cancer", "MONDO_0005575"),
    ("core", "melanoma", "MONDO_0005105"),
    ("core", "lung cancer", "MONDO_0005138"),
    ("core", "glaucoma", "MONDO_0005041"),
    ("core", "age-related macular degeneration", "MONDO_0005150"),
    ("core", "osteoporosis", "MONDO_0005298"),
    ("core", "gout", "MONDO_0005393"),
    ("core", "type 1 diabetes", "MONDO_0005147"),
    ("core", "coeliac disease", "MONDO_0005130"),
    ("core", "hypothyroidism", "MONDO_0005420"),
    # --- extended ---
    ("extended", "rheumatoid arthritis", "MONDO_0008383"),
    ("extended", "inflammatory bowel disease", "MONDO_0005265"),
    ("extended", "psoriasis", "MONDO_0005083"),
    ("extended", "multiple sclerosis", "MONDO_0005301"),
    ("extended", "systemic lupus erythematosus", "MONDO_0007915"),
    ("extended", "ankylosing spondylitis", "MONDO_0005306"),
    ("extended", "asthma", "MONDO_0004979"),
    ("extended", "chronic obstructive pulmonary disease", "MONDO_0005002"),
    ("extended", "Parkinson disease", "MONDO_0005180"),
    ("extended", "migraine", "MONDO_0005277"),
    ("extended", "Alzheimer disease", "MONDO_0004975"),
    ("extended", "atopic dermatitis", "MONDO_0004980"),
    ("extended", "ovarian cancer", "MONDO_0005140"),
    ("extended", "pancreatic cancer", "MONDO_0005192"),
    ("extended", "bladder cancer", "MONDO_0004986"),
    ("extended", "testicular cancer", "MONDO_0005447"),
    ("extended", "thyroid cancer", "MONDO_0015075"),
    ("extended", "endometrial cancer", "MONDO_0002447"),
    ("extended", "basal cell carcinoma", "MONDO_0020804"),
    ("extended", "osteoarthritis", "MONDO_0006629"),
    ("extended", "hyperthyroidism", "MONDO_0010138"),
    # PCOS dropped: no PGS in the Catalog under any current ontology id
    # (MONDO_0008559 / EFO_0000660 both return 0). Replaced with multiple myeloma
    # (grade A, ~1.05M-sample GWAS) — fills the hematologic-cancer gap.
    ("extended", "multiple myeloma", "MONDO_0009693"),
    ("extended", "cataract", "MONDO_0005129"),
    ("extended", "fatty liver disease", "MONDO_0013209"),
    ("extended", "gallstones", "MONDO_0012672"),
    ("extended", "abdominal aortic aneurysm", "MONDO_0005350"),
    ("extended", "C-reactive protein", "EFO_0004458"),
    ("extended", "HbA1c", "EFO_0004541"),
    ("extended", "height", "OBA_VT0001253"),
    ("extended", "waist-hip ratio", "EFO_0004343"),
    # --- next tranche (2026-07): well-powered disease/cancer + biomarker (see docs/candidate-traits.md) ---
    ("extended", "peripheral arterial disease", "MONDO_0005386"),
    ("extended", "endometriosis", "MONDO_0005133"),
    ("extended", "obstructive sleep apnea", "MONDO_0005296"),
    ("extended", "esophageal cancer", "MONDO_0007576"),
    ("extended", "gastric cancer", "MONDO_0001056"),
    ("extended", "cervical cancer", "MONDO_0005131"),
    ("extended", "bone mineral density", "EFO_0009270"),
    ("extended", "estimated glomerular filtration rate", "OBA_0003747"),
    ("extended", "systolic blood pressure", "EFO_0006335"),
    # --- OmniGen additions (2026-07): quantitative / behavioral traits OmniGen renders ---
    ("extended", "chronotype", "EFO_0008328"),          # benign "fun" trait; morningness/eveningness
    # --- gated (opt-in; sensitive / low-actionability) ---
    ("gated", "schizophrenia", "MONDO_0005090"),
    ("gated", "major depressive disorder", "MONDO_0002009"),
    ("gated", "bipolar disorder", "MONDO_0004985"),
    ("gated", "ADHD", "MONDO_0007743"),
    ("gated", "autism spectrum disorder", "MONDO_0005258"),
    ("gated", "intelligence", "EFO_0004337"),
    # behavioral / psych traits carry disclosure + misuse risk -> gated (opt-in) like intelligence
    ("gated", "educational attainment", "EFO_0011015"),
    ("gated", "neuroticism", "EFO_0007660"),
    ("gated", "loneliness", "EFO_0007865"),
]

# Pinned best-evidence scores for the OmniGen-additions traits. select_pgs() normally
# picks the top-graded Catalog score per trait, but for these we pin the specific,
# verified, deployable PGS ID (resolved live against the PGS Catalog REST API,
# 2026-07-11). Where the requested landmark study was NOT deposited with usable
# weights, we pin the documented genome-wide LDpred2 substitute and LABEL it.
# main() forces the pinned id to rank 1 and records the note in the metadata.
PINNED = {
    "height": dict(pgs_id="PGS002804", author="Yengo", year=2022, pmid="36224396",
                   n_variants=1099005, training_ancestry="European", base_grade="A",
                   multi_ancestry_eval=True, gwas_n=5380080, method="C+T (GIANT-UKB EUR)",
                   note="Yengo 2022 GIANT height, EUR weights (PGS002802 = multi-ancestry ALL)"),
    "bone mineral density": dict(pgs_id="PGS000657", author="Forgetta", year=2020,
                   pmid="32614825", n_variants=21716, training_ancestry="European",
                   base_grade="A", multi_ancestry_eval=False, gwas_n=426824,
                   method="ML gradient-boosted (gSOS)",
                   note="Forgetta 2020 gSOS heel-BMD proxy (Morris 2019 eBMD not deposited as a weighted score)"),
    "educational attainment": dict(pgs_id="PGS002231", author="Prive", year=2022,
                   pmid="34995502", n_variants=950845, training_ancestry="European",
                   base_grade="A", multi_ancestry_eval=False, gwas_n=391124,
                   method="LDpred2 (UKB)",
                   note="SUBSTITUTE: EA4/Okbay 2022 (PMID 35361970) not deposited with usable weights; "
                        "PGS002231 is the genome-wide UKB LDpred2 reconstruction"),
    "chronotype": dict(pgs_id="PGS002209", author="Prive", year=2022, pmid="34995502",
                   n_variants=955439, training_ancestry="European", base_grade="A",
                   multi_ancestry_eval=False, gwas_n=449734, method="LDpred2 (UKB)",
                   note="SUBSTITUTE: Jones 2019 (PMID 30696823) not deposited; PGS002209 scores "
                        "'more_evening' — higher percentile = MORE evening (INVERT for morningness)"),
    "neuroticism": dict(pgs_id="PGS002213", author="Prive", year=2022, pmid="34995502",
                   n_variants=950183, training_ancestry="European", base_grade="A",
                   multi_ancestry_eval=False, gwas_n=380506, method="LDpred2 (UKB)",
                   note="Prive 2022 UKB LDpred2 neuroticism"),
    "loneliness": dict(pgs_id="PGS001091", author="Tanigawa", year=2022, pmid="35324888",
                   n_variants=660, training_ancestry="European", base_grade="C",
                   multi_ancestry_eval=False, gwas_n=337000, method="Sparse GBE (snpnet)",
                   note="only Catalog loneliness score; thin (660 variants) -> expect grade C/D self-downgrade"),
}

# Traits requested for OmniGen but NOT addable via the Catalog: the landmark scores are
# 23andMe-held and not deposited with weights. Left as documented stubs — do NOT fabricate
# scorefiles. Add later via the --scorefile-cache path only if licensed weights are obtained.
DEFERRED = {
    "beat/rhythm synchronization": dict(author="Niarchou", year=2022, pmid="35618881",
        reason="23andMe-held; no weighted score in the PGS Catalog"),
    "general risk tolerance": dict(author="Karlsson Linner", year=2019, pmid="30643258",
        reason="23andMe-held; general-risk-tolerance weights not deposited in the PGS Catalog"),
}
_TIERS = {"core": {"core"}, "extended": {"core", "extended"},
          "all": {"core", "extended", "gated"}, "gated": {"gated"}}
TRAITS = {label: tid for tier, label, tid in LAUNCH_SET
          if tier in _TIERS.get(TIER, {"core", "extended"})}
EUR_KEYS = {"EUR", "European"}


def get(path):
    for attempt in range(3):
        try:
            req = urllib.request.Request(f"{BASE}/{path}",
                                         headers={"User-Agent": "poly-suite/0.1"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt == 2:
                raise
            time.sleep(2)


def scores_for_trait(efo):
    """All score objects for a trait (paginated)."""
    out, offset = [], 0
    while True:
        j = get(f"score/search?trait_id={efo}&limit=50&offset={offset}")
        res = j.get("results", [])
        out.extend(res)
        if not j.get("next"):
            break
        offset += 50
    return out


def gwas_n(s):
    return sum((x.get("sample_number") or 0) for x in (s.get("samples_variants") or []))


def ancestry_profile(s):
    """(training_ancestry_majority, has_multiancestry_eval)."""
    dist = s.get("ancestry_distribution") or {}
    train = dist.get("gwas") or dist.get("dev") or {}
    tdist = train.get("dist") or {}
    train_major = max(tdist, key=tdist.get) if tdist else "unknown"
    ev = (dist.get("eval") or {}).get("dist") or {}
    multi_eval = len([k for k in ev if k]) >= 2 or any(k not in EUR_KEYS for k in ev)
    return train_major, multi_eval, bool(ev)


def grade(s):
    n = gwas_n(s)
    train_major, multi_eval, has_eval = ancestry_profile(s)
    if n >= 100_000 and has_eval and (multi_eval or n >= 250_000):
        g = "A"
    elif n >= 50_000 and has_eval:
        g = "B"
    elif n >= 10_000:
        g = "C"
    else:
        g = "D"
    rationale = (f"GWAS N~{n:,}; train~{train_major}; "
                 f"eval={'multi-ancestry' if multi_eval else ('present' if has_eval else 'none')}")
    return g, rationale, n, train_major, multi_eval


GRADE_RANK = {"A": 0, "B": 1, "C": 2, "D": 3}


def main():
    meta, summary = {}, []
    for label, efo in TRAITS.items():
        try:
            cand = scores_for_trait(efo)
        except Exception as e:
            summary.append((label, None, f"API error: {e}"))
            continue
        if not cand:
            summary.append((label, None, "no scores in Catalog for this trait_id"))
            continue
        graded = []
        for s in cand:
            g, rat, n, train, multi = grade(s)
            graded.append((s, g, rat, n, train, multi))
        # best: grade, then GWAS N, then variant count
        graded.sort(key=lambda t: (GRADE_RANK[t[1]], -t[3],
                                    -(t[0].get("variants_number") or 0)))
        # pinned verified/substitute score -> force to rank 1 (see PINNED above)
        pin = PINNED.get(label)
        if pin:
            pid = pin["pgs_id"]
            idx = next((i for i, t in enumerate(graded) if t[0].get("id") == pid), None)
            if idx is not None:
                graded.insert(0, graded.pop(idx))          # promote the real Catalog record
            else:                                          # not returned (offline/paginated) -> synthesize
                synth = ({"id": pid, "method_name": pin.get("method"),
                          "variants_number": pin.get("n_variants"),
                          "publication": {"PMID": pin.get("pmid"), "doi": pin.get("doi"),
                                          "firstauthor": pin.get("author"),
                                          "pub_year": pin.get("year")}},
                         pin.get("base_grade", "B"), "pinned (static metadata)",
                         pin.get("gwas_n", 0), pin.get("training_ancestry", "European"),
                         pin.get("multi_ancestry_eval", False))
                graded.insert(0, synth)
        picks = []
        for rank, chosen in enumerate(graded[:TOP_N], 1):    # top-N per trait -> consensus
            s = chosen[0]
            pgs_id = s.get("id")
            pub = s.get("publication") or {}
            meta[pgs_id] = {
                "trait": label, "efo_id": efo,
                "n_variants": s.get("variants_number"),
                "method": s.get("method_name"),
                "training_ancestry": chosen[4],
                "multi_ancestry_eval": chosen[5],
                "gwas_n": chosen[3],
                # PMID may be null for preprints (e.g. medRxiv/bioRxiv scores
                # not yet in PubMed). Keep the DOI + journal so downstream can
                # still cite honest provenance instead of an empty cell.
                "pmid": pub.get("PMID"), "doi": pub.get("doi"),
                "journal": pub.get("journal"),
                "author": pub.get("firstauthor"),
                "year": pub.get("pub_year"),
                "base_grade": chosen[1], "rationale": chosen[2],
                "trait_rank": rank, "n_candidates": len(cand),
            }
            if pin and rank == 1 and pgs_id == pin["pgs_id"]:
                meta[pgs_id]["pinned"] = True
                meta[pgs_id]["pin_note"] = pin["note"]
                meta[pgs_id]["rationale"] = f"{chosen[2]} | PINNED: {pin['note']}"
            picks.append((pgs_id, chosen[1]))
        summary.append((label, ",".join(p for p, _ in picks),
                        f"top-{len(picks)} of {len(cand)} | grades {','.join(g for _, g in picks)}"))
        time.sleep(0.3)  # be polite to the API

    with open(OUT, "w") as fh:
        json.dump(meta, fh, indent=2)

    print(f"\n=== evidence-graded PGS selection ===  -> {OUT}\n")
    for label, pgs_id, note in summary:
        pid = pgs_id or "—"
        print(f"  {label:<26} {pid:<12} {note}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
