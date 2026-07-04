#!/usr/bin/env python3
"""Resolve a curated trait list to PGS Catalog trait_ids + score counts (the
'enough data' check), for building the tiered launch set. Prints TRAITS-dict-ready
lines; scores<MIN flagged so thin traits move to the exploratory/opt-in tier.

Tiers: core (grade-A actionable), extended (grade>=B informative), gated (opt-in,
sensitive). The evidence grade in select_pgs is the final data-sufficiency filter;
this just resolves ids + surfaces availability.
"""
import json, urllib.request, urllib.parse, time, sys

B = "https://www.pgscatalog.org/rest"

CANDIDATES = {
"core": [
    "coronary artery disease","type 2 diabetes","atrial fibrillation","ischemic stroke",
    "LDL cholesterol","HDL cholesterol","triglycerides","lipoprotein A","body mass index",
    "chronic kidney disease","heart failure","venous thromboembolism","hypertension",
    "breast cancer","prostate cancer","colorectal cancer","melanoma","lung carcinoma",
    "glaucoma","age-related macular degeneration","osteoporosis","gout","type 1 diabetes",
    "coeliac disease","hypothyroidism",
],
"extended": [
    "rheumatoid arthritis","Crohn disease","ulcerative colitis","psoriasis","multiple sclerosis",
    "systemic lupus erythematosus","ankylosing spondylitis","asthma","chronic obstructive pulmonary disease",
    "Parkinson disease","migraine disorder","epilepsy","Alzheimer disease","atopic dermatitis",
    "ovarian cancer","pancreatic cancer","renal cell carcinoma","bladder carcinoma","testicular cancer",
    "thyroid carcinoma","endometrial cancer","esophageal carcinoma","basal cell carcinoma",
    "osteoarthritis","hyperthyroidism","polycystic ovary syndrome","cataract","myopia",
    "non-alcoholic fatty liver disease","gallstones","irritable bowel syndrome","kidney stone",
    "benign prostatic hyperplasia","hemoglobin A1c measurement","fasting blood glucose measurement",
    "peripheral arterial disease","abdominal aortic aneurysm","C-reactive protein measurement",
    "estimated glomerular filtration rate","height","waist-hip ratio",
],
"gated": [
    "schizophrenia","bipolar disorder","major depressive disorder","attention deficit hyperactivity disorder",
    "autism spectrum disorder","alcohol dependence","cognitive ability",
],
}
MIN_SCORES = 3


def g(p):
    try:
        with urllib.request.urlopen(f"{B}/{p}", timeout=60) as r:
            return json.load(r)
    except Exception as e:
        return {"_err": str(e)}


def best_id(term):
    j = g(f"trait/search?term={urllib.parse.quote(term)}")
    res = j.get("results") or []
    if not res:
        return None, 0, ""
    best = max(res, key=lambda t: len(t.get("associated_pgs_ids") or []))
    return best.get("id"), len(best.get("associated_pgs_ids") or []), best.get("label", "")


def main():
    total = thin = 0
    for tier, terms in CANDIDATES.items():
        print(f"\n# --- tier: {tier} ---")
        for term in terms:
            tid, n, label = best_id(term)
            total += 1
            flag = "" if n >= MIN_SCORES else "   # THIN (opt-in)"
            if n < MIN_SCORES:
                thin += 1
            key = label or term
            print(f'    ("{tier}", "{key}", "{tid}"),   # {n} scores{flag}'
                  if tid else f'    # ("{tier}", "{term}", None),  # no match')
            time.sleep(0.25)
    sys.stderr.write(f"\n{total} traits resolved; {thin} thin (<{MIN_SCORES} scores)\n")


if __name__ == "__main__":
    main()
