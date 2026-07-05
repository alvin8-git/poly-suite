#!/usr/bin/env python3
"""Render the poly-suite graded contract to a self-contained HTML report
(spec §2 standalone deliverable). No external dependencies — inline CSS + an
embedded heading font. The only JS is two inline onclick handlers for
expand/collapse-all; every row is a native <details> that works without it.

Written for a NON-scientific reader, tuned for scanning 60+ traits:
 - a top summary of the few high-confidence elevated findings,
 - traits GROUPED BY EVIDENCE CONFIDENCE (High -> Insufficient), so a
   high-percentile-but-low-grade score (e.g. grade-D hyperthyroidism) sinks
   to the bottom instead of leading on a number you shouldn't trust,
 - each trait is a compact 2-line row: [organ icon] name + likelihood bar,
   then a one-line plain-language interpretation. A disclosure triangle
   reveals the clinical layer (risk translation, what it means, per-score
   table, coverage) — collapsed by default so the whole report stays scannable.

Importable (render) and runnable: python3 bin/report_html.py [results_dir]
"""
import os, sys, csv, json, html
from collections import OrderedDict, Counter
import absolute_risk

SEX_WORD = {"female": "women", "male": "men"}

GRADE_COLOR = {"A": "#1a7f37", "B": "#9a6700", "C": "#bc4c00", "D": "#82071e"}
GRADE_RANK = {"A": 4, "B": 3, "C": 2, "D": 1}
# grade -> (plain confidence label, css class)
CONFIDENCE = {"A": ("High confidence", "hi"), "B": ("Good confidence", "ok"),
              "C": ("Limited confidence", "low"), "D": ("Insufficient evidence", "none")}
TIERS = [("A", "High confidence", "Large, replicated, ancestry-checked scores."),
         ("B", "Good confidence", "Solid scores; treat as informative."),
         ("C", "Limited confidence", "Interpret with caution."),
         ("D", "Insufficient evidence", "Shown for completeness — not reliable enough to interpret.")]
ELEV, PROT, AVGC, WEAKHI, WEAK, UNCAL = (
    "elevated", "protective", "average", "weak-high", "weak", "uncalibrated")

# --- organ / body-system icon per trait (ponytail: keyword heuristic; inline
# SVG glyphs so the file stays one offline asset AND renders identically on
# every OS, themeable via currentColor). chip tint per system. -----------------
SYS = {
    "neuro": "#6d5ae0", "cardio": "#d6455c", "resp": "#4aa3c7", "msk": "#b07a3c",
    "metab": "#c68a2e", "endo": "#2fa6a0", "gut": "#4c8c5a", "eye": "#3c6db0",
    "renal": "#8a5a3c", "skin": "#c0745a", "cancer": "#7a4fb0", "anthro": "#6b7a8c",
    "default": "#64748b",
}
# 24x24 FILLED, cartoon-ish organ glyphs (fill:currentColor + white highlights)
# so they read as little illustrations, not line icons. The gut (tube) and DNA
# fall back to thick strokes where a silhouette wouldn't read.
SVG = {
    "cardio": '<path d="M12 20.6 4.4 13C2.5 11 2.5 8 4.4 6.1a4.4 4.4 0 0 1 6.3 0l1.3 1.3 1.3-1.3a4.4 4.4 0 0 1 6.3 0c1.9 1.9 1.9 4.9 0 6.9L12 20.6z"/><path fill="none" stroke="#fff" stroke-opacity=".45" stroke-width="1.1" stroke-linecap="round" d="M6.9 7.6a2.4 2.4 0 0 0-.5 2.7"/>',
    "neuro": '<path d="M9.2 4.3A2.8 2.8 0 0 0 6.4 6 2.5 2.5 0 0 0 4.9 8.4 2.6 2.6 0 0 0 4.7 12a2.6 2.6 0 0 0 .7 3.5A2.5 2.5 0 0 0 7.3 18.7 2.7 2.7 0 0 0 9.9 20c.5.2 1.3.3 2.1.3s1.6-.1 2.1-.3a2.7 2.7 0 0 0 2.6-1.3 2.5 2.5 0 0 0 1.9-2.7 2.6 2.6 0 0 0 .7-3.5 2.6 2.6 0 0 0-.2-3.6A2.5 2.5 0 0 0 17.6 6 2.8 2.8 0 0 0 14.8 4.3 2.9 2.9 0 0 0 12 3.9a2.9 2.9 0 0 0-2.8.4z"/><path fill="none" stroke="#fff" stroke-opacity=".55" stroke-width="1" stroke-linecap="round" d="M12 4.6v14.8M9 8.2c-1 .4-1.1 2 0 2.6M15 8.2c1 .4 1.1 2 0 2.6M9.6 13c-1 .5-1 2.1 0 2.7M14.4 13c1 .5 1 2.1 0 2.7"/>',
    "resp": '<path d="M11.25 4.2a.75.75 0 0 1 1.5 0v5.3c.7-.5 1.6-.6 2.4-.2 1.7.9 2.8 3.4 3.05 5.8.2 1.9-.3 3.6-2 4.1-1.8.5-3.9-.6-4.7-2.5-.15-.4-.25-.8-.25-1.2V4.2z"/><path d="M12.75 4.2a.75.75 0 0 0-1.5 0v5.3c-.7-.5-1.6-.6-2.4-.2-1.7.9-2.8 3.4-3.05 5.8-.2 1.9.3 3.6 2 4.1 1.8.5 3.9-.6 4.7-2.5.15-.4.25-.8.25-1.2V4.2z"/>',
    "msk": '<path d="M8.2 4.6a2 2 0 0 1 2.6 3l5.6 5.6a2 2 0 1 1 .6 3.4 2 2 0 1 1-3.4.6l-5.6-5.6a2 2 0 0 1-3-.6 2 2 0 0 1 .6-2.6 2 2 0 0 1 .6-3.4 2 2 0 0 1 2 .6z"/>',
    "metab": '<path d="M12 3.4c-.3 0-.5.1-.7.3C10 5.3 6 10 6 13.2A6 6 0 0 0 18 13.2c0-3.2-4-7.9-5.3-9.5-.2-.2-.4-.3-.7-.3z"/><path fill="none" stroke="#fff" stroke-opacity=".45" stroke-width="1.1" stroke-linecap="round" d="M9.3 13.5a2.7 2.7 0 0 0 2.5 2.7"/>',
    "endo": '<path d="M11.3 9.3C10.4 7.3 8.6 6 6.8 6.3 5.1 6.6 4.2 8.2 4.5 10.1c.3 2.2 2.1 4.1 4.3 4.1 1.6 0 2.6-1.3 3-3.2v3.6c0 .3.4.5.9.5s.9-.2.9-.5v-3.6c.4 1.9 1.4 3.2 3 3.2 2.2 0 4-1.9 4.3-4.1.3-1.9-.6-3.5-2.3-3.8-1.8-.3-3.6 1-4.5 3z"/><path fill="none" stroke="#fff" stroke-opacity=".5" stroke-width=".9" stroke-linecap="round" d="M7.6 8.6c-.5.6-.7 1.4-.5 2.2M16.4 8.6c.5.6.7 1.4.5 2.2"/>',
    "gut": '<path fill="none" stroke="currentColor" stroke-width="3.4" stroke-linecap="round" stroke-linejoin="round" d="M7 18V9.5a3 3 0 0 1 6 0v5a3 3 0 0 0 6 0V8"/>',
    "eye": '<path d="M12 5.6C6.6 5.6 3.2 11.4 3 11.8a.5.5 0 0 0 0 .4C3.2 12.6 6.6 18.4 12 18.4s8.8-5.8 9-6.2a.5.5 0 0 0 0-.4C20.8 11.4 17.4 5.6 12 5.6z"/><circle cx="12" cy="12" r="3.3" fill="#fff"/><circle cx="12" cy="12" r="1.6" fill="currentColor"/><circle cx="13.2" cy="10.9" r=".6" fill="#fff"/>',
    "renal": '<path d="M14.5 5C10.4 5 7 8.4 7 12.5S10 20 13 20c1.9 0 3.1-1.2 3.1-3 0-1.4-.9-2.1-.9-3.5S16.4 11 17.5 10s2-2 2-3.4C19.5 5.5 16.6 5 14.5 5z"/>',
    "skin": '<path d="M8.5 3.8 3.8 8.5a3 3 0 0 0 0 4.2l6.5 6.5a3 3 0 0 0 4.2 0l4.7-4.7a3 3 0 0 0 0-4.2l-6.5-6.5a3 3 0 0 0-4.2 0z"/><g fill="#fff" fill-opacity=".55"><circle cx="10.2" cy="12" r=".7"/><circle cx="12" cy="10.2" r=".7"/><circle cx="12" cy="13.8" r=".7"/><circle cx="13.8" cy="12" r=".7"/></g>',
    "cancer": '<path fill-rule="evenodd" d="M12 3.5a3.6 3.6 0 0 0-3 5.6l1.9 2.8-3.6 6.8a.6.6 0 0 0 .8.8l3.9-2.1 3.9 2.1a.6.6 0 0 0 .8-.8l-3.6-6.8 1.9-2.8A3.6 3.6 0 0 0 12 3.5zm0 2a1.6 1.6 0 1 1 0 3.2 1.6 1.6 0 0 1 0-3.2z"/>',
    "anthro": '<path d="M3.8 14.4 14.4 3.8a1 1 0 0 1 1.4 0l4.4 4.4a1 1 0 0 1 0 1.4L9.6 20.2a1 1 0 0 1-1.4 0l-4.4-4.4a1 1 0 0 1 0-1.4z"/><path fill="none" stroke="#fff" stroke-opacity=".6" stroke-width="1" stroke-linecap="round" d="M7.3 11.1 8.7 12.5M9.8 8.6l1.4 1.4M12.3 6.1l1.4 1.4"/>',
    "default": '<path fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" d="M8 4c0 4 8 4 8 8s-8 4-8 8M16 4c0 4-8 4-8 8s8 4 8 8"/><path fill="none" stroke="currentColor" stroke-width="1.3" stroke-linecap="round" d="M9.2 7h5.6M9.2 17h5.6M8.5 9.5h7M8.5 14.5h7"/>',
}
# first substring hit wins — cancer routes ahead of organ so "thyroid cancer" -> cancer
_KW = [
    ("cancer", "cancer"), ("carcinoma", "cancer"), ("melanoma", "cancer"),
    ("thyroid", "endo"),
    ("diabetes", "metab"), ("hba1c", "metab"), ("cholesterol", "metab"),
    ("lipoprotein", "metab"), ("hdl", "metab"), ("ldl", "metab"), ("glucose", "metab"),
    ("gout", "metab"), ("triglyceride", "metab"), ("c-reactive", "metab"), ("liver", "metab"),
    ("parkinson", "neuro"), ("alzheimer", "neuro"), ("dementia", "neuro"),
    ("multiple sclerosis", "neuro"), ("migraine", "neuro"), ("depress", "neuro"),
    ("bipolar", "neuro"), ("schizo", "neuro"), ("autism", "neuro"), ("adhd", "neuro"),
    ("intelligence", "neuro"), ("cognitive", "neuro"),
    ("coronary", "cardio"), ("atrial", "cardio"), ("hypertension", "cardio"),
    ("heart", "cardio"), ("venous thrombo", "cardio"), ("stroke", "cardio"),
    ("aortic", "cardio"), ("aneurysm", "cardio"),
    ("chronic obstructive", "resp"), ("pulmonary", "resp"), ("copd", "resp"), ("asthma", "resp"),
    ("rheumatoid", "msk"), ("osteoarthritis", "msk"), ("osteoporosis", "msk"),
    ("ankylosing", "msk"), ("arthritis", "msk"),
    ("inflammatory bowel", "gut"), ("coeliac", "gut"), ("celiac", "gut"),
    ("crohn", "gut"), ("colitis", "gut"), ("gallstone", "gut"), ("gallbladder", "gut"),
    ("lupus", "msk"),
    ("kidney", "renal"), ("renal", "renal"),
    ("macular", "eye"), ("cataract", "eye"), ("glaucoma", "eye"),
    ("dermatitis", "skin"), ("atopic", "skin"), ("psoriasis", "skin"), ("eczema", "skin"),
    ("height", "anthro"), ("body mass", "anthro"), ("waist", "anthro"),
]
_SYSNAME = {"neuro": "brain / nervous system", "cardio": "heart / circulation",
            "resp": "lungs", "msk": "bones & joints", "metab": "blood sugar / lipids",
            "endo": "thyroid / hormones", "gut": "digestive tract", "eye": "eyes",
            "renal": "kidneys", "skin": "skin", "cancer": "cancer", "anthro": "body measure",
            "default": "general"}


def _sys(trait):
    t = (trait or "").lower()
    for kw, s in _KW:
        if kw in t:
            return s
    return "default"


def _f(v):
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def _fmt(v, pct=False, default="—"):
    if v in (None, "", "NA"):
        return default
    x = _f(v)
    if x is None:
        return html.escape(str(v))
    return f"{x*100:.0f}%" if pct else f"{x:g}"


def _is_boiler(s):
    return (s.startswith(("explains", "screening context"))
            or "reduced accuracy outside" in s)


def _trait_caveat(rows):
    segs = []
    for r in rows:
        for s in (r.get("allowed_statement") or "").split(" | "):
            s = s.strip()
            if (s and not s.startswith("absolute risk") and not _is_boiler(s)
                    and s not in segs):
                segs.append(s)
    return " | ".join(segs)


def _best_pct(rows):
    vals = [p for p in (_f(r.get("percentile")) for r in rows) if p is not None]
    return max(vals) if vals else -1.0


def _font_css():
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "resources", "heading-font.css")
    try:
        with open(p) as f:
            return f.read()
    except OSError:
        return ""


def _rep(trows):
    """Representative score for the headline: best evidence grade, then percentile."""
    return max(trows, key=lambda r: (GRADE_RANK.get(r.get("evidence_grade", "D"), 0),
                                     _f(r.get("percentile")) or -1))


def _rank_phrase(p):
    if p >= 50:
        top = 100 - p
        return "top 1%" if top < 1 else f"top {top:.0f}%"
    return "bottom 1%" if p < 1 else f"bottom {p:.0f}%"


def _nat_n(x):
    v = _f(x)
    if v is None or v <= 0 or v >= 1:
        return None
    return round(1 / v)


def _dot_array(you, avg):
    out = []
    for i in range(100):
        c = "base" if i < avg else ("on" if i < you else "off")
        out.append(f'<i class="d {c}"></i>')
    return "".join(out)


def _category(rep):
    p, g = _f(rep.get("percentile")), rep.get("evidence_grade", "D")
    if p is None:
        return UNCAL, None
    if g in ("A", "B"):
        return (ELEV if p >= 90 else PROT if p <= 10 else AVGC), p
    return (WEAKHI if p >= 90 else WEAK), p


# generic (trait-agnostic) explanations keyed on the interpretation category
MEANS = {
    ELEV: ("A raised <em>baseline</em> genetic likelihood over your lifetime. Genetics is "
           "one factor among many — age, sex and lifestyle also matter.",
           "Not a diagnosis and not a certainty. Most people at this rank never develop it."),
    PROT: ("A lower-than-average genetic likelihood — reassuring, but not protective on its own.",
           "Low genetic likelihood does not rule it out; other risk factors still apply."),
    WEAKHI: ("Nothing actionable yet. If larger, replicated studies appear, this can be re-scored.",
             "The high rank is noise from a low-quality score, not evidence of raised risk."),
    WEAK: ("The evidence behind this score is too thin to interpret.",
           "A low grade is about the science, not your result — treat it as no information."),
}

THYROID_NOTE = (
    "Both hypo- and hyperthyroidism score elevated here — that's expected, not "
    "contradictory. Graves' (over-active) and Hashimoto's (under-active) are both "
    "<b>autoimmune thyroid disease</b> and share risk genes, so a high score for one "
    "often means a high score for the other. One person can have both across a lifetime: "
    "Hashimoto's can begin with a brief over-active phase, and treated Graves' frequently "
    "ends up under-active.")


def render(results_dir, out=None):
    contract = os.path.join(results_dir, "pgs_scores.tsv")
    if not os.path.exists(contract):
        raise SystemExit(f"no contract at {contract} — run grade_pgs first")
    with open(contract) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    prov = {}
    pp = os.path.join(results_dir, "provenance.json")
    if os.path.exists(pp):
        prov = json.load(open(pp))
    out = out or os.path.join(results_dir, "report.html")

    by_trait = OrderedDict()
    for r in rows:
        by_trait.setdefault(r.get("trait", "?"), []).append(r)

    tnames = {t.strip().lower() for t in by_trait}
    thyroid_pair = {"hypothyroidism", "hyperthyroidism"} <= tnames

    # sample sex (for labelling sex-specific risk figures) + which traits carry a
    # sex-specific baseline (so we only annotate "in women/men" where it matters).
    sample_sex = prov.get("sample_sex")
    if not sample_sex:
        sp = os.path.join(results_dir, "sample_sex.txt")
        if os.path.exists(sp):
            sample_sex = open(sp).read().strip() or None
    base = absolute_risk.load_baselines()
    sexspec = ({t for (t, sx, a) in base if sx in ("male", "female")}
               - {t for (t, sx, a) in base if sx == "overall"})
    sexinfo = (SEX_WORD.get((sample_sex or "").lower()), sexspec)

    # ---- summary -------------------------------------------------------------------
    # one grade per trait (its representative score) so this reconciles with the
    # per-trait tier counts below; the raw per-score count lives in the header.
    tally = Counter(_rep(trows).get("evidence_grade", "D") for trows in by_trait.values())
    tally_html = " ".join(
        f'<span class="tg"><i class="dot" style="background:{GRADE_COLOR.get(g, "#57606a")}">'
        f'</i>{tally[g]}&nbsp;{g}</span>'
        for g in ("A", "B", "C", "D") if tally.get(g))
    chips, n_uninterp = [], 0
    for trait, trows in by_trait.items():
        cat, p = _category(_rep(trows))
        if cat == ELEV:
            s = _sys(trait)
            chips.append((p, f'<span class="chip"><span class="ci" style="color:'
                          f'{SYS.get(s, SYS["default"])}">{_sys_icon(s)}</span>'
                          f'{html.escape(trait)} <b>{_rank_phrase(p)}</b></span>'))
        elif cat in (WEAKHI, WEAK):
            n_uninterp += 1
    chips.sort(key=lambda x: x[0], reverse=True)
    n_notable = len(chips)
    n_other = len(by_trait) - n_notable - n_uninterp
    notable_html = ("".join(c for _, c in chips[:12]) if chips
                    else '<span class="muted">none at high confidence</span>')

    # ---- trait rows, grouped by confidence tier ------------------------------------
    groups = {g: [] for g, _, _ in TIERS}
    for trait, trows in by_trait.items():
        rep = _rep(trows)
        groups.setdefault(rep.get("evidence_grade", "D"), []).append((trait, trows, rep))
    for g in groups:
        groups[g].sort(key=lambda x: _f(x[2].get("percentile")) or -1, reverse=True)

    sections = []
    for g, label, blurb in TIERS:
        items = groups.get(g) or []
        if not items:
            continue
        rows_html = "".join(
            _row(trait, trows, rep, thyroid_pair, sexinfo) for trait, trows, rep in items)
        n = len(items)
        sections.append(
            f'<h2 class="tier"><span class="tdot" style="background:{GRADE_COLOR[g]}"></span>'
            f'{label}<span class="tn">{n} {"trait" if n == 1 else "traits"}</span></h2>'
            f'<p class="tblurb">{blurb}</p><div class="tbody">{rows_html}</div>')

    # ---- document ------------------------------------------------------------------
    sample = html.escape(str(rows[0].get("sample", "?"))) if rows else "?"
    pops = sorted({str(r.get("most_similar_pop", "")).strip() for r in rows
                   if str(r.get("most_similar_pop", "")).strip() not in ("", "—")})
    anc_txt = f" · {html.escape(', '.join(pops))} ancestry" if pops else ""
    ptxt = (f"poly-suite v{prov.get('poly_suite_version','?')} · "
            f"{'calibrated' if prov.get('calibrated') else 'uncalibrated'} · "
            f"build {prov.get('target_build','?')} · "
            f"generated {html.escape(str(prov.get('generated_at','?')))} · "
            f"contract sha {str(prov.get('contract_sha256',''))[:12]}")
    font_face = _font_css()

    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>poly-suite PGS report — {sample}</title>
<style>
{font_face}
:root{{--bg:#f6f8fb;--panel:#fff;--raise:#fbfcfe;--ink:#182430;--soft:#516170;--faint:#7c8b99;
  --line:#dbe3ec;--line2:#eef2f6;--accent:#2f6f8f;--up:#b45309;--dn:#0f766e;--muted:#64748b;
  --hi:#0f766e;--ok:#2f6f8f;--low:#b45309;--none:#64748b;
  --dot-on:#c2703a;--dot-ref:#d9c3ab;--dot-off:#e4eaf0;--shadow:0 1px 2px rgba(20,40,60,.04),0 3px 12px rgba(20,40,60,.05)}}
*{{box-sizing:border-box}}
body{{font:15.5px/1.55 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
  background:var(--bg);color:var(--ink);max-width:820px;margin:0 auto;padding:2rem 1.1rem 4rem;-webkit-font-smoothing:antialiased}}
h1,h2,h3,.eyebrow,.big{{font-family:"Report Display",ui-sans-serif,system-ui,sans-serif}}
h1{{font-size:1.6rem;font-weight:600;letter-spacing:-.015em;margin:0 0 .25rem}}
.sub{{color:var(--soft);margin:0 0 1.1rem;font-size:.95rem}}
.disc{{background:#eef6ff;border:1px solid #bcd9f5;border-radius:8px;padding:.6rem .85rem;font-size:.86rem;margin:0 0 1.1rem;color:#1d3b57}}
.summary{{background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:1.1rem 1.25rem;margin:0 0 1.6rem;box-shadow:var(--shadow)}}
.summary h2{{font-size:1rem;font-weight:600;margin:0 0 .5rem}}
.chips{{display:flex;flex-wrap:wrap;gap:.4rem;margin:.2rem 0 .2rem}}
.chip{{display:inline-flex;align-items:center;gap:.32rem;background:var(--raise);border:1px solid var(--line);border-radius:2rem;padding:.14rem .6rem .14rem .4rem;font-size:.82rem;text-transform:capitalize}}
.chip .ci{{display:inline-flex;align-items:center;line-height:1}} .chip .ci .ig{{width:1rem;height:1rem}}
.chip b{{color:var(--up)}} .muted{{color:var(--soft)}}
.glbl{{color:var(--soft);font-weight:500}}
.tally{{display:flex;gap:1.1rem;flex-wrap:wrap;font-size:.82rem;color:var(--soft);margin-top:.85rem;padding-top:.75rem;border-top:1px solid var(--line2)}}
.tally b{{color:var(--ink)}} .grades{{display:flex;gap:1rem;flex-wrap:wrap;font-weight:600;margin-top:.5rem;font-size:.82rem}}
.tg{{display:inline-flex;align-items:center;gap:.35rem}} .dot{{width:.58rem;height:.58rem;border-radius:50%}}
/* ---- confidence tiers ---- */
.tier{{display:flex;align-items:center;gap:.5rem;font-size:1.05rem;font-weight:600;margin:1.7rem 0 .1rem;letter-spacing:-.01em}}
.tier .tdot{{width:.62rem;height:.62rem;border-radius:50%}}
.tier .tn{{margin-left:auto;font-family:ui-sans-serif,system-ui;font-size:.8rem;font-weight:500;color:var(--faint)}}
.tblurb{{margin:.1rem 0 .5rem 1.12rem;font-size:.8rem;color:var(--faint)}}
.tbody{{display:flex;flex-direction:column;gap:.5rem}}
.tools{{display:flex;gap:.5rem;margin:.7rem 0 0}}
.tbtn{{font:inherit;font-size:.78rem;font-weight:600;color:var(--accent);background:var(--panel);
  border:1px solid var(--line);border-radius:2rem;padding:.24rem .8rem;cursor:pointer}}
.tbtn:hover{{background:var(--raise);border-color:var(--accent)}}
/* ---- compact trait row (details/summary, no JS) ---- */
.t{{background:var(--panel);border:1px solid var(--line);border-radius:12px;box-shadow:var(--shadow);overflow:hidden;--stripe:var(--line)}}
.t.elevated{{--stripe:var(--up)}} .t.protective{{--stripe:var(--dn)}}
.t.weak-high,.t.weak,.t.uncalibrated{{--stripe:var(--muted)}} .t.average{{--stripe:var(--line)}}
.t>summary{{list-style:none;cursor:pointer;display:grid;grid-template-columns:auto 1fr auto;
  grid-template-rows:auto auto;column-gap:.7rem;row-gap:.12rem;align-items:center;
  padding:.62rem .85rem .62rem 0;border-left:3px solid var(--stripe)}}
.t>summary::-webkit-details-marker{{display:none}}
.t>summary:hover{{background:var(--raise)}}
.ico{{grid-row:1/3;justify-self:center;width:2.35rem;height:2.35rem;margin-left:.7rem;border-radius:10px;
  display:flex;align-items:center;justify-content:center;color:var(--sys);
  background:color-mix(in srgb,var(--sys) 14%,#fff);border:1px solid color-mix(in srgb,var(--sys) 30%,transparent)}}
.ico .ig{{width:1.4rem;height:1.4rem}}
.r1{{grid-column:2;grid-row:1;display:flex;align-items:center;gap:.5rem;flex-wrap:wrap;min-width:0}}
.tname{{font-weight:600;font-size:1.02rem;text-transform:capitalize;letter-spacing:-.01em}}
.lbar{{position:relative;flex:1;min-width:80px;max-width:150px;height:7px;border-radius:5px;
  background:linear-gradient(90deg,#e7edf3,#eef2f6);margin-left:auto}}
.lbar .avg{{position:absolute;left:50%;top:-2px;bottom:-2px;width:1px;background:var(--faint);opacity:.45}}
.lbar .mk{{position:absolute;top:50%;width:12px;height:12px;border-radius:50%;transform:translate(-50%,-50%);border:2px solid var(--panel);box-shadow:0 0 0 1px rgba(0,0,0,.12)}}
.conf{{display:inline-flex;align-items:center;gap:.3rem;font-size:.72rem;font-weight:600;padding:.1rem .45rem;border-radius:2rem;border:1px solid var(--line);white-space:nowrap}}
.conf .cdot{{width:.44rem;height:.44rem;border-radius:50%}}
.conf.hi{{color:var(--hi);border-color:color-mix(in srgb,var(--hi) 40%,transparent)}} .conf.hi .cdot{{background:var(--hi)}}
.conf.ok{{color:var(--ok);border-color:color-mix(in srgb,var(--ok) 40%,transparent)}} .conf.ok .cdot{{background:var(--ok)}}
.conf.low{{color:var(--low);border-color:color-mix(in srgb,var(--low) 40%,transparent)}} .conf.low .cdot{{background:var(--low)}}
.conf.none{{color:var(--none)}} .conf.none .cdot{{background:var(--none)}}
.interp{{grid-column:2;grid-row:2;font-size:.86rem;color:var(--soft);min-width:0}}
.interp b.up{{color:var(--up)}} .interp b.dn{{color:var(--dn)}} .interp .vs{{color:var(--faint)}}
.tri{{grid-column:3;grid-row:1/3;justify-self:center;margin-right:.2rem;color:var(--faint);font-size:.7rem;
  transition:transform .15s ease;transform:rotate(0deg)}}
.t[open]>summary .tri{{transform:rotate(90deg)}}
.t[open]>summary{{border-bottom:1px solid var(--line2)}}
.cbadge{{font-size:.68rem;font-weight:600;padding:.1rem .45rem;border-radius:2rem}}
.cbadge.ok{{background:#dafbe1;color:#1a7f37}} .cbadge.low{{background:#ffebe9;color:#cf222e}}
/* ---- expanded body ---- */
.body{{padding:.35rem 1rem 1rem 1.15rem}}
.tnote{{background:#eef6ff;border:1px solid #cfe4f7;border-radius:8px;padding:.55rem .8rem;font-size:.84rem;color:#1d3b57;margin:.7rem 0 0}}
.tnote b{{color:#0d2b45}}
.risk{{text-align:center;margin:1rem auto 0;max-width:280px}}
.array{{display:grid;grid-template-columns:repeat(20,1fr);gap:2px;max-width:230px;margin:0 auto}}
.array .d{{aspect-ratio:1;border-radius:50%}} .array .on{{background:var(--dot-on)}} .array .base{{background:var(--dot-ref)}} .array .off{{background:var(--dot-off)}}
.big{{font-size:1.25rem;font-weight:600;margin:.5rem 0 0}} .big b{{color:var(--up)}}
.rsub{{font-size:.82rem;color:var(--soft)}}
.lgd{{display:inline-flex;flex-wrap:wrap;justify-content:center;gap:.3rem;align-items:center;font-size:.75rem;color:var(--faint);margin-top:.3rem}}
.lgd i{{width:.58rem;height:.58rem;border-radius:50%}}
.notinterp{{background:var(--raise);border:1px dashed var(--line);border-radius:10px;padding:.7rem .85rem;color:var(--soft);font-size:.86rem;margin:.8rem 0 0}}
.notinterp b{{color:var(--ink)}}
.means{{display:grid;grid-template-columns:1fr 1fr;gap:1.1rem;margin:1rem 0 0;padding-top:.9rem;border-top:1px solid var(--line2)}}
.means h3{{font-size:.72rem;font-weight:700;letter-spacing:.05em;text-transform:uppercase;color:var(--faint);margin:0 0 .25rem}}
.means p{{margin:0;font-size:.87rem;color:var(--soft)}}
.todo{{background:#f0f9f4;border:1px solid #bfe6cd;border-radius:8px;padding:.55rem .8rem;font-size:.85rem;color:#123a25;margin:.8rem 0 0}}
.todo b{{color:#0b5132}}
.cdt{{margin:1rem 0 0;padding-top:.85rem;border-top:1px solid var(--line2)}}
.cdt h3{{font-size:.72rem;font-weight:700;letter-spacing:.05em;text-transform:uppercase;color:var(--faint);margin:0 0 .35rem}}
.ctab{{width:100%;border-collapse:collapse;font-size:.82rem;font-variant-numeric:tabular-nums}}
.ctab th{{text-align:left;color:var(--faint);font-weight:500;font-size:.76rem;padding:.2rem .5rem;border-bottom:1px solid var(--line)}}
.ctab td{{padding:.26rem .5rem;border-bottom:1px solid var(--line2)}}
.ctab .mono{{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.78rem;color:var(--soft)}}
.ctab a.mono{{color:var(--accent);text-decoration:none}} .ctab a.mono:hover{{text-decoration:underline}}
.sources{{font-size:.8rem;color:var(--soft);margin:.7rem 0 0}}
.sources b{{color:var(--soft)}} .sources a{{color:var(--accent);text-decoration:none}}
.sources a:hover{{text-decoration:underline}} .sources .lmute{{color:var(--faint)}}
.gb{{color:#fff;font-weight:600;font-size:.72rem;padding:.03rem .38rem;border-radius:2rem}}
.cav{{background:#fff8c5;border:1px solid #eac54f;border-radius:6px;padding:.45rem .65rem;font-size:.82rem;color:#4d3800;margin:.55rem 0 0}}
.tested{{font-size:.8rem;color:var(--faint);margin:.7rem 0 0}}
.tested b{{color:var(--soft)}}
.key{{background:var(--raise);border:1px solid var(--line);border-radius:12px;padding:1rem 1.2rem;margin:1.8rem 0 0}}
.key h2{{font-size:.92rem;font-weight:600;margin:0 0 .6rem}}
.key dl{{display:grid;grid-template-columns:auto 1fr;gap:.3rem 1rem;margin:0;font-size:.86rem}}
.key dt{{font-weight:600;color:var(--accent)}} .key dd{{margin:0;color:var(--soft)}}
footer{{color:var(--faint);font-size:.78rem;margin-top:1.6rem;border-top:1px solid var(--line);padding-top:.8rem;font-family:ui-monospace,monospace}}
@media (max-width:560px){{.means{{grid-template-columns:1fr}} .lbar{{max-width:110px}}}}
</style></head><body>
<h1>What your genetic results suggest — {sample}</h1>
<p class="sub">poly-suite · {len(by_trait)} traits · {len(rows)} scores{anc_txt} · educational / research use only</p>
<p class="disc"><b>Not a diagnostic test.</b> Polygenic scores estimate <em>genetic likelihood</em>
from common variants; they explain only a fraction of trait risk and are less accurate outside the
training ancestry. Genetics is one factor among many. Discuss anything actionable with a clinician
or genetic counselor.</p>
<section class="summary">
  <h2>The results worth your attention</h2>
  <div class="chips">{notable_html}</div>
  <div class="tally"><span><b>{n_notable}</b> worth attention (high confidence &amp; elevated)</span>
    <span><b>{n_other}</b> around average or not notable</span>
    <span><b>{n_uninterp}</b> shown but not reliable enough to interpret</span></div>
  <div class="grades"><span class="glbl">Evidence grade — one per trait:</span>{tally_html}</div>
</section>
<p class="sub" style="margin:.2rem 0 0">Traits are grouped by how much to trust the score — most trustworthy first.
Tap any row for the full clinical detail.</p>
<div class="tools"><button type="button" class="tbtn" onclick="document.querySelectorAll('details.t').forEach(function(d){{d.open=true}})">Expand all</button><button type="button" class="tbtn" onclick="document.querySelectorAll('details.t').forEach(function(d){{d.open=false}})">Collapse all</button></div>
{''.join(sections)}
<section class="key">
  <h2>How to read every row</h2>
  <dl>
    <dt>Icon</dt><dd>The body system the trait affects (brain, heart, thyroid, …) — a visual anchor, not a result.</dd>
    <dt>Likelihood bar</dt><dd>Where your score sits vs others of the same ancestry (left = lower, right = higher).
      The tick is the population average.</dd>
    <dt>Confidence</dt><dd>How much to trust the number. <b>High</b> = large, replicated, ancestry-checked
      studies. <b>Insufficient</b> = thin science; we don't translate it into a risk.</dd>
    <dt>Lifetime chance</dt><dd>The percentile turned into an actual probability using the disease's
      effect size and how common it is — shown only when confidence allows.</dd>
    <dt>Coverage</dt><dd>How much of the score we could measure in your data. Completeness, not quality.</dd>
  </dl>
</section>
<footer>{ptxt}</footer>
</body></html>"""
    with open(out, "w") as fh:
        fh.write(doc)
    return out


# OBO prefix -> EBI OLS4 ontology slug, for a direct term page. MONDO is handled
# separately (Monarch). A malformed id (non-numeric local part) falls to search.
_OBO_OLS = {"EFO": "efo", "HP": "hp", "ORPHANET": "ordo", "ORPHA": "ordo",
            "DOID": "doid", "NCIT": "ncit", "GO": "go", "OBA": "oba"}


def _onto_url(eid):
    """(url, label) for a trait ontology id -> the term's own page: Monarch for
    MONDO, the EBI OLS4 class page for other OBO ontologies; OLS search as a last
    resort for anything malformed/unknown. None if absent."""
    if not eid or eid in ("NA", "—"):
        return None, None
    lab = eid.replace("_", ":")
    if eid.startswith("MONDO_"):
        return f"https://monarchinitiative.org/{lab}", lab
    pre, _, loc = eid.partition("_")
    slug = _OBO_OLS.get(pre.upper())
    if slug and loc.isdigit():
        return f"https://www.ebi.ac.uk/ols4/ontologies/{slug}/classes?obo_id={lab}", lab
    return f"https://www.ebi.ac.uk/ols4/search?q={html.escape(eid)}", lab


def _sys_icon(s, cls="ig"):
    """Inline filled SVG glyph for a body system (currentColor, sized by CSS class)."""
    return (f'<svg class="{cls}" viewBox="0 0 24 24" fill="currentColor" '
            f'aria-hidden="true">{SVG.get(s, SVG["default"])}</svg>')


def _row(trait, trows, rep, thyroid_pair, sexinfo=(None, frozenset())):
    """One compact, collapsible trait row (summary = 2 lines; body = clinical layer)."""
    cat, p = _category(rep)
    sex_word, sexspec = sexinfo
    sexlab = f' in {sex_word}' if (sex_word and trait.strip().lower() in sexspec) else ''
    g = rep.get("evidence_grade", "D")
    clabel, ccls = CONFIDENCE.get(g, ("—", "none"))
    anc = html.escape(str(rep.get("most_similar_pop") or rep.get("training_ancestry") or "—"))
    s = _sys(trait)
    syscol = SYS.get(s, SYS["default"])

    # consensus badge (>=2 scores)
    n = str(trows[0].get("robustness_n_scores", "1"))
    conc = trows[0].get("robustness_concordance", "NA")
    low = any("robustness LOW" in (r.get("allowed_statement") or "") for r in trows)
    cbadge = ""
    if n not in ("1", "", "NA") and conc != "NA":
        cbadge = ('<span class="cbadge low">2 disagree</span>' if low
                  else f'<span class="cbadge ok">{n} agree</span>')

    # likelihood bar
    lbar = ""
    if p is not None:
        barcol = {ELEV: "var(--up)", PROT: "var(--dn)", AVGC: "var(--accent)"}.get(cat, "var(--muted)")
        lbar = (f'<span class="lbar" title="{p:.0f}th percentile"><span class="avg"></span>'
                f'<span class="mk" style="left:{max(2,min(98,p)):.1f}%;background:{barcol}"></span></span>')

    # one-line interpretation
    ar_n = _nat_n(rep.get("absolute_risk"))
    base_n = _nat_n(rep.get("baseline_incidence"))
    if cat == ELEV:
        extra = ""
        if ar_n:
            extra = (f' · about <b>1 in {ar_n}</b> lifetime{sexlab}'
                     + (f' <span class="vs">(vs 1 in {base_n} typical)</span>' if base_n else ''))
        interp = (f'Genetic likelihood <b class="up">higher than average</b> — '
                  f'{_rank_phrase(p)} for {anc} ancestry{extra}')
    elif cat == PROT:
        interp = (f'Genetic likelihood <b class="dn">lower than average</b> — '
                  f'{_rank_phrase(p)} for {anc} ancestry')
    elif cat == AVGC:
        interp = f'Around <b>typical</b> genetic likelihood — {p:.0f}th percentile'
    elif cat == WEAKHI:
        interp = f'Ranks high ({_rank_phrase(p)}) but <b>evidence too weak</b> to interpret'
    elif cat == WEAK:
        interp = 'Evidence <b>insufficient</b> to interpret'
    else:
        interp = 'Not calibrated — needs an ancestry reference panel'

    # ---- expanded body ----
    thy = ""
    if thyroid_pair and trait.strip().lower() in ("hypothyroidism", "hyperthyroidism"):
        thy = f'<div class="tnote">{THYROID_NOTE}</div>'

    risk = ""
    if cat == ELEV and ar_n:
        you = min(100, max(1, round(1 / ar_n * 100)))
        avg = min(you, round(1 / base_n * 100)) if base_n else 0
        base_txt = f' · vs about 1 in {base_n} on average' if base_n else ''
        legend = (f'<p class="lgd"><i style="background:var(--dot-ref)"></i>'
                  f'<span>average</span><i style="background:var(--dot-on)"></i>'
                  f'<span>your added likelihood</span></p>' if base_n else '')
        risk = (f'<div class="risk"><div class="array">{_dot_array(you, avg)}</div>'
                f'<p class="big">about <b>1 in {ar_n}</b></p>'
                f'<p class="rsub">estimated lifetime chance{sexlab}{base_txt}</p>{legend}</div>')
    elif cat == WEAKHI:
        risk = ('<div class="notinterp"><b>Why no number?</b> This score comes from small or '
                'unreplicated studies (grade D). Turning a confident-looking rank from weak '
                'science into "1 in X" would be false precision, so we don\'t. Raw values are '
                'in the table below.</div>')

    means = ""
    if cat in MEANS:
        mt, mn = MEANS[cat]
        means = (f'<div class="means"><div><h3>What this means</h3><p>{mt}</p></div>'
                 f'<div><h3>What it doesn\'t mean</h3><p>{mn}</p></div></div>')

    crows = []
    for r in sorted(trows, key=lambda r: _f(r.get("percentile")) or -1, reverse=True):
        rg = r.get("evidence_grade", "D")
        sid = html.escape(r.get("pgs_id", "?"))
        sid_html = (f"<a class='mono' href='https://www.pgscatalog.org/score/{sid}/' "
                    f"target='_blank' rel='noopener noreferrer'>{sid}</a>"
                    if sid.startswith("PGS") else f"<span class='mono'>{sid}</span>")
        crows.append(
            f"<tr><td>{sid_html}</td>"
            f"<td><span class='gb' style='background:{GRADE_COLOR.get(rg,'#57606a')}'>{rg}</span></td>"
            f"<td>{_fmt(r.get('percentile'))}</td><td>{_fmt(r.get('z_score'))}</td>"
            f"<td>{_fmt(r.get('risk_ratio'))}</td>"
            f"<td>{_fmt(r.get('absolute_risk'), pct=True)}</td>"
            f"<td>{_fmt(r.get('match_rate'), pct=True)}</td></tr>")
    cav = _trait_caveat(trows)
    cav_row = (f'<p class="cav">{html.escape(cav)}</p>') if cav else ""
    clinical = (f'<div class="cdt"><h3>Clinical detail</h3><table class="ctab">'
                f'<tr><th>score</th><th>grade</th><th>percentile</th><th>Z</th>'
                f'<th>risk ratio</th><th>abs. risk</th><th>coverage</th></tr>'
                f'{"".join(crows)}</table>{cav_row}</div>')

    # authoritative links: the MONDO disease class + the study behind the top score.
    # Per-score cohorts are reachable via the linked PGS Catalog score IDs above.
    ou, olab = _onto_url(rep.get("efo_id") or "")
    pmid = str(rep.get("source_pmid") or "").strip()
    lparts = []
    if ou:
        lparts.append(f'<a href="{ou}" target="_blank" rel="noopener noreferrer">{olab}</a>'
                      f' <span class="lmute">disease class</span>')
    if pmid and pmid != "NA":
        lparts.append(f'<a href="https://pubmed.ncbi.nlm.nih.gov/{html.escape(pmid)}/" '
                      f'target="_blank" rel="noopener noreferrer">study PMID&nbsp;{html.escape(pmid)}</a>')
    sources = (f'<p class="sources"><b>Learn more:</b> {" · ".join(lparts)}'
               f' · each score ID above opens its PGS Catalog page.</p>' if lparts else '')

    nsc = len(trows)
    cov = _fmt(rep.get("match_rate"), pct=True)
    tested = (f'<p class="tested"><b>Tested:</b> {nsc} published '
              f'{"scores" if nsc != 1 else "score"} for {html.escape(trait)}, '
              f'{anc}-ancestry GWAS; {cov} of variants measured.</p>')

    # one-line, non-prescriptive next step for high-confidence elevated findings
    todo = ('<p class="todo"><b>What to do:</b> nothing urgent — this is background '
            'likelihood, not a diagnosis. Worth mentioning to your doctor or a genetic '
            'counsellor, especially if it runs in your family.</p>' if cat == ELEV else '')

    rich = cat in (ELEV, PROT, WEAKHI)
    return (f'<details class="t {cat}" style="--sys:{syscol}">'
            f'<summary><span class="ico" role="img" aria-label="{_SYSNAME.get(s,s)} icon">{_sys_icon(s)}</span>'
            f'<span class="r1"><span class="tname">{html.escape(trait)}</span>{cbadge}'
            f'{lbar}<span class="conf {ccls}"><span class="cdot"></span>{clabel}</span></span>'
            f'<span class="interp">{interp}</span>'
            f'<span class="tri">▶</span></summary>'
            f'<div class="body">{thy}{risk if rich else ""}{means if rich else ""}{todo}'
            f'{clinical}{sources}{tested}</div></details>')


if __name__ == "__main__":
    rd = sys.argv[1] if len(sys.argv) > 1 else "results"
    print(f"html report -> {render(rd)}")
