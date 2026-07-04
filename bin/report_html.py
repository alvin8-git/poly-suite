#!/usr/bin/env python3
"""Render the poly-suite graded contract to a self-contained HTML report
(spec §2 standalone deliverable). No external dependencies — inline CSS.

Grouped BY TRAIT: one section per trait, a consensus badge when >=2 scores were
run (concordant / robustness-LOW), and each score's grade + percentile + absolute
risk + coverage in a compact table. Provenance footer for reproducibility.

Importable (render) and runnable: python3 bin/report_html.py [results_dir]
"""
import os, sys, csv, json, html
from collections import OrderedDict

GRADE_COLOR = {"A": "#1a7f37", "B": "#9a6700", "C": "#bc4c00", "D": "#82071e"}


def _fmt(v, pct=False, default="—"):
    if v in (None, "", "NA"):
        return default
    try:
        f = float(v)
        return f"{f*100:.0f}%" if pct else (f"{f:g}")
    except (ValueError, TypeError):
        return html.escape(str(v))


def _trait_caveat(rows):
    """Shared caveat for a trait: dedup segments, drop per-score absolute-risk
    lines (shown in the table), keep variance/portability/coverage/consensus."""
    segs = []
    for r in rows:
        for s in (r.get("allowed_statement") or "").split(" | "):
            s = s.strip()
            if s and not s.startswith("absolute risk") and s not in segs:
                segs.append(s)
    return " | ".join(segs)


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

    cards = []
    for trait, trows in by_trait.items():
        n = str(trows[0].get("robustness_n_scores", "1"))
        conc = trows[0].get("robustness_concordance", "NA")
        low = any("robustness LOW" in (r.get("allowed_statement") or "") for r in trows)
        cbadge = ""
        if n not in ("1", "", "NA") and conc != "NA":
            if low:
                cbadge = '<span class="cbadge low">consensus: robustness LOW</span>'
            else:
                cbadge = f'<span class="cbadge ok">consensus: concordant ({conc})</span>'

        score_rows = []
        for r in trows:
            g = r.get("evidence_grade", "D")
            color = GRADE_COLOR.get(g, "#57606a")
            pctl = _fmt(r.get("percentile"))
            ar = _fmt(r.get("absolute_risk"), pct=True)
            base = _fmt(r.get("baseline_incidence"), pct=True)
            risk = f"{ar} (vs {base})" if ar != "—" else "—"
            score_rows.append(
                f"<tr><td class='pgs'>{html.escape(r.get('pgs_id','?'))}</td>"
                f"<td><span class='badge sm' style='background:{color}'>{html.escape(g)}</span></td>"
                f"<td>{pctl if pctl != '—' else 'uncalibrated'}</td>"
                f"<td>{risk}</td>"
                f"<td>{_fmt(r.get('match_rate'), pct=True)}</td></tr>")

        pop = html.escape(str(trows[0].get("most_similar_pop", "—")))
        cards.append(f"""
        <section class="card">
          <div class="head"><h2>{html.escape(trait)}</h2>{cbadge}
            <span class="pop">pop {pop}</span></div>
          <table class="scores">
            <tr><th>score</th><th>grade</th><th>percentile</th><th>absolute risk</th><th>coverage</th></tr>
            {''.join(score_rows)}
          </table>
          <p class="caveat">{html.escape(_trait_caveat(trows))}</p>
        </section>""")

    ptxt = (f"poly-suite v{prov.get('poly_suite_version','?')} · "
            f"{'calibrated' if prov.get('calibrated') else 'uncalibrated'} · "
            f"build {prov.get('target_build','?')} · "
            f"generated {html.escape(str(prov.get('generated_at','?')))} · "
            f"contract sha {str(prov.get('contract_sha256',''))[:12]}")
    sample = html.escape(str(rows[0].get("sample", "?"))) if rows else "?"

    doc = f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>poly-suite PGS report — {sample}</title>
<style>
 body{{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:880px;margin:2rem auto;padding:0 1rem;color:#1f2328}}
 h1{{font-size:1.4rem;margin:0 0 .2rem}} .sub{{color:#57606a;margin:0 0 1.5rem}}
 .card{{border:1px solid #d0d7de;border-radius:10px;padding:1rem 1.2rem;margin:1rem 0}}
 .head{{display:flex;align-items:center;gap:.6rem;flex-wrap:wrap}}
 .head h2{{font-size:1.1rem;margin:0;flex:1;text-transform:capitalize}}
 .pop{{color:#57606a;font-family:ui-monospace,monospace;font-size:.8rem}}
 .cbadge{{font-size:.75rem;font-weight:600;padding:.15rem .55rem;border-radius:2rem}}
 .cbadge.ok{{background:#dafbe1;color:#1a7f37;border:1px solid #4ac26b}}
 .cbadge.low{{background:#ffebe9;color:#cf222e;border:1px solid #ff8182}}
 .badge.sm{{color:#fff;font-weight:600;font-size:.72rem;padding:.05rem .4rem;border-radius:2rem}}
 table.scores{{border-collapse:collapse;margin:.7rem 0 0;width:100%;font-size:.9rem}}
 table.scores th{{text-align:left;color:#57606a;font-weight:500;font-size:.78rem;border-bottom:1px solid #d0d7de;padding:.2rem .5rem}}
 table.scores td{{padding:.28rem .5rem;border-bottom:1px solid #eaeef2}}
 td.pgs{{font-family:ui-monospace,monospace;font-size:.82rem;color:#57606a}}
 .caveat{{background:#fff8c5;border:1px solid #eac54f;border-radius:6px;padding:.5rem .7rem;font-size:.84rem;color:#4d3800;margin:.6rem 0 0}}
 footer{{color:#57606a;font-size:.8rem;margin-top:2rem;border-top:1px solid #d0d7de;padding-top:.8rem;font-family:ui-monospace,monospace}}
 .disc{{background:#ddf4ff;border:1px solid #54aeff;border-radius:6px;padding:.6rem .8rem;font-size:.85rem;margin:0 0 1.2rem}}
</style></head><body>
<h1>Polygenic risk report — {sample}</h1>
<p class="sub">poly-suite · {len(by_trait)} traits · educational / research use only</p>
<p class="disc"><b>Not a diagnostic test.</b> Polygenic scores estimate relative predisposition
from common variants; they explain only a fraction of trait risk and are less accurate outside
the training ancestry. Where two scores were run per trait, the <b>consensus</b> badge shows
whether they agree. Discuss anything actionable with a clinician or genetic counselor.</p>
{''.join(cards)}
<footer>{ptxt}</footer>
</body></html>"""
    with open(out, "w") as fh:
        fh.write(doc)
    return out


if __name__ == "__main__":
    rd = sys.argv[1] if len(sys.argv) > 1 else "results"
    print(f"html report -> {render(rd)}")
