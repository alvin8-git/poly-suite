#!/usr/bin/env python3
"""Render the poly-suite graded contract to a self-contained HTML report
(spec §2 standalone deliverable). No external dependencies — inline CSS.

Grouped BY TRAIT: one section per trait, a consensus badge when >=2 scores were
run (concordant / robustness-LOW), and each score's grade + percentile + absolute
risk + coverage in a compact table. Provenance footer for reproducibility.

Importable (render) and runnable: python3 bin/report_html.py [results_dir]
"""
import os, sys, csv, json, html
from collections import OrderedDict, Counter

GRADE_COLOR = {"A": "#1a7f37", "B": "#9a6700", "C": "#bc4c00", "D": "#82071e"}


def _fmt(v, pct=False, default="—"):
    if v in (None, "", "NA"):
        return default
    try:
        f = float(v)
        return f"{f*100:.0f}%" if pct else (f"{f:g}")
    except (ValueError, TypeError):
        return html.escape(str(v))


def _is_boiler(s):
    """Universal disclaimer text (variance / cross-ancestry accuracy / not-diagnostic)
    shown once at the top — not worth repeating in every card's caveat box."""
    return (s.startswith(("explains", "screening context"))
            or "reduced accuracy outside" in s)


def _trait_caveat(rows):
    """Trait-SPECIFIC caveats only: drop per-score absolute-risk lines (shown in the
    table) and the universal boilerplate (hoisted to the top disclaimer), leaving the
    real warnings — robustness-LOW, low coverage, uncalibrated. Empty -> no caveat box."""
    segs = []
    for r in rows:
        for s in (r.get("allowed_statement") or "").split(" | "):
            s = s.strip()
            if (s and not s.startswith("absolute risk") and not _is_boiler(s)
                    and s not in segs):
                segs.append(s)
    return " | ".join(segs)


def _best_pct(rows):
    """Highest percentile across a trait's scores, or -1 if none are calibrated."""
    vals = []
    for r in rows:
        try:
            vals.append(float(r.get("percentile") or ""))
        except (ValueError, TypeError):
            pass
    return max(vals) if vals else -1.0


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
    # surface the actionable findings: most-elevated traits first, uncalibrated last
    by_trait = OrderedDict(
        sorted(by_trait.items(), key=lambda kv: _best_pct(kv[1]), reverse=True))

    # summary: grade tally (doubles as the colour legend) + notable elevated findings
    tally = Counter(r.get("evidence_grade", "D") for r in rows)
    tally_html = " ".join(
        f'<span class="tg"><i class="dot" style="background:{GRADE_COLOR.get(g, "#57606a")}"></i>'
        f'{tally[g]}&nbsp;{g}</span>'
        for g in ("A", "B", "C", "D") if tally.get(g))
    chips = []
    for trait, trows in by_trait.items():
        best = max(trows, key=lambda r: _best_pct([r]))
        p, g = _best_pct([best]), best.get("evidence_grade", "D")
        if p >= 90 and g in ("A", "B"):
            chips.append(
                f'<span class="chip">{html.escape(trait)} <b>{p:.0f}th</b>'
                f'<span class="badge sm" style="background:{GRADE_COLOR[g]}">{g}</span></span>')
    notable_html = ("".join(chips[:10]) if chips
                    else '<span class="muted">none ≥90th percentile at grade A/B</span>')

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

        cav = _trait_caveat(trows)
        cav_html = f'<p class="caveat">{html.escape(cav)}</p>' if cav else ""
        cards.append(f"""
        <section class="card">
          <div class="head"><h2>{html.escape(trait)}</h2>{cbadge}</div>
          <table class="scores">
            <tr><th>score</th><th>grade</th><th>percentile</th><th>absolute risk</th><th>coverage</th></tr>
            {''.join(score_rows)}
          </table>
          {cav_html}
        </section>""")

    ptxt = (f"poly-suite v{prov.get('poly_suite_version','?')} · "
            f"{'calibrated' if prov.get('calibrated') else 'uncalibrated'} · "
            f"build {prov.get('target_build','?')} · "
            f"generated {html.escape(str(prov.get('generated_at','?')))} · "
            f"contract sha {str(prov.get('contract_sha256',''))[:12]}")
    sample = html.escape(str(rows[0].get("sample", "?"))) if rows else "?"
    pops = sorted({str(r.get("most_similar_pop", "")).strip() for r in rows
                   if str(r.get("most_similar_pop", "")).strip() not in ("", "—")})
    anc = f" · ancestry {html.escape(', '.join(pops))}" if pops else ""

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
 .disc{{background:#ddf4ff;border:1px solid #54aeff;border-radius:6px;padding:.6rem .8rem;font-size:.85rem;margin:0 0 1rem}}
 .summary{{background:#f6f8fa;border:1px solid #d0d7de;border-radius:8px;padding:.8rem 1rem;margin:0 0 1.5rem;display:flex;flex-direction:column;gap:.6rem}}
 .tally{{display:flex;gap:1.1rem;flex-wrap:wrap;font-size:.85rem;font-weight:600}}
 .tg{{display:inline-flex;align-items:center;gap:.35rem}}
 .dot{{width:.6rem;height:.6rem;border-radius:50%}}
 .notable{{display:flex;flex-wrap:wrap;gap:.4rem;align-items:baseline;font-size:.85rem}}
 .notable .lbl{{color:#57606a;font-weight:600}}
 .chip{{display:inline-flex;align-items:center;gap:.3rem;background:#fff;border:1px solid #d0d7de;border-radius:2rem;padding:.12rem .5rem;font-size:.8rem;text-transform:capitalize}}
 .muted{{color:#57606a}}
</style></head><body>
<h1>Polygenic risk report — {sample}</h1>
<p class="sub">poly-suite · {len(by_trait)} traits{anc} · educational / research use only</p>
<p class="disc"><b>Not a diagnostic test.</b> Polygenic scores estimate relative predisposition
from common variants; they explain only a fraction of trait risk and are less accurate outside
the training ancestry. Where two scores were run per trait, the <b>consensus</b> badge shows
whether they agree. Discuss anything actionable with a clinician or genetic counselor.</p>
<div class="summary">
  <div class="tally">{tally_html}</div>
  <div class="notable"><span class="lbl">Notable elevated:</span> {notable_html}</div>
</div>
{''.join(cards)}
<footer>{ptxt}</footer>
</body></html>"""
    with open(out, "w") as fh:
        fh.write(doc)
    return out


if __name__ == "__main__":
    rd = sys.argv[1] if len(sys.argv) > 1 else "results"
    print(f"html report -> {render(rd)}")
