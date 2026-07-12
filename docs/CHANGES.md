# Changelog

Notable changes to poly-suite. The shipped-feature list lives in
[TODO.md](../TODO.md); this file records dated releases and the runs behind them.

## 2026-07-12 — 74-trait launch set + full HG002 run

**Trait expansion (70 → 74).** `bin/select_pgs.py` `LAUNCH_SET` now defines **74
traits** — 25 core + 40 extended + 9 gated — which resolve to **142 PGS scorefiles**
at `TOP_N=2`. The five newly wired / re-pinned traits:

| Trait | Tier | PGS ID (pinned) | Note |
|---|---|---|---|
| multiple myeloma | extended | (Catalog, grade A) | replaces the dropped PCOS (no Catalog PGS); SEER-baselined |
| chronotype / morningness | extended | `PGS002209` | Privé 2022 LDpred2 substitute; direction "more evening" → invert |
| educational attainment | gated | `PGS002231` | Privé 2022 UKB LDpred2 reconstruction (EA4/Okbay not deposited) |
| neuroticism | gated | `PGS002213` | Privé 2022 UKB LDpred2 |
| loneliness | gated | `PGS001091` | Tanigawa 2022; thin (660 variants) → expect grade C/D self-downgrade |

height and bone mineral density were already in the set (both scored in the prior
launch70 run). The four behavioural/quantitative additions render in **predisposition
mode** (direction + percentile, no "1 in N"), flowing through the unchanged 24-column
`pgs_scores.tsv` contract with the absolute-risk columns left `NA`.

**Deep ancestry — Neanderthal card.** `bin/neanderthal.py` +
`resources/neanderthal_panel.tsv` emit `score/neanderthal.tsv`
(`sample, neanderthal_pct, method`), surfaced as a "Deep ancestry" card in
`report.html`. The shipped panel is a documented SEED, so percentages are flagged
PROVISIONAL.

**HG002 launch run.** Scored end-to-end on GIAB HG002 via
`bin/run.sh --sample HG002 … --outdir results/launch74/HG002 --batch` (pgsc_calc,
batch mode: 5 bins + FRAPOSA ancestry). Outputs:
`results/launch74/HG002/score/pgs_scores.tsv` (+ `.json`, `provenance.json`,
`report.html`). The demo report is regenerated at `docs/demo/HG002_report.html`.

Doc/count updates: README launch-set badge 70 → 74; README + Documentation.md §6
counts (extended 30 → 40, gated 6 → 9, full set → 74; 119 → 142 scorefiles); §6 trait
tables extended with chronotype, educational attainment, neuroticism, loneliness.

### HG002 scored results (verified from `pgs_scores.tsv`)

**74 unique traits, 141 score rows** of the 142 requested. Grade distribution across
the 141 scores: **87 A · 18 B · 19 C · 17 D**.

**One score dropped: `PGS005285`** (the *second* heart-failure score). Only ~1.6% of its
~994k variants matched the HG002 target (49.7% unmatched, 48.6% excluded), far below the
`min_overlap=0.1` gate, so pgsc_calc marked it `score_pass=false` and it was excluded.
Heart failure is **still covered** by its other score `PGS005097` (90.9% match, grade A),
so all 74 traits retain at least one score — hence 74 traits / 141 scores, not 73.

**All five new traits scored** (none failed):

| Trait | PGS | Grade | Match | Percentile |
|---|---|---|---|---|
| multiple myeloma | PGS000653 / PGS000654 | C / C | 0.68 / 0.67 | 94.2 / 92.8 |
| chronotype | PGS002209 | A | 0.90 | 68.8 |
| educational attainment | PGS002231 | A | 0.90 | 70.9 |
| neuroticism | PGS002213 | A | 0.90 | 23.1 |
| loneliness | PGS001091 | D | 0.71 | 60.7 (grade D → interpret with caution) |

loneliness self-downgrades to grade D as expected (thin, 660-variant score). The demo
report `docs/demo/HG002_report.html` renders all 74 traits with the Deep ancestry —
Neanderthal card.
