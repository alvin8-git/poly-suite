# Changelog

Notable changes to poly-suite. The shipped-feature list lives in
[TODO.md](../TODO.md); this file records dated releases and the runs behind them.

## 2026-07-24 — Class-A runtime optimizations (rigor-preserving; outputs unchanged)

Faster runs **without changing any score or output**. Baseline: THAL1/THAL2
(2026-07-23), ~3.5 h each, from an intact Nextflow trace. Nothing here touches the
scorefile set/count, scoring parameters, or the graded contract — only where time
is spent. Validated statically (`bash -n bin/run.sh`, `--dry-run` plan, and a
`nextflow config` parse of `conf/rootless.config`); **no scoring run was launched.**

**1. `DOWNLOAD_SCOREFILES` hang fixed + scorefile cache is now the default path.**
- *Symptom.* pgsc_calc's `DOWNLOAD_SCOREFILES` fetches scorefiles over the
  PGS-Catalog FTP with no per-file timeout. A stalled connection makes the process
  hang with **no exit status** (observed ~66 min idle on 2026-07-23). Nextflow's
  process `time` directive is **not enforced by the local executor**, so it cannot
  kill the hang; and `errorStrategy 'retry'` only fires on non-zero exits, so a
  silent stall never retries either.
- *(b) Cache-as-default (`bin/run.sh`).* A persistent local cache of harmonized
  hmPOS scorefiles already existed behind `--scorefile-cache`; it is now the
  **default**. When the flag is omitted, run.sh falls back to `$POLY_SCOREFILE_CACHE`,
  then to the first existing known cache dir (`results/launch70/scorefile_cache`,
  136 ids, or `cache/scorefiles`). The cache lookup is **all-or-nothing**: if it does
  not fully cover the requested `--pgs`, the run transparently reverts to `--pgs_id`
  (download), so **outputs are identical** either way. On a full hit,
  `DOWNLOAD_SCOREFILES` is skipped entirely (both the harmonize and scoring runs).
- *(a) Bounded timeout + retry (`bin/run.sh`).* The residual download path is the
  cache-miss harmonize run. It is now wrapped in `timeout -k 30 ${HM_TIMEOUT:-1500}s`
  with `${HM_RETRIES:-2}` retries. A hung fetch is killed after 25 min/attempt and
  **fails loud** (`rc=124`) instead of hanging indefinitely.
- *Config net (`conf/rootless.config`).* `withName:'.*DOWNLOAD_SCOREFILES'` gets
  `errorStrategy 'retry'`, `maxRetries 2` so *erroring* (non-zero) fetches retry.
- *Expected saving.* **−13 min** (skips a healthy download) **+ removes the hang**.

**2. Threaded the reference-panel PLINK2 scoring (`conf/rootless.config`).**
- `PLINK2_SCORE` (used for both the `--run_ancestry` reference panel and the single
  sample) runs at `process_low` (cpus=2) upstream; the reference-panel task
  (HGDP+1kGP, ~4k genomes) is the long pole (~26 min). Bumped to `cpus = 8` via
  `withName:'.*PLINK2_SCORE'` (bounded by `params.max_cpus = 16`).
- *Score-safe.* `plink2 --score --threads N --seed 31` is **deterministic in N** —
  identical `.sscore` values at any thread count — so this parallelizes without
  changing results. Per-sample scoring is already fast; extra threads are harmless.
- *Expected saving.* **−10 to −15 min**, only when ancestry/`--panel` is used.

**3. Deliberately NOT touched: `FORMAT_SCOREFILES` / `MATCH_COMBINE`.**
These are single-process and memory-bound (89–108 GB RSS in the THAL trace).
`cpus`/`maxForks` parallelism directives do **nothing** for them — do not add them
here; raising concurrency only risks the OOM already documented (2026-07-12). The
memory fix for those is `--batch` (see `docs/memory-optimization-plan.md`), not cores.

## 2026-07-15 — honest PMID/DOI provenance (`source_doi` added to the contract)

**Problem.** A handful of scored PGS carried an **empty `source_pmid`** in
`pgs_scores.tsv` — PGS004923 (type 2 diabetes), PGS005258/PGS005259 (thyroid
cancer), PGS005267/PGS005268 (hypothyroidism). An empty cell was ambiguous: it
could mean "the Catalog has no PMID" *or* "we failed to fetch metadata".

**Root cause (not a bug in the fetch).** These five scores are **medRxiv
preprints** (Ritchie SC 2024; White SL 2025). Their PGS Catalog `publication`
record genuinely has `PMID: null` — no PubMed ID exists upstream yet — but each
*does* carry a **DOI** (e.g. `10.1101/2024.08.22.24312440`), plus first author,
journal and PGP id. `bin/select_pgs.py` only captured `publication.PMID`
(dropping the DOI), and `bin/grade_pgs.py` had no fallback, so the DOI was lost
and the cell went blank.

**Fix.**
- `bin/select_pgs.py` — now also captures `doi` and `journal` from the Catalog
  `publication` object (previously discarded).
- `bin/grade_pgs.py` — new **`source_doi`** column in `CONTRACT_COLS` (24 → 25
  cols; scores append, no other column moved). `load_catalog_meta()` threads the
  DOI through. Provenance semantics are now explicit and **distinguishable**:
  - `source_pmid` non-empty → peer-reviewed publication (PubMed ID).
  - `source_pmid` `""` **+ `source_doi` set** → provenance *is* known but the
    source is a preprint/DOI-only record with no PMID upstream. **Not** a failure.
  - both `""` → provenance could not be resolved (unknown score / fetch failure).
- `bin/report_html.py` — when a score has no PMID it now renders a **`study DOI`**
  link (`https://doi.org/…`, labelled "preprint — no PMID yet") instead of an
  empty "Learn more" line.
- No PMIDs were fabricated. Scores with neither PMID nor DOI upstream stay
  explicitly blank in both columns (the unresolved case).

Verified against the live Catalog for all five IDs and with a new unit test
(`test_provenance_pmid_doi_fallback`). Regenerating `pgs_scores.tsv` requires a
scoring run; existing `pgs_catalog_meta.json` files predate the DOI capture, so
the `source_doi` column will populate on the next `select_pgs.py` run.

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
