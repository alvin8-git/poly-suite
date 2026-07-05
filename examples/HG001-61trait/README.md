# HG001 — 61-trait launch-set card

The full launch set scored end-to-end on GIAB **HG001 / NA12878** (a public reference
sample, WGS BAM), ancestry-calibrated against the HGDP+1kGP panel.

- **118 scores across 60 traits** (`select_pgs.py … all`, top-2/trait → 119 ids, 1 dropped
  on coverage), force-genotyped at the union of **12.6M autosomal loci** (~99.9% genotyped).
- Ancestry: **EUR** (FRAPOSA most-similar-pop).
- Evidence grades: **73 A / 15 B / 16 C / 14 D** — the honesty gates and multi-score
  robustness flags fire per trait.
- Notable elevated calls: hypothyroidism 99th (A), Parkinson 98th (A), rheumatoid
  arthritis 97th (A), type 1 diabetes 97th (B), inflammatory bowel disease 97th (B).

Files: [`report.html`](report.html) (standalone) · `pgs_scores.tsv` / `.json` (the contract)
· `provenance.json` (versions + sha256) · `pop_summary.csv` (ancestry call).

**The report** is patient-first: 60 traits grouped by evidence **confidence** (High → Insufficient,
so grade-D hyperthyroidism sinks to the bottom), each a compact row with a filled organ icon, a
likelihood bar, and a plain-language verdict that expands to the clinical detail. Two elevated
grade-A/B traits with a verified effect show an absolute risk as a natural frequency — rheumatoid
arthritis **~1 in 14 (women)** and IBD **~1 in 31**; the 97th–99th-percentile T1D / Parkinson /
hypothyroidism cards stay number-less on purpose (no standard per-SD effect or no defensible
baseline). Each score's clinical row also carries its publication-reported AUROC/C-index/R² beside
the grade, and links out to the PGS Catalog score, the MONDO disease class, and the source study.

Reproduce with `bin/run.sh … --pgs-meta <select_pgs output>` (see [../../README.md](../../README.md)
and [docs/Documentation.md §6](../../docs/Documentation.md#6-the-tiered-launch-set)).
Research/educational output — not a diagnostic test.
