# HG001 — 70-trait launch-set card

The full launch set scored end-to-end on GIAB **HG001 / NA12878** (a public reference
sample, WGS BAM), ancestry-calibrated against the HGDP+1kGP panel. This supersedes the
earlier [61-trait card](../HG001-61trait/) with the next tranche of traits (see
[docs/candidate-traits.md](../../docs/candidate-traits.md)).

- **135 scores across 69 traits.** The launch set requests **70** traits, but polycystic
  ovary syndrome has **no PGS in the Catalog** under its ontology id, so it drops — the
  report counts what was actually scored (69), not what was requested. Force-genotyped at
  the union of **~12.6M autosomal loci** (~99.9% genotyped).
- Ancestry: **EUR** (FRAPOSA most-similar-pop).
- Evidence grades (one per trait): **48 A / 9 B / 7 C / 5 D**.

**New in this tranche:**
- **Biomarker predisposition mode** — bone mineral density, eGFR, and systolic blood
  pressure render as *genetic predisposition* (direction + percentile, neutral colour, no
  "1 in N"), not disease risk. HG001's **systolic BP is 97th percentile (grade A)** →
  "Higher genetic predisposition — top 3%".
- **New disease/cancer traits with curated absolute risk** — peripheral arterial disease,
  endometriosis, esophageal and gastric cancer get a natural-frequency figure where elevated.
- **Percentile-only by design** — obstructive sleep apnea (no standard per-SD effect) and
  cervical cancer (HPV-driven; carries an on-card HPV note) show a rank but no risk number.

Files: [`report.html`](report.html) (standalone) · `pgs_scores.tsv` / `.json` (the contract)
· `provenance.json` (versions + sha256) · `pgs_performance.tsv` (reported AUROC/C-index/R²)
· `pop_summary.csv` (ancestry call).

Reproduce with `bin/run.sh … --pgs-meta <select_pgs output>` (see [../../README.md](../../README.md)
and [docs/Documentation.md §6](../../docs/Documentation.md#6-the-tiered-launch-set)).
Research/educational output — not a diagnostic test.
