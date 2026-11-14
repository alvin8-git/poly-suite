# Changelog

S1 and S2 stand for two non-GIAB samples (names withheld); GIAB samples are named as is.

Notable changes to poly-suite. The shipped-feature list lives in
[TODO.md](../TODO.md); this file records dated releases and the runs behind them.

## 2026-11-08 - grade_pgs: portability downgrade recognises every European-only label
- Defect: `bin/grade_pgs.py` compared `training_ancestry == "European"`, so contract labels `EUR`, `EUR100%` and `EUR (UKB-trained...` never triggered the portability downgrade for a non-EUR most-similar population.
- Fix: helper `is_european_only()` counts `EUR`, `EUR100%`, `European` and strings starting `EUR (` (case-insensitive) as European; mixed / multi-ancestry labels (`EUR, EAS`, `EUR/AFR`, `multi`, `mixed`) stay non-European and are not downgraded.
- Test: `test_is_european_only_labels`. Docs: Documentation.md section 5 gate table. Re-graded (no rescoring): results/launch78/S1 (`pgs_scores.tsv.pre_portability` keeps the old file).

## 2026-11-08 - `--vcf` ALT fill: multi-record positions (plink2 duplicate-ID abort)
- Failure: launch77 (S1/S2) aborted in pgsc_calc `PLINK2_ORIENT`: `--alt1-allele variant ID '9:135681662:G:T' appears multiple times`. Source (S1 allsites input): `bcftools call -m` all-sites output has, at some positions, a variant record plus a separate hom-ref record with ALT `.` (e.g. chr9:135681662 `G T 1/1` + `G . 0/0`; chr3:64572064 `A AAAGTAGG 0/1` + `A . 0/0`). Filling the hom-ref twin gave two records with the same CHROM:POS:REF:ALT. Measured on the restricted S1 file: 32,870 positions with >1 record (23,057 hom-ref+variant, 8,423 variant+variant, 994 variant+hom-ref, 396 hom-ref+hom-ref) of which 990 had identical REF/ALT pairs after the fill.
- Fix (`fill_homref_alt()`): records are buffered by (CHROM,POS). A position with >1 record is never filled; its hom-ref (ALT `.`) records are dropped when a variant record is present; a record repeating an earlier REF/ALT is dropped. No two output records share CHROM:POS:REF:ALT. Log `[vcf] alt_filled=N of M; skipped_multirecord=K; dropped_homref_dup=D; dropped_dup=E` and the same fields in `vcf_input.json` / `provenance.json` `parameters.vcf_input_alt_fill`; `--check-input` also reports `duplicate_allele_ids=` (must be 0). Effect on scoring: the dropped hom-ref twin of a variant is not needed (the variant record carries the alleles); unfilled hom-ref+hom-ref positions stay ALT `.` as before.
- Test: `test_vcf_homref_alt_filled_from_targets` extended (SNP + hom-ref twin, indel + exact repeat + hom-ref, two hom-ref). Rescoring: results/launch78/{S1,S2} (launch77 aborted).

## 2026-11-08 - `--vcf`: fill ALT on hom-ref target rows (match-rate defect)
- Defect: the restriction kept `bcftools call -m` hom-ref rows with `ALT "."` (e.g. `chr1 68082 T . GT:DP:AD 0/0:47:47`); pgsc_calc MATCH needs REF/ALT = scorefile alleles, so every hom-ref target site was unmatched. launch76 S1: median `match_rate` 0.479 (151/152 scores < 0.6; pgsc_calc summary ~45-53 % unmatched per score) vs median ~0.90 on the BAM path (launch75 HG002), where `call -C alleles` carries the target ALT. Explicit-genotype fraction was fine (0.9997), so the input itself was correct. The `results/track_d/HG001allsitesB69` ~0.90 figure was a prepped (allele-aware) VCF, not a raw one; the raw path had never been measured at that rate.
- Fix: `bin/scoring_targets.py` output is now also kept as an allele table (`vcf_targets_<hash>.alleles.tsv.gz`+`.tbi`, CHROM POS REF first-ALT; old caches rebuild once). Each per-contig worker pipes `bcftools view -r chrN -T targets -Ov` through `python3 bin/vcf_input.py fill ALLELES CHROM COUNTS` (`fill_homref_alt()`, `load_alleles()`, `alleles_path()`) and `bgzip`; ALT set only when ALT is `.` and REF equals the table REF, GT stays `0/0`, `AD/ADF/ADR` get a trailing `,0`. Variant rows untouched. Multi-ALT sites: first table ALT only (documented limitation). Streaming, parallel-safe.
- Reporting: log line `[vcf] alt_filled=N of M hom-ref target rows`; `--check-input` adds `alt_filled=`/`hom_ref_target_rows=`; `run.sh` writes `<outdir>/vcf_input.json` -> `score/vcf_input.json` -> `grade_pgs.py` puts it in `provenance.json` `parameters.vcf_input_alt_fill`.
- Test: `test_vcf_homref_alt_filled_from_targets`. Rescoring: results/launch77/{S1,S2}.

## 2026-11-07 - launch71: base-69 + 13 CircleDNA-gap polygenic scores
- Why: the OmniGen DTC report needs polygenic cards for CircleDNA-equivalent traits the base set lacks. One PGS per trait, chosen from the PGS Catalog REST API.
- Criteria: GRCh38 harmonised scorefile available; largest/most recent development GWAS; at least one independent evaluation in the catalog performance records; multi-ancestry preferred. `grade_suggestion` A only with two or more evaluation cohorts outside the UKB development data and AUC/R2/OR reported, else B.
- Added (trait: id): lean body mass PGS003531; sleep duration PGS003480; insomnia PGS000908; waist circumference PGS003555; male testosterone PGS003559; neuroticism PGS004608; educational attainment PGS003573; allergic disease PGS000800; idiopathic pulmonary fibrosis PGS001791; psoriatic arthritis PGS000198; vitiligo PGS000738; fasting insulin adj. BMI (HOMA-IR) PGS005276; acne PGS019637 (East Asian development only: weak transfer to other ancestries). Androgenetic alopecia: no score exists in the catalog.
- Machine-readable record: a local working file (not in the repo) (PGS id, n_variants, development ancestry/N, evaluation, scorefile URL). Metadata: `results/launch71/pgs_catalog_meta.json` = launch70 file + 13 entries (`trait_rank` 1, extra `category` key); launch70 entries unchanged. Ids: `results/launch71/pgs_ids.txt` (82 = 69 + 13).
- Scorefiles added to `cache/scorefiles` (82/82 present). Scoring parameters unchanged. Driver `scratch/run_launch71.sh` runs HG001 then HG002 sequentially from the BAM (as launch70 did) so the new scores' loci are genotyped; the launch70 prepped VCFs cover only launch70 targets and were not reused. Output `results/launch71/<S>`, timing in `results/launch71/timing.log`.
- Unverified: Ding 2023 scores (PGS003531/3555/3559/3573) have no development N in the catalog; evaluation is UKB held-out partitions only (grade B).

## 2026-11-07 - `--vcf` restricted to PGS target positions (all-sites input)
- Why: pgsc_calc's `PLINK2_VCF` runs out of memory on a whole-genome all-sites VCF (HG001, 20 GB, 2.88 G rows) after the autosome staging had copied the whole file for 69 min. A VCF restricted to the PGS loci scores in ~1 h 45 min at ~90 % median scorefile match; a variant-only VCF halves the match rate.
- `run.sh --vcf`: index the input if needed (`bcftools index -t`, logged) -> build the target list for `--pgs` (`bin/scoring_targets.py`, same as `genotype_prep.sh`; cached per PGS-set hash as `<scorefile-cache>/vcf_targets_<hash>.tsv.gz`+`.tbi`) -> per-contig parallel `bcftools view -r chrN -T targets` (`NPROC` workers, default 24; single `-R` over ~10 M points was >67 min, killed) + `concat --naive` -> `<outdir>/<S>.targeted.vcf.gz` (HG001 20 GB: 3 min with 24 workers, 4 min 19 s with 8) -> sex-chromosome exclusion on that small file -> `<S>.autosomes.vcf.gz`. Requires the scorefile cache to cover every `--pgs` id and chr-prefixed contigs. One log line: input rows, targets, restricted rows, explicit-genotype fraction, wall time.
- Records are kept whole, so hom-ref `0/0` rows of an all-sites input survive (unit-tested); a variant-only input yields a sparse file. The README input requirement is an all-sites VCF.
- `--check-input` now needs `--pgs` and reports index, input/restricted record counts, target positions, explicit-genotype fraction (denominator: targets inside the input's own per-contig footprint, so region slices are judged on what they cover), inferred sex; `CHECK_INPUT_FAILED: input is not an all-sites VCF` when the fraction is < 0.5. A normal run only warns.
- `bin/vcf_input.py`: `targets`, `afrac` subcommands. Tests: `test_vcf_target_list_builder_and_cache`, `test_vcf_restriction_keeps_homref_drops_offtarget`.
- Note: base-69 gives 9,786,436 autosomal target positions here (not 12.6 M, which was the HG001 prepped set's size).

## 2026-11-07 - raw-VCF input fix + scoring watchdog
- `run.sh --vcf`: stage autosomes-only (chrX/chrY/chrM excluded) before pgsc_calc import; plink2 aborted on chrX with no sex. Sex from `--sex` or VCF chrX non-PAR het (`bin/vcf_input.py`). 4,592 chrX rows lost across base-69 (4,576 in PGS004799).
- `--check-input`: staging + plink2 import test, then exit.
- `run_nf` watchdog: `--nf-timeout`/`NF_TIMEOUT` (8 h) and kill after `NF_FAIL_GRACE` (10 min) once the pipeline logs an error; failed scoring now exits non-zero.

## 2026-07-24 — Class-A runtime optimizations (rigor-preserving; outputs unchanged)

Faster runs **without changing any score or output**. Baseline: S1/S2
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
These are single-process and memory-bound (89–108 GB RSS in the S1/S2 trace).
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
