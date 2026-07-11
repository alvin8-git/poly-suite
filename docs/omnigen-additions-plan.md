# poly-suite additions for OmniGen — implementation plan

**Status:** IMPLEMENTED (PGS expansion + Neanderthal %), 2026-07-11 — code in the working tree, in-session tests passing; real percentiles/percentages still need a full pgsc_calc run (command below). Haplogroups (item 3) remain a routing recommendation only — implemented in pgx-suite per that recommendation, NOT here. See the "Implementation status" section at the end.
**Scope:** three additions OmniGen will consume — (1) PGS expansion, (2) Neanderthal ancestry %, (3) routing assessment for mtDNA / Y haplogroups.

Investigated against the real repos:
- poly-suite: `bin/run.sh` (orchestrator), `bin/select_pgs.py` (launch-set + Catalog evidence grading), `bin/genotype_prep.sh` + `bin/scoring_targets.py` (force-genotype at scoring loci), `bin/grade_pgs.py` (contract writer), `bin/report_html.py` (standalone report), `bin/cache_scorefiles.py` / `cache_scorefile_paths()` (scorefile cache), `bin/partition_scores.py`, `conf/rootless.config`.
- OmniGen consumers: `prototype/report_pgs.py` (reads `**/pgs_scores.tsv`), `prototype/report_ancestry.py` (reads `*_popsimilarity.txt.gz`), `docs/implementation-placement.md`, `docs/blueprint-expansion.md`.
- pgx-suite: `docker/Snakefile` + `docker/genes.tsv` (gene×tool DAG), `mutserve` (MT-RNR1, baked-in jar), `docker/test_parsers.py::parse_mutserve`.

Key structural facts that constrain every option below:
- **The contract OmniGen reads is `results/launch70/<S>/score/pgs_scores.tsv`**, 24 columns, one row per (sample, trait): `sample, trait, efo_id, pgs_id, source_pmid, n_variants, training_ancestry, n_matched, match_rate, inferred_ancestry, most_similar_pop, ancestry_distance, percentile, z_score, ci_low, ci_high, absolute_risk, baseline_incidence, risk_ratio, robustness_n_scores, robustness_concordance, evidence_grade, portability_flag, allowed_statement`. `grade_pgs.py` reads `score/aggregated_scores.txt.gz` (pgsc_calc long format) generically per PGS ID — **new PGS IDs auto-appear as new rows; quantitative traits leave `absolute_risk/baseline_incidence/risk_ratio` = NA, which the contract already supports.**
- **`bin/scoring_targets.py` is hard-restricted to autosomes** (`AUTOSOMES = {chr1..chr22}`; line ~24/78 explicitly `exclude chrX/Y/M`). The prepped VCF has **0 chrM / 0 chrY / 0 chrX** records. This is decisive for additions 2 and 3.

---

## 1. PGS expansion

### 1a. Where it slots in
Two mechanisms already exist; both are used, no new plumbing:

- **Launch-set membership (recommended for the durable ones):** append `(tier, trait_label, trait_id)` rows to `LAUNCH_SET` in `bin/select_pgs.py`. `select_pgs.py` fetches candidate scores per trait from the PGS Catalog REST API, computes an evidence grade + data-sufficiency filter, and writes `results/pgs_catalog_meta.json`, which `grade_pgs.py` loads. Thin traits self-downgrade (grade D) rather than being excluded — so low-power scores are safe to add.
- **Direct IDs / scorefile cache (for a specific pinned score, or a non-Catalog custom score):** pass the PGS IDs in `--pgs "PGS...,PGS..."` to `bin/run.sh`, and pre-place the harmonized scorefile `<ID>_hmPOS_GRCh38.txt.gz` in the scorefile cache dir (`--scorefile-cache cache/scorefiles`). `cache_scorefile_paths()` resolves them, symlinks into `$OUTDIR/scorefiles/`, and pgsc_calc is invoked with `--scorefile "$SFDIR/*.txt.gz"` instead of `--pgs_id`, skipping the download step. `bin/cache_scorefiles.py` fetches+caches Catalog scorefiles.

Flow through the pgsc_calc invocation is **identical** to existing scores — same `nextflow run pgscatalog/pgsc_calc … --min_overlap 0.1 -c conf/rootless.config`. No change to `genotype_prep` (new loci are force-genotyped automatically because the harmonized scorefiles are the target source), no change to FRAPOSA (ancestry is score-independent).

### 1b. Tool + container
None new. pgsc_calc (already pinned via `pgscatalog/pgsc_calc` Nextflow + Docker `conf/rootless.config`) scores every added ID. `select_pgs.py` / `cache_scorefiles.py` hit the Catalog API via `urllib` (curl/wget blocked on this box).

### 1c. Contract output — NO schema change (confirmed)
Added scores produce additional (sample, trait) rows in the existing `pgs_scores.tsv`. Quantitative traits (height, EA, chronotype, neuroticism, BMD, loneliness) populate `percentile` + `z_score` and leave `absolute_risk/baseline_incidence/risk_ratio = NA` — exactly how `grade_pgs.py` already handles non-disease traits (it only fills absolute risk from `resources/baseline_incidence.tsv`, which has no entry for these). `evidence_grade`, `portability_flag`, `allowed_statement` are computed generically. **OmniGen's `report_pgs.py` reads the same file with zero changes.**

### 1d. Specific PGS IDs (resolved live against the Catalog REST API, 2026-07-11)

| Trait (requested) | Deployable PGS ID | Provenance | nv | Notes |
|---|---|---|---|---|
| Height (Yengo 2022) | **PGS002804** (EUR) or **PGS002802** (multi-anc "ALL") | GIANT, Yengo 2022, PMID 36224396 | ~1.1M | EUR matches the EUR calibration panel; use PGS002802 if multi-ancestry portability preferred. Trait id `OBA_VT0001253` (101 scores). |
| Educational attainment (EA4) | **PGS002231** (Privé LDpred2 `years_of_edu`) | Privé 2022, PMID 34995502 | 950,845 | **EA4/Okbay 2022 (PMID 35361970) is NOT deposited with usable weights** under `EFO_0011015`; PGS002231 is the genome-wide UKB LDpred2 reconstruction — the defensible substitute. Flag substitution in meta rationale. |
| Chronotype / morningness (Jones 2019) | **PGS002209** (Privé LDpred2 `more_evening`) | Privé 2022, PMID 34995502 | 955,439 | Jones 2019 (PMID 30696823) not deposited as a weighted score under `EFO_0008328`. **Direction = "more evening"; OmniGen renders morningness → invert the percentile interpretation** (document in meta). |
| Beat / rhythm synchronization (Niarchou 2022) | **not in Catalog** | Niarchou 2022, PMID 35618881 (23andMe) | — | 23andMe restricts weight deposition; no `--pgs_id`. **Deferred** unless supplemental weights are obtained → then add via the custom-scorefile path (§1a, cache a hand-built `<ID>_hmPOS_GRCh38.txt.gz`). |
| Bone mineral density (Morris 2019) | **PGS000657** (gSOS, Forgetta) | Forgetta 2020, PMID 32614825 | 21,716 | Morris 2019 eBMD (PMID 30598549) not a downloadable weighted score; gSOS is the grade-A heel-BMD proxy under `EFO_0009270`. Alternatively add the trait id and let `select_pgs.py` pick top-2 by evidence. |
| Neuroticism | **PGS002213** (Privé LDpred2) | Privé 2022, PMID 34995502 | 950,183 | Trait `EFO_0007660`. |
| Risk tolerance | **not in Catalog (clean)** | Karlsson Linnér 2019, PMID 30643258 | — | General-risk-tolerance weights not deposited (23andMe). Defer / custom-scorefile path. |
| Loneliness | **PGS001091** (Tanigawa GBE) | Tanigawa 2022, PMID 35324888 | 660 | Only score under `EFO_0007865`; thin → expect self-downgrade to grade C/D (acceptable, honest). |

**Recommended action:** add durable trait ids to `LAUNCH_SET` in `select_pgs.py` — height (`OBA_VT0001253`), educational attainment (`EFO_0011015`), chronotype (`EFO_0008328`), bone mineral density (`EFO_0009270` — already candidate in docs), neuroticism (`EFO_0007660`), loneliness (`EFO_0007865`) — tier `extended` (or `gated` for EA/risk-adjacent behavioral traits per the existing sensitivity convention: intelligence/ADHD/autism are already `gated`). Pin height/EA/chronotype/neuroticism to the specific IDs above via meta if `select_pgs.py`'s top-2 picks drift. Beat-sync + risk-tolerance: leave a documented stub (custom-scorefile path) until weights are sourced.

### 1e. Report integration (poly's own report.html)
Automatic. `report_html.py::render()` groups `by_trait` from the contract and renders every trait row. Quantitative traits already render as percentile/z rows without an absolute-risk block. One optional polish: `report_html._sys(trait)` maps trait→organ icon; add icons for height/EA/chronotype/neuroticism/BMD/loneliness (cosmetic, non-blocking).

### 1f. Testing
- **In-session (no full run):** (i) run `python3 bin/cache_scorefiles.py` for a new ID and confirm it writes `<ID>_hmPOS_GRCh38.txt.gz` in Catalog scorefile format; (ii) confirm `cache_scorefile_paths cache/scorefiles "<ID>"` resolves it; (iii) confirm `bin/partition_scores.py` bins it; (iv) confirm `scoring_targets.py` accepts its loci (autosomal filter). (v) Add a synthetic row to a copy of `aggregated_scores.txt.gz` and re-run `grade_pgs.py` to confirm a new contract row with NA risk columns and a sane grade. (vi) `tests/validate_contract.py` must still pass.
- **Needs a full run** (BAM + panel): real `percentile`/`z_score` for the new scores (requires the FRAPOSA-calibrated scoring pass). Fold into the next `launch70` re-run.

---

## 2. Neanderthal ancestry %

### 2a. Where it slots in
**Cannot reuse the prepped VCF** — it is autosomal PGS loci only, and archaic-introgression-informative SNPs are largely not at PGS loci. Neanderthal % needs its **own target panel**. Two options; recommend Option A (max infra reuse):

- **Option A — custom scorefile through the existing scoring path (recommended).** Build the archaic panel as a PGS-Catalog-format scorefile (`effect_allele` = archaic-derived allele, `effect_weight` = 1 per allele) and feed it through the *same* `genotype_prep` + pgsc_calc pipeline as an extra "score". Its raw `SUM` = count of archaic alleles the sample carries at panel loci; `DENOM` = alleles genotyped. A new small step `bin/neanderthal.py` reads that score's row from `aggregated_scores.txt.gz`, converts dosage fraction → % via a reference calibration (1kGP EUR mean ≈ 1.8–2.0%), and writes `neanderthal.tsv`. **Note:** the panel loci must be included in `genotype_prep` targets — since Option A's panel is genome-wide (not autosome-restricted for this purpose but archaic SNPs are autosomal anyway), it flows through `scoring_targets.py` cleanly.
- **Option B — standalone step.** `bin/neanderthal.sh`: reuse `genotype_prep.sh`'s mpileup machinery with the archaic panel as its sole "scorefile" to force-genotype the BAM at panel loci, then `neanderthal.py` computes archaic-allele fraction directly from the resulting VCF. Cleaner separation, one extra bcftools pass; does not touch the PGS score set. Recommend B if keeping the PGS score list clean matters.

Slots into `bin/run.sh` as **step 3b** (after scoring, before grade), guarded so it is best-effort (never fails a completed run), mirroring the `cache_performance.py` pattern.

### 2b. Reference panel + tool + container
- **Panel (`resources/neanderthal_panel.tsv`, built once, committed):** archaic-introgression-informative SNPs. Construction: intersect the public high-coverage Neanderthal genome (Prüfer 2017 Vindija/Altai VCFs, Max Planck) with 1000 Genomes to select sites where the archaic-derived allele is present in Neanderthal, near-absent in African (YRI/LWK) samples, and present in non-Africans — i.e. Sankararaman 2014 / Vernot & Akey 2016 (`S*`) / Browning 2018 (Sprime) introgression tag SNPs, or Chen 2020 (IBDmix) segments. Store `chrom, pos, ref, archaic_allele, weight(=1)`; target a few-thousand high-|Δaf| tag SNPs. Format also emitted as a `_hmPOS_GRCh38.txt.gz` scorefile for Option A.
- **Calibration:** compute the panel-dosage→% mapping once by scoring 1kGP super-populations (or reuse published per-population Neanderthal fractions: EUR≈1.8–2.4%, EAS≈2.3–2.6%, AFR≈0.3%). Store slope/intercept in `resources/neanderthal_calibration.json`.
- **Tool/container:** none new — bcftools (already `quay.io/biocontainers/bcftools:1.21`) for genotyping; plain Python for `neanderthal.py`. No archaic reference genome ships at runtime (panel is pre-built).

### 2c. Contract output + schema
`results/launch70/<S>/score/neanderthal.tsv`:
```
sample	neanderthal_pct	method
HG002	1.94	dosage-over-archaic-tag-panel-v1 (calibrated to 1kGP)
```
Required columns exactly as OmniGen expects: `sample, neanderthal_pct, method`. Optionally append `n_snps_genotyped`, `raw_archaic_dosage` for provenance (OmniGen ignores extras). Written next to `pgs_scores.tsv` in the same score dir.

### 2d. Report integration
Add `report_html._deep_ancestry_section(results_dir)`: reads `neanderthal.tsv` (and `haplogroups.tsv` if present, §3) and renders a compact "Deep ancestry" card near the ancestry summary. Non-fatal if the file is absent (the section simply omits). Small addition to the `sections` assembly in `render()`.

### 2e. Testing
- **In-session:** panel construction + calibration math is fully testable — feed `neanderthal.py` a synthetic genotype vector with a known archaic-allele fraction and assert the % output; unit-test the calibration transform. Validate the panel scorefile parses through `scoring_targets.py`.
- **Needs a full run:** real per-sample % requires the mpileup pass over the BAM at panel loci (Option A folds into the scoring run; Option B is one extra bcftools pass). Sanity target: HG001/NA12878 (CEU) should land ≈1.8–2.2%.

---

## 3. Routing assessment — mtDNA (haplogrep2) + Y (yhaplo): poly vs pgx

**Recommendation: put BOTH uniparental haplogroups in pgx-suite, NOT poly-suite.** (This is a decision the user makes at approval; the poly-suite contract is sketched below as the fallback.)

### Why pgx-suite wins
1. **pgx-suite already runs a mitochondrial caller.** `mutserve` (baked-in jar) runs per sample for MT-RNR1, wired through `docker/genes.tsv` + `docker/Snakefile` with a status-sentinel rule and a `parse_mutserve` parser (`docker/test_parsers.py`). Extending it to whole-chrM calls + haplogrep2/haplocheck classification is a small delta on infra that already touches chrM. poly-suite has **no** mito caller.
2. **poly-suite's "already mpileups the whole BAM" advantage is illusory.** `scoring_targets.py` is hard-restricted to autosomes (`AUTOSOMES = chr1..chr22`, explicit `exclude chrX/Y/M`); the prepped VCF has **0 chrM and 0 chrY** records. Adding uniparental calling in poly means a **brand-new mpileup pass over a brand-new target set** (all of chrM + Y-informative SNPs) that reuses none of the PGS-loci prep. poly gains nothing from its existing genotype-prep.
3. **pgx-suite's architecture is the pattern-matched home.** Its `genes.tsv`-driven Snakefile adds a BAM→caller as a new rule exactly the way `HLA→optitype` and `MT-RNR1→mutserve` are wired (non-fatal status sentinels, `--cores` scheduling). A `chrY→yhaplo` rule and a `chrM→haplogrep2` rule drop in naturally.
4. **Keep both maternal (mtDNA) and paternal (Y) lineages in one module.** pgx already owns half of it. Splitting mtDNA→pgx and Y→poly is the worst outcome.
5. **Thematic adjacency to ancestry lives in OmniGen's rendering, not the producer.** OmniGen already merges contracts from both pipelines (`report_ancestry.py`, `report_pgs.py`); it renders haplogroups next to FRAPOSA ancestry regardless of which pipeline produced them. Co-location with FRAPOSA in poly is not required.

### Tradeoffs / counter-arguments (honest)
- **Against pgx:** haplogroups aren't pharmacogenomics — mild thematic scope creep for a PGx pipeline. Y requires a new tool (**yhaplo**, Python) + container; mtDNA requires adding **haplogrep2** (Java jar) or **haplocheck** on top of mutserve's chrM output. Y must be sex-gated (male-only; use pgx's existing sex handling).
- **For poly (the only real pull):** poly owns FRAPOSA ancestry, and haplogroups are deep-ancestry — a tidy single "ancestry & deep-ancestry" producer. But this is outweighed by items 1–2 (production cost + zero reuse in poly).

### If the user chooses poly-suite anyway — contract sketch
`results/launch70/<S>/score/haplogroups.tsv`:
```
sample	mt_haplogroup	y_haplogroup
HG002	H1	R1b1a2
```
Optional provenance columns: `mt_quality` (haplogrep2 quality score), `y_confidence`, `y_call_rate`; `y_haplogroup = NA` when sex≠male. Implementation in poly: new `bin/haplogroups.sh` after `genotype_prep` as a parallel branch off the same BAM — (i) mpileup+call whole chrM → VCF → haplogrep2; (ii) mpileup Y-informative SNP panel → yhaplo (gated by `infer_sex.py` = male). New containers: haplogrep2 + yhaplo. Rendered by the same `_deep_ancestry_section` helper (§2d). Note this re-implements in poly the chrM calling that pgx already does — reinforcing the recommendation against it.

### Testing (either home)
- **In-session:** classifier parsing is testable with canned haplogrep2 / yhaplo output (pgx already does this in `test_parsers.py::parse_mutserve` — mirror it). Y sex-gating logic testable with a stubbed sex call.
- **Needs a full run:** chrM/chrY mpileup over a real BAM. Sanity: HG002 mt≈H, Y≈J/R depending on reference truth; NA12878 (female) → y_haplogroup=NA.

---

## Summary of files touched (when implemented — not now)
- **PGS:** `bin/select_pgs.py` (LAUNCH_SET rows), optional pins in `results/pgs_catalog_meta.json`, `bin/report_html.py::_sys` icons. No schema change.
- **Neanderthal:** new `resources/neanderthal_panel.tsv` + `resources/neanderthal_calibration.json`, new `bin/neanderthal.py` (+ optional `bin/neanderthal.sh`), `bin/run.sh` step 3b, `bin/report_html.py::_deep_ancestry_section`. New contract `score/neanderthal.tsv`.
- **Haplogroups:** recommended in **pgx-suite** (`docker/genes.tsv` + `Snakefile` new rules, haplogrep2 + yhaplo containers, parser). Fallback poly-suite: `bin/haplogroups.sh` + contract `score/haplogroups.tsv` + report section.

---

## Implementation status (2026-07-11)

**Shipped in the working tree (not committed — user handles git):**
- `bin/select_pgs.py` — added 4 launch-set traits (`chronotype` extended; `educational attainment`, `neuroticism`, `loneliness` gated). Added `PINNED` (6 verified/substitute IDs, forced to rank 1 in `main()` with the rationale recorded in `pgs_catalog_meta.json`) and `DEFERRED` (beat-sync, risk-tolerance — documented, NOT fabricated).
- `bin/neanderthal.py` — new module: seed-panel loader, pure `pct_from_counts` calibration math, `archaic_counts_from_vcf` (pysam), and the `sample/neanderthal_pct/method` contract writer/reader + CLI.
- `resources/neanderthal_panel.tsv` (SEED + construction recipe) and `resources/neanderthal_calibration.json` (PROVISIONAL slope/intercept).
- `bin/report_html.py` — new `_deep_ancestry_section()` (Neanderthal card), new quantitative traits added to `BIOMARKERS` (predisposition framing) + icon keywords.
- `bin/run.sh` — best-effort, non-fatal Neanderthal step before grading (so the render picks up the card).
- `tests/tests.py` — +7 tests (launch-set wiring, pins, deferred, no-schema-change on a novel PGS, Neanderthal math, contract round-trip, seed-panel load). **23/23 pass in-session** (`tests/selftest.sh` also green).

**No schema change:** verified in-session — a novel PGS grades into the same 24-column `pgs_scores.tsv`; quantitative traits leave the absolute-risk columns `NA`.

**Ran in-session (deterministic):** unit tests, selftest, and a synthetic end-to-end (`grade_pgs.py` + `report_html.py`) confirming the new pinned traits render and the provisional Neanderthal card appears.

**Needs a full run (real scores — download scorefiles + pgsc_calc + FRAPOSA; ~hours):**
```bash
# 1. resolve the launch set incl. the new traits + pinned IDs (writes pgs_catalog_meta.json)
python3 bin/select_pgs.py results/pgs_catalog_meta.json 2 all
# 2. score a sample end-to-end (adds the pinned IDs to --pgs; --pgs-meta carries the pins)
bin/run.sh --sample HG001 --pgs "PGS002804,PGS000657,PGS002231,PGS002209,PGS002213,PGS001091,<existing ids>" \
  --outdir results/launch70/HG001 \
  --bam /data/alvin/ref/GIAB/HG001.bwa.sortdup.bqsr.bam \
  --bootstrap-vcf data/HG001_GRCh38_benchmark.vcf.gz \
  --panel /data/alvin/ref/pgsc/pgsc_HGDP+1kGP_v1.tar.zst \
  --pgs-meta results/pgs_catalog_meta.json --scorefile-cache cache/scorefiles --work-cache work
```
**Real Neanderthal % (separate, needs the panel force-genotyped off the BAM — the PGS-only prepped VCF does not cover archaic loci):** build the production panel per the `resources/neanderthal_panel.tsv` recipe, re-fit `neanderthal_calibration.json` against 1000G, then force-genotype and score:
```bash
# force-genotype the archaic panel (panel expressed as a scorefile), then estimate:
bin/genotype_prep.sh <bam> <ref> results/launch70/HG001/neanderthal.panel.vcf.gz <panel_scorefile.txt.gz>
python3 bin/neanderthal.py results/launch70/HG001/score --targets results/launch70/HG001/neanderthal.panel.vcf.gz --sample HG001
```
Until the seed panel is replaced, `neanderthal.tsv`/the report card are flagged **PROVISIONAL-SEED — not for interpretation**.
