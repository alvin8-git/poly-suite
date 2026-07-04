# poly-suite — Standalone PGS Pipeline Spec

Sibling to pgx-suite and SVcaller: a self-contained Nextflow pipeline that turns a
sample (BAM or VCF) into a graded, referenced, ancestry-calibrated polygenic-risk
deliverable. OmniGen's PGS tier is *only* a consumer of this pipeline's output —
nothing PGS-related is computed inside OmniGen except final rendering.

## 0. Three hard requirements (these will be interrogated during the build)

1. **Self-contained deliverable.** Everything OmniGen needs for PGS is produced here
   (§2 output contract). OmniGen adds no computation — it reads the contract and renders,
   exactly as it consumes `all_genes_summary.tsv` (pgx-suite) and `*.filtered.tsv` (SVcaller).
2. **Runs independently + adds value over existing PGS tools.** Usable standalone
   (`nextflow run`), with its own HTML report + validation, AND it must beat "just run
   pgsc_calc" on defensible grounds (§1). If a step's only justification is "wrap the tool",
   it's cut.
3. **Robust output.** Input-format robustness, QC gates, ancestry-validity, absolute risk,
   provenance, and validation — so the deliverable to OmniGen is trustworthy, not a raw sum.

**Honest scope boundary (state this when interrogated):** the poly-suite does NOT reinvent
scoring or ancestry PCA — **pgsc_calc is the scoring core** (plink2 scoring + PGS harmonization
+ FRAPOSA ancestry), the same way pgx-suite wraps Cyrius/StellarPGx and SVcaller wraps Manta/
GATK. The value-add is the pipeline around the core. Claiming novel scoring would be
indefensible and violate OmniGen's "orchestrate best-in-class, don't reinvent" principle.

## 1. Value-add vs the PGS space (the defensibility case)

Existing tools are **scoring engines**: pgsc_calc (scoring + PCA + normalization), PRSice-2,
PLINK `--score`, LDpred2, the PGS Catalog web Calculator. What they leave to the user — and
what poly-suite delivers — is the pipeline from *raw sample* to *interpreted, graded finding*:

| Gap in raw tools | poly-suite value-add | Novel? |
|---|---|---|
| Assume you already have a full-genotype/imputed VCF; a BAM or variant-only/callset VCF breaks scoring (the exact bug we hit) | **Input-robustness / genotype-prep stage**: BAM → genotype at scoring + PCA loci (force-call), or impute a partial VCF, → scoring-ready set | Integration, not method |
| You pick PGS IDs blindly | **Evidence-graded PGS selection** from PGS Catalog metadata (sample size, ancestry-eval, reported performance) per trait | Yes — most tools don't |
| One score per trait, no robustness signal | **Multi-score consensus**: run ≥N published PGS per trait, report agreement/robustness (the pgx multi-caller idea, for PRS) | Yes |
| Ancestry info reported, but no decision | **Ancestry-validity gating**: suppress/flag when the sample's FRAPOSA distance from the PGS training ancestry is too large (portability) | Interpretation layer |
| Output stops at percentile | **Absolute-risk conversion**: percentile → absolute lifetime risk via baseline incidence + published OR/HR (ancestry-specific where available) | Yes — high user value |
| Output is a scores table | **QC + evidence grading + controlled caveats + provenance** (referenced, reproducible) | The OmniGen honesty layer |
| No self-validation | **Validation harness**: reproduce expected per-population percentile distributions (1000G/HGDP), regression on pgsc_calc's worked example | Defensibility |

One-line thesis: **existing PGS tools score; poly-suite delivers a graded, referenced,
ancestry-valid, absolute-risk PGS finding from a raw sample.** The moat is the same as
OmniGen's: robustness + interpretation + honesty, not the math.

## 2. Output contract — the deliverable OmniGen consumes (STABLE SCHEMA)

Per sample × per trait, one row (`pgs_scores.tsv` + `pgs_scores.json`), plus a standalone
HTML report and `provenance.json` (mirrors pgx-suite/SVcaller). Columns:

```
sample, trait, efo_id,
pgs_id, pgs_version, source_pmid, n_variants, training_ancestry,      # provenance
n_matched, match_rate, mean_imputed_frac,                              # QC (the min_overlap lesson)
inferred_ancestry, most_similar_pop, ancestry_distance,               # FRAPOSA
percentile, z_score, ci_low, ci_high,                                 # ancestry-adjusted
absolute_risk, baseline_incidence, risk_ratio,                        # actionable layer
robustness_n_scores, robustness_concordance,                          # multi-score
evidence_grade, portability_flag,                                     # grading
allowed_statement                                                     # controlled-vocab caveat
```

Gates baked into the contract (OmniGen enforces on render, but the pipeline SETS them):
- `match_rate` below threshold → `evidence_grade = D` + suppress from headline (QC gate).
- `ancestry_distance` beyond the PGS's eval range → `portability_flag = true`, cap grade, caveat.
- `absolute_risk` null where baseline incidence unavailable → report percentile only, say so.
- `allowed_statement` is controlled vocabulary (no "you will get X"; "explains ~Y% of variance;
  reduced accuracy outside training ancestry; not diagnostic").

## 3. Pipeline stages (Nextflow, `-profile docker`, like the siblings)

1. **Ingest & QC** — samplesheet (`bam` | `vcf`); coverage/callable QC; detect variant-only
   input and route to genotype-prep.
2. **Genotype preparation** (the input-robustness value-add) — BAM → force-genotype at the
   union of scoring-file + PCA loci (bcftools/GATK), or impute a partial VCF; emit a
   scoring-ready genotype set. This is the stage that fixes the variant-only-VCF failure.
3. **Scoring core = pgsc_calc** (orchestrated) — harmonize scores to target build, plink2
   scoring, FRAPOSA ancestry PCA + normalization → raw sums + ancestry-adjusted percentiles.
4. **Evidence selection + multi-score** — resolve traits → best-evidence PGS set (PGS Catalog
   API metadata); score ≥N per trait; compute robustness/concordance.
5. **Calibration & absolute risk** — percentile → absolute risk (baseline-incidence table +
   published effect); ancestry-validity gating.
6. **Grading & caveats** — evidence grade, QC gates (§2), controlled caveats.
7. **Report + provenance** — structured contract (§2) + HTML + `provenance.json` + validation summary.

## 4. Independent operation

```
nextflow run <poly-suite> -profile docker \
  --input samplesheet.csv --target_build GRCh38 \
  --traits "coronary_artery_disease,type_2_diabetes,breast_cancer" \  # or --pgs_id
  --run_ancestry <panel> --baseline_incidence resources/incidence.tsv \
  --outdir results_pgs
```
Standalone HTML + `pgs_scores.tsv`; usable with no OmniGen present. Ships validation
samplesheets (1000G/HGDP samples) as SVcaller ships GIAB ones.

## 5. Robustness & validation (defensibility)

- **Inputs handled:** BAM, full-genotype VCF, partial/benchmark VCF (flagged + imputed),
  array/imputed. No silent acceptance of an input that can't support a valid PRS.
- **QC gates:** match_rate, coverage, ancestry_distance — each can downgrade/suppress.
- **Validation:** (a) per-population percentile distributions in 1000G/HGDP match expectation;
  (b) multi-score concordance per trait; (c) regression on the pgsc_calc worked example;
  (d) score AUC/OR reproduces the source publication on a labelled cohort where available.
- **Reproducibility:** pinned PGS Catalog + panel versions, containerized tools, `provenance.json`,
  pipeline version — a report is reproducible to the release.

## 6. Build tasks (concrete, ordered)

1. Scoring core: wrap pgsc_calc (DONE — pulled rev 72ee54f). Prove end-to-end on a
   full-genotype VCF (not the variant-only benchmark — see tier-c-tooling.md §f).
2. Genotype-prep stage (BAM/partial-VCF → scoring-ready). **This is the flagship value-add
   and the fix for the input problem; build it first after the core.**
3. Evidence-selection + multi-score modules (PGS Catalog API).
4. Absolute-risk module + baseline-incidence resource table.
5. Grading + caveat engine → the §2 output contract.
6. Validation harness (1000G/HGDP).
7. OmniGen PGS-tier consumer: reads `pgs_scores.tsv` → report-card row (analogous to the
   pgx/SV tiers in `prototype/report_hg002.py`). This is the last mile into OmniGen.

## 7. Open dependencies to resolve during build

- **Full-genotype input for validation samples** (benchmark VCFs are variant-only) — genotype
  from BAM or use 1000G callsets.
- **Baseline-incidence data** per trait × ancestry × sex/age (for absolute risk) — source
  (e.g., SEER for cancers, national registries for CVD); a real data-gathering task.
- **16 GB ancestry panel** (`pgsc_HGDP+1kGP_v1.tar.zst`) — one-time download.
- **Which traits/PGS** are in the launch set, and the evidence-grade thresholds.
