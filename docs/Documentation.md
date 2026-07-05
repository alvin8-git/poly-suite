# poly-suite — Documentation (the 101)

This is the "why and how" for poly-suite: the science of polygenic scores, the
decisions the pipeline makes, and how the pieces fit. If you want to run it, start
with the [README](../README.md). If you want to understand it, read this.

---

## 1. What poly-suite is, in one paragraph

poly-suite turns a raw whole-genome sample (a BAM or VCF) into a graded,
ancestry-calibrated, absolute-risk polygenic report across a curated set of
traits. It does **not** reinvent scoring: [`pgsc_calc`](https://github.com/pgscatalog/pgsc_calc)
is the scoring core (plink2 scoring + PGS harmonization + FRAPOSA ancestry).
poly-suite is the pipeline **around** that core: the parts raw PGS tools leave to
the user. Input-robustness, evidence-graded score selection, ancestry-validity
gating, absolute-risk conversion, multi-score consensus, QC and honesty grading,
and a reproducible, referenced report.

The moat is the same as its sibling OmniGen: **robustness + interpretation +
honesty, not the math.**

---

## 2. The science 101: what a polygenic score actually is

A polygenic score (PGS, also PRS) is a weighted sum:

```
score = sum over scoring variants of ( dosage_i * weight_i )
```

`dosage_i` is how many copies of the effect allele you carry at variant i (0, 1,
or 2). `weight_i` is that variant's effect size, estimated by a genome-wide
association study (GWAS). The arithmetic is trivial. Every hard problem is in the
five steps **around** the sum:

1. **Get dosages** at the scoring variants (genotype or impute).
2. **Match** the scoring file to the sample (genome build, strand, effect-allele,
   rsID remapping).
3. **Sum** (plink2 `--score`).
4. **Ancestry-adjust.** Raw sums live on different scales across ancestries because
   allele frequencies differ. A raw sum is uninterpretable without a reference
   population to compare against.
5. **Interpret.** Percentile then relative risk then absolute risk then "so what."

Tools like pgsc_calc, PRSice, and LDpred2 own step 3 and nibble at 2 and 4. Nobody
owns the honest version of 4 and 5 for a real consumer sample. That gap is what
poly-suite fills.

### Why ancestry is the central problem

Most published scores are trained on European-ancestry GWAS (UK Biobank especially).
They lose 2 to 5 times accuracy in African-ancestry individuals, and less but still
meaningfully in East and South Asian ancestries. For admixed individuals (common in
the real consumer market) it is worse. A percentile computed against the wrong
reference distribution is not just imprecise, it is misleading.

**The distinction that matters most: calibration is not accuracy.** poly-suite (via
pgsc_calc's FRAPOSA) normalizes your percentile so it is *interpretable* for your
ancestry. But a European-trained score is still *less accurate* in a non-European
person, because the effect sizes and linkage patterns do not transfer. Fixing the
scale does not fix the prediction. That is why the pipeline emits a
`portability_flag` and a "reduced accuracy outside training ancestry" caveat rather
than pretending a score works equally everywhere. See §6.

### FRAPOSA, explained

FRAPOSA ("Fast and Robust Ancestry Prediction by using Online SVD and Shrinkage
Adjustment") is the ancestry engine `pgsc_calc --run_ancestry` calls. Its job:
place your one sample in genetic-ancestry space so the score can be normalized
against the right population.

```
reference panel (HGDP + 1000 Genomes, ~3,300 people, precomputed PCA)
        │
        │  1. PROJECT your sample onto the panel's principal components
        │     (online SVD — no need to recompute PCA on panel + you)
        ▼
   raw projected PCs  ──▶ biased toward the origin (a new sample that did not
        │                 contribute to the PCA lands too central)
        │  2. SHRINKAGE ADJUSTMENT: analytically correct that bias so your PCs
        │     sit on the same scale as the reference samples
        ▼
   calibrated PCs
        │
        ├─▶ 3a. assign MOST-SIMILAR POPULATION (nearest reference cluster)
        └─▶ 3b. NORMALIZE the score: regress PGS on the PCs across the panel,
                standardize your score against that model → percentile
```

The shrinkage step is the load-bearing trick. Without it, projected samples cluster
too tightly near the center and ancestry is mis-assigned. The `Z_norm1` (mean-
adjusted) and `Z_norm2` (mean + variance-adjusted) columns in the output are this
PC-based normalization. Its correctness is exactly what
[`tests/validate_calibration.py`](../tests/validate_calibration.py) checks:
reference-panel percentiles must come out uniform on 0 to 100, and they do
(mean 50.1, max decile deviation 0.001 on real data).

---

## 3. Architecture: the pipeline

```
 FASTQ ─┐
 BAM  ──┼─▶ (align once, if FASTQ)             one canonical BAM
 VCF  ──┘                                            │
                                                     ▼
  [1] GENOTYPE PREP  (scoring_targets.py + genotype_prep.sh)
      force-genotype the BAM at the scoring loci so hom-ref sites become explicit
      0/0 (dosage 0) instead of being absent. Fixes the variant-only-VCF coverage gap.
                                                     │
                                                     ▼
  [2] SCORING CORE = pgsc_calc  (orchestrated, not reinvented)
      harmonize scores to GRCh38 → plink2 --score → raw sums + match rate
                                                     │
                                                     ▼
  [3] ANCESTRY = FRAPOSA  (pgsc_calc --run_ancestry)
      project onto HGDP+1kGP PCs → most-similar-pop + ancestry-normalized percentile
                                                     │
                                                     ▼
  [4] SELECTION + GRADING + INTERPRETATION  (the poly-suite layer)
      select_pgs (evidence-graded, tiered) · grade_pgs (QC/ancestry/grade gates)
      absolute_risk (odds-scale) · consensus (multi-score robustness) · infer_sex
                                                     │
                                                     ▼
  [5] OUTPUT CONTRACT + REPORT
      pgs_scores.tsv/.json + provenance.json + report.html (patient-first, grouped by confidence)
```

`bin/run.sh` is the single orchestrator that composes all of this.

### Producer / consumer boundary

poly-suite is a **standalone deliverable** (sibling to pgx-suite and SVcaller). It
produces `results/.../pgs_scores.tsv` (a stable 24-column contract). OmniGen, the
wider genomics-report project, is *only a consumer* of that contract. Nothing PGS
is computed in OmniGen except final rendering. This keeps poly-suite independently
testable and runnable.

---

## 4. Key decisions (the "why")

**Wrap pgsc_calc, do not reinvent scoring.** Claiming novel scoring would be
indefensible and duplicate a mature, validated tool. pgsc_calc is the scoring core
the same way pgx-suite wraps Cyrius and SVcaller wraps Manta. The value-add is the
pipeline around it.

**Force-genotype at scoring loci (the variant-only-VCF fix).** A GIAB benchmark VCF
lists only sites where the sample differs from the reference. So a scoring variant
where the sample is hom-ref is *absent*, not dosage 0. pgsc_calc then treats it as
missing and mean-imputes, biasing the score, and coverage collapses (CAD dropped to
27.5%). The fix is to force-call the BAM at the scoring alleles
(`bcftools mpileup | call -A -C alleles`) so hom-ref becomes an explicit 0/0. This
lifted CAD coverage 27.5% → **99.8%**. The `-A` (keep-alts) flag is load-bearing:
without it, hom-ref sites drop their ALT allele and coverage gets *worse* — a bug
only visible by re-scoring, not by eyeballing the genotypes.

**Autosomes only.** chrX needs sex-aware (haploid-male) dosage that pgsc_calc/plink2
reject without sample sex, and chrX is a small fraction of scoring variants. A
documented limitation, not an oversight.

**The evidence grade is the data-sufficiency filter.** Rather than hard-coding which
traits have "enough data," `select_pgs` grades every candidate score from PGS
Catalog metadata (GWAS sample size, ancestry evaluation, method). A trait with a
large multi-ancestry GWAS lands grade A; a thin one (AMD's 6 scores) self-downgrades
to C and renders with a loud caveat. You include broadly and let the grade disclose
confidence.

**Absolute risk on the odds scale, not linear.** `percentile → Z → relative risk
(via published per-SD OR/HR) → absolute risk`. The relative risk is calibrated so
the population *mean* RR = 1 (a median-percentile person sits slightly below
baseline; the high-PGS tail carries the excess). Critically, the RR is applied on
the **odds scale**, not by multiplying the baseline: linear `baseline × RR`
overshoots for common diseases (a 2× RR on a 49% CAD baseline gave an impossible
99.9%; the odds-scale model gives 66.8%).

**The report is patient-first, with a clinical layer underneath.** The same contract
renders for a lay reader: traits grouped by evidence confidence (a high-percentile
grade-D score sinks to the bottom instead of leading), each a compact row (organ icon +
likelihood bar + plain verdict) that expands to the detail. Absolute risk is shown as a
natural frequency ("about 1 in N", sex-labelled) *only* for grade A/B traits with a
verified effect + baseline; a grade C/D "elevated" score gets an explicit refusal, not a
number. The clinical layer links out to the PGS Catalog score, the MONDO disease class,
and the source study.

**Evidence grade is not predictive accuracy.** The A–D grade rates the *strength of the
evidence* (GWAS size + ancestry evaluation, from Catalog metadata — poly-suite's own
call, not the Catalog's). It is distinct from *discrimination*: the report shows each
score's publication-reported AUROC / C-index / R² next to its grade so the two aren't
conflated (a grade-A score can still have a modest AUROC).

**Raw caller output is never a reportable finding.** Every tier has a gate. See §6.

**Bash orchestrator, not a Nextflow rewrite (for now).** pgsc_calc is already a
Nextflow pipeline. Building a second one to orchestrate it is symmetry-driven
over-engineering until the genotype-prep DAG (real imputation) justifies it. The
value-adds are a thin prep front-end plus a Python post-processor. `run.sh` composes
the proven stages; a Nextflow scaffold is the productionization step, deferred.

---

## 5. The honesty layer (per-trait gates)

The whole point: **a raw PGS number is not a reportable finding.** `grade_pgs.py`
applies gates before anything is shown:

| Gate | Fires when | Effect |
|---|---|---|
| Coverage / match-rate | matched scoring variants below the overlap floor | downgrade to D, "not a valid PRS" caveat, suppress headline |
| Calibration | no ancestry percentile available | grade D, "uncalibrated — needs `--run_ancestry`" |
| Ancestry portability | sample outside the score's training ancestry | `portability_flag`, cap grade, "reduced accuracy" caveat |
| Absolute risk | percentile + baseline + effect available | convert (odds-scale, sex-resolved); else percentile-only, say so |
| Multi-score robustness | ≥2 scores per trait disagree (different risk tertiles) | "robustness LOW — treat as uncertain" |
| Evidence grade | always | A (large multi-ancestry GWAS) … D (thin/unreplicated) |
| Controlled vocabulary | always | `allowed_statement`: no "you will get X"; "explains ~Y% of variance; not diagnostic" |

Every claim traces to a source, a grade, and a confidence. The report also surfaces each
score's publication-reported discrimination (AUROC / C-index / R²) beside its grade, so
evidence *strength* and predictive *accuracy* are never conflated.

---

## 6. The tiered launch set

Traits are organized by tier so the default view is confident and actionable, and
sensitive/thin traits are opt-in. `select_pgs.py OUT TOP_N {core|extended|all|gated}`.

| Tier | Count | Posture |
|---|---|---|
| **core** | 25 | grade-A-capable, actionable (headline) |
| **extended** | 30 | grade ≥ B, informative |
| **= default suite** | **55** | core + extended |
| **gated** | 6 | opt-in (schizophrenia, depression, bipolar, ADHD, autism, cognitive) |

On a real 55-trait selection: 83 of 108 scores graded A, 11 B, and the 14 C/D are
exactly the thin-evidence traits (AMD, ankylosing spondylitis, coeliac, migraine…),
auto-caveated rather than hidden or hard-excluded.

At `TOP_N=2` (two scores per trait for consensus), the full 61-trait set resolves to
**119 PGS scorefiles** (some thin traits have only one usable score). The full-set run
force-genotypes at the union of their loci (**12.6M autosomal sites** on HG001).

### The full trait list

Every trait in `LAUNCH_SET` (`bin/select_pgs.py`), with the PGS Catalog trait ontology
ID each is resolved against. Grades are per-score and data-driven (§4), so a trait's
tier is its *posture*, not a fixed grade.

**Core — 25 (grade-A-capable, headline):**

| Trait | Ontology ID | | Trait | Ontology ID |
|---|---|---|---|---|
| coronary artery disease | MONDO_0005010 | | age-related macular degeneration | MONDO_0005150 |
| type 2 diabetes | MONDO_0005148 | | osteoporosis | MONDO_0005298 |
| atrial fibrillation | MONDO_0004981 | | gout | MONDO_0005393 |
| ischemic stroke | HP_0002140 | | type 1 diabetes | MONDO_0005147 |
| LDL cholesterol | EFO_0004611 | | coeliac disease | MONDO_0005130 |
| HDL cholesterol | EFO_0004612 | | hypothyroidism | MONDO_0005420 |
| triglycerides | EFO_0004530 | | breast cancer | MONDO_0004989 |
| lipoprotein(a) | EFO_0006925 | | prostate cancer | MONDO_0005159 |
| body mass index | EFO_0004340 | | colorectal cancer | MONDO_0005575 |
| chronic kidney disease | MONDO_0005300 | | melanoma | MONDO_0005105 |
| heart failure | MONDO_0005252 | | lung cancer | MONDO_0005138 |
| venous thromboembolism | MONDO_0005399 | | glaucoma | MONDO_0005041 |
| hypertension | MONDO_0001134 | | | |

**Extended — 30 (grade ≥ B, informative):**

| Trait | Ontology ID | | Trait | Ontology ID |
|---|---|---|---|---|
| rheumatoid arthritis | MONDO_0008383 | | thyroid cancer | MONDO_0015075 |
| inflammatory bowel disease | MONDO_0005265 | | endometrial cancer | MONDO_0002447 |
| psoriasis | MONDO_0005083 | | basal cell carcinoma | MONDO_0020804 |
| multiple sclerosis | MONDO_0005301 | | osteoarthritis | MONDO_0006629 |
| systemic lupus erythematosus | MONDO_0007915 | | hyperthyroidism | MONDO_0010138 |
| ankylosing spondylitis | MONDO_0005306 | | multiple myeloma | MONDO_0009693 |
| asthma | MONDO_0004979 | | cataract | MONDO_0005129 |
| chronic obstructive pulmonary disease | MONDO_0005002 | | fatty liver disease | MONDO_0013209 |
| Parkinson disease | MONDO_0005180 | | gallstones | MONDO_0012672 |
| migraine | MONDO_0005277 | | abdominal aortic aneurysm | MONDO_0005350 |
| Alzheimer disease | MONDO_0004975 | | C-reactive protein | EFO_0004458 |
| atopic dermatitis | MONDO_0004980 | | HbA1c | EFO_0004541 |
| ovarian cancer | MONDO_0005140 | | height | OBA_VT0001253 |
| pancreatic cancer | MONDO_0005192 | | waist-hip ratio | EFO_0004343 |
| bladder cancer | MONDO_0004986 | | | |
| testicular cancer | MONDO_0005447 | | | |

**Next tranche (added 2026-07; see [candidate-traits.md](candidate-traits.md)):** peripheral
arterial disease (`MONDO_0005386`), endometriosis (`MONDO_0005133`), obstructive sleep apnea
(`MONDO_0005296`), esophageal cancer (`MONDO_0007576`), gastric cancer (`MONDO_0001056`),
cervical cancer (`MONDO_0005131`) — plus three **biomarkers** rendered in a *predisposition*
mode (direction + percentile, no "1 in N"): bone mineral density (`EFO_0009270`), estimated
glomerular filtration rate (`OBA_0003747`), systolic blood pressure (`EFO_0006335`).
Absolute-risk inputs were curated for PAD / endometriosis / esophageal / gastric; OSA and
cervical stay percentile-only (no standard per-SD effect; cervical is HPV-driven, with an
on-card note). Polycystic ovary syndrome was later **dropped** (no Catalog PGS under any
ontology id) and replaced with **multiple myeloma** (`MONDO_0009693`, grade A, SEER-baselined)
to fill the hematologic-cancer gap.

**Gated — 6 (opt-in; sensitive / behavioural — off by default):**

| Trait | Ontology ID |
|---|---|
| schizophrenia | MONDO_0005090 |
| major depressive disorder | MONDO_0002009 |
| bipolar disorder | MONDO_0004985 |
| ADHD | MONDO_0007743 |
| autism spectrum disorder | MONDO_0005258 |
| intelligence | EFO_0004337 |

### On ancestry sources beyond European

Non-European GWAS exist and are growing: Biobank Japan (East Asian), the Million
Veteran Program (~20% African), All of Us (~50% non-European), H3Africa, and the
Global Biobank Meta-analysis Initiative (multi-ancestry meta-analysis). The PGS
Catalog carries some multi-ancestry and non-European scores, just far fewer. The
right response, and a planned enhancement, is **ancestry-aware selection**: use the
sample's inferred ancestry (already available from FRAPOSA) to prefer ancestry-
matched or multi-ancestry-trained scores, instead of always the best-evidence
(usually European) one.

---

## 7. Performance and caching

The expensive stage is genotype-prep (force-genotyping millions of scattered loci,
seek-bound on the BAM). It is chunk-parallel: the coord-sorted targets are split into
`NPROC * 4` contiguous genomic blocks and pulled dynamically by `NPROC` workers
(`xargs -P`), so the deepest-region chunk — the mpileup straggler — is ~1/4 the size
and no worker sits idle at the tail. On the 61-trait set (12.6M loci, 60 workers) prep
runs in ~15 min at ~99.9% genotyped.

Three cheap prep-stage wins keep that number down:

- **Memoized target orientation** (`scoring_targets.py`): overlapping scorefiles share
  most loci, so each `(chrom,pos,allele)` is oriented against the FASTA once, not once
  per scorefile that references it. On 119 heavily-overlapping scorefiles this collapses
  tens of millions of FASTA fetches to the unique-locus count.
- **No unused FORMAT fields**: mpileup drops `FORMAT/DP,FORMAT/AD` — plink2 scores off
  `GT`, nothing downstream reads them, so computing them is pure waste.
- **ALT-allele cap** (`MAX_ALTS=100`): plink2's VCF import aborts above 254 ALTs, and at
  scale a single indel hotspot can accumulate hundreds of distinct ALTs across scorefiles.
  Capping well under the limit keeps the score pass from crashing (~10 of 12.6M sites).

Caching, in order of real impact:

- **`--reuse-prep`** (the real win): a cached full-union prepped VCF skips harmonize
  and the ~1.5–2 h genotype-prep. A re-run or a new trait whose loci are already
  prepped becomes score-only (~90–110s of compute).
- **`--scorefile-cache`**: locally cached PGS Catalog scoring files feed pgsc_calc
  `--scorefile`, skipping the ~40–58s per-run download — **and** on a BAM run with every
  score cached, the harmonize nextflow pass is skipped entirely (the cached hmPOS files
  go straight into genotype-prep), removing one of two pgsc_calc cold-starts.
- **`--work-cache`**: persistent nextflow work dir + `-resume` reuses the 16 GB
  panel extraction across calibrated runs.
- **`make_pgen.sh` / `format=pfile`**: marginal (plink2's VCF→pgen is already ~13s;
  kept for very large inputs, not a headline).

---

## 8. Repository map

```
bin/     pipeline (run once, orchestrated by run.sh)
  run.sh              end-to-end orchestrator
  select_pgs.py       evidence-graded, tiered PGS selection (PGS Catalog API)
  scoring_targets.py  build force-call targets (REF/ALT oriented, autosomes, SNP+indel)
  genotype_prep.sh    force-genotype the BAM at scoring loci (chunk-parallel)
  infer_sex.py        X/Y-coverage sex inference (for sex-dimorphic absolute risk)
  grade_pgs.py        QC / ancestry / grade gates -> the 24-column output contract
  absolute_risk.py    percentile -> odds-scale absolute risk (sex- and baseline-aware)
  consensus.py        multi-score per-trait concordance (robustness flag)
  ensemble.py         meta-PGS (linear combination of standardized scores)
  provenance.py       reproducibility record (versions, ref, sha256)
  report_html.py      self-contained patient-first HTML report (confidence tiers, links, perf)
  cache_performance.py  cache reported AUROC/C-index/R² per score (PGS Catalog perf API)
  make_pgen.sh · cache_scorefiles.py · fetch_panel.py · resolve_traits.py  (helpers/caching)
tests/   validation (no network/BAM/panel needed unless noted)
  tests.py            16 unit tests · selftest.sh  suite runner
  validate_contract.py     schema + gate-invariant checks on pgs_scores.tsv
  validate_calibration.py  calibration self-consistency (reference-panel uniformity)
conf/    rootless.config (Docker rootless fix)
resources/  pgs_effect.tsv, baseline_incidence.tsv (absolute-risk inputs, sourced)
            pgs_performance.tsv (reported AUROC/C-index/R² per score, PGS Catalog)
docs/    this file + pgs-pipeline-spec.md
```

---

## 9. Gotchas (learned on real hardware)

- **Rootless Docker + pgsc_calc = permission denied.** pgsc_calc passes
  `-u $(id -u):$(id -g)`; under rootless Docker container-root already maps to the
  host user, so forcing a UID breaks workdir writes. Fix: `conf/rootless.config`
  sets `docker.runOptions = ""`; always pass `-c conf/rootless.config`.
- **pgsc_calc samplesheet `path_prefix` must be absolute.** It is resolved relative
  to pgsc_calc's own base, not your cwd. `run.sh` runs every path through
  `readlink -f`. Two real bugs came from this (relative bootstrap VCF, relative outdir).
- **`--run_ancestry` output is named differently:** `*_pgs.txt.gz` (not
  `aggregated_scores.txt.gz`), it includes ~3,300 reference-panel rows, and the
  ancestry call is in `pop_summary.csv`. The grader handles all three.
- **chrX in the scoring set breaks plink2 without sample sex.** poly-suite restricts
  targets to autosomes.
- **plink2 VCF import caps at 254 ALT alleles per site.** Across many scorefiles one
  indel hotspot can pile up hundreds of distinct ALTs; the score pass aborts with
  "variant with N ALT alleles; limited to 254". `scoring_targets.py` caps at
  `MAX_ALTS=100`. (Surfaced only at 61-trait scale — 18 traits never hit it.)
- **Harmonized scorefiles come in two shapes.** Raw PGS-Catalog `*_hmPOS_GRCh38.txt.gz`
  carry a `#` comment header and put GRCh38 coords in `hm_chr`/`hm_pos` (the plain
  `chr_name`/`chr_position` are the *original* build), and some omit `other_allele`
  entirely (only `hm_inferOtherAllele`). pgsc_calc's own `normalised_*` files are
  headerless with GRCh38 already in `chr_name`/`chr_position`. `scoring_targets.py`
  reads both, so the harmonize-skip fast path (feeding cached raw files straight to
  genotype-prep) produces identical targets.
- **pgsc_calc `FORMAT_SCOREFILES` OOMs under parallel runs.** The step (`pgscatalog-format`)
  loads *all* selected scorefiles into memory at once; the large biomarker scores (systolic
  BP / eGFR / BMD, ~135 MB gzipped each) make it **~35 GB per run**. Running the 7 GIAB
  samples 4-wide OOM-killed three (exit 137). Cap concurrent runs at **≤2 on a 128 GB box**
  (serial is safest). It's a concurrency/memory limit, not a bug — HG001 solo and genotype-prep
  are fine. Note `--reuse-prep` does NOT avoid it: `FORMAT_SCOREFILES` is in the *score* step,
  so a reuse-prep re-run still hits it (reuse-prep only skips harmonize + prep). The failed
  samples' `*.prepped.vcf.gz` survive, so re-run them with `--reuse-prep … --pgs-meta …` at
  low concurrency.

---

## 10. Validation

- **Unit tests** (`tests/selftest.sh`): 16 tests over the pure logic (absolute-risk
  monotonicity + odds-scale bounding, sex precedence, consensus concordance, grade
  gates, launch-set integrity, contract schema).
- **Calibration self-consistency** (`tests/validate_calibration.py`): the reference
  panel's ancestry-normalized percentiles must be ~uniform 0–100, proving the FRAPOSA
  normalization is sound. Verified on real data (all scores mean 50.1, max decile
  deviation 0.001).
- **Contract validation** (`tests/validate_contract.py`): the 24-column output plus
  gate invariants (grade range, no absolute-risk-without-a-percentile, low-coverage
  caveat present).

---

## 11. Known limitations

- **European-ancestry bias** in the underlying scores (a field-wide data problem, not
  a pipeline bug). poly-suite calibrates the percentile and flags portability, but it
  cannot manufacture accuracy that the training data lacks.
- **chrX excluded** (sex-aware dosage deferred).
- **Array / low-pass input** needs an imputation branch (GLIMPSE/Beagle), not built;
  the WGS force-genotype path is the current focus.
- **Absolute-risk inputs are illustrative launch values** (per-SD effects and
  baseline incidences) that need per-score verification against source publications.
- **Reported performance metrics are publication-reported, not computed here.** The
  AUROC / C-index / R² shown come from each score's PGS Catalog performance records
  (median across mixed cohorts/ancestries, EUR-preferred), not evaluated on the sample.
- **Not a diagnostic test.** Educational / research use. Discuss anything actionable
  with a clinician or genetic counselor.
