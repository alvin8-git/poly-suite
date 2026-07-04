# poly-suite

![tests](https://img.shields.io/badge/tests-16%20passing-brightgreen)
![python](https://img.shields.io/badge/python-3.8%2B-blue)
![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A5%2025.10-brightgreen)
![scoring core](https://img.shields.io/badge/scoring%20core-pgsc__calc-informational)
![ancestry](https://img.shields.io/badge/ancestry-FRAPOSA-9cf)
![launch set](https://img.shields.io/badge/launch%20set-61%20traits-orange)
![license](https://img.shields.io/badge/license-MIT-green)

Standalone polygenic-score pipeline: a raw sample (WGS BAM/VCF) → a graded,
referenced, ancestry-calibrated PGS deliverable. Sibling to pgx-suite and SVcaller.
Independent of OmniGen — OmniGen only ever consumes this pipeline's output contract
(`results/pgs_scores.tsv`); nothing PGS-related is computed downstream.

Scoring core is **pgsc_calc** (plink2 scoring + PGS harmonization + FRAPOSA ancestry) —
not reinvented. The value-add is the pipeline around it: input-robustness, evidence-graded
score selection, ancestry-validity gating, absolute-risk conversion, QC + honesty grading,
reproducible/referenced output.

**Docs:** [Documentation.md](docs/Documentation.md) — the 101 (science, decisions, architecture) · [pgs-pipeline-spec.md](docs/pgs-pipeline-spec.md) — build spec.

**Example output:** [`examples/HG001-18trait/`](examples/HG001-18trait/) — the calibrated 18-trait card for GIAB HG001 (report.html + pgs_scores.tsv/.json + provenance).

## Layout

```
data/       samplesheet.csv + localized test input (GIAB HG001)
conf/       rootless.config — Docker rootless fix (see Gotchas)
resources/  pgs_effect.tsv, baseline_incidence.tsv — absolute-risk inputs (sourced)
bin/        run.sh            — end-to-end orchestrator (one entry point)
            select_pgs.py     — evidence-graded, TIERED PGS selection (PGS Catalog API)
            resolve_traits.py — resolve trait names -> Catalog ids + score counts
            make_pgen.sh      — pre-convert cached prep to plink2 .pgen (pfile fast path)
            cache_scorefiles.py — cache PGS Catalog harmonized scorefiles (--scorefile, skip download)
            scoring_targets.py + genotype_prep.sh — force-genotype at scoring loci
            infer_sex.py      — X/Y-coverage sex inference
            grade_pgs.py      — QC/ancestry/grade gates -> the output contract
            absolute_risk.py · consensus.py · ensemble.py — actionable + robustness + meta-PGS
            provenance.py · report_html.py  — reproducibility + standalone HTML
            score_selected.sh — score the auto-selected set (calls the stages)
tests/      tests.py + selftest.sh — 16 unit tests + self-check runner
            validate_contract.py — schema + gate-invariant validation of the output
            validate_calibration.py — calibration self-consistency (reference-panel uniformity)
results/    pgsc_calc output + deliverable: pgs_scores.tsv/.json + provenance.json + report.html
work/       Nextflow work dir
docs/       Documentation.md (the 101) + pgs-pipeline-spec.md
```

## Run

One command (orchestrator) — BAM → prep → score → ancestry → grade → report:
```bash
bin/run.sh --sample HG001 --pgs "PGS000018,PGS000004" --outdir results/hg001 \
  --bam /data/alvin/ref/GIAB/HG001.bwa.sortdup.bqsr.bam \
  --bootstrap-vcf data/HG001_GRCh38_benchmark.vcf.gz \
  --panel /data/alvin/ref/pgsc/pgsc_HGDP+1kGP_v1.tar.zst \    # optional; omit = uncalibrated
  --scorefile-cache cache/scorefiles --work-cache work        # optional; reuse downloads + panel
# add --dry-run to print the plan. Or run the stages manually:
```

Full launch set (61 traits, 119 scorefiles at top-2): resolve with
`python3 bin/select_pgs.py results/meta.json 2 all`, then pass the resolved ids as `--pgs`.
See [Documentation.md §6](docs/Documentation.md#6-the-tiered-launch-set) for the trait list.

Manual stages:
```bash
# 0. (BAM input) genotype-prep: force-genotype at scoring loci -> full-genotype VCF.
#    Needs the harmonized scorefiles pgsc_calc produces (work/**/normalised_*.txt.gz).
NPROC=32 bin/genotype_prep.sh \
  /data/alvin/ref/GIAB/HG001.bwa.sortdup.bqsr.bam \
  /data/alvin/ref/GRCh38/hg38.canonical.fa \
  results/prep/HG001.vcf.gz  work/**/normalised_PGS000018_hmPOS_GRCh38.txt.gz ...
#    then point the samplesheet at results/prep/HG001 (path_prefix, format vcf)

# 1. score (pgsc_calc, containerized). Add --run_ancestry <panel> for calibrated percentiles.
nextflow run pgscatalog/pgsc_calc -profile docker \
  --input data/samplesheet.csv --target_build GRCh38 \
  --pgs_id PGS000018,PGS000004 --min_overlap 0.1 \
  -c conf/rootless.config --outdir results -work-dir work

# 2. grade + render the contract
python3 bin/grade_pgs.py results
```

Note the ordering wrinkle: genotype-prep needs the *harmonized* scorefiles, which
pgsc_calc only writes on a first pass. So the real flow is score-once (discover +
harmonize scorefiles) → genotype-prep at those loci → score again on the prepped VCF.
A future Nextflow scaffold folds this into one DAG.

`grade_pgs.py` reads `results/score/aggregated_scores.txt.gz` and writes
`results/pgs_scores.tsv` (+ `.json`) — the stable schema OmniGen consumes — applying
the QC gate (coverage/match-rate), the ancestry/portability gate, evidence grade, and
controlled-vocab caveats.

## Gotchas (learned on this box)

- **Rootless Docker + pgsc_calc = permission denied.** pgsc_calc's config passes
  `-u $(id -u):$(id -g)`; under rootless Docker container-root already maps to the host
  user, so forcing a UID breaks workdir writes (`.command.trace: Permission denied`).
  Fix: `conf/rootless.config` sets `docker.runOptions = ""`. Always pass `-c conf/rootless.config`.
- **`curl`/`wget` blocked on this box.** pgsc_calc downloads scorefiles inside its
  containers (unaffected); manual panel downloads need `python3`+requests, not wget.
- **Variant-only VCFs cannot support a valid PRS.** A GIAB benchmark VCF lists only
  non-ref sites, so scoring loci where the sample is hom-ref are *missing*, not dosage 0.
  Coverage comes out low and the score is biased. The grader correctly grades these **D**
  and refuses a percentile. The genotype-prep stage (`bin/genotype_prep.sh`) fixes this by
  force-calling the BAM at the scoring alleles — CAD coverage 27.5% → 99.8% (below).

## Resources on this box (`/data/alvin/ref`)

- **`ref/GIAB/`** — analysis-ready **full BAMs** HG001–HG007 (`HG00N.bwa.sortdup.bqsr.bam`
  + `.bai`), plus HG002 truth/SV VCFs. These are the **full-genotype inputs** the
  genotype-prep stage needs: `HG001.bwa.sortdup.bqsr.bam` → call genome-wide (or
  force-call at scoring loci) → an accurate PRS, replacing the variant-only benchmark VCF
  that caps coverage at ~28% (below).
- **`ref/GRCh38/`** — hg38 FASTA, **bwa-mem2 indexed** (`.0123/.bwt.2bit.64/...`) +
  `GRCh38.dict` — everything alignment + calling needs. `hg38.alphabetical.fa` is the
  canonical reference.
- **`ref/pgsc/`** — the 16 GB HGDP+1kGP ancestry panel (`pgsc_HGDP+1kGP_v1.tar.zst`) for
  `--run_ancestry`; present and wired. `ref/vep_cache/`, `ref/annotsv/` — annotation.

## Validated (HG001, variant-only benchmark VCF)

End-to-end run proves the honesty layer on real data — and empirically confirms the
input-robustness gap:

| Score | Trait | Coverage | Grade | Why |
|---|---|---|---|---|
| PGS000018 | coronary artery disease | **27.5%** (481,004 / 1,745,179) | **D** | variant-only VCF → 1.26M scoring variants unmatched; uncalibrated |
| PGS000004 | breast cancer | 40.3% (126 / 313) | **D** | same |

The gate correctly **refuses a percentile** and returns grade D ("not a valid PRS")
instead of a falsely-precise number. Reproduce: `python3 bin/grade_pgs.py results`.

**Genotype-prep fixes the coverage cap.** Force-genotyping the `ref/GIAB` HG001 BAM at
the scoring alleles (`bin/genotype_prep.sh`) turns absent hom-ref loci into explicit
`0/0` (dosage-0) calls, then re-scoring lifts coverage:

| Score | Coverage: benchmark VCF | Coverage: after genotype-prep |
|---|---|---|
| PGS000004 (breast) | 40.3% (126/313) | **85.6% (268/313)** — v2 w/ indels (was 74% SNP-only) |
| PGS000018 (CAD) | 27.5% (481k/1.74M) | **99.8% (1,742,125/1,745,179)** |

Breast ceiling is ~85.6%: the remaining 31 (9.9%) are strand-ambiguous SNPs (A/T, C/G)
that pgsc_calc excludes by design — genotype-prep can't recover those, and shouldn't.

Remaining breast gap: 48 indels/multiallelics (v1 force-calls SNPs only) + 31
strand-ambiguous SNPs pgsc_calc excludes. Including indels is the next increment.

## Status

- [x] **Full 61-trait launch run (HG001)** — `select_pgs.py … all` → 119 scorefiles at top-2,
      force-genotyped at the **union of 12.6M autosomal loci** (99.9% genotyped), scored +
      ancestry-calibrated + graded. Runs on the fully-optimized path (scorefile-cache →
      harmonize-skip, reused panel extraction). See `results/launch61/`.
- [x] **Genotype-prep perf pass** — memoized target orientation (each locus oriented vs the
      FASTA once across overlapping scorefiles, not once per scorefile), over-split chunks
      (`NPROC*4`, dynamic `xargs -P` pull) to kill the mpileup straggler, dropped the unused
      `FORMAT/DP,FORMAT/AD`, and record-count from the index (no full-decompress).
- [x] **Scale fixes exposed at 61 traits** — `scoring_targets.py` reads both harmonized
      formats (raw `*_hmPOS_*` with `hm_chr`/`hm_pos` + `#` header, and `normalised_*`),
      tolerates a missing `other_allele` column (falls back to `hm_inferOtherAllele`), and
      caps ALTs at 100 (plink2's VCF import aborts above 254). `run.sh --scorefile-cache` now
      also skips the harmonize nextflow run on a full cache hit.
- [x] **Prep-once + caching speedup** (`run.sh --reuse-prep`, `--work-cache`) — reuse a cached
      full-union prepped VCF (skips harmonize + the multi-hour genotype-prep; a re-run/new-trait
      drops from hours to minutes) and share a nextflow work dir so the 16GB panel extraction is
      reused across runs instead of re-unpacked. **The real win is --reuse-prep** (skips the
      ~1.5-2h genotype-prep -> a re-run is score-only, ~90-110s compute). `bin/make_pgen.sh` +
      `format=pfile` is MARGINAL (measured: plink2 VCF->pgen is only ~13s, so pgen saves ~3s) —
      kept for large-input edge cases, not a headline. Next real lever: cache scorefile downloads
      (~40-58s/run): DONE — `bin/cache_scorefiles.py` + `run.sh --scorefile-cache DIR` feeds pgsc_calc
      local `--scorefile` (one-time fetch, fixed per PGS id) instead of re-downloading each run.
- [x] **Tiered launch set** (`select_pgs.py`, `bin/resolve_traits.py`) — 25 core (grade-A actionable)
      + 30 extended (grade>=B) = **55 default suite**, +6 gated (opt-in sensitive) = 61 traits.
      `select_pgs.py OUT TOP_N {core|extended|all|gated}`. IDs+counts API-resolved; grade is the
      final data-sufficiency filter (thin traits self-downgrade).

- [x] **Full top-2 calibrated run (16 scores / 8 traits, HG001)** — all 8 traits' two scores
      **concordant** (0.80–0.99), all calibrations self-consistent, contract validates. Notable: T1D 82/97th (elevated).
      Ran via `bin/run.sh` end-to-end (fixed: samplesheet path_prefix must be absolute — bootstrap, outdir).

- [x] Scoring core (pgsc_calc) wired, rootless fix, rooted in this repo
- [x] Grading + contract + report (`bin/grade_pgs.py`) — QC + ancestry/portability gates
- [x] **Genotype-prep stage** (`bin/genotype_prep.sh` + `scoring_targets.py`) — the
      input-robustness value-add: force-genotype a BAM at scoring loci (bcftools
      `mpileup | call -A -C alleles`, chunk-parallel), REF/ALT oriented vs the FASTA,
      hom-ref → `0/0`. Proven: breast coverage 40%→74% on the HG001 BAM.
- [x] **Genotype-prep v2: indels + multiallelics** — anchored-indel REF/ALT oriented vs
      the FASTA (longest ref-match wins), multiple ALTs grouped per site. Breast 74%→85.6%.
- [ ] Genotype-prep v3: imputation branch for array/low-pass input (GLIMPSE/Beagle)
- [x] **Evidence-graded PGS selection** (`bin/select_pgs.py`) — queries the PGS Catalog
      API per trait, grades candidates by GWAS N + ancestry eval + method, picks the
      best-evidence score with rationale → `results/pgs_catalog_meta.json`, loaded by
      `grade_pgs.py`. Launch set: 8 traits, grade-A auto-selections (CAD/T2D/breast/
      prostate/AF/T1D/Alzheimer/LDL). Catalog taxonomy quirk handled: scores tag specific
      MONDO terms (breast→`MONDO_0004989`, prostate→`MONDO_0005159`), not the obvious EFO.
- [x] **Absolute-risk conversion** (`bin/absolute_risk.py` + `resources/pgs_effect.tsv`,
      `baseline_incidence.tsv`) — percentile → Z → per-SD OR/HR → **odds-scale** absolute
      risk × baseline incidence (odds scale, not linear — linear overshoots >90% on common
      diseases). Wired into grade_pgs (fills absolute_risk/baseline_incidence/risk_ratio).
      Demoed: breast 87th-pct→18%, prostate 42nd-pct→9%. Sex-specific traits inferred
      (breast/prostate); sex-dimorphic (CAD) needs sample sex plumbed.
- [x] **Provenance + HTML report** (`bin/provenance.py`, `bin/report_html.py`) — grade_pgs
      now emits the full standalone deliverable: `pgs_scores.tsv`/`.json` + `provenance.json`
      (versions, ref, scores, calibrated flag, contract sha256) + self-contained `report.html`
      (grade badges, absolute risk, caveats, disclaimer). Captures nextflow/samtools/pgsc_calc versions.
- [x] **Multi-score consensus** (`bin/consensus.py`, `select_pgs.py --top-N`) — select_pgs now picks TOP-2 scores per trait (16 across 8), so grade_pgs computes per-trait concordance + a "robustness LOW" flag on tertile disagreement. Demonstrated on real HG001 CAD: PGS000018 16th vs PGS004941 43rd pct → flagged uncertain.
      scores → robustness flag; wired into the contract (`robustness_n_scores`/`_concordance`)
      + a "robustness LOW" caveat when scores disagree.
- [x] **End-to-end orchestrator** (`bin/run.sh`) — one entry: BAM→harmonize→prep→score→grade→report
      or VCF→score→grade→report; optional `--panel` (ancestry) and `--sex`/auto-infer. `--dry-run` plans.
- [x] **Test suite** (`tests/tests.py` + `tests/selftest.sh`) — 16 unit tests + module self-checks +
      grader end-to-end; no network/BAM/panel needed. `tests/selftest.sh` is green.
- [x] **Scored the auto-selected 8-score set** (`bin/score_selected.sh`) — harmonize → genotype-prep
      at 8.1M autosomal loci → score. Coverage: CAD/breast/T1D/LDL ~100%, T2D 93%, AF 85%, prostate 77%.
      Then calibrated: **8-trait card, all grade A, EUR** (T1D 82nd, AF 4th, LDL 8th, CAD 43rd...). chrX fix noted.
      autosome-only prep.)
- [x] **Sample-sex inference** (`bin/infer_sex.py`) — X/Y coverage via samtools idxstats;
      HG001→female (x=0.99, y=0.06). Wired into grade_pgs (`{results}/sample_sex.txt` or
      `POLY_SUITE_SEX`). CAD 90th-pct→~45% (female baseline); unknown sex → percentile-only, no guess.
- [x] **Ensemble / meta-PGS** (`bin/ensemble.py`) — score-level linear combination of
      standardized PRSs (the metaGRS approach): `meta_z = Σ w·z / Σ w`, evidence-weighted
      (GWAS N) or equal → a combined predictor per trait. The safest "novel score" route:
      validated inputs, cheap to validate, submittable to the Catalog. Real output needs
      ≥2 calibrated scores/trait; percentile flagged approximate pending reference re-standardization.
- [x] **Ensemble v2** (`ensemble.py:meta_from_calibrated`) — exact meta-PGS percentile by ranking
      the sample's meta_z against the reference panel (removes the v1 'approximate' caveat);
      activates on a run with >=2 scores/trait (the 16-score set).
- [ ] PRS-CSx multi-ancestry derivation (deferred, research)
- [x] **`--run_ancestry` calibration** (16 GB HGDP+1kGP panel) → real percentiles. First
      calibrated card (HG001, EUR): **CAD 15.7th pct, grade A, ~19% risk (female)**;
      **breast 31.2nd pct, grade C** (73% coverage). Grader adapted to the ancestry output
      schema (`*_pgs.txt.gz` + `pop_summary.csv`, panel rows filtered).
- [x] **Contract validation** (`tests/validate_contract.py`) — schema (24 cols) + gate invariants
      (grade range, portability boolean, no absolute-risk-without-percentile, low-coverage caveat
      present, ranges) — in `selftest.sh`.
- [x] **Validation harness** (`tests/validate_calibration.py`) — calibration self-consistency:
      the ~3,330 reference-panel samples' ancestry-normalized percentiles must be ~uniform 0-100.
      Verified on real data: all 8 scores mean 50.1, max decile dev 0.001 (FRAPOSA normalization sound).
- [ ] Nextflow scaffolding — productionization once the flow is validated (bash orchestrator covers it now)
