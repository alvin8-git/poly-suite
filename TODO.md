# poly-suite — roadmap & changelog

Extracted from the README `## Status` section. Open work first, shipped features below.

## Open / roadmap

- [ ] Genotype-prep v3: imputation branch for array/low-pass input (GLIMPSE/Beagle)

- [ ] PRS-CSx multi-ancestry derivation (deferred, research)

- [ ] Nextflow scaffolding — productionization once the flow is validated (bash orchestrator covers it now)


## Shipped

- [x] **Next-tranche traits (61 → 70)** — added the [candidate-traits.md](docs/candidate-traits.md)
      recommended set to `select_pgs.py`: peripheral arterial disease, endometriosis, obstructive
      sleep apnea, esophageal / gastric / cervical cancer, and three **biomarkers** (bone mineral
      density, eGFR, systolic blood pressure). Biomarkers render in a new **predisposition mode**
      (direction + percentile, neutral colour, "what this measures" — no "1 in N" / "what to do",
      excluded from the disease "worth attention" chips). Absolute-risk inputs curated (sourced,
      PMIDs) for PAD / endometriosis / esophageal / gastric; OSA + cervical stay percentile-only
      (no standard per-SD effect; cervical is HPV-driven → on-card HPV note via a general
      `TRAIT_NOTES` mechanism). New icons wired (arterial→heart, apnea→lungs, bone→bone,
      glomerular→kidney, systolic→heart).

- [x] **Patient-first report redesign** (`bin/report_html.py`) — a lay-reader-first
      `report.html` that scales to 60+ traits. Traits grouped by evidence **confidence**
      (High/Good/Limited/Insufficient), so a high-percentile grade-D score no longer leads;
      each trait is a compact 2-line row (filled SVG **organ icon** + likelihood bar +
      plain-language verdict) that expands via a native `<details>` to a clinical layer.
      Notable-trait extras: absolute risk as a **natural frequency** ("about 1 in N" + a
      100-person dot array, sex-labelled "in women/men", gated to grade A/B *with* a verified
      effect — grade C/D elevated get a "why no number" refusal instead of false precision),
      "what this means / doesn't mean", a non-prescriptive next step, and a thyroid
      co-occurrence note. Expand/collapse-all (tiny inline JS; rows work without it). Embedded
      Space Grotesk heading font. Self-contained, CSS-only styling.

- [x] **Reported performance + evidence-grade legend** (`bin/cache_performance.py`,
      `resources/pgs_performance.tsv`) — each score's publication-reported discrimination
      (AUROC → C-index → R², EUR-preferred, median + range; physically out-of-range values
      dropped) shown beside its grade in the clinical table, plus a legend spelling out A–D
      (thresholds + downgrades) and that the grade is poly-suite's evidence-strength call, NOT
      the same as accuracy. `run.sh` auto-fetches after grading (`timeout 180`, best-effort:
      offline/slow Catalog just leaves "—") to a run-local cache, then re-renders.

- [x] **Outbound provenance links** (`report_html.py`) — each score ID links to its PGS
      Catalog score page (cohort), a "Learn more" line links the **MONDO disease class**
      (Monarch; EFO/HP → EBI OLS term page) and the source study (PubMed). Links only, so the
      report stays a self-contained offline file.

- [x] **Absolute-risk resources expanded + sex-aware risk** — per-SD effects + baselines added
      for rheumatoid arthritis, IBD, atrial fibrillation, colorectal cancer, Alzheimer
      (`resources/pgs_effect.tsv`, `baseline_incidence.tsv`, sourced with PMIDs), so real
      elevated traits get a figure (HG001: RA ~1 in 14, IBD ~1 in 31). T1D/Parkinson/
      hypothyroidism deliberately left number-less (no standard per-SD OR or no defensible
      baseline). Sample sex now flows through `provenance.py` (`sample_sex`) → the report labels
      sex-specific figures.

- [x] **`run.sh --pgs-meta FILE`** — copies the run's `select_pgs` metadata into the score
      dir so `grade_pgs` resolves trait names + evidence grades for custom (non-starter) score
      sets. Fixes the full 61-trait card rendering `trait == pgs_id` (it had fallen back to the
      36-entry starter meta). Replaces the manual `cp meta results/<run>/score/pgs_catalog_meta.json`.

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
