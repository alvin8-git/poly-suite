# poly-suite — roadmap & changelog

Extracted from the README `## Status` section. Open work first, shipped features below.

## Open / roadmap

- [ ] Genotype-prep v3: imputation branch for array/low-pass input (GLIMPSE/Beagle)

- [ ] PRS-CSx multi-ancestry derivation (deferred, research)

- [ ] Nextflow scaffolding — productionization once the flow is validated (bash orchestrator covers it now)


## Shipped

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
