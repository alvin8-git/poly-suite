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

**Example output:** [`examples/HG001-61trait/`](examples/HG001-61trait/) — the full calibrated 61-trait launch-set card for GIAB HG001 (118 scores, EUR-calibrated, 73 A / 15 B / 16 C / 14 D) · [`examples/HG001-18trait/`](examples/HG001-18trait/) — the earlier 18-trait card. Each: report.html + pgs_scores.tsv/.json + provenance.

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

Roadmap and the shipped-feature changelog live in **[TODO.md](TODO.md)**.
