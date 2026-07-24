# poly-suite on a whole-exome (WES) BAM — limitations

**Scope.** This documents what poly-suite (pgsc_calc scoring + FRAPOSA ancestry)
does when handed a **whole-exome capture (WES) BAM** instead of a WGS BAM.
It is an assessment of existing behaviour — **no pipeline logic was changed.**

## Verdict (one line)

poly-suite on a WES BAM **runs and is fail-safe, but produces no valid
polygenic output.** Essentially every PGS collapses to grade **C/D** carrying
the hard caveat *"NOT a valid PRS"*, and the ancestry call is unreliable.
**Safe, but not useful on WES.**

## Why WES breaks polygenic scoring (premise)

Exome capture targets only exons — roughly **1–2 % of the genome**. Published
polygenic scores use genome-wide SNPs that are **overwhelmingly intergenic and
intronic**; only a small minority sit in exons. On a WES BAM the vast majority
of PGS scoring loci are **off-target** — no reads, no genotype. Ancestry
markers (the reference-panel SNPs FRAPOSA projects onto) are likewise
genome-wide and mostly off-target.

## How poly-suite genotypes — and why the loci fall off-target

poly-suite does **not** call variants and intersect; it **force-genotypes the
BAM at the PGS scoring loci themselves**:

- `bin/genotype_prep.sh:3-5` — force-genotype at PGS loci so covered
  non-variant sites become `0/0`, **uncovered sites become `./.` (honest
  missing)**. The targets come from the scorefiles, not from a capture region.
- `bin/scoring_targets.py:2-3,22-24` — targets are built from the PGS
  harmonized scorefiles (genome-wide autosomal loci; chrX/Y/M excluded).

On WES, a PGS locus that lies outside a capture target has zero depth →
`bcftools call` emits `./.` → the site is **missing**, not matched. Since PGS
loci are overwhelmingly non-exonic, **almost all become `./.`**.

## The fail-safe: three independent gates, none of which a WES score can pass as valid

### Gate 1 — pgsc_calc `--min_overlap 0.1` (score is dropped entirely)

Every scoring invocation passes `--min_overlap 0.1`:
`bin/run.sh:135`, `bin/score_selected.sh:20,45`. A scorefile whose matched
fraction is below **10 %** is marked `score_pass=false` by pgsc_calc and
**excluded from output** (see the min_overlap lesson, `docs/CHANGES.md:131`).
On WES most scores fall here and are never scored at all — they simply do not
appear.

### Gate 2 — poly-suite QC / coverage gate, `MATCH_RATE_MIN = 0.75` (two-grade downgrade + hard caveat)

`bin/grade_pgs.py:38` sets `MATCH_RATE_MIN = 0.75` ("pgsc_calc's own default
overlap floor"). The gate, `bin/grade_pgs.py:168-173`:

```python
if match_rate is not None and match_rate < MATCH_RATE_MIN:
    grade = downgrade(grade, 2)
    notes.append(f"low variant coverage ({match_rate:.0%} of {nvar:,} scoring "
                 "variants matched) — input is not full-genotype; raw score is "
                 "biased and NOT a valid PRS")
```

- `match_rate` = matched / total scoring variants, from pgsc_calc's own match
  summary (`bin/grade_pgs.py:117-142`).
- `downgrade(g, 2)` drops **two** grades and is **capped at "D"**
  (`bin/grade_pgs.py:68-69`: `GRADES = ["A","B","C","D"]`). So base-grade A → **C**,
  base-grade B → **D**, base-grade C → **D**. **No score can remain grade A or B.**
- The note text is unconditional: **"input is not full-genotype; raw score is
  biased and NOT a valid PRS"** — attached verbatim to `allowed_statement`
  in the contract OmniGen consumes.

A score that slipped past Gate 1 (10 % < match_rate < 75 %) still lands here.
On WES match_rate is a few percent, so scores that appear at all are
downgraded and flagged.

### Gate 3 — calibration / ancestry gate → hard grade "D"

`bin/grade_pgs.py:175-181`: if no percentile is available (`pct is None`),
grade is set to a hard **"D"** with note *"UNCALIBRATED — raw sum only …"*.
Any score run without a working `--run_ancestry` calibration cannot be graded
above D.

**Combined effect.** A WES score that produces any output at all is graded
**C at best (base-grade-A scores) and D otherwise**, and always carries the
*"NOT a valid PRS"* caveat. There is **no code path by which a WES-derived score
is reported as grade A/B or as a valid PRS.** The gates are fail-safe.

## Expected match rate on WES (order of magnitude)

Exons are ~1–2 % of the genome and PGS variants are largely non-exonic, so the
fraction of a genome-wide PGS's loci that fall within exome-capture targets is
**of order a few percent** (single digits; PGS with more coding SNPs sit at the
high end, sparse genome-wide scores at the low end). Concretely for the launch
set (`bin/grade_pgs.py:43-56`):

- PGS000018 CAD — 1,745,179 variants (genome-wide) → single-digit % matched.
- PGS000014 T2D — 6,917,436 variants → single-digit % matched.
- PGS000004 breast (313 variants) / PGS000662 prostate (269 variants) — small
  hand-curated sets, still overwhelmingly non-exonic.

Every one of these is **far below the 0.75 QC floor**, and typically below the
**0.10 min_overlap floor** — i.e. essentially all scores are either dropped
(Gate 1) or downgraded-and-flagged (Gate 2). **Essentially all WES scores →
grade C/D, none valid.**

## Ancestry (FRAPOSA): produces a result, but unreliable and largely unflagged

- Ancestry runs only when a panel is supplied: `bin/run.sh:61`
  (`--run_ancestry "$PANEL"`), projected via FRAPOSA (`bin/run.sh:172-177`).
- FRAPOSA projects the sample onto the reference-panel PCA space using the
  **overlapping** panel SNPs. On WES the overlap is the same tiny off-target
  fraction, so the projection is computed from very few markers and the
  most-similar-population call is **statistically unreliable**.
- **Nothing in poly-suite flags ancestry as unreliable on low-coverage / WES
  input.** FRAPOSA still emits a `MostSimilarPop` and percentile; the grader
  consumes it (`bin/grade_pgs.py:175-186`) and only applies a generic
  *portability* note when an EUR-trained score meets a non-EUR call. There is
  **no WES-specific ancestry-reliability gate** — a known gap. Treat any WES
  ancestry inference as indicative only.

## Does poly-suite accept a target / exome BED? — No

Confirmed: `bin/run.sh:20-33` enumerates every accepted flag —
`--sample --pgs --outdir --vcf --bam --panel --sex --ref --bootstrap-vcf
--reuse-prep --work-cache --scorefile-cache --pgs-meta --batch
--max-batch-bytes --dry-run`. There is **no `--bed`, `--target-bed`, `--exome`,
`--capture`, `--interval`, or `--region` parameter** anywhere in `bin/`
(grep for `--(bed|target-bed|exome|capture|interval|region)` returns nothing).
The pipeline cannot be told "this is a WES capture; restrict to on-target
loci." It always force-genotypes the full genome-wide PGS locus set, and
off-target loci simply return `./.`.

## Bottom line

On a WES BAM poly-suite **runs to completion and is fail-safe**: the
`min_overlap 0.1` drop (Gate 1), the `MATCH_RATE_MIN = 0.75` two-grade
downgrade with the *"NOT a valid PRS"* caveat (Gate 2), and the uncalibrated
hard-D gate (Gate 3) together guarantee that **no WES-derived polygenic score
is ever presented as a valid, grade-A/B PRS.** Essentially every score
collapses to grade C/D ("insufficient / not a valid PRS"), and ancestry is
unreliable and unflagged. **The pipeline is safe on WES, but produces no valid
polygenic output — WES is not a supported input for polygenic scoring.**

## Evidence index (file:line)

| Claim | Evidence |
|---|---|
| QC floor is 0.75 | `bin/grade_pgs.py:38` |
| <0.75 → downgrade 2 + "NOT a valid PRS" | `bin/grade_pgs.py:168-173` |
| downgrade capped at D; grades A–D | `bin/grade_pgs.py:68-69` |
| match_rate = matched/total from match summary | `bin/grade_pgs.py:117-142` |
| uncalibrated → hard grade D | `bin/grade_pgs.py:175-181` |
| force-genotype PGS loci; uncovered → `./.` | `bin/genotype_prep.sh:3-5` |
| targets from scorefiles, autosomes only | `bin/scoring_targets.py:2-3,22-24` |
| pgsc_calc `--min_overlap 0.1` (Gate 1) | `bin/run.sh:135`, `bin/score_selected.sh:20,45` |
| min_overlap excludes low-overlap scores | `docs/CHANGES.md:131` |
| ancestry via `--run_ancestry`/FRAPOSA | `bin/run.sh:61,172-177` |
| no BED/exome/target/interval flag | `bin/run.sh:20-33` (full flag list) |
