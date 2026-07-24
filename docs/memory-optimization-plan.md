# Memory-optimization plan — averting the pgsc_calc OOM

Status: **VALIDATED** (2026-07-08) · **independently re-assessed → GO** (2026-07-24, see §8).
Behind `--batch`, old path unchanged. Owner: poly-suite.
Trigger: the 70-trait launch set OOM-hung the box under parallel runs (2026-07-06).

> **Integration defect found 2026-07-24:** `bin/run.sh` (`--batch` path, committed in
> 5c69e06) calls `bin/partition_scores.py` and `bin/merge_scores.py`, but those two
> scripts — plus `fit_budget.py`, `measure_mem.sh`, `test_batch_equiv.sh` and this doc —
> were left **untracked** on `main`. So `main`'s `--batch` mode is half-committed and
> breaks on a fresh clone. Fix: commit the helper scripts (done on branch `perf/mem-batch`).

Built: `bin/partition_scores.py` (size-aware FFD bin packing, self-test), `bin/merge_scores.py`
(disjoint row-concat + ancestry-drift assertion, self-test), `run.sh --batch`/`--max-batch-bytes`
(default 525 MB → 5 bins on the launch set, shared per-sample work dir + `-resume` for ancestry
reuse), and `bin/test_batch_equiv.sh` (single vs batched equivalence gate). Unit self-tests pass.

Equivalence gate (§6.1) PASSED on HG001 (6 scores, single vs `--batch` into 2 bins): the
**matched-variant set is identical per score** (the exact invariant — batching does not change
which variants match) and the **ancestry percentile is byte-identical** (→ identical grades).
The raw `SUM` differs only by floating-point summation-order noise (~5e-6 relative): PLINK2
reduces 1e5–1e6 dosage×weight terms and a different scorefile grouping reorders the sum. So
batching is **grade-identical and matched-set-exact**, not bit-exact on SUM — which is expected
and clinically irrelevant (rank/percentile/grade are unchanged).

The 525 MB default is calibrated from the §4 measurement pass (`bin/measure_mem.sh` +
`bin/fit_budget.py`): `peak_rss(GB) = 36.33 × scorefile_GB + 15.8` → ~34.9 GB worst bin, safe
2–3 wide. Tune with `--max-batch-bytes` (1200 MB → 3 bins / ~59 GB for serial-fast).

---

## 1. The problem, measured

From HG001's own `pipeline_info/execution_trace` (a run that *succeeded* solo):

```
 process (pgsc_calc)        peak_rss   peak_vmem   realtime   scope
 ─────────────────────────  ────────   ─────────   ────────   ───────────────
 MATCH_COMBINE   (match)     108 GB      269 GB     1h31m     per-SAMPLE
 FORMAT_SCOREFILES (format)   88 GB       92 GB       57m     per-SCORE-SET
 MATCH_VARIANTS  (match)      71 GB      126 GB        8m     per-SAMPLE
 PLINK2_SCORE ×4               7 GB        8 GB     ~25m ea   per-sample
 everything else              <4 GB                            —
```

Box = 125 GB. **A single run peaks at 108 GB (86% of RAM) for ~90 min.** HG001
solo survived on luck. Any second run overlapping that window → global OOM →
kernel thrash → server hang (confirmed via dmesg: 9 `global_oom` kills, every
victim `pgscatalog-match`/`-format`).

The memory note's "cap at 2" was wrong: at a 108 GB single-run peak, **2-wide is
already unsafe** — two runs sitting near 90–108 GB for 90 min collide almost
surely, which is exactly what killed the batch.

## 2. Root cause

`pgscatalog-format` (FORMAT_SCOREFILES) and `pgscatalog-match` (MATCH_VARIANTS /
MATCH_COMBINE) load their whole working set into one in-memory frame. Peak RSS
scales with **total variants across the score set**. The launch set is 136
scores, but ~8 of them are genome-wide (73–135 MB gz, millions of variants each):

```
 scorefile bytes (gz), cache/scorefiles/  — the memory driver
 PGS002867  135 MB ┓
 PGS003724  135 MB ┃
 PGS004693  127 MB ┃  ~8 genome-wide scores
 PGS012532  122 MB ┣━ dominate FORMAT + MATCH memory
 PGS012533  122 MB ┃
 PGS003212   99 MB ┃
 PGS004798   73 MB ┃
 PGS004799   73 MB ┛
 …128 others  <20 MB each (long tail, cheap)
```

Everything runs in one `nextflow run` over the full `--scorefile *.txt.gz` glob
(run.sh:133), so all 136 land in one frame → the 108 GB peak.

**Levers we do NOT control:** `pgscatalog-*` are upstream tools; rewriting them
for streaming is an innovation-token spend we should avoid. The lever we DO
control is the wrapper: how many scores we hand to each invocation.

## 3. Design — score-set batching + merge

PGS scores are independent: the final `aggregated_scores` has one column per
(score × sample). Splitting the score set across invocations and concatenating
the resulting columns is **exact, not approximate**. So:

```
                          bin/run.sh  (score step, batched)
 ┌────────────────────────────────────────────────────────────────────────┐
 │ partition 136 scorefiles  ──►  size-aware bins (byte budget B)           │
 │                                                                          │
 │   bin_small : 128 tiny scores           ~15–25 GB peak                   │
 │   bin_big_1 : 2 genome-wide scores      ~30–40 GB peak   } sequential    │
 │   bin_big_2 : 2 genome-wide scores      ~30–40 GB peak   } per sample,   │
 │   bin_big_3 : 2 genome-wide scores      ~30–40 GB peak   } shared work   │
 │   bin_big_4 : 2 genome-wide scores      ~30–40 GB peak   } dir + -resume │
 │        │                                                                 │
 │        ▼   pgsc_calc per bin  (ancestry computed ONCE in bin 1,          │
 │            reused by bins 2..K via -resume on the per-sample work dir)   │
 │        │                                                                 │
 │        ▼   K × aggregated_scores.txt.gz  (disjoint score columns)        │
 │            bin/merge_scores.py  ── outer join on (sampleset, IID) ──►     │
 │            score/aggregated_scores.txt.gz  (all 136 columns, as before)  │
 └────────────────────────────────────────────────────────────────────────┘
        │
        ▼  grade_pgs.py + report_html.py  ── UNCHANGED (read the merged file)
```

### 3.1 Partition (size-aware bin packing)
- Sort the resolved scorefiles by gz size desc.
- Greedy first-fit into bins with a total-bytes budget `B` (tunable; §4).
- The ~8 big scores each land near-alone; the 128 small ones fill one bin.
- Deterministic (size then id) so `-resume` cache keys are stable across reruns.

### 3.2 Ancestry reuse (the one non-obvious bit)
`--run_ancestry` (16 GB panel extract + FRAPOSA projection) is **score-independent**
— identical across a sample's bins. Running it per bin would K× the panel work.
Fix: give each *sample* one shared work dir and run its bins **sequentially with
`-resume`** (this is exactly what `--work-cache` already does for the panel).
Ancestry computes in bin 1, the rest cache-hit it. Bin order is fixed so bin 1
is always the same → stable cache.

Across *different* samples: separate work dirs, so samples could run parallel
once per-bin peak is low enough (§4 decides how wide). Start serial-across-
samples, widen after measuring.

### 3.3 Merge — `bin/merge_scores.py` (new, ~40 lines)
- Input: K per-bin `aggregated_scores.txt.gz`.
- Outer join on `(sampleset, IID)`; columns are disjoint per bin, so no
  collisions. Assert row-key sets are identical across bins (same sample →
  same IIDs); fail loud if not.
- Emit the combined file where `grade_pgs.py` expects it. Nothing downstream
  changes.

### 3.4 run.sh changes
- New `--max-batch-bytes B` (default from §4).
- Replace the single score invocation (run.sh:133) with: partition → per-bin
  loop (shared per-sample work dir, `-resume`) → `merge_scores.py`.
- `--reuse-prep` path unchanged; batching is orthogonal to prep reuse.
- Behind a flag (`--batch`); no-flag path stays the current single invocation,
  so we can A/B and roll back trivially.

## 4. Memory model — MEASURED (2026-07-08)

`bin/measure_mem.sh` ran FORMAT+MATCH on HG001 over 5 subsets (no ancestry — the
hogs are score-driven), `bin/fit_budget.py` fit the curve:

```
 subset   scorefile GB  files   overall peak   FORMAT  MATCH_VARIANTS  MATCH_COMBINE
 b1              0.14      1        18.5 GB      18.5        9.7            7.0
 b2              0.27      2        29.2 GB      29.2       13.7           11.4
 b4              0.52      4        35.3 GB      35.3       20.2           19.9
 b8              0.88      8        45.0 GB      45.0       30.9           34.6
 small           1.56    128        73.5 GB      53.0       49.1           73.5
```

**Fit: `peak_rss(GB) = 36.33 × scorefile_GB + 15.8`** (5 points; predicts small
72.5 vs actual 73.5). The 15.8 GB intercept is a fixed per-invocation floor
(nextflow + base working set), so many tiny bins each cost ~18 GB — but they run
serially per sample sharing one work dir, so it's wall-clock, not peak. FORMAT
dominates the peak up to ~1 GB/bin; MATCH_COMBINE overtakes it above that.

Chosen budgets (invert the fit):
- **2–3 wide parallel (default):** `B = 525 MB` → worst bin **~34.9 GB** (3 × 35 ≈
  105 GB < 125). Launch set (2.44 GB, 136 scores) → **5 bins**. This is now the
  `run.sh` default.
- **Serial-fast:** `B = 1200 MB` → worst bin ~59.4 GB, launch set → **3 bins**
  (fewer cold-starts; safe only at MAXJ=1). Pass `--max-batch-bytes 1200000000`.

Caveat: measured without `--run_ancestry`. In a real batched run ancestry runs
once in bin 1 (~16 GB panel extract) and cache-hits the rest; the solo trace
(§1) peaked at 108 GB on MATCH_COMBINE with ancestry present, so ancestry does
not stack above the score-driven peak. Re-measure if that assumption is ever in
doubt. Recalibrate the curve with `bin/measure_mem.sh` if pgsc_calc changes.

## 5. Tradeoffs

| | Batching (this plan) | Serial + no code (fallback) | Drop big scores |
|---|---|---|---|
| Per-run peak | ~35 GB (tunable) | 108 GB | ~25 GB |
| Parallelism | 2–3 wide | 1 (serial) | 3–4 wide |
| 7-sample wall-clock | ~½ of serial | ~17 h | fast, fewer traits |
| Code change | run.sh + merge (~2 files) | none | select_pgs list |
| Traits reported | all 70 | all 70 | 70 − biomarkers |
| Risk | merge correctness (tested) | none | clinical scope call |

Costs of batching, stated honestly:
- **K× nextflow cold-start** (~1–2 min each) + K× MATCH/FORMAT fixed overhead
  per sample. With ~5 bins that's ~10 min extra/sample — cheap vs the 90-min
  MATCH_COMBINE we're shrinking.
- **Merge is a new correctness surface.** Mitigated by §6.
- Sequential bins per sample cap intra-sample parallelism; the win is
  cross-sample. Acceptable.

Rejected alternatives:
- **Format-once cache reuse:** FORMAT is only #2. The #1 hog MATCH_COMBINE is
  per-sample, so caching FORMAT can't dodge it. Weak.
- **cgroup / systemd-run memory cap:** not available here (`systemd --user`
  bus absent, `/sys/fs/cgroup` is tmpfs). Docker `--memory` on the pgsc_calc
  containers is also unenforced for the same reason (that's *why* it went
  global-OOM instead of a contained container kill). So no clean cap exists in
  this env — serial is the only hard guarantee without batching.
- **ulimit -v:** MATCH_COMBINE reserves 269 GB vmem to use 108 GB RSS; a vmem
  cap kills it far too early. Unusable.

## 5b. Genome parallelization — considered, does NOT solve the peak

Natural question: split the single run by chromosome to parallelize? It doesn't
help the binding constraint, and it trades against memory. The pipeline is a
map→reduce over the genome:

```
 step                chrom-shardable?   memory driver           genome-parallel helps peak?
 FORMAT_SCOREFILES    no (per score-set) #scores × variants      no
 MATCH_VARIANTS 71GB  yes  (MAP)         one chrom's variants    yes — but only 8 min
 MATCH_COMBINE 108GB  NO   (REDUCE)      #scores × ALL variants  NO  ← the actual peak
 PLINK2_SCORE   7GB   yes  (MAP)         already tiny            already fine
 SCORE_AGGREGATE      reduce (tiny)      —                       n/a
```

The 108 GB hog `MATCH_COMBINE` is the **reduce**: it must see every matched
variant genome-wide at once to compute each score's overlap and build the
combined scoring file. Sharding splits the map steps; the reduce recombines
them, so its peak is unchanged. FORMAT (#2) isn't genomic at all.

And the memory/parallelism tension:

```
 chrom shards in PARALLEL → same aggregate peak, spread over concurrent tasks   faster, not lighter
 chrom shards SERIAL      → lower peak (one chrom working set at a time)         lighter, not faster
 fewer SCORES per task    → lower peak AND parallelizable                        lighter AND faster ← §3
```

Sharding only lowers peak if shards run serially (slower). The working set is
dominated by the *score* dimension (the ~8 genome-wide scores), not chromosome,
so score-batching (§3) is the lever. Per-chromosome target split is a possible
*complement inside a batch* to cap FORMAT/MATCH_VARIANTS, run serially — but only
if the measurement pass (§4) shows the map steps, not COMBINE, are binding after
score-batching. Don't add the orchestration blind.

## 6. Test coverage (non-negotiable)

The merge and partition are money-path logic → they get checks:

1. **Equivalence test (the important one) — PASSED 2026-07-08:** batched vs
   unbatched on HG001. The assertion evolved during validation: raw `SUM` is
   *not* bit-identical (float summation-order noise, ~5e-6 rel), so the gate
   asserts the **exact** invariants instead — identical matched-variant set per
   score (from the match summary) and byte-identical ancestry percentile (→ same
   grade) — and allows SUM to drift within float tolerance. This is the
   correctness proof for the whole plan. `bin/test_batch_equiv.sh`.
2. **Merge self-check** (`merge_scores.py __main__`): synthetic K=3 frames with
   disjoint columns + one deliberately mismatched IID set → assert join is
   exact and the mismatch raises. No framework, one `assert`-based `demo()`.
3. **Partition invariant:** every input scorefile lands in exactly one bin;
   no bin exceeds `B` unless it's a single oversized score. Assert in-code.

## 7. Rollout (incremental, reversible)

1. Measurement pass on HG001 → bytes→RSS curve → pick `B`.
2. Build `merge_scores.py` + its self-check.
3. Add `--batch` path to run.sh (old path stays default).
4. Equivalence test on HG001 (§6.1). Gate: identical scores.
5. Flip the launch batch to `--batch`, run 2-wide, watch `peak_rss`/`free`.
6. Update the OOM memory note + `docs/Documentation.md` §9 with the real
   108 GB figure and the batching mechanism.

Nothing here is a big-bang: the old path is the fallback at every step.

## 8. Independent re-assessment — 2026-07-24 (verdict: **GO**)

Re-evaluated by an independent in-repo pass (feasibility + correctness + payoff), not
re-running the full 136-scorefile pipeline. Evidence below is measured on this host
(125 GB RAM, 64 cores) unless noted.

### 8.1 Partitionability — CLEAN (this is what decides everything)
The batch split is safe because the heavy steps are per-PGS independent:
- **MATCH_COMBINE (the 108 GB reduce)** combines match *candidates* within its working
  set, but each scoring file's variants are matched against the target independently —
  a variant shared by two PGS is scored in both (long format, one row per PGS). There is
  **no cross-PGS dedup or global normalization** that a batch split would perturb.
- `--min_overlap 0.1` is a **per-scorefile** threshold, so it is unaffected by grouping.
- Ancestry (`pop_summary.csv`) is **score-independent** — computed from the reference-panel
  projection, identical across a sample's bins; the merge asserts byte-identity.
Verdict: partitioning by scorefile-batch does **not** change per-PGS matched sets, ancestry,
percentile, or grade. Merge is a disjoint-PGS row concatenation → exact, not approximate.

### 8.2 Correctness — re-verified (measured, against the 2026-07-08 `scratch/batch_equiv` outputs)
Independent diff of the retained single-vs-batched HG001 outputs (6 scores):
- **matched-variant counts identical** — 873,227 both paths (exact set invariant holds).
- **max |percentile diff| = 0.000e+00** → grades are byte-identical.
- **`pop_summary.csv` byte-identical** across bins (ancestry unperturbed).
- **max relative SUM drift = 4.7e-6** — float summation-order noise (PLINK2 reorders the
  1e5–1e6-term reduction under a different grouping); clinically irrelevant, rank/grade unchanged.
All three helper self-tests (`partition_scores.py`, `merge_scores.py`, `fit_budget.py`
`--self-test`) pass. Equivalence is **proven-by-run** for matched-set + grade identity.

### 8.3 Payoff — real, but it is *throughput*, not single-sample latency
Partition of the real 2.3 GB / 136-scorefile cache (measured):
- `--max-batch-bytes 525M` → **5 bins**, worst bin ~34.9 GB (fit predicts, ≈ measured).
- `--max-batch-bytes 1200M` → **3 bins**, worst bin ~59.4 GB.
- Monolithic (2.4 GB) → fit predicts 103 GB, matches the observed 108 GB peak.

Safe concurrency `K = floor(~110 GB usable / per-bin peak)`:
| mode | per-run peak | K (safe concurrent samples) |
|---|---|---|
| monolithic | 108 GB | **1** |
| batch 1200M | 59 GB | 1 |
| batch 525M | **35 GB** | **3** |

- **Serial batches (single sample):** same-or-slightly-slower wall-clock than monolithic
  (total FORMAT+MATCH work is unchanged; +K× nextflow cold-start + ~15.8 GB fixed floor per
  bin). Its only value alone is bounding peak so the run cannot OOM the box.
- **Parallel batches (the prize):** at 35 GB/bin the host safely runs **3 samples at once**
  vs 1 for monolithic. For a 7-sample cohort: monolithic ≈ 7 × 2.5 h serial ≈ 17 h;
  batched ≈ ceil(7/3) waves × ~3 h ≈ **~9 h (~2× cohort throughput)** — plus it removes the
  global-OOM/server-hang failure mode entirely (no cgroup/ulimit cap exists in this env, §5).

### 8.4 Overhead — does not erode the win
K× cold-start (~1–2 min) + 15.8 GB fixed floor per bin ≈ ~10 min extra/sample at 5 bins —
small against the ~90 min MATCH_COMBINE. `fit_budget.py` inverts the measured curve to pick
`B`; 525 MB (5 bins / 35 GB, 3-wide) is the sweet spot for this host. Going finer buys no
extra concurrency (already RAM-bound at K=3) and only adds cold-starts.

### 8.5 Recommendation
**GO.** The design is correct (partitionable, equivalence proven-by-run), the payoff is real
for multi-sample cohorts (~2× throughput + OOM safety), and effort/risk is low (helper scripts
+ a flagged run.sh path already written and self-tested; old path is the untouched default).
The **one required fix** is committing the untracked helper scripts so `main`'s already-committed
`--batch` path is not broken on a fresh clone. Prototype + this assessment live on branch
`perf/mem-batch`; not merged to `main`.
