# Fastpath: fast per-sample polygenic scoring

`bin/fastpath.py` scores a sample with **the same scores as the pgsc_calc path** while doing the
sample-independent work once. It is selected with `bin/run.sh --fastpath` (or `POLY_FASTPATH=1`);
**pgsc_calc stays the default** until the fastpath is accepted (status at the bottom).

## Why

On the server a 154-score HG002 run spent most of its 4 h in steps that do not depend on the sample:
`FORMAT_SCOREFILES` 30-60 min (95-150 GB), `INTERSECT_VARIANTS` 18-30 min (one core, reads the 84 M
variant panel), `MATCH_VARIANTS` 4-9 min (77-176 GB), `MATCH_COMBINE` 1.7-2.4 h (120 GB), plus the
`genotype_prep.sh` target-list build (8.5 min, rebuilt per sample from the same scoring files). The
memory peaks limit a 503 GB node to 2-3 samples at a time.

## What is cached, what is per sample

| Step (pgsc_calc process) | Fastpath |
|---|---|
| `scoring_targets.py` target site list | build, once (`cache/fastpath/<id>/target/targets.tsv`) |
| `FORMAT_SCOREFILES` | build, once: the PGS Catalog's own `pgscatalog-format`, same flags |
| `EXTRACT_DATABASE` (16 GB panel) | build, once, into `<id>/ref/` (+ a tabix index of the panel pvar) |
| `MATCH_VARIANTS` | build, once, against the **canonical** target site set (`pgscatalog-match --only_match`) -> `match/candidates.ipc.zst` |
| `INTERSECT_VARIANTS` | build, once, on the canonical target (`pgscatalog-intersect`) |
| `PLINK2_VCF` | per sample, identical command (17 s) |
| `MATCH_COMBINE` | per sample: **restrict** the cached candidates to the sample's own pvar, then the PGS Catalog's `pgscatalog-matchmerge` entry point (see "Per-sample missing sites") |
| `INTERSECT_VARIANTS` | per sample: derived from the cached table (seconds) |
| ancestry projection chain, `RELABEL_*`, `PLINK2_SCORE` (sample **and** reference panel), `SCORE_AGGREGATE`, `ANCESTRY_ANALYSIS`, `SCORE_REPORT` | per sample, the same commands/flags as pgsc_calc 2.3.0 (they depend on the sample's matched variant set) |

The reference-panel scores are **not** cacheable across samples: `plink2 --score` sums in single
precision, so the sum over a *different* scoring-file variant list differs by ~1e-7 relative. Because
the match depends on the sample's own variant list (next section), the scoring file differs between
samples, so the panel is re-scored per sample with exactly the scoring file pgsc_calc would have used.

## Per-sample missing sites (how pgsc_calc treats them, and how the fastpath reproduces it)

pgsc_calc matches the scoring files against the sample's **own** pvar (`PLINK2_VCF` output). A target
site contributes only if it is in that pvar. Sites vanish per sample because

* `plink2 --max-alleles 2` drops every VCF record with more than two alleles ("Too many alleles"); and
* `bcftools call -C alleles -A` sometimes emits a different indel representation than the one in the
  target list (HG002: 15.8 k sample-only IDs, 19.6 k target IDs absent), so the sample pvar is not a
  subset of the canonical site set.

What pgsc_calc does with a scoring variant whose site is absent: no match candidate exists, the row is
`unmatched`, counted as a miss in `match_rate` (`fail_rate = unmatched / all rows`), and a score is
dropped when `match_rate < min_overlap` (0.1). A site that is **present but with a missing genotype**
(`./.`) is *kept* in the scoring file: `INTERSECT_VARIANTS` does not filter on missingness (only the PCA
eligibility column does) and plink2 mean-imputes with the reference allele frequency.

The match candidates are one row per (scoring row, target pvar row), so a sample's candidates are
exactly the canonical candidates whose target ID is in the sample pvar, plus the candidates of the
sample-only rows. The fastpath therefore:

1. restricts the cached canonical candidates to the sample pvar;
2. runs the real `pgscatalog-match --only_match` on just the sample-only rows (a few thousand) and appends;
3. runs the real `pgscatalog-matchmerge` (`merge_cli.run_merge`) with `--filter_IDs` (the canonical
   reference-biallelic ID list: a superset, and every candidate is already in the sample pvar) and the
   same `--min_overlap 0.1`, so labelling, duplicate handling, overlap fractions, `_summary.csv`
   (matched counts / coverage) and the plink scoring files are the library's own;
4. derives the sample's `matched_variants.txt.gz` from the cached canonical one (AAF, F_MISS_DOSAGE and
   PCA_ELIGIBLE recomputed with the tool's formula) and runs the real `pgscatalog-intersect` on the
   sample-only rows against only the reference rows at their positions.

It refuses (exit 2/3) rather than guess when an allele key collides between a sample-only and a cached
variant (the tool's positional pairing could not be reproduced).

Two behaviour-preserving patches are applied inside the container (checked by
`bin/test_fastpath_unit.py`): the per-row python UDF in `label._encode_match_priority` is replaced by
the identical vectorised expression, and the multi-GB per-variant `<S>_log.csv.gz` (read by no
consumer) is only written with `--full-log`.

## Cache

`bin/fastpath.py build --pgs "PGS...,..." --scorefile-cache cache/scorefiles --panel panel.tar.zst --ref hg38.fa`
writes `cache/fastpath/fp-<id>/` (gitignored with `cache/`): `target/`, `formatted/`, `match/`, `ref/`,
`intersect/`, `manifest.json`. The scoring files themselves come from the existing scorefile cache
(`bin/cache_scorefiles.py`); nothing is duplicated.

The **key** (`manifest.json: key`) is the sha256 of: the sorted PGS list, the content sha256 of every
scoring file, the sha256 of the target site list, the reference panel (size + head/tail sha256), the
tool images and pgsc_calc revision, `min_overlap`, target build and the build stage list. The PGS
metadata (traits, grades) is not an input of any cached product and is not part of the key.
`score` recomputes the key on every run and **fails loud (exit 2)** if the cache is missing, incomplete
or does not match; there is no fallback to pgsc_calc and no partial score is ever written.

## Outputs

`<outdir>/score/<S>/match/` (`<S>_ALL_*.scorefile.gz`, `<S>_summary.csv`), `<outdir>/score/<S>/score/`
(`<S>_pgs.txt.gz`, `_popsimilarity.txt.gz`, `_info.json.gz`, `pop_summary.csv`, `report.html`) and
`<outdir>/score/pipeline_info/fastpath_timings.json` (wall time and peak RSS per step). `run.sh` then
runs the unchanged `grade_pgs.py` etc., so `pgs_scores.tsv/.json`, `provenance.json` and the report
are produced exactly as for pgsc_calc.

## Check

`bin/compare_fastpath.py TRUTH_SCORE_DIR NEW_SCORE_DIR SAMPLE` compares every column of the match
summary, scoring files, `_pgs.txt.gz` (SUM/AVG, Z, percentiles; sample and reference rows),
popsimilarity, info.json and the graded `pgs_scores.tsv/.json`; exit 1 on any difference beyond
1e-9 relative. `bin/test_fastpath_unit.py` holds the unit checks.

## Acceptance status

**Status: ACCEPTED as an opt-in path (`--fastpath`); pgsc_calc stays the default.** Exact on 3 samples end to end through
`bin/run.sh --fastpath` (HG002, HG005, HG006): 0 differences in 4,830,368 values each, graded files included; 23-32 min per
sample (see "v3" below). Known limits: the reference panel is identified by size plus head/tail sha256 (not a full hash);
the delta match step's RSS is about 174 GB including memory-mapped arrow files (container capped at 90 GB real);
`--ref` (reference FASTA) is required by `build` and `score` and must match the cache (it is part of the cache id).

(History of the earlier, pre-acceptance designs follows.)

Design now in the tree (v2): the one-off build also runs the canonical labelling (`bin/fastpath_label.py`:
the library's own `label_matches`, in 16 accession bins, 11.7 min, <= 100 GB) and caches the labelled
candidates plus the pre-merge summary counts. Per sample `bin/fastpath_delta.py` deletes the candidates of
the sample's absent variants, adds the candidates of sample-only rows (real `pgscatalog-match`), re-labels only
the closure of affected (accession, row) groups with `label_matches` and rebuilds summary counts additively;
the plink scoring files come from the library's `PlinkFrames` on the sample's matched rows in canonical
candidate order. The reference panel is re-scored per sample with the sample's own scoring file (needed for
bit-identical single-precision sums).

* HG002 from scratch (results/fastpath/HG002_v2, truth results/concord/HG002/score, 154 scores):
  `bin/compare_fastpath.py` -> 4,830,368 values compared, **0 differences** (max relative 0): match summary,
  both plink scoring files, `_pgs.txt.gz` (raw SUM, Z, percentiles, sample and reference rows), popsimilarity,
  info.json, and the graded `pgs_scores.tsv/.json` (run with the same `pgs_catalog_meta.json`; `grade_pgs.py`).
  (An earlier full run on the first design was also exact; its output was overwritten.)
* HG006 (results/fastpath/HG006_v1, the server's prepped VCF, truth = server batch/HG006/score): same comparison,
  4,830,368 values, **0 differences** (14,856 sample-only pvar rows, 30,365 affected scoring rows). Total 2146 s.
  HG005 not compared (its server score dir was incomplete when checked).
* HG002 per-sample timings, v2 (wall / peak RSS): plink2 import 17 s / 1.4 GB; derive intersect 3.4 min;
  delta match+merge 14.8 min (affected rows 30,268 of 154 scores; peak RSS 160 GB counts mmap'd arrow files, the
  container is capped at 90 GB; 10 min of it is the plink wide pivot of 130 M rows); ancestry chain 1.3 min /
  8.3 GB; relabel 7.5-8 min x2 (parallel) / 3.8 GB; sample scoring 35 s / 2.2 GB; panel scoring 5-6 min x2 /
  9.7 GB; aggregate 3.1 min / 2.5 GB; ancestry-analysis 2 min; report 9 s. Total 38 min (2278 s).
  pgsc_calc on the same sample: about 4 h (concord trace: FORMAT 63 min, MATCH_COMBINE 106 min, ...).
* Build (once, 154 scores): targets 8.5 min, canonical import 22 s, format 23 min (peak 140 GB RSS incl. mmap,
  capped 105 GB), match 5 min (175 GB RSS incl. mmap), panel extract 1.7 min, intersect 18 min (0.3 GB),
  ref index 1 min, label 11.7 min (108 GB RSS incl. mmap, capped 100 GB). About 70 min.
* Not at the "few minutes, < 20 GB" target: remaining per-sample cost is the wide pivot (10 min, large),
  relabel (8 min), panel scoring (6 min), aggregate (3 min). Ideas: write the sample's wide files by
  filtering a cached canonical wide file when no row is re-labelled (pure deletion, as for HG002); relabel with
  a vectorised join (verify equal to `pgscatalog-relabel`); a first aggregate on cached panel columns.

## v3 (current): cached wide files + awk relabel, validated end to end

* Build adds stage `wide`: the canonical wide plink scoring files (`label/wide/CANON_ALL_*.scorefile.gz`) plus
  split statistics (`idsplit.ipc`, `split_counts.ipc`, `passing.txt`); build total about 80 min once.
* Per sample, when re-labelling leaves no candidate and the sample adds none (a pure deletion of whole variant
  IDs; HG002 and HG006), the wide scoring files are the cached canonical ones with the deleted IDs' rows removed
  (streaming awk). The plink split structure (file order, strictly ordered sizes, accession columns, passing
  scores, removed-row counts per file) is checked against the cached statistics, and any doubt falls back to the
  library `PlinkFrames` pivot (exact, ~10 min). Samples that ADD candidates (HG005) take the library path for the
  re-labelled groups and the pivot.
* Relabel is a streaming awk map, `fastpath.relabel_awk`: byte-identical (decompressed) to `pgscatalog-relabel` on
  the real HG002 file (3.74 GB) and in `bin/test_fastpath_unit.py`; an unmapped ID is a hard error. 62 s, 1.4 GB.
* End to end with `bin/run.sh --fastpath` (new dirs), `compare_fastpath.py` against the pgsc_calc truth including
  the graded files, 4,830,368 values each, **0 differences**: HG002 (HG002_v3, 23.4 min wall), HG006 (HG006_v2,
  24.0 min), HG005 (HG005_v2, 32.3 min; the candidate-adding case).
* HG002_v3 steps (wall s / peak RSS GB): plink2 import 17/1.3; derive intersect 204; delta match+merge 358/174
  (RSS counts the mmap'd labelled arrow files; container capped at 90 GB; this is the remaining memory hog,
  HG005 832/175); ancestry chain: filter 174/7.3, thinned 9.5/0.8, pca 26/1.1, project 4; relabel 61 and 79/1.4;
  relabel afreq 89/2.5; sample scoring 20 and 25/1.8; panel scoring 299 and 336/9.6; aggregate 186/2.5;
  ancestry analysis 118/0.3; report 9. Total 1326 s. HG006: 1352 s, same profile.
* Still above the "under 20 GB" goal: only the delta match step (mmap'd arrow scans). Panel scoring (~6 min),
  aggregate (3 min) and derive intersect (3.4 min) remain per sample.

## Resource knobs and server numbers

`bin/run.sh --fastpath` passes `--threads`, `--score-threads` and `--mem-gb` to `fastpath.py score` from two env vars:

* `POLY_FASTPATH_THREADS` (default: `NPROC`, else 8) - threads for plink2, polars and the containers.
* `POLY_FASTPATH_MEM_GB` (default 100, unchanged on dev) - docker memory cap of the heavy containers (delta match/merge).

Server (CentOS 7, kernel 3.10, root docker): its cgroup accounting needs more headroom than dev. The one-off
`build` was killed in MATCH_VARIANTS at `--mem-gb 110` and succeeded with **`--mem-gb 220`**; a per-sample `score` was
killed in the delta match step (rc=137, 161 GB RSS) at the default 100 GB. Set `POLY_FASTPATH_MEM_GB` (e.g. 200+)
on the server, and pass `--mem-gb 220` to `fastpath.py build` there.
