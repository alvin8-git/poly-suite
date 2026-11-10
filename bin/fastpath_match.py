#!/usr/bin/env python3
"""fastpath per-sample MATCH_COMBINE: runs INSIDE the pygscatalog container (polars, pgscatalog.*).

pgsc_calc's MATCH_VARIANTS joins every scoring variant against the sample's own pvar, then
MATCH_COMBINE (pgscatalog-matchmerge) labels / filters / writes the PLINK scoring files.  Each
match candidate is one (scoring row, target pvar row) pair, so the sample's candidates are
exactly the cached CANONICAL candidates (matched once against the full target site set)
restricted to the rows whose target ID exists in the sample's pvar.  This script does that
restriction and then runs the PGS Catalog's own `pgscatalog-matchmerge` entry point
(run_merge) on the restricted file, with the same flags pgsc_calc passes.

Two behaviour-preserving performance patches (verified in bin/test_fastpath_unit.py):
  * label._encode_match_priority: the python-UDF (map_elements) per candidate row is replaced by
    the identical vectorised when/then mapping;
  * _write.write_log: the per-variant `<S>_log.csv.gz` (GBs, read by no consumer) is only written
    with FASTPATH_FULL_LOG=1; the small `<S>_summary.csv` that grade_pgs/the report read is
    always written.

usage: fastpath_match.py CAND.ipc.zst SAMPLE.pvar.zst FILTER_IDS.gz OUTDIR DATASET MIN_OVERLAP EXTRA.pvar|- SCOREFILE...
"""
import os, sys, time, pathlib
import polars as pl
from xopen import xopen


def _vector_encode_match_priority(df):
    """== label._encode_match_priority without the per-row python lambda."""
    mp = (pl.when(pl.col("best_match").is_null()).then(0)
          .when(pl.col("best_match")).then(1).otherwise(3)).cast(pl.UInt8)
    ex = (pl.when(pl.col("exclude").is_null()).then(0)
          .when(pl.col("exclude")).then(2).otherwise(0)).cast(pl.UInt8)
    mx = pl.max_horizontal(mp, ex)
    status = (pl.when(pl.col("best_match").is_null() & pl.col("exclude").is_null()).then(pl.lit(None, dtype=pl.Utf8))
              .when(mx == 0).then(pl.lit("unmatched")).when(mx == 1).then(pl.lit("matched"))
              .when(mx == 2).then(pl.lit("excluded")).otherwise(pl.lit("not_best")))
    return df.with_columns(match_status=status.cast(pl.Categorical))


def patch():
    import pgscatalog.match.lib._match.label as label
    label._encode_match_priority = _vector_encode_match_priority
    import pgscatalog.match.cli._write as w
    _orig = w.write_log

    def write_log(matchresults, score_df):
        if os.environ.get("FASTPATH_FULL_LOG") == "1":
            return _orig(matchresults, score_df)
        matchresults.summary_log.write_csv(w.Config.OUTDIR / f"{w.Config.DATASET}_summary.csv")
    w.write_log = write_log


def restrict(cand, pvar, out, extra_ipc=None):
    t0 = time.time()
    ids = set()
    with xopen(pvar, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            ids.add(line.split("\t", 3)[2])
    ids.discard(".")
    print(f"[fastpath_match] sample pvar IDs: {len(ids):,}", flush=True)
    with pl.StringCache():
        lf = pl.scan_ipc(cand)
        total = lf.select(pl.len()).collect().item()
        keep = lf.filter(pl.col("ID").cast(pl.Utf8).is_in(pl.Series(list(ids), dtype=pl.Utf8)))
        kept = keep.collect()
        if extra_ipc is not None:
            # sample-only rows: candidates computed by the real pgscatalog-match on just those rows
            kept = pl.concat([kept, pl.scan_ipc(extra_ipc).collect()], how="vertical")
        kept.write_ipc(out)
        n = pl.scan_ipc(out).select(pl.len()).collect().item()
    print(f"[fastpath_match] candidates {total:,} -> {n:,} restricted to sample pvar ({time.time()-t0:.0f}s)",
          flush=True)


def main():
    cand, pvar, filt, outdir, dataset, min_overlap, extra_pvar, *scorefiles = sys.argv[1:]
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    restricted = outdir / f"{dataset}_candidates.ipc.zst"
    extra_ipc = None
    if extra_pvar != "-":
        import subprocess, glob
        dd = outdir / "delta"
        dd.mkdir()
        t0 = time.time()
        subprocess.run(["pgscatalog-match", "--dataset", dataset, "--scorefile", *scorefiles, "--target", extra_pvar,
                        "--only_match", "--outdir", str(dd), "-v"], check=True)
        hits = glob.glob(str(dd / "*.ipc.zst"))
        assert len(hits) == 1, hits
        extra_ipc = hits[0]
        print(f"[fastpath_match] delta MATCH_VARIANTS on sample-only rows {time.time()-t0:.0f}s", flush=True)
    restrict(cand, pvar, restricted, extra_ipc)
    patch()
    from pgscatalog.match.cli import merge_cli
    sys.argv = ["pgscatalog-matchmerge", "--dataset", dataset, "--scorefile", *scorefiles,
                "--matches", str(restricted), "--min_overlap", min_overlap, "--filter_IDs", filt,
                "--outdir", str(outdir), "--combined", "-v"]
    t0 = time.time()
    merge_cli.run_merge()
    print(f"[fastpath_match] matchmerge {time.time()-t0:.0f}s", flush=True)
    restricted.unlink()


if __name__ == "__main__":
    main()
