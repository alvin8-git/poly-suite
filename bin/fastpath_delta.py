#!/usr/bin/env python3
"""fastpath per-sample match/merge (runs inside the pygscatalog container).

Equivalent of pgsc_calc's MATCH_VARIANTS + MATCH_COMBINE for one sample, computed as a DELTA against the
cached labelled canonical candidates (bin/fastpath_label.py):

 * ABSENT = canonical target IDs missing from the sample pvar; their candidates disappear.
 * EXTRA  = candidates of sample-only pvar rows (real `pgscatalog-match` on just those rows).
 * Every label (best match, duplicate best match, duplicate ID, ambiguity, filter) is keyed by
   (accession, row_nr) except duplicate_ID, which couples rows that share a candidate ID inside one
   accession.  AFFECTED = the closure of scoring rows that touch an ABSENT/EXTRA ID under that coupling.
 * AFFECTED rows are re-labelled with the library's own `label_matches`; all other rows keep their cached
   labels, so the sample's labelled table equals what matching against its own pvar would give.
 * summary counts are additive over rows -> cached counts - old(AFFECTED) + new(AFFECTED); the pass/fail
   rate is the library formula; the plink scoring files come from `PlinkFrames` on the sample's matched rows
   in canonical candidate order, exactly as pgscatalog-matchmerge writes them.

usage: fastpath_delta.py LABELDIR ABSENT_IDS.txt EXTRA.ipc|- FILTER_IDS OUTDIR DATASET MIN_OVERLAP [canon]

  canon   build mode: no absent IDs; also emits the canonical wide scoring files + split statistics into
          OUTDIR (cached as LABELDIR/wide) used by the per-sample fast path below.

Fast path (per sample): when re-labelling leaves no candidate and the sample adds none (a pure deletion of
whole variant IDs), the sample's wide scoring files are the cached canonical ones with the rows of the
deleted IDs removed; the plink split structure (file order, sizes, accession columns) is checked against the
cached split statistics.  Any doubt -> the library's PlinkFrames path (exact, slower).
"""
import os, sys, time, glob, pathlib
import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fastpath_match

RAW = ["row_nr", "chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight", "effect_type",
       "accession", "effect_allele_FLIP", "other_allele_FLIP", "ID", "REF", "ALT", "is_multiallelic",
       "matched_effect_allele", "match_type"]
KEYS = ["dataset", "accession", "ambiguous", "is_multiallelic", "match_flipped", "duplicate_best_match",
        "duplicate_ID", "match_IDs", "match_status"]
PLINK_COLS = ["row_nr", "accession", "effect_type", "chr_name", "ID", "matched_effect_allele", "effect_weight"]
SUMMARY_COLS = ["dataset", "accession", "score_pass", "match_status", "ambiguous", "is_multiallelic",
                "duplicate_best_match", "duplicate_ID", "match_flipped", "match_IDs", "count", "percent"]


def lap(t0, msg):
    print(f"[fastpath_delta] {msg} ({time.time()-t0:.0f}s)", flush=True)


def main():
    labeldir, absent_f, extra_ipc, filt, outdir, dataset, min_overlap = sys.argv[1:8]
    canon_mode = len(sys.argv) > 8 and sys.argv[8] == "canon"
    min_overlap = float(min_overlap)
    out = pathlib.Path(outdir); out.mkdir(parents=True, exist_ok=True)
    fastpath_match.patch()
    from pgscatalog.match.lib._match.label import label_matches
    from pgscatalog.match.lib._plinkframe import PlinkFrames
    t0 = time.time()
    with pl.StringCache():
        absent = pl.Series([l.strip() for l in open(absent_f) if l.strip()], dtype=pl.Utf8)
        files = sorted(glob.glob(os.path.join(labeldir, "L_*.ipc")))
        L = pl.concat([pl.scan_ipc(f) for f in files])
        EMPTY = pl.read_ipc(files[0], n_rows=0)
        extra = None
        if extra_ipc != "-":
            extra = pl.scan_ipc(extra_ipc).collect()
        # ---- affected closure -----------------------------------------------------------------
        K = ["accession", "row_nr"]
        A = L.filter(pl.col("ID").cast(pl.Utf8).is_in(absent)).select(K).unique().collect()
        if extra is not None and extra.height:
            A = pl.concat([A, extra.select(K)]).unique()
        n_abs_rows = A.height
        while True:
            ids = L.join(A.lazy(), on=K, how="semi").select(["accession", "ID"]).unique().collect()
            if extra is not None and extra.height:
                ids = pl.concat([ids, extra.select(["accession", "ID"])]).unique()
            A2 = pl.concat([A, L.join(ids.lazy(), on=["accession", "ID"], how="semi").select(K).unique().collect()]).unique()
            if A2.height == A.height:
                break
            A = A2
        lap(t0, f"affected scoring rows {A.height:,} (touching absent/extra IDs: {n_abs_rows:,})")
        # ---- re-label affected rows ----------------------------------------------------------------
        if A.height == 0:
            sub = old_all = EMPTY.clone()
        else:
            sub = L.join(A.lazy(), on=K, how="semi").filter(~pl.col("ID").cast(pl.Utf8).is_in(absent)).collect()
            old_all = L.join(A.lazy(), on=K, how="semi").collect()
        sub_raw = sub.select(RAW + ["pos"])
        if extra is not None and extra.height:
            e = (extra.select(RAW).with_columns([pl.col(c).cast(sub_raw.schema[c]) for c in RAW])
                 .with_columns(pos=pl.int_range(pl.len(), dtype=pl.UInt32) + 2**31))
            sub_raw = pl.concat([sub_raw.with_columns(pl.col("pos").cast(pl.UInt32)), e], how="vertical")
        sub_raw = sub_raw.sort("pos")
        if sub_raw.height == 0:
            newL = EMPTY.clone()
        else:
            newL = label_matches(sub_raw.lazy(), keep_first_match=False, remove_ambiguous=True, remove_multiallelic=True,
                                 skip_flip=False, filter_IDs=filt).collect()
        lap(t0, f"re-labelled {newL.height:,} candidates")

        # ---- summary: cached counts - old(A) + new(A) ----------------------------------------------
        def contrib(Akeys, labelled):
            best = labelled.lazy().filter(pl.col("best_match"))
            return (Akeys.lazy().join(best, on=K, how="left")
                    .with_columns(pl.col("match_status").fill_null(value="unmatched"), dataset=pl.lit("x"))
                    .group_by(KEYS).len().rename({"len": "count"}).collect())
        S = pl.read_ipc(os.path.join(labeldir, "summary_counts.ipc"))
        i64 = lambda d: d.select(S.columns).with_columns(pl.col("count").cast(pl.Int64))
        if A.height == 0:
            old_c = new_c = i64(S.head(0))
        else:
            old_c = i64(contrib(A, old_all)).with_columns(pl.col("count") * -1)
            new_c = i64(contrib(A, newL))
        tot = (pl.concat([i64(S), old_c, new_c])
               .group_by(KEYS).agg(pl.col("count").sum()).filter(pl.col("count") != 0))
        N = pl.read_ipc(os.path.join(labeldir, "scoring_rows.ipc")).rename({"len": "n_rows"})
        matched = (tot.filter(pl.col("match_status").cast(pl.Utf8) == "matched")
                   .group_by(pl.col("accession").cast(pl.Utf8)).agg(pl.col("count").sum().alias("n_matched")))
        rate = (N.join(matched, on="accession", how="left").with_columns(pl.col("n_matched").fill_null(0))
                .with_columns(fail_rate=(pl.col("n_rows") - pl.col("n_matched")) / pl.col("n_rows")))
        rate = rate.with_columns(score_pass=pl.col("fail_rate") <= (1 - min_overlap))
        passing = set(rate.filter(pl.col("score_pass"))["accession"].to_list())
        if not passing:
            raise SystemExit("All scores fail to meet the match threshold")
        sp = rate.select(["accession", "score_pass"]).with_columns(pl.col("accession").cast(pl.Categorical))
        summ = (tot.join(sp, on="accession", how="left")
                .with_columns((pl.col("count") / pl.sum("count") * 100).over(["dataset", "accession"]).alias("percent"))
                .with_columns(dataset=pl.lit(dataset)).select(SUMMARY_COLS))
        summ.write_csv(out / f"{dataset}_summary.csv")
        lap(t0, "summary written")

        # ---- scoring files ---------------------------------------------------------------------------
        wide = os.path.join(labeldir, "wide")
        no_new = newL.height == 0 and (extra is None or extra.height == 0)
        if not canon_mode and no_new and os.path.isdir(wide):
            if fast_wide(L, A, passing, wide, out, dataset, pl, K, lap, t0):
                return
            print("[fastpath_delta] fast wide path declined -> library PlinkFrames path", flush=True)
        matched_cols = PLINK_COLS + ["pos"]
        base = (L.filter(pl.col("match_status").cast(pl.Utf8) == "matched").select(matched_cols)
                .join(A.lazy(), on=K, how="anti").collect())
        newm = newL.filter(pl.col("match_status").cast(pl.Utf8) == "matched").select(matched_cols)
        F = (pl.concat([base, newm.with_columns(pl.col("pos").cast(base.schema["pos"]))]).sort("pos")
             .filter(pl.col("accession").cast(pl.Utf8).is_in(list(passing))))
        lap(t0, f"matched rows {F.height:,}")
        outfs, idsplit, counts = [], [], []
        for frame in PlinkFrames.from_matchresult(F.lazy()):
            outfs.append(frame.write(directory=out, split=False, dataset=dataset))
            if canon_mode:
                d = frame.df.collect()
                et = str(frame.effect_type)
                idsplit.append(d.select(["ID", "matched_effect_allele"]).unique()
                               .with_columns(effect_type=pl.lit(et), n=pl.lit(int(frame.n))))
                counts.append(d.group_by(pl.col("accession").cast(pl.Utf8)).len()
                              .with_columns(effect_type=pl.lit(et), n=pl.lit(int(frame.n))))
        lap(t0, f"scoring files written: {[str(x[0].name) for x in outfs]}")
        if canon_mode:
            pl.concat(idsplit).with_columns(pl.col("matched_effect_allele").cast(pl.Utf8)).write_ipc(out / "idsplit.ipc")
            pl.concat(counts).rename({"len": "rows"}).write_ipc(out / "split_counts.ipc")
            (out / "passing.txt").write_text("\n".join(sorted(passing)) + "\n")


def fast_wide(L, A, passing, wide, out, dataset, pl, K, lap, t0):
    """Per-sample wide scoring files by deleting whole variant IDs from the cached canonical files."""
    import subprocess, gzip
    if set(open(os.path.join(wide, "passing.txt")).read().split()) != set(passing):
        print("[fastpath_delta] set of passing scores differs from canonical", flush=True)
        return False
    idsplit = pl.read_ipc(os.path.join(wide, "idsplit.ipc"))
    sc = pl.read_ipc(os.path.join(wide, "split_counts.ipc"))
    dele = (L.filter(pl.col("match_status").cast(pl.Utf8) == "matched").join(A.lazy(), on=K, how="semi")
            .select(["ID", "matched_effect_allele", "effect_type", "accession"]).collect()
            .with_columns(pl.col("matched_effect_allele").cast(pl.Utf8), pl.col("effect_type").cast(pl.Utf8),
                          pl.col("accession").cast(pl.Utf8), pl.col("ID").cast(pl.Utf8))
            .filter(pl.col("accession").is_in(list(passing))))
    j = dele.join(idsplit, on=["ID", "matched_effect_allele", "effect_type"], how="left")
    if j["n"].null_count():
        print("[fastpath_delta] deleted row without canonical split entry", flush=True); return False
    lens = sc.group_by(["effect_type", "n"]).agg(pl.col("rows").sum().alias("len"))
    dl = j.group_by(["effect_type", "n"]).len().rename({"len": "del"})
    dac = j.group_by(["effect_type", "n", "accession"]).len().rename({"len": "delac"})
    nl = lens.join(dl, on=["effect_type", "n"], how="left").with_columns(pl.col("del").fill_null(0)) \
        .with_columns(new=pl.col("len") - pl.col("del"))
    for et in nl["effect_type"].unique().to_list():
        v = nl.filter(pl.col("effect_type") == et).sort("n")
        o, nw = v["len"].to_list(), v["new"].to_list()
        if min(nw) <= 0 or any(o[i] <= o[i + 1] for i in range(len(o) - 1)) or any(nw[i] <= nw[i + 1] for i in range(len(nw) - 1)):
            print(f"[fastpath_delta] split sizes {o} -> {nw} not strictly ordered", flush=True); return False
    cc = sc.join(dac, on=["effect_type", "n", "accession"], how="left").with_columns(pl.col("delac").fill_null(0))
    if (cc["rows"] - cc["delac"]).min() <= 0:
        print("[fastpath_delta] an accession column would vanish from a split", flush=True); return False
    uniq = j.select(["ID", "matched_effect_allele", "effect_type", "n"]).unique()
    exp = {(r[0], r[1]): r[2] for r in uniq.group_by(["effect_type", "n"]).len().iter_rows()}
    ids = out / "deleted_ids.txt"
    ids.write_text("".join(i + "\n" for i in sorted(set(j["ID"].to_list()))))
    for et, n in sorted({(r[0], r[1]) for r in lens.iter_rows()}):
        src = os.path.join(wide, f"CANON_ALL_{et}_{n}.scorefile.gz")
        dst = out / f"{dataset}_ALL_{et}_{n}.scorefile.gz"
        cnt = out / f"removed_{et}_{n}.count"
        subprocess.run(["bash", "-c", "set -euo pipefail; zcat \"$1\" | awk -F'\\t' -v c=\"$4\" 'NR==FNR{d[$1];next} FNR==1||!($1 in d){print;next} {r++} END{print r+0 > c}' \"$2\" - | gzip -c > \"$3\"",
                        "_", src, str(ids), str(dst), str(cnt)], check=True)
        removed = int(cnt.read_text())
        want = exp.get((et, n), 0)
        if removed != want:
            print(f"[fastpath_delta] removed {removed} wide rows from {et}_{n}, expected {want}", flush=True)
            dst.unlink(); return False
    lap(t0, "scoring files by deleting rows of the cached canonical wide files")
    return True


if __name__ == "__main__":
    main()
