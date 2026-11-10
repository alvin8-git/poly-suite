#!/usr/bin/env python3
"""fastpath BUILD step "label" (runs inside the pygscatalog container).

Labels the cached canonical match candidates with the PGS Catalog's own `label_matches` (the work
`pgscatalog-matchmerge` does per sample in pgsc_calc), in disjoint accession bins so peak memory stays
bounded (every label operation is keyed by accession), and caches
  label/L_<bin>.ipc      labelled candidates, with `pos` = row index in the canonical candidate file
  label/summary_counts.ipc   counts of (scoring row x its best candidate state), the pre-merge match summary
The per-sample step (fastpath_delta.py) then only re-labels the groups a sample's variant list changes.

usage: fastpath_label.py CAND.ipc FILTER_IDS OUTDIR NBINS SCOREFILE...
"""
import os, sys, time, pathlib, tempfile
import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fastpath_match

KEYS = ["dataset", "accession", "ambiguous", "is_multiallelic", "match_flipped", "duplicate_best_match",
        "duplicate_ID", "match_IDs", "match_status"]


def main():
    cand, filt, outdir, nbins, *scorefiles = sys.argv[1:]
    nbins = int(nbins)
    out = pathlib.Path(outdir); out.mkdir(parents=True, exist_ok=True)
    fastpath_match.patch()
    from pgscatalog.match.lib import ScoringFileFrame
    from pgscatalog.match.lib._match.label import label_matches
    tmp = pathlib.Path(tempfile.mkdtemp(dir=outdir))
    with ScoringFileFrame(paths=[pathlib.Path(f) for f in scorefiles], chrom=None, cleanup=True, tmpdir=tmp) as score_df, \
            pl.StringCache():
        t0 = time.time()
        accs = sorted(score_df.select(pl.col("accession").cast(pl.Utf8)).unique().collect()["accession"].to_list())
        sizes = (score_df.group_by(pl.col("accession").cast(pl.Utf8)).len().collect())
        size = dict(zip(sizes["accession"].to_list(), sizes["len"].to_list()))
        bins, load = [[] for _ in range(nbins)], [0] * nbins
        for a in sorted(accs, key=lambda a: -size[a]):          # greedy balance by scoring rows
            i = load.index(min(load)); bins[i].append(a); load[i] += size[a]
        counts = []
        for i, b in enumerate(bins):
            lf = (pl.scan_ipc(cand).with_row_index("pos")
                  .filter(pl.col("accession").cast(pl.Utf8).is_in(b)))
            lab = label_matches(lf, keep_first_match=False, remove_ambiguous=True, remove_multiallelic=True,
                                skip_flip=False, filter_IDs=filt).collect()
            lab.write_ipc(out / f"L_{i:02d}.ipc")
            best = lab.lazy().filter(pl.col("best_match"))
            sf = score_df.filter(pl.col("accession").cast(pl.Utf8).is_in(b))
            c = (sf.join(best, on=["row_nr", "accession"], how="full")
                 .with_columns(pl.col("match_status").fill_null(value="unmatched"), dataset=pl.lit("x"))
                 .group_by(KEYS).len().rename({"len": "count"}).collect())
            counts.append(c)
            print(f"[label] bin {i+1}/{nbins}: {lab.height:,} candidates, {time.time()-t0:.0f}s", flush=True)
            del lab, best
        pl.concat(counts).write_ipc(out / "summary_counts.ipc")
        # total scoring rows per accession (for the pass/fail rate)
        sizes.write_ipc(out / "scoring_rows.ipc")


if __name__ == "__main__":
    main()
