#!/usr/bin/env python3
"""Merge per-bin pgsc_calc outputs back into one score dir (the --batch merge step).

Score-set batching runs pgsc_calc K times over disjoint subsets of the scorefiles,
then this reconstructs the single-run layout grade_pgs.py expects:

  <out>/<sample>/score/<sample>_pgs.txt.gz     merged score rows (all PGS)
  <out>/<sample>/score/pop_summary.csv         ancestry (identical across bins)
  <out>/<sample>/match/<sample>_summary.csv    merged per-score coverage rows

The score file and the match summary are LONG format (one row per PGS, or per
PGS x individual), and bins hold DISJOINT PGS, so merging is a row concatenation
with a shared header. It is exact, not approximate. We assert:
  - every bin shares the same header,
  - the PGS sets across bins are disjoint (no score double-counted),
  - pop_summary.csv is byte-identical across bins (ancestry is score-independent;
    if batching perturbed it, this fires) — the core correctness check.

Usage: merge_scores.py <out_score_dir> <sample> <bin_dir_1> <bin_dir_2> ...
       merge_scores.py --self-test
"""
import sys, os, glob, gzip, csv


def _find(bin_dir, pats):
    for pat in pats:
        hits = sorted(glob.glob(os.path.join(bin_dir, pat), recursive=True))
        if hits:
            return hits[0]
    return None


def find_score_file(bin_dir):
    # same precedence grade_pgs.py uses: ancestry-calibrated *_pgs.txt.gz, else aggregated
    return _find(bin_dir, ["**/*_pgs.txt.gz", "**/aggregated_scores*.txt.gz"])


def find_match_summary(bin_dir):
    return _find(bin_dir, ["**/match/*_summary.csv"])


def find_pop_summary(bin_dir):
    return _find(bin_dir, ["**/pop_summary.csv"])


def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def _pgs_col(header, sep):
    cols = header.rstrip("\n").split(sep)
    for name in ("PGS", "accession"):
        if name in cols:
            return cols.index(name)
    raise ValueError(f"no PGS/accession column in header: {cols[:12]}")


def concat_rows(paths, sep, out_path):
    """Concatenate data rows across bin files, one shared header, assert PGS disjoint.
    Returns the set of PGS accessions written (for cross-file sanity)."""
    header = None
    pgs_idx = None
    seen = set()
    opener = gzip.open if out_path.endswith(".gz") else open
    with opener(out_path, "wt") as out:
        for i, p in enumerate(paths):
            with _open(p) as fh:
                h = fh.readline()
                if header is None:
                    header = h
                    pgs_idx = _pgs_col(header, sep)
                    out.write(header)
                elif h != header:
                    raise ValueError(f"header mismatch in {p}:\n  {h!r}\n vs {header!r}")
                bin_pgs = set()
                for line in fh:
                    if not line.strip():
                        continue
                    pgs = line.rstrip("\n").split(sep)[pgs_idx]
                    bin_pgs.add(pgs)
                    out.write(line)
                dup = seen & bin_pgs
                if dup:
                    raise ValueError(f"PGS in >1 bin (not disjoint): {sorted(dup)[:5]}")
                seen |= bin_pgs
    return seen


def merge(out_dir, sample, bin_dirs):
    assert len(bin_dirs) >= 1, "need >=1 bin dir"
    score_dir = os.path.join(out_dir, sample, "score")
    match_dir = os.path.join(out_dir, sample, "match")
    os.makedirs(score_dir, exist_ok=True)
    os.makedirs(match_dir, exist_ok=True)

    score_files = [find_score_file(b) for b in bin_dirs]
    match_files = [find_match_summary(b) for b in bin_dirs]
    pop_files = [find_pop_summary(b) for b in bin_dirs]
    for b, s in zip(bin_dirs, score_files):
        if s is None:
            raise FileNotFoundError(f"no score file (*_pgs.txt.gz / aggregated_scores) in {b}")

    # ancestry must be score-independent -> identical across bins. Load once, verify all.
    pop_ref = None
    for p in pop_files:
        if p is None:
            continue
        data = open(p, "rb").read()
        if pop_ref is None:
            pop_ref = data
        elif data != pop_ref:
            raise ValueError("pop_summary.csv differs across bins — batching perturbed ancestry")

    score_out = os.path.join(score_dir, f"{sample}_pgs.txt.gz")
    pgs_score = concat_rows(score_files, "\t", score_out)

    if all(m is not None for m in match_files):
        match_out = os.path.join(match_dir, f"{sample}_summary.csv")
        pgs_match = concat_rows(match_files, ",", match_out)
        # match summary lists accessions; every scored PGS should appear there
        missing = pgs_score - pgs_match
        if missing:
            print(f"warn: {len(missing)} PGS in scores but not match summary (e.g. {sorted(missing)[:3]})",
                  file=sys.stderr)

    if pop_ref is not None:
        with open(os.path.join(score_dir, "pop_summary.csv"), "wb") as fh:
            fh.write(pop_ref)

    print(f"merged {len(bin_dirs)} bins -> {score_out}  ({len(pgs_score)} PGS)")
    return pgs_score


def _self_test():
    import tempfile, shutil
    d = tempfile.mkdtemp()
    try:
        # 3 bins, disjoint PGS, long format (2 individuals x N PGS), same pop_summary
        hdr = "sampleset\tFID\tIID\tPGS\tSUM\tpercentile_MostSimilarPop\n"
        pop = "sampleset,MostSimilarPop,p\nS,EUR,0.9\n"
        bins = []
        pgs_by_bin = [["PGS0001", "PGS0002"], ["PGS0003"], ["PGS0004", "PGS0005"]]
        for i, pgs_list in enumerate(pgs_by_bin):
            b = os.path.join(d, f"bin{i}", "S", "score")
            m = os.path.join(d, f"bin{i}", "S", "match")
            os.makedirs(b); os.makedirs(m)
            with gzip.open(os.path.join(b, "S_pgs.txt.gz"), "wt") as fh:
                fh.write(hdr)
                for pgs in pgs_list:
                    for iid in ("SAMPLE", "REF1"):  # sample + one reference individual
                        fh.write(f"S\t{iid}\t{iid}\t{pgs}\t1.0\t50.0\n")
            with open(os.path.join(m, "S_summary.csv"), "w") as fh:
                fh.write("dataset,accession,count\n")
                for pgs in pgs_list:
                    fh.write(f"S,{pgs},10\n")
            open(os.path.join(b, "pop_summary.csv"), "w").write(pop)
            bins.append(os.path.join(d, f"bin{i}"))

        out = os.path.join(d, "score")
        got = merge(out, "S", bins)
        assert got == {"PGS0001", "PGS0002", "PGS0003", "PGS0004", "PGS0005"}, got
        # merged score file has 5 PGS x 2 individuals = 10 data rows + 1 header
        with gzip.open(os.path.join(out, "S", "score", "S_pgs.txt.gz"), "rt") as fh:
            rows = [r for r in fh]
        assert len(rows) == 11, f"expected 11 lines, got {len(rows)}"
        assert open(os.path.join(out, "S", "score", "pop_summary.csv")).read() == pop

        # disjointness violation must raise
        b_dup = os.path.join(d, "bindup", "S", "score")
        os.makedirs(b_dup)
        with gzip.open(os.path.join(b_dup, "S_pgs.txt.gz"), "wt") as fh:
            fh.write(hdr + "S\tSAMPLE\tSAMPLE\tPGS0001\t1.0\t50.0\n")  # PGS0001 already in bin0
        open(os.path.join(d, "bindup", "S", "score", "pop_summary.csv"), "w").write(pop)
        try:
            merge(os.path.join(d, "score2"), "S", [bins[0], os.path.join(d, "bindup")])
            raise AssertionError("expected disjointness violation to raise")
        except ValueError as e:
            assert "not disjoint" in str(e), e

        # ancestry drift must raise
        b_drift = os.path.join(d, "bindrift", "S", "score")
        os.makedirs(b_drift)
        with gzip.open(os.path.join(b_drift, "S_pgs.txt.gz"), "wt") as fh:
            fh.write(hdr + "S\tSAMPLE\tSAMPLE\tPGS9999\t1.0\t50.0\n")
        open(os.path.join(d, "bindrift", "S", "score", "pop_summary.csv"), "w").write(
            "sampleset,MostSimilarPop,p\nS,EAS,0.9\n")  # different ancestry
        try:
            merge(os.path.join(d, "score3"), "S", [bins[0], os.path.join(d, "bindrift")])
            raise AssertionError("expected ancestry-drift to raise")
        except ValueError as e:
            assert "perturbed ancestry" in str(e), e

        print("merge_scores self-test PASSED")
    finally:
        shutil.rmtree(d)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
        _self_test()
    elif len(sys.argv) >= 4:
        merge(sys.argv[1], sys.argv[2], sys.argv[3:])
    else:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
