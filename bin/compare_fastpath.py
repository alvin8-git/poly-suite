#!/usr/bin/env python3
"""Compare two poly-suite score directories (a pgsc_calc run vs a fastpath run, or any two runs).

Exits 1 on ANY difference beyond floating-point formatting (relative 1e-9), 0 when identical.

usage: compare_fastpath.py TRUTH_SCORE_DIR NEW_SCORE_DIR SAMPLE [--rtol 1e-9] [--id-map NEWNAME=SAMPLE]

--id-map declares that the new run's sample is called NEWNAME (file names and values) where the truth says SAMPLE;
only that exact ID is mapped, every other value must still be identical.

Compares, for every score and every column:
  * match/<S>_summary.csv        matched counts / coverage rows (exact, row for row)
  * match/<S>_ALL_*.scorefile.gz the scoring files plink2 reads (same lines, order-insensitive)
  * score/<S>_pgs.txt.gz         raw SUM/AVG + ancestry-adjusted Z / percentiles, sample AND reference rows
  * score/<S>_popsimilarity.txt.gz, <S>_info.json.gz, pop_summary.csv
  * pgs_scores.tsv / pgs_scores.json (the graded output OmniGen reads), when present in both dirs
Numbers equal after float() within --rtol count as "formatting only"; everything else is listed.
"""
import csv, glob, gzip, json, math, os, sys


def opener(p):
    return gzip.open(p, "rt", newline="") if p.endswith(".gz") else open(p, newline="")


def find(d, pat):
    h = sorted(glob.glob(os.path.join(d, "**", pat), recursive=True))
    return h[0] if h else None


class Rep:
    def __init__(self, rtol):
        self.rtol, self.diffs, self.stats = rtol, [], {}

    def num(self, what, a, b):
        s = self.stats.setdefault(what, {"exact": 0, "fmt_only": 0, "max_rel": 0.0, "n": 0})
        s["n"] += 1
        if a == b:
            s["exact"] += 1
            return
        try:
            fa, fb = float(a), float(b)
        except ValueError:
            self.diffs.append(f"{what}: {a!r} != {b!r}")
            return
        if math.isnan(fa) and math.isnan(fb):
            s["exact"] += 1
            return
        rel = abs(fa - fb) / max(abs(fa), abs(fb), 1e-300)
        s["max_rel"] = max(s["max_rel"], rel)
        if rel <= self.rtol:
            s["fmt_only"] += 1
        else:
            self.diffs.append(f"{what}: {a} vs {b} (rel {rel:.3e})")

    def eq(self, what, a, b):
        s = self.stats.setdefault(what, {"exact": 0, "fmt_only": 0, "max_rel": 0.0, "n": 0})
        s["n"] += 1
        if a == b:
            s["exact"] += 1
        else:
            self.diffs.append(f"{what}: {a!r} != {b!r}")


ID_MAP = {}          # declared sample-ID mapping NEWNAME -> SAMPLE (--id-map); only that exact string is mapped


def _m(x):
    return ID_MAP.get(x, x) if isinstance(x, str) else x


def rows(path, delim):
    with opener(path) as f:
        r = csv.reader(f, delimiter=delim)
        hdr = [_m(x) for x in next(r)]
        return hdr, [[_m(x) for x in row] for row in r]


def cmp_table(rep, tag, pa, pb, delim, key_cols):
    ha, ra = rows(pa, delim)
    hb, rb = rows(pb, delim)
    if ha != hb:
        rep.diffs.append(f"{tag}: header differs {ha} vs {hb}"); return
    ki = [ha.index(k) for k in key_cols if k in ha]
    da = {tuple(r[i] for i in ki) + (n,): r for n, r in enumerate(ra)} if not ki else None
    if ki:
        from collections import defaultdict
        ga, gb = defaultdict(list), defaultdict(list)
        for r in ra: ga[tuple(r[i] for i in ki)].append(r)
        for r in rb: gb[tuple(r[i] for i in ki)].append(r)
        if set(ga) != set(gb):
            rep.diffs.append(f"{tag}: key sets differ (only-truth {len(set(ga)-set(gb))}, only-new {len(set(gb)-set(ga))})")
        for k in sorted(set(ga) & set(gb)):
            if len(ga[k]) != len(gb[k]):
                rep.diffs.append(f"{tag} {k}: {len(ga[k])} vs {len(gb[k])} rows"); continue
            for x, y in zip(sorted(ga[k]), sorted(gb[k])):
                for c, (u, v) in zip(ha, zip(x, y)):
                    rep.num(f"{tag}:{c}", u, v)
    else:
        if len(ra) != len(rb):
            rep.diffs.append(f"{tag}: {len(ra)} vs {len(rb)} rows"); return
        for x, y in zip(ra, rb):
            for c, (u, v) in zip(ha, zip(x, y)):
                rep.num(f"{tag}:{c}", u, v)


def cmp_json(rep, tag, a, b, path=""):
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            rep.diffs.append(f"{tag}{path}: keys differ"); return
        for k in a:
            cmp_json(rep, tag, a[k], b[k], f"{path}/{k}")
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            rep.diffs.append(f"{tag}{path}: list length {len(a)} vs {len(b)}"); return
        for i, (x, y) in enumerate(zip(a, b)):
            cmp_json(rep, tag, x, y, f"{path}[{i}]")
    elif isinstance(a, str) and isinstance(b, str) and ID_MAP:
        rep.eq(f"{tag}:json-value", _m(a), _m(b))
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        rep.num(f"{tag}:json-number", repr(a), repr(b))
    else:
        rep.eq(f"{tag}:json-value", a, b)


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    rtol = 1e-9
    if "--rtol" in sys.argv:
        rtol = float(sys.argv[sys.argv.index("--rtol") + 1]); args = [x for x in args if x != str(rtol) and x != sys.argv[sys.argv.index("--rtol") + 1]]
    if "--id-map" in sys.argv:
        new, truth = sys.argv[sys.argv.index("--id-map") + 1].split("=")
        ID_MAP[new] = truth
        args = [x for x in args if x != sys.argv[sys.argv.index("--id-map") + 1]]
    if len(args) != 3:
        sys.exit(__doc__)
    A, B, S = args
    NEW = next(iter(ID_MAP), S)
    rep = Rep(rtol)
    pairs = [("match-summary", f"{S}_summary.csv", ",", ["accession", "match_status", "score_pass", "ambiguous",
                                                         "is_multiallelic", "duplicate_best_match", "duplicate_ID", "match_flipped"]),
             ("pgs", f"{S}_pgs.txt.gz", "\t", ["sampleset", "IID", "PGS"]),
             ("popsimilarity", f"{S}_popsimilarity.txt.gz", "\t", ["sampleset", "IID"]),
             ("pop_summary", "pop_summary.csv", ",", ["Population"]),
             ("graded-tsv", "pgs_scores.tsv", "\t", ["sample", "pgs_id"])]
    for tag, fname, delim, keys in pairs:
        pa, pb = find(A, fname), find(B, fname.replace(f"{S}_", f"{NEW}_", 1) if fname.startswith(f"{S}_") else fname)
        if pa is None and pb is None:
            continue
        if pa is None or pb is None:
            if tag == "graded-tsv":
                continue
            rep.diffs.append(f"{tag}: {fname} present in only one dir ({pa} / {pb})"); continue
        cmp_table(rep, tag, pa, pb, delim, keys)
    pa, pb = find(A, f"{S}_info.json.gz"), find(B, f"{NEW}_info.json.gz")
    if pa and pb:
        cmp_json(rep, "info.json", json.load(opener(pa)), json.load(opener(pb)))
    # scoring files: same set, same lines
    fa = {os.path.basename(p): p for p in glob.glob(os.path.join(A, "**", f"{S}_ALL_*.scorefile.gz"), recursive=True)}
    fb = {os.path.basename(p).replace(f"{NEW}_", f"{S}_", 1): p
          for p in glob.glob(os.path.join(B, "**", f"{NEW}_ALL_*.scorefile.gz"), recursive=True)}
    if set(fa) != set(fb):
        rep.diffs.append(f"scorefiles: names differ {sorted(fa)} vs {sorted(fb)}")
    for n in sorted(set(fa) & set(fb)):
        la, lb = sorted(opener(fa[n]).read().split("\n")), sorted(opener(fb[n]).read().split("\n"))
        rep.eq(f"scorefile-lines:{n}", len(la), len(lb))
        rep.eq(f"scorefile-content:{n}", la == lb, True)
    ga, gb = find(A, "pgs_scores.json"), find(B, "pgs_scores.json")
    if ga and gb:
        cmp_json(rep, "graded-json", json.load(open(ga)), json.load(open(gb)))
    print(f"{'comparison':44s} {'n':>10s} {'exact':>10s} {'fmt-only':>9s} {'max rel diff':>13s}")
    tot = {"n": 0, "exact": 0, "fmt_only": 0}
    for k in sorted(rep.stats):
        s = rep.stats[k]
        if s["n"] == 0:
            continue
        print(f"{k:44s} {s['n']:>10d} {s['exact']:>10d} {s['fmt_only']:>9d} {s['max_rel']:>13.3e}")
        for z in tot: tot[z] += s[z]
    print(f"{'TOTAL':44s} {tot['n']:>10d} {tot['exact']:>10d} {tot['fmt_only']:>9d}")
    if rep.diffs:
        print(f"\nDIFFERENCES ({len(rep.diffs)}):")
        for d in rep.diffs[:60]:
            print("  " + d)
        if len(rep.diffs) > 60:
            print(f"  ... {len(rep.diffs) - 60} more")
        sys.exit(1)
    print("\nIDENTICAL (differences, if any, are float formatting <= rtol)")


if __name__ == "__main__":
    main()
