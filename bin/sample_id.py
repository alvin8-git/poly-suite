#!/usr/bin/env python3
"""Internal sampleset IDs for pgsc_calc and the fastpath.

pgsc_calc's samplesheet schema allows only [A-Za-z0-9] in a sampleset name (and its ancestry step splits file
names on '_'), but production sample names contain '_' / '.' / '-'.  Internally the pipeline therefore uses
the sample name with every non-alphanumeric character removed; after scoring, `rename` maps the score
directory back to the real name so everything downstream (graded tsv/json, provenance, report) sees it.

  sample_id.py id NAME                  print the internal ID (exit 2 if empty / reserved)
  sample_id.py check NAME [--others N1,N2,...]   same, plus collision check against other sample names
  sample_id.py rename SCORE_DIR SID NAME   rewrite SCORE_DIR/<SID> outputs to NAME (no-op if SID == NAME)
"""
import glob, gzip, json, os, re, sys

RESERVED = {"reference"}


def sid(name):
    s = re.sub(r"[^A-Za-z0-9]", "", name)
    if not s:
        raise SystemExit(f"sample name {name!r} has no alphanumeric character: cannot form a sampleset id")
    if s.lower() in RESERVED:
        raise SystemExit(f"sample name {name!r} maps to the reserved sampleset id {s!r}")
    return s


def check(name, others=()):
    s = sid(name)
    for o in others:
        if o != name and sid(o) == s:
            raise SystemExit(f"sample names {name!r} and {o!r} both map to internal id {s!r}: collision")
    return s


def _op(p, mode, gz=None):
    gz = p.endswith(".gz") if gz is None else gz
    return gzip.open(p, mode + "t", newline="") if gz else open(p, mode, newline="")


def _swap_fields(path, sep, old, new):
    tmp = path + ".tmp"
    n = 0
    with _op(path, "r") as fi, _op(tmp, "w", gz=path.endswith(".gz")) as fo:
        for line in fi:
            body = line.rstrip("\r\n")
            end = line[len(body):]
            f = body.split(sep)
            if old in f:
                f = [new if x == old else x for x in f]; n += 1
            fo.write(sep.join(f) + end)
    os.replace(tmp, path)
    return n


def _swap_json(o, old, new):
    if isinstance(o, dict):
        return {(new if k == old else k): _swap_json(v, old, new) for k, v in o.items()}
    if isinstance(o, list):
        return [_swap_json(v, old, new) for v in o]
    return new if o == old else o


def rename(scoredir, old, new):
    if old == new:
        return
    d_old, d_new = os.path.join(scoredir, old), os.path.join(scoredir, new)
    if not os.path.isdir(d_old):
        raise SystemExit(f"rename: expected {d_old} (internal sampleset id {old}); the score step wrote something else")
    pgs = glob.glob(os.path.join(d_old, "**", f"{old}_pgs.txt.gz"), recursive=True)
    if len(pgs) != 1:
        raise SystemExit(f"rename: no unique {old}_pgs.txt.gz under {d_old}: {pgs}")
    for p in glob.glob(os.path.join(d_old, "**", "*"), recursive=True):
        if not os.path.isfile(p):
            continue
        b = os.path.basename(p)
        if b.endswith("_pgs.txt.gz") or b.endswith("_popsimilarity.txt.gz"):
            n = _swap_fields(p, "\t", old, new)
            if n == 0:
                raise SystemExit(f"rename: {p} does not mention sample id {old}")
        elif b.endswith("_summary.csv"):
            _swap_fields(p, ",", old, new)
        elif b == "pop_summary.csv":
            _swap_fields(p, ",", old, new)
        elif b.endswith("_info.json.gz"):
            with gzip.open(p, "rt") as f:
                j = json.load(f)
            with gzip.open(p, "wt") as f:
                json.dump(_swap_json(j, old, new), f)
        elif b == "report.html":
            t = open(p, encoding="utf-8", errors="surrogateescape").read()
            open(p, "w", encoding="utf-8", errors="surrogateescape").write(re.sub(rf"\b{re.escape(old)}\b", new, t))
        if b.startswith(old + "_"):
            os.replace(p, os.path.join(os.path.dirname(p), new + b[len(old):]))
    os.replace(d_old, d_new)
    # consistency: the internal id must not survive as a sampleset / IID in the graded input
    p = glob.glob(os.path.join(d_new, "**", f"{new}_pgs.txt.gz"), recursive=True)[0]
    with gzip.open(p, "rt") as f:
        h = next(f).rstrip("\n").split("\t")
        i = h.index("sampleset")
        for line in f:
            if line.split("\t", i + 1)[i] == old:
                raise SystemExit(f"rename: internal id {old} still present in {p}")


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) >= 2 and a[0] == "id":
        print(sid(a[1]))
    elif len(a) >= 2 and a[0] == "check":
        others = a[3].split(",") if len(a) > 3 and a[2] == "--others" else []
        print(check(a[1], others))
    elif len(a) == 4 and a[0] == "rename":
        rename(a[1], a[2], a[3])
    else:
        sys.exit(__doc__)
