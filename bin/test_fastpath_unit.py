#!/usr/bin/env python3
"""Unit checks for the fastpath's two behaviour-preserving patches (run inside the pygscatalog image):

  docker run --rm -v /data:/data -w $PWD ghcr.io/pgscatalog/pygscatalog:pgscatalog-utils-2.0.4 \
      python bin/test_fastpath_unit.py

  1. fastpath_match._vector_encode_match_priority == pgscatalog label._encode_match_priority
     on every (best_match, exclude) combination including nulls.
  2. derive_intersect(): restricting the canonical intersect table to a sample's pvar and recomputing
     AAF / F_MISS_DOSAGE / PCA_ELIGIBLE reproduces what pgscatalog-intersect would write
     (tiny synthetic example with a missing genotype, an absent site and an ambiguous SNP).
Exit 1 on any failure.
"""
import gzip, json, os, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import itertools


def test_encode():
    import polars as pl
    import fastpath_match
    from pgscatalog.match.lib._match import label
    vals_b, vals_e = [True, False, None], [True, False, None]
    rows = list(itertools.product(vals_b, vals_e)) * 3
    df = pl.DataFrame({"best_match": [r[0] for r in rows], "exclude": [r[1] for r in rows]},
                      schema={"best_match": pl.Boolean, "exclude": pl.Boolean}).lazy()
    with pl.StringCache():
        old = label._encode_match_priority(df).collect()
        new = fastpath_match._vector_encode_match_priority(df).collect()
        assert old.columns == new.columns, (old.columns, new.columns)
        assert old.schema == new.schema, (old.schema, new.schema)
        assert old["match_status"].cast(pl.Utf8).to_list() == new["match_status"].cast(pl.Utf8).to_list()
    print("ok  encode_match_priority identical on", len(rows), "rows")


def test_derive_intersect():
    import shutil
    if not shutil.which("zstd"):
        raise ImportError("zstd binary not on PATH (run this check on the host)")
    import fastpath
    t = tempfile.mkdtemp()
    hdr = ["CHR:POS:A0:A1", "ID_REF", "REF_REF", "IS_INDEL", "STRANDAMB", "IS_MA_REF", "ID_TARGET", "REF_TARGET",
           "IS_MA_TARGET", "AAF", "F_MISS_DOSAGE", "SAME_REF", "PCA_ELIGIBLE"]
    canon = [
        ["1:10:A:C", "1:10:A:C", "A", "False", "False", "False", "1:10:A:C", "A", "False", "0.0", "0", "True", "True"],
        ["1:20:A:T", "1:20:A:T", "A", "False", "True", "False", "1:20:A:T", "A", "False", "0.0", "0", "True", "False"],
        ["1:30:C:G", "1:30:G:C", "G", "False", "False", "False", "1:30:C:G", "C", "False", "0.0", "0", "True", "True"],
        ["1:40:A:G", "1:40:A:G", "A", "False", "False", "True", "1:40:A:G", "A", "False", "0.0", "0", "True", "False"],
    ]
    cm = os.path.join(t, "c.txt.gz")
    with gzip.open(cm, "wt", newline="") as f:
        f.write("\r\n".join("\t".join(r) for r in [hdr] + canon) + "\r\n")
    open(os.path.join(t, "cc.txt"), "w").write("\n".join(["4", "99", "4"]))
    # sample: site 1:30 present with missing genotype (fmiss 1), 1:20 absent, 1:10 and 1:40 present
    pv = os.path.join(t, "s.pvar")
    open(pv, "w").write("#CHROM\tPOS\tID\tREF\tALT\n1\t10\t1:10:A:C\tA\tC\n1\t30\t1:30:C:G\tC\tG\n1\t40\t1:40:A:G\tA\tG\n")
    subprocess.run(["zstd", "-q", "-f", pv, "-o", pv + ".zst"], check=True)
    with gzip.open(os.path.join(t, "s.afreq.gz"), "wt") as f:
        f.write("#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n1\t1:10:A:C\tA\tC\t0.5\t2\n1\t1:30:C:G\tC\tG\t1\t2\n1\t1:40:A:G\tA\tG\t0\t2\n")
    with gzip.open(os.path.join(t, "s.vmiss.gz"), "wt") as f:
        f.write("#ID\tF_MISS_DOSAGE\tF_MISS\n1:10:A:C\t0\t0\n1:30:C:G\t1\t1\n1:40:A:G\t0\t0\n")
    cp = os.path.join(t, "canon.pvar")
    open(cp, "w").write("#CHROM\tPOS\tID\tREF\tALT\n" + "".join(
        f"1\t{p}\t{i}\t{i.split(':')[2]}\t{i.split(':')[3]}\n" for p, i in
        [(10, "1:10:A:C"), (20, "1:20:A:T"), (30, "1:30:C:G"), (40, "1:40:A:G")]))
    subprocess.run(["zstd", "-q", "-f", cp, "-o", cp + ".zst"], check=True)
    out = os.path.join(t, "o.txt.gz")
    st = fastpath.derive_intersect(cm, os.path.join(t, "cc.txt"), cp + ".zst", pv + ".zst", os.path.join(t, "s.afreq.gz"),
                                   os.path.join(t, "s.vmiss.gz"), out, os.path.join(t, "oc.txt"), os.path.join(t, "x"))
    got = [l.split("\t") for l in gzip.open(out, "rt", newline="").read().split("\r\n") if l]
    assert [g[1] for g in got[1:]] == ["1:10:A:C", "1:30:G:C", "1:40:A:G"], got      # 1:20 dropped (absent)
    byid = {g[6]: g for g in got[1:]}
    assert byid["1:10:A:C"][9:] == ["0.5", "0", "True", "True"], byid["1:10:A:C"]
    assert byid["1:30:C:G"][9:] == ["1.0", "1", "True", "False"], byid["1:30:C:G"]  # missing genotype -> not PCA eligible
    assert byid["1:40:A:G"][12] == "False"                                           # IS_MA_REF stays excluded
    assert open(os.path.join(t, "oc.txt")).read() == "3\n99\n3" and st["matched"] == 3
    # sample-only row (not in the canonical site set) is handed to the delta runner and merged in key order
    open(pv, "a").write("1\t25\t1:25:AT:A\tAT\tA\n")
    subprocess.run(["zstd", "-q", "-f", pv, "-o", pv + ".zst"], check=True)
    with gzip.open(os.path.join(t, "s.afreq.gz"), "at") as f:
        f.write("1\t1:25:AT:A\tAT\tA\t0.5\t2\n")
    with gzip.open(os.path.join(t, "s.vmiss.gz"), "at") as f:
        f.write("1:25:AT:A\t0\t0\n")
    drow = ["1:25:A:AT", "1:25:AT:A", "AT", "True", "False", "False", "1:25:AT:A", "AT", "False", "0.5", "0", "True", "False"]
    got_pos = []
    st = fastpath.derive_intersect(cm, os.path.join(t, "cc.txt"), cp + ".zst", pv + ".zst", os.path.join(t, "s.afreq.gz"),
                                   os.path.join(t, "s.vmiss.gz"), out, os.path.join(t, "oc.txt"), os.path.join(t, "x"),
                                   delta_runner=lambda pre, pos: (got_pos.extend(pos), [drow])[1])
    got = [l.split("\t") for l in gzip.open(out, "rt", newline="").read().split("\r\n") if l]
    assert got_pos == [("1", 25)] and st["sample_only_rows"] == 1 and st["matched"] == 4, (got_pos, st)
    assert [g[0] for g in got[1:]] == ["1:10:A:C", "1:25:A:AT", "1:30:C:G", "1:40:A:G"], [g[0] for g in got[1:]]
    print("ok  derive_intersect restriction + recomputed columns + delta merge")


def test_relabel_awk_matches_tool():
    """relabel_awk == pgscatalog-relabel byte for byte (decompressed), on a synthetic map + scorefile."""
    import shutil
    if not shutil.which("docker") or not shutil.which("zstd"):
        raise ImportError("needs docker and run on the host")
    import fastpath
    t = tempfile.mkdtemp()
    hdr = ["CHR:POS:A0:A1", "ID_REF", "REF_REF", "IS_INDEL", "STRANDAMB", "IS_MA_REF", "ID_TARGET", "REF_TARGET",
           "IS_MA_TARGET", "AAF", "F_MISS_DOSAGE", "SAME_REF", "PCA_ELIGIBLE"]
    rows = [["k", "1:10:A:C", "A", "False", "False", "False", "1:10:C:A", "C", "False", "0.1", "0", "True", "True"],
            ["k", "1:20:A:G", "A", "False", "False", "False", "1:20:A:G", "A", "False", "0.2", "0", "True", "True"]]
    with gzip.open(os.path.join(t, "m.gz"), "wt", newline="") as f:
        f.write("\r\n".join("\t".join(r) for r in [hdr] + rows) + "\r\n")
    with gzip.open(os.path.join(t, "s.gz"), "wt") as f:
        f.write("ID\teffect_allele\tPGS1\tPGS2\n1:10:C:A\tC\t0.5\t0\n1:20:A:G\tG\t-1e-05\t2.25\n")
    fastpath.relabel_awk(os.path.join(t, "m.gz"), os.path.join(t, "s.gz"), os.path.join(t, "awk.gz"))
    subprocess.run(fastpath.docker_base(t) + [fastpath.IMG_UTILS, "bash", "-c",
                    f"pgscatalog-relabel --maps {t}/m.gz --col_from ID_TARGET --col_to ID_REF --target_file {t}/s.gz "
                    f"--target_col ID --dataset reference --combined --outdir . && mv reference_ALL_relabelled.gz tool.gz"],
                   check=True, capture_output=True)
    a = gzip.open(os.path.join(t, "awk.gz"), "rb").read(); b = gzip.open(os.path.join(t, "tool.gz"), "rb").read()
    assert a == b, (a, b)
    print("ok  relabel_awk byte-identical to pgscatalog-relabel")


def test_sample_id_roundtrip():
    """AB_12.x -> internal AB12x -> score dir mapped back; collisions / empty / reserved fail loud."""
    import sample_id as sid
    assert sid.sid("AB_12.x") == "AB12x" and sid.sid("Proband01") == "Proband01"
    for bad in ("_.-", "reference", "Reference"):
        try:
            sid.sid(bad); raise AssertionError(bad)
        except SystemExit:
            pass
    try:
        sid.check("AB_12", ["AB.12"]); raise AssertionError("collision")
    except SystemExit:
        pass
    assert sid.check("AB_12", ["AB_12", "CD_3"]) == "AB12"
    t = tempfile.mkdtemp(); d = os.path.join(t, "AB12x"); os.makedirs(d + "/match"); os.makedirs(d + "/score")
    with gzip.open(d + "/score/AB12x_pgs.txt.gz", "wt") as f:
        f.write("sampleset\tFID\tIID\tPGS\tSUM\nAB12x\tAB12x\tAB12x\tPGS1\t1.5\nreference\tHG1\tHG1\tPGS1\t2\n")
    with gzip.open(d + "/score/AB12x_popsimilarity.txt.gz", "wt") as f:
        f.write("sampleset\tFID\tIID\tPC1\nAB12x\tAB12x\tAB12x\t0.1\n")
    with gzip.open(d + "/score/AB12x_info.json.gz", "wt") as f:
        f.write('{"AB12x": {"x": "AB12x"}}')
    open(d + "/score/pop_summary.csv", "w").write("Most similar population,AB12x,reference\nEUR,1 (100%),2 (50%)\n")
    open(d + "/match/AB12x_summary.csv", "w").write("dataset,accession,count\nAB12x,PGS1,3\n")
    open(d + "/score/report.html", "w").write("<h1>AB12x</h1> AB12xy")
    sid.rename(t, "AB12x", "AB_12.x")
    n = os.path.join(t, "AB_12.x")
    assert os.path.isdir(n) and not os.path.exists(d)
    rows = gzip.open(n + "/score/AB_12.x_pgs.txt.gz", "rt").read().split("\n")
    assert rows[1] == "AB_12.x\tAB_12.x\tAB_12.x\tPGS1\t1.5" and rows[2].startswith("reference\tHG1"), rows
    assert "AB_12.x" in gzip.open(n + "/score/AB_12.x_popsimilarity.txt.gz", "rt").read()
    assert json.load(gzip.open(n + "/score/AB_12.x_info.json.gz", "rt")) == {"AB_12.x": {"x": "AB_12.x"}}
    assert open(n + "/score/pop_summary.csv").readline().strip() == "Most similar population,AB_12.x,reference"
    assert open(n + "/match/AB_12.x_summary.csv").read().startswith("dataset,accession,count\nAB_12.x,")
    assert open(n + "/score/report.html").read() == "<h1>AB_12.x</h1> AB12xy"
    print("ok  sample id mapping (AB_12.x <-> AB12x)")


if __name__ == "__main__":
    for fn in (test_sample_id_roundtrip, test_encode, test_derive_intersect, test_relabel_awk_matches_tool):
        try:
            fn()
        except ImportError as e:
            print(f"SKIP {fn.__name__}: {e}")
        except AssertionError as e:
            print(f"FAIL {fn.__name__}: {e}"); sys.exit(1)
    print("fastpath unit checks passed")
