#!/usr/bin/env python3
"""poly-suite unit tests — pure-logic checks for the modules (no pytest, no
network, no BAM needed). Run: python3 bin/tests.py  (or bin/selftest.sh).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import absolute_risk as AR
import consensus as C
import ensemble as E
import validate_contract as V
import validate_calibration as VC
import select_pgs as SP
import grade_pgs as G


def test_absolute_risk_monotonic_bounded():
    eff, base = AR.load_effects(), AR.load_baselines()
    f = lambda p: AR.estimate("breast cancer", p, eff=eff, base=base)["absolute_risk"]
    assert f(20) < f(50) < f(80) < f(95) < 1.0, "risk must increase with percentile, stay <1"


def test_absolute_risk_odds_scale_caps_common_disease():
    # CAD baseline 49% (male) at extreme percentile must stay <1 (the odds-scale fix)
    eff, base = AR.load_effects(), AR.load_baselines()
    hi = AR.estimate("coronary artery disease", 99.9, pgs_id="PGS000018",
                     sex="male", eff=eff, base=base)["absolute_risk"]
    assert 0.49 < hi < 0.95, f"CAD 99.9th should be high-but-bounded, got {hi}"


def test_median_below_baseline():
    # mean-RR=1 calibration => median (50th) sits just below baseline
    eff, base = AR.load_effects(), AR.load_baselines()
    b = base[("breast cancer", "female", "overall")][0]
    assert AR.estimate("breast cancer", 50, eff=eff, base=base)["absolute_risk"] < b


def test_sex_precedence():
    eff, base = AR.load_effects(), AR.load_baselines()
    # trait-mandated sex wins even if sample sex differs
    assert AR.estimate("prostate cancer", 90, sex="female", eff=eff, base=base)["sex"] == "male"
    # sex-dimorphic trait uses sample sex
    assert AR.estimate("coronary artery disease", 90, pgs_id="PGS000018",
                       sex="female", eff=eff, base=base)["sex"] == "female"
    # sex-dimorphic trait with unknown sex -> no estimate (honest)
    assert AR.estimate("coronary artery disease", 90, pgs_id="PGS000018",
                       eff=eff, base=base) is None


def test_absolute_risk_needs_percentile():
    assert AR.estimate("breast cancer", None) is None


def test_consensus_concordance():
    c = C.consensus([
        {"trait": "CAD", "pgs_id": "A", "percentile": 82},
        {"trait": "CAD", "pgs_id": "B", "percentile": 88},
        {"trait": "T2D", "pgs_id": "C", "percentile": 15},
        {"trait": "T2D", "pgs_id": "D", "percentile": 71},
        {"trait": "LDL", "pgs_id": "E", "percentile": 60},
    ])
    assert c["CAD"]["concordance"] > 0.9 and c["CAD"]["tertile_concordant"]
    assert c["T2D"]["concordance"] < 0.6 and not c["T2D"]["tertile_concordant"]
    assert c["LDL"]["concordance"] is None and c["LDL"]["n_scores"] == 1


def test_consensus_ignores_uncalibrated():
    c = C.consensus([{"trait": "X", "pgs_id": "A", "percentile": None},
                     {"trait": "X", "pgs_id": "B", "percentile": None}])
    assert c["X"]["n_scores"] == 2 and c["X"]["concordance"] is None


def test_ensemble_meta():
    rows = [{"trait": "CAD", "pgs_id": "A", "z": 1.5},
            {"trait": "CAD", "pgs_id": "B", "z": 2.5},   # equal meta_z = 2.0
            {"trait": "T2D", "pgs_id": "C", "z": 0.3}]   # single -> no meta
    m = E.meta(rows)
    assert abs(m["CAD"]["meta_z"] - 2.0) < 1e-9
    assert m["CAD"]["meta_percentile"] > 95            # z=2 -> ~97.7th
    assert "T2D" not in m
    mw = E.meta(rows, weights={"A": 1.0, "B": 3.0})    # (1.5*1+2.5*3)/4 = 2.25
    assert abs(mw["CAD"]["meta_z"] - 2.25) < 1e-9 and mw["CAD"]["weighting"] == "evidence"


def test_ensemble_ignores_uncalibrated():
    assert E.meta([{"trait": "X", "pgs_id": "A", "z": None},
                   {"trait": "X", "pgs_id": "B", "z": None}]) == {}


def test_validate_contract():
    header = list(G.CONTRACT_COLS)
    good = {c: "NA" for c in header}
    good.update(evidence_grade="A", portability_flag="false", allowed_statement="ok",
                percentile="90", absolute_risk="0.2", match_rate="0.99")
    assert V.validate([good], header) == []
    # absolute risk without a percentile must be caught
    bad = dict(good); bad["percentile"] = "NA"
    assert V.validate([bad], header)
    # bad grade + empty caveat -> >=2 problems
    bad2 = dict(good); bad2["evidence_grade"] = "Z"; bad2["allowed_statement"] = ""
    assert len(V.validate([bad2], header)) >= 2
    # low coverage without a coverage caveat must be caught
    bad3 = dict(good); bad3["match_rate"] = "0.3"; bad3["allowed_statement"] = "no note"
    assert V.validate([bad3], header)
    # missing a required column
    assert V.validate([good], header[:-1])


def test_calibration_uniformity():
    uniform = [i % 100 + 0.5 for i in range(2000)]     # even across 0-100
    dev, _ = VC.uniformity(uniform)
    assert dev < 0.02
    skewed = [5.0] * 1000 + [15.0] * 1000              # all in bottom 2 deciles
    dev2, _ = VC.uniformity(skewed)
    assert dev2 > 0.1
    assert VC.uniformity([1, 2, 3]) is None            # too few


def test_ensemble_from_calibrated():
    rows = []
    for iid, z in [("r1", -2), ("r2", -1), ("r3", 0), ("r4", 1)]:
        rows += [{"sampleset": "reference", "IID": iid, "PGS": "A", "z": z},
                 {"sampleset": "reference", "IID": iid, "PGS": "B", "z": z}]
    rows += [{"sampleset": "S", "IID": "t", "PGS": "A", "z": 0.5},
             {"sampleset": "S", "IID": "t", "PGS": "B", "z": 0.5}]
    m = E.meta_from_calibrated(rows, {"A": "CAD", "B": "CAD"})
    assert m["CAD"]["exact_percentile"] == 75.0        # target 0.5 > 3 of 4 ref meta_z
    assert m["CAD"]["n_scores"] == 2 and m["CAD"]["n_reference"] == 4
    assert E.meta_from_calibrated(rows, {"A": "CAD"}) == {}   # only 1 pgs for the trait


def test_launch_set_tiers():
    ls = SP.LAUNCH_SET
    assert all(t in ("core", "extended", "gated") for t, _, _ in ls)
    assert all(tid and "_" in tid for _, _, tid in ls)          # ontology-id shaped
    labels = [l for _, l, _ in ls]
    assert len(labels) == len(set(labels))                      # no duplicate traits
    core = [l for t, l, _ in ls if t == "core"]
    assert 20 <= len(core) <= 30                                # ~25 core
    # tier filter is monotone: core ⊂ extended ⊂ all
    assert SP._TIERS["core"] < SP._TIERS["extended"] < SP._TIERS["all"]


def test_grade_downgrade():
    assert G.downgrade("A") == "B"
    assert G.downgrade("A", 2) == "C"
    assert G.downgrade("D") == "D"          # floor
    assert G.downgrade("C", 5) == "D"       # clamps


def test_norm_pgs():
    assert G.norm_pgs("PGS000018_hmPOS_GRCh38") == "PGS000018"
    assert G.norm_pgs("PGS000004") == "PGS000004"


def test_sample_sex_parsing(tmp="/tmp/poly_sex_test"):
    os.makedirs(tmp, exist_ok=True)
    open(os.path.join(tmp, "sample_sex.txt"), "w").write("female\n")
    assert G.sample_sex(tmp) == "female"
    open(os.path.join(tmp, "sample_sex.txt"), "w").write("garbage\n")
    assert G.sample_sex(tmp) is None        # rejects non male/female


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print("PASS", t.__name__)
    print(f"\n{len(tests)} unit tests passed")
