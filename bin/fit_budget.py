#!/usr/bin/env python3
"""Fit the scorefile-bytes -> peak_rss curve and recommend a --max-batch-bytes B.

Input: a results dir from measure_mem.sh holding
  subsets.tsv          name<TAB>total_bytes<TAB>n_files   (one row per measured subset)
  <name>.trace.txt     that subset's pgsc_calc execution_trace (peak_rss per process)

The memory driver is total scorefile bytes handled in one pgsc_calc invocation
(FORMAT_SCOREFILES / MATCH_VARIANTS / MATCH_COMBINE all load their whole working
set). We fit overall_peak_GB = a*bytes_GB + b by least squares across the subsets,
then invert to the byte budget whose worst bin peaks at a memory target.

Usage: fit_budget.py <results_dir>
       fit_budget.py --self-test
"""
import sys, os, glob

MEM_STEPS = ("FORMAT_SCOREFILES", "MATCH_VARIANTS", "MATCH_COMBINE")
# memory targets (GB) -> parallelism story, on a 125 GB box
TARGETS = {"serial-safe (1 run)": 60.0, "2-3 wide parallel": 35.0}


def to_gb(s):
    """'108 GB' / '512 MB' / '-' -> float GB."""
    s = (s or "").strip()
    if not s or s == "-":
        return 0.0
    parts = s.split()
    num = float(parts[0])
    unit = parts[1].upper() if len(parts) > 1 else "B"
    return num * {"B": 1e-9, "KB": 1e-6, "MB": 1e-3, "GB": 1.0, "TB": 1e3}[unit]


def peak_from_trace(path):
    """Return (overall_peak_gb, {step: peak_gb}) from one execution_trace file."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ni, pi = header.index("name"), header.index("peak_rss")
        by_step = {}
        overall = 0.0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ni, pi):
                continue
            gb = to_gb(f[pi])
            overall = max(overall, gb)
            name = f[ni]
            for step in MEM_STEPS:
                if step in name:
                    by_step[step] = max(by_step.get(step, 0.0), gb)
    return overall, by_step


def lsq(xs, ys):
    """Least-squares slope,intercept for y = a*x + b."""
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    denom = n * sxx - sx * sx
    if denom == 0:
        raise ValueError("degenerate x (all subsets same size?) — cannot fit")
    a = (n * sxy - sx * sy) / denom
    b = (sy - a * sx) / n
    return a, b


def recommend(points):
    """points: list of (bytes, overall_peak_gb). Returns (a,b,budgets{label:bytes})."""
    xs = [p[0] / 1e9 for p in points]   # scorefile GB
    ys = [p[1] for p in points]         # peak RSS GB
    a, b = lsq(xs, ys)
    budgets = {}
    for label, target in TARGETS.items():
        gb_budget = (target - b) / a if a > 0 else float("inf")
        budgets[label] = max(0.0, gb_budget) * 1e9
    return a, b, budgets


def main(results_dir):
    tsv = os.path.join(results_dir, "subsets.tsv")
    rows = [l.split("\t") for l in open(tsv).read().splitlines() if l.strip()]
    points = []
    print(f"{'subset':<12} {'GB':>6} {'files':>6} {'overall':>8}  " + "  ".join(f"{s[:6]:>6}" for s in MEM_STEPS))
    for name, bytes_, nfiles in rows:
        tp = os.path.join(results_dir, f"{name}.trace.txt")
        if not os.path.exists(tp):
            print(f"{name:<12} (no trace — subset did not complete)")
            continue
        overall, by = peak_from_trace(tp)
        points.append((int(bytes_), overall))
        cells = "  ".join(f"{by.get(s, 0):6.1f}" for s in MEM_STEPS)
        print(f"{name:<12} {int(bytes_)/1e9:6.2f} {nfiles:>6} {overall:8.1f}  {cells}")
    if len(points) < 2:
        print("\nneed >=2 completed subsets to fit a curve", file=sys.stderr)
        sys.exit(1)
    a, b, budgets = recommend(points)
    print(f"\nfit: peak_rss(GB) = {a:.2f} * scorefile_GB + {b:.1f}   ({len(points)} points)")
    print("recommended --max-batch-bytes:")
    for label, B in budgets.items():
        print(f"  {label:<22} target -> B = {B/1e6:7.0f} MB  ({B/1e9:.2f} GB of scorefiles/bin)")
    # sanity: what the current 600 MB proxy implies
    proxy_peak = a * 0.6 + b
    print(f"\ncurrent 600 MB proxy -> predicted worst-bin peak {proxy_peak:.1f} GB")


def _self_test():
    import tempfile, shutil
    d = tempfile.mkdtemp()
    try:
        # synthetic: peak = 40*GB + 5, exactly. 3 subsets.
        subsets = [("s1", 0.5e9, 10), ("s2", 1.0e9, 20), ("s3", 2.0e9, 40)]
        with open(os.path.join(d, "subsets.tsv"), "w") as fh:
            for name, by, n in subsets:
                fh.write(f"{name}\t{int(by)}\t{n}\n")
                peak = 40 * (by / 1e9) + 5
                with open(os.path.join(d, f"{name}.trace.txt"), "w") as tr:
                    tr.write("task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tstart\trealtime\t%cpu\tpeak_rss\n")
                    # MATCH_COMBINE carries the overall peak; FORMAT a bit less
                    tr.write(f"1\th\t1\tFOO:MATCH_COMBINE (S)\tOK\t0\t-\t-\t1h\t100%\t{peak:.1f} GB\n")
                    tr.write(f"2\th\t2\tFOO:FORMAT_SCOREFILES (S)\tOK\t0\t-\t-\t1h\t100%\t{peak*0.8:.1f} GB\n")
                    tr.write("3\th\t3\tFOO:PLINK2_SCORE (S)\tOK\t0\t-\t-\t1h\t100%\t512 MB\n")
        pts = []
        for name, by, n in subsets:
            overall, by_step = peak_from_trace(os.path.join(d, f"{name}.trace.txt"))
            assert abs(overall - (40 * (by / 1e9) + 5)) < 0.11, overall
            assert "MATCH_COMBINE" in by_step and "FORMAT_SCOREFILES" in by_step
            pts.append((by, overall))
        a, b, budgets = recommend(pts)
        assert abs(a - 40) < 0.5 and abs(b - 5) < 0.5, (a, b)
        # target 35 GB -> (35-5)/40 = 0.75 GB = 750 MB
        assert abs(budgets["2-3 wide parallel"] - 0.75e9) < 5e6, budgets
        assert abs(to_gb("512 MB") - 0.512) < 1e-9 and to_gb("-") == 0.0 and abs(to_gb("2 TB") - 2000) < 1e-6
        print("fit_budget self-test PASSED")
    finally:
        shutil.rmtree(d)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
        _self_test()
    elif len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
