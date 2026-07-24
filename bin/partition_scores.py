#!/usr/bin/env python3
"""Partition scorefiles into size-aware bins for the --batch score step.

Peak RSS of pgscatalog-format/-match scales with total scorefile bytes handled in
one pgsc_calc invocation (the launch set's ~2.4 GB of scorefiles drove a 108 GB
peak). Splitting the score set across invocations caps the per-invocation working
set. This does first-fit-decreasing bin packing under a byte budget B: biggest
scores first so the ~8 genome-wide ones land near-alone; the small tail fills bins.

Output: one line per bin, tab-separated absolute scorefile paths. Deterministic
(sort by size desc, then name) so -resume cache keys are stable across reruns.
A single score larger than B gets its own bin (never dropped).

Usage: partition_scores.py <scorefile_dir> <max_bytes>   [glob defaults to *.txt.gz]
       partition_scores.py --self-test
"""
import sys, os, glob


def partition(paths, max_bytes):
    items = sorted(((os.path.getsize(p), p) for p in paths), key=lambda t: (-t[0], t[1]))
    bins = []  # each: [total_bytes, [paths]]
    for size, p in items:
        for b in bins:
            if b[0] + size <= max_bytes:
                b[0] += size
                b[1].append(p)
                break
        else:
            bins.append([size, [p]])  # new bin (also the oversized-single case)
    # invariant: every input lands in exactly one bin
    placed = [p for b in bins for p in b[1]]
    assert sorted(placed) == sorted(paths), "partition lost or duplicated a scorefile"
    return bins


def main(dir_, max_bytes, pattern="*.txt.gz"):
    paths = [os.path.abspath(p) for p in glob.glob(os.path.join(dir_, pattern))]
    if not paths:
        print(f"no scorefiles matching {pattern} in {dir_}", file=sys.stderr)
        sys.exit(1)
    bins = partition(paths, max_bytes)
    for b in bins:
        print("\t".join(sorted(b[1])))
    print(f"partitioned {len(paths)} scorefiles into {len(bins)} bins "
          f"(budget {max_bytes/1e6:.0f} MB, worst bin {max(b[0] for b in bins)/1e6:.0f} MB)",
          file=sys.stderr)


def _self_test():
    # sizes in bytes; FFD with B=100 should pack without exceeding B (except oversized singles)
    import tempfile, shutil
    d = tempfile.mkdtemp()
    try:
        sizes = {"a": 90, "b": 60, "c": 60, "d": 30, "e": 150, "f": 10}
        for name, n in sizes.items():
            open(os.path.join(d, f"{name}.txt.gz"), "wb").write(b"x" * n)
        bins = partition([os.path.join(d, f"{n}.txt.gz") for n in sizes], 100)
        # every file placed exactly once
        placed = [os.path.basename(p) for b in bins for p in b[1]]
        assert sorted(placed) == sorted(f"{n}.txt.gz" for n in sizes), placed
        # no bin exceeds B unless it is a single oversized score
        for total, ps in bins:
            assert total <= 100 or len(ps) == 1, f"bin {ps} = {total} > B and not a single"
        # the oversized 'e' (150) is alone
        assert any(len(ps) == 1 and os.path.basename(ps[0]) == "e.txt.gz" for _, ps in bins)
        # FFD: a=90 alone (next smallest 60 won't fit), b=60+d=30 together, c=60+f=10
        print(f"partition self-test PASSED ({len(bins)} bins)")
    finally:
        shutil.rmtree(d)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
        _self_test()
    elif len(sys.argv) >= 3:
        main(sys.argv[1], int(sys.argv[2]), *(sys.argv[3:4]))
    else:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
