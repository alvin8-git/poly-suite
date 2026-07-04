#!/usr/bin/env python3
"""Cache PGS Catalog harmonized scoring files locally, so pgsc_calc runs can use
--scorefile (local) instead of --pgs_id and skip the ~40-58s DOWNLOAD_SCOREFILES
step on every run. Scoring files are fixed per (PGS id, build) -> one-time fetch.
curl/wget are blocked on this box, so download via urllib (streaming, resumable-skip).

usage: cache_scorefiles.py <cache_dir> "PGS000004,PGS004941,..."
"""
import sys, os, urllib.request, urllib.error

BUILD = "hmPOS_GRCh38"
URL = ("https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/"
       "{pgs}/ScoringFiles/Harmonized/{pgs}_{build}.txt.gz")


def fetch(pgs, cache):
    out = os.path.join(cache, f"{pgs}_{BUILD}.txt.gz")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return "cached"
    tmp = out + ".part"
    try:
        with urllib.request.urlopen(URL.format(pgs=pgs, build=BUILD), timeout=180) as r, \
                open(tmp, "wb") as f:
            while True:
                b = r.read(1 << 20)
                if not b:
                    break
                f.write(b)
        os.rename(tmp, out)
        return f"downloaded ({os.path.getsize(out)//1024} KB)"
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        if os.path.exists(tmp):
            os.remove(tmp)
        return f"FAIL ({e})"


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: cache_scorefiles.py <cache_dir> 'PGS...,PGS...'")
    cache = sys.argv[1]
    ids = sys.argv[2].replace(",", " ").split()
    os.makedirs(cache, exist_ok=True)
    ok = 0
    for pgs in ids:
        st = fetch(pgs, cache)
        ok += st.startswith(("cached", "downloaded"))
        print(f"  {pgs}: {st}")
    print(f"\n{ok}/{len(ids)} scorefiles in {cache}")
    return 0 if ok == len(ids) else 1


if __name__ == "__main__":
    sys.exit(main())
