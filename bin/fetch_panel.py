#!/usr/bin/env python3
"""Stream-download the pgsc_calc ancestry reference panel (curl/wget are blocked
on this box). Resumable via HTTP Range; verifies size; prints progress.

usage: fetch_panel.py [URL] [OUT]
default: HGDP+1kGP v1 panel -> /data/alvin/ref/pgsc/pgsc_HGDP+1kGP_v1.tar.zst
"""
import sys, os, urllib.request

URL = sys.argv[1] if len(sys.argv) > 1 else \
    "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources/pgsc_HGDP+1kGP_v1.tar.zst"
OUT = sys.argv[2] if len(sys.argv) > 2 else \
    "/data/alvin/ref/pgsc/pgsc_HGDP+1kGP_v1.tar.zst"

os.makedirs(os.path.dirname(OUT), exist_ok=True)
tmp = OUT + ".part"
pos = os.path.getsize(tmp) if os.path.exists(tmp) else 0
req = urllib.request.Request(URL, headers={"Range": f"bytes={pos}-"} if pos else {})

CHUNK = 8 << 20  # 8 MB
with urllib.request.urlopen(req, timeout=120) as r:
    partial = (r.status == 206)
    if pos and not partial:          # server ignored Range -> restart clean
        pos = 0
    total = int(r.headers.get("Content-Length", 0)) + (pos if partial else 0)
    mode = "ab" if (pos and partial) else "wb"
    done = pos if partial else 0
    mark = 0
    with open(tmp, mode) as f:
        while True:
            buf = r.read(CHUNK)
            if not buf:
                break
            f.write(buf)
            done += len(buf)
            if done - mark >= (512 << 20):   # every ~512 MB
                mark = done
                pct = f"{done/total*100:.0f}%" if total else "?"
                print(f"  {done/1e9:.1f} / {total/1e9:.1f} GB ({pct})", flush=True)

sz = os.path.getsize(tmp)
if total and sz < total:
    print(f"INCOMPLETE: {sz}/{total} bytes — rerun to resume", flush=True)
    sys.exit(1)
os.rename(tmp, OUT)
print(f"DONE {OUT} ({sz/1e9:.2f} GB)", flush=True)
