#!/usr/bin/env python3
"""poly-suite fast per-sample polygenic scoring path ("fastpath").

Gives the SAME scores as the pgsc_calc path (bin/run.sh default) while doing the
sample-independent work once.  Two subcommands:

  build  one-off, per (PGS list, score files, target site set, reference panel):
         formats the scoring files, imports the canonical target site set, matches
         every scoring variant against it, intersects it with the reference panel.
         Result: cache/fastpath/<id>/ (see docs/fastpath.md).
  score  per sample: PLINK2 import of the prepped VCF -> restrict the cached match
         candidates to the sample's own variant list -> pgscatalog label/filter/write
         -> relabel -> reference-panel + sample scoring -> ancestry -> aggregate.
         Writes the pgsc_calc output layout so grade_pgs.py / OmniGen need no change.

Every heavy step is the PGS Catalog's own command (pgscatalog-format / -match /
-intersect / -matchmerge / -relabel / -aggregate / -ancestry-adjust, plink2, fraposa)
from the SAME container images and with the SAME flags pgsc_calc 2.3.0 uses; this file
only sequences them and swaps the sample-independent parts for cached results.

A missing or mismatching cache is a hard error (exit 2); there is no fallback.
usage:
  fastpath.py build --pgs "PGS...,PGS..." --scorefile-cache DIR --panel panel.tar.zst
                    --ref hg38.fa [--cache-root cache/fastpath] [--threads 8] [--mem-gb 100]
  fastpath.py score --sample S --prepped PREFIX_OR_VCF --outdir DIR --pgs "..." ...
"""
import argparse, glob, gzip, hashlib, json, os, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
IMG_UTILS = os.environ.get("PGSC_UTILS_IMG", "ghcr.io/pgscatalog/pygscatalog:pgscatalog-utils-2.0.4")
IMG_PLINK2 = os.environ.get("PLINK2_IMG", "ghcr.io/pgscatalog/plink2:2.00a5.10")
IMG_FRAPOSA = os.environ.get("FRAPOSA_IMG", "ghcr.io/pgscatalog/fraposa_pgsc:v1.0.2")
IMG_ZSTD = os.environ.get("ZSTD_IMG", "ghcr.io/pgscatalog/zstd:2-beta")
PGSC_ASSETS = os.environ.get("PGSC_ASSETS", os.path.expanduser("~/.nextflow/assets/pgscatalog/pgsc_calc"))
FASTPATH_VERSION = "1"
TARGET_BUILD = "GRCh38"
MIN_OVERLAP = "0.1"        # run.sh passes --min_overlap 0.1 to pgsc_calc
STAGES = ["targets", "canon_import", "format", "match", "ref_extract", "intersect", "ref_index", "intersect_stats", "label", "wide"]


def die(msg, code=2):
    sys.stderr.write(f"[fastpath] FATAL: {msg}\n")
    sys.exit(code)


def log(msg):
    sys.stderr.write(f"[fastpath {time.strftime('%H:%M:%S')}] {msg}\n")
    sys.stderr.flush()


def sha256_file(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def docker_base(workdir, name=None):
    """`docker run` prefix honouring the server patch's env vars (same defaults as apply_server_patches.sh P1):
    OG_MOUNT (default /data) is mounted at the same path; OG_DOCKER_USER (e.g. "--user 1006:1007") is split into args."""
    mount = os.environ.get("OG_MOUNT", "/data")
    cmd = ["docker", "run", "--rm"] + (["--name", name] if name else []) + os.environ.get("OG_DOCKER_USER", "").split()
    return cmd + ["-v", f"{mount}:{mount}", "-w", workdir]


def dkr(image, script, workdir, mem_gb=None, env=None, cpus=None, name=None):
    """Run a bash script in a container with $OG_MOUNT (default /data) mounted; raises on non-zero."""
    cmd = docker_base(workdir)
    if mem_gb:
        cmd += ["--memory", f"{mem_gb}g", "--memory-swap", f"{mem_gb}g"]
    if cpus:
        cmd += ["--cpus", str(cpus)]
    for k, v in (env or {}).items():
        cmd += ["-e", f"{k}={v}"]
    cmd += [image, "bash", "-euo", "pipefail", "-c", script]
    t0 = time.time()
    r = subprocess.run(cmd)
    if r.returncode != 0:
        die(f"container step failed (rc={r.returncode}): {name or script[:80]}", 3)
    return time.time() - t0


def tool_versions():
    return {"fastpath": FASTPATH_VERSION, "pgsc_calc_assets": PGSC_ASSETS, "utils_image": IMG_UTILS,
            "plink2_image": IMG_PLINK2, "fraposa_image": IMG_FRAPOSA,
            "pgsc_calc_revision": _git_rev(PGSC_ASSETS)}


def _git_rev(d):
    try:
        return subprocess.check_output(["git", "-C", d, "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def panel_identity(panel):
    st = os.stat(panel)
    h = hashlib.sha256()
    with open(panel, "rb") as f:                       # head+tail+size: cheap, catches a swapped panel
        h.update(f.read(1 << 26)); f.seek(max(0, st.st_size - (1 << 26))); h.update(f.read(1 << 26))
    return {"path": os.path.abspath(panel), "size": st.st_size, "mtime_ns": st.st_mtime_ns,
            "head_tail_sha256": h.hexdigest()}


def score_paths(cache_dir, pgs):
    out = []
    for p in pgs:
        f = os.path.join(cache_dir, f"{p}_hmPOS_GRCh38.txt.gz")
        if not os.path.isfile(f):
            die(f"scorefile cache does not cover {p}: {f} (run bin/cache_scorefiles.py)")
        out.append(os.path.realpath(f))
    return out


def compute_key(pgs, sf_sha, targets_sha, panel_id, versions):
    blob = json.dumps({"pgs": pgs, "scorefiles": sf_sha, "targets": targets_sha,
                       "panel": {"size": panel_id["size"], "head_tail_sha256": panel_id["head_tail_sha256"]},
                       "tools": {k: v for k, v in versions.items() if k != "pgsc_calc_assets"},
                       "min_overlap": MIN_OVERLAP, "build": TARGET_BUILD, "stages": STAGES}, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()


def save_manifest(d, m):
    tmp = os.path.join(d, "manifest.json.tmp")
    json.dump(m, open(tmp, "w"), indent=1, sort_keys=True)
    os.replace(tmp, os.path.join(d, "manifest.json"))


# ------------------------------------------------------------------------------ build
def cmd_build(a):
    pgs = sorted(x for x in a.pgs.replace(",", " ").split() if x)
    sfs = score_paths(a.scorefile_cache, pgs)
    log(f"hashing {len(sfs)} scoring files")
    sf_sha = {p: sha256_file(f) for p, f in zip(pgs, sfs)}
    pre = hashlib.sha256(json.dumps([sf_sha, os.path.realpath(a.ref)], sort_keys=True).encode()).hexdigest()[:16]
    root = os.path.abspath(a.cache_root)
    cdir = os.path.join(root, f"fp-{pre}")
    os.makedirs(cdir, exist_ok=True)
    mpath = os.path.join(cdir, "manifest.json")
    versions = tool_versions()
    panel_id = panel_identity(a.panel)
    m = json.load(open(mpath)) if os.path.exists(mpath) else {"stages": {}, "complete": False}
    m.update({"pgs": pgs, "scorefile_sha256": sf_sha, "panel": panel_id, "tools": versions,
              "min_overlap": MIN_OVERLAP, "target_build": TARGET_BUILD, "reference_fasta": os.path.realpath(a.ref),
              "pgs_meta_note": "score metadata (traits/grades) is not an input of any cached product; not part of the key"})
    m.setdefault("stages", {})
    save_manifest(cdir, m)
    thr, mem = str(a.threads), a.mem_gb

    def stage(name, fn):
        if m["stages"].get(name, {}).get("done"):
            log(f"stage {name}: cached"); return
        log(f"stage {name}: start")
        t0 = time.time()
        extra = fn() or {}
        m["stages"][name] = {"done": True, "wall_s": round(time.time() - t0, 1), **extra}
        save_manifest(cdir, m)
        log(f"stage {name}: done in {m['stages'][name]['wall_s']} s")

    tdir = os.path.join(cdir, "target"); os.makedirs(tdir, exist_ok=True)
    targets = os.path.join(tdir, "targets.tsv")

    def s_targets():
        env = dict(os.environ)
        subprocess.run([sys.executable, os.path.join(HERE, "scoring_targets.py"), targets, a.ref] + sfs,
                       check=True, env=env)
        m["target_sites_sha256"] = sha256_file(targets)
        m["n_target_sites"] = sum(1 for _ in open(targets))
    stage("targets", s_targets)
    key = compute_key(pgs, sf_sha, m["target_sites_sha256"], panel_id, versions)
    m["key"] = key; save_manifest(cdir, m)
    log(f"cache id fp-{pre}  key {key[:16]}")

    canon = os.path.join(tdir, f"{TARGET_BUILD}_CANON_ALL")

    def s_canon():
        vcf = os.path.join(tdir, "canon.vcf")
        with open(targets) as f, open(vcf, "w") as o:
            o.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCANON\n")
            for line in f:
                c, p, al = line.rstrip("\n").split("\t")
                ref, alts = al.split(",", 1)
                o.write(f"{c}\t{p}\t.\t{ref}\t{alts}\t.\t.\t.\tGT\t0/0\n")
        subprocess.run(["bgzip", "-f", vcf], check=True)
        dkr(IMG_PLINK2, plink2_vcf_cmd(os.path.join(tdir, "canon.vcf.gz"), canon, thr, mem), tdir,
            name="canonical PLINK2_VCF")
        os.remove(os.path.join(tdir, "canon.vcf.gz"))
    stage("canon_import", s_canon)

    fdir = os.path.join(cdir, "formatted"); os.makedirs(fdir, exist_ok=True)

    def s_format():
        # FORMAT_SCOREFILES (pgsc_calc modules/local/format_scorefiles.nf, liftover=false)
        sdir = os.path.join(cdir, "raw"); shutil.rmtree(sdir, ignore_errors=True); os.makedirs(sdir)
        for f in sfs:
            os.symlink(f, os.path.join(sdir, os.path.basename(f)))
        shutil.rmtree(os.path.join(cdir, "fmt_work"), ignore_errors=True)
        os.makedirs(os.path.join(cdir, "fmt_work"))
        dkr(IMG_UTILS, f"mkdir formatted; pgscatalog-format -s {sdir}/*.txt.gz -t {TARGET_BUILD} --threads {thr} "
                       f"-o formatted/ -l log_scorefiles.json -v", os.path.join(cdir, "fmt_work"), mem_gb=mem,
            name="FORMAT_SCOREFILES")
        for f in glob.glob(os.path.join(cdir, "fmt_work", "formatted", "*")):
            shutil.move(f, os.path.join(fdir, os.path.basename(f)))
        shutil.rmtree(os.path.join(cdir, "fmt_work")); shutil.rmtree(sdir)
        n = len(glob.glob(os.path.join(fdir, "normalised_*")))
        if n != len(pgs):
            die(f"FORMAT_SCOREFILES wrote {n} files for {len(pgs)} scores")
        return {"n_formatted": n}
    stage("format", s_format)

    mdir = os.path.join(cdir, "match"); os.makedirs(mdir, exist_ok=True)

    def s_match():
        # MATCH_VARIANTS (match_variants.nf) against the canonical target; candidates only (--only_match)
        formatted = sorted(glob.glob(os.path.join(fdir, "normalised_*")))
        shutil.rmtree(os.path.join(mdir, "w"), ignore_errors=True); os.makedirs(os.path.join(mdir, "w"))
        dkr(IMG_UTILS, f"export POLARS_MAX_THREADS={thr}; pgscatalog-match --dataset CANON "
                       f"--scorefile {' '.join(formatted)} --target {canon}.pvar.zst --only_match "
                       f"--outdir $PWD -v", os.path.join(mdir, "w"), mem_gb=mem, name="MATCH_VARIANTS")
        ipcs = glob.glob(os.path.join(mdir, "w", "*.ipc.zst"))
        if len(ipcs) != 1:
            die(f"expected one candidate file, got {ipcs}")
        shutil.move(ipcs[0], os.path.join(mdir, "candidates.ipc.zst"))
        shutil.rmtree(os.path.join(mdir, "w"))
        return {"candidates_bytes": os.path.getsize(os.path.join(mdir, "candidates.ipc.zst"))}
    stage("match", s_match)

    rdir = os.path.join(cdir, "ref"); os.makedirs(rdir, exist_ok=True)

    def s_ref():
        # EXTRACT_DATABASE once (the per-run 55 s / 16 GB extraction of pgsc_calc)
        subprocess.run(f"tar --use-compress-program=zstd -xf '{os.path.abspath(a.panel)}' --wildcards "
                       f"'{TARGET_BUILD}*' meta.txt 2>/dev/null", shell=True, cwd=rdir)
        if open(os.path.join(rdir, "meta.txt")).read().strip() != "v0.1":
            die("reference panel format version is not v0.1")
        for f in glob.glob(os.path.join(rdir, "._*")):
            os.remove(f)
        return {"files": sorted(os.listdir(rdir))}
    stage("ref_extract", s_ref)

    idir = os.path.join(cdir, "intersect"); os.makedirs(idir, exist_ok=True)

    def s_intersect():
        # INTERSECT_VARIANTS (ancestry/intersect_variants.nf) on the canonical target (once)
        w = os.path.join(idir, "w"); shutil.rmtree(w, ignore_errors=True); os.makedirs(w)
        refvar = glob.glob(os.path.join(rdir, f"{TARGET_BUILD}_*_ALL.pvar.zst"))[0]
        dkr(IMG_UTILS, f"pgscatalog-intersect --ref {refvar} --target {canon}.pvar.zst --chrom ALL "
                       f"--maf_target 0 --geno_miss 0.1 --outdir . -v", w, name="INTERSECT_VARIANTS")
        shutil.move(os.path.join(w, "matched_variants.txt.gz"), os.path.join(idir, "matched_variants.txt.gz"))
        shutil.copy(os.path.join(w, "intersect_counts_ALL.txt"), os.path.join(idir, "intersect_counts_ALL.txt"))
        shutil.rmtree(w)
        # filter_ids exactly as match_combine.nf derives it (ref-biallelic target IDs); superset for every sample
        subprocess.run("awk 'NR > 1 && $6 == \"False\" { print $7 } ' <(zcat %s) | gzip > %s" % (
            os.path.join(idir, "matched_variants.txt.gz"), os.path.join(idir, "filter_ids.txt.gz")),
            shell=True, check=True, executable="/bin/bash")
    stage("intersect", s_intersect)

    def s_ref_index():
        # random access to the reference pvar (per-sample delta intersect reads only a few positions)
        rv = os.path.join(rdir, f"{TARGET_BUILD}_HGDP+1kGP_ALL.pvar.zst")
        gz = os.path.join(rdir, f"{TARGET_BUILD}_HGDP+1kGP_ALL.pvar.gz")
        subprocess.run(f"zstd -dc '{rv}' | bgzip -@ 4 -c > '{gz}'", shell=True, check=True, executable="/bin/bash")
        subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", gz], check=True)
        hdr = subprocess.run(f"zcat '{gz}' | head -n 5000 | grep -E '^#' ", shell=True, capture_output=True, text=True).stdout
        open(os.path.join(rdir, "pvar.header"), "w").write(hdr)
    stage("ref_index", s_ref_index)

    def s_intersect_stats():
        import csv
        seen, dup = set(), 0
        with gzip.open(os.path.join(idir, "matched_variants.txt.gz"), "rt", newline="") as f:
            r = csv.reader(f, delimiter="\t"); next(r)
            for row in r:
                if row[0] in seen:
                    dup += 1
                seen.add(row[0])
        return {"canonical_matched_rows": len(seen) + dup, "canonical_duplicate_keys": dup}
    stage("intersect_stats", s_intersect_stats)

    ldir = os.path.join(cdir, "label")

    def s_label():
        # canonical labelling (the heavy half of MATCH_COMBINE), once, in memory-bounded accession bins
        shutil.rmtree(ldir, ignore_errors=True); os.makedirs(ldir)
        formatted = sorted(glob.glob(os.path.join(fdir, "normalised_*")))
        dkr(IMG_UTILS, f"export POLARS_MAX_THREADS={thr}; python {os.path.join(HERE, 'fastpath_label.py')} "
                       f"{os.path.join(mdir, 'candidates.ipc.zst')} {os.path.join(idir, 'filter_ids.txt.gz')} {ldir} 16 "
                       + " ".join(formatted), cdir, mem_gb=mem, name="canonical label")
        return {"label_bytes": sum(os.path.getsize(f) for f in glob.glob(os.path.join(ldir, "*")))}
    stage("label", s_label)

    def s_wide():
        # canonical wide plink scoring files + split statistics (the per-sample fast path deletes rows from these)
        wdir = os.path.join(ldir, "wide"); shutil.rmtree(wdir, ignore_errors=True); os.makedirs(wdir)
        open(os.path.join(wdir, "absent.txt"), "w").close()
        dkr(IMG_UTILS, f"export POLARS_MAX_THREADS={thr}; python {os.path.join(HERE, 'fastpath_delta.py')} {ldir} "
                       f"{wdir}/absent.txt - {os.path.join(idir, 'filter_ids.txt.gz')} {wdir} CANON {MIN_OVERLAP} canon",
            cdir, mem_gb=mem, name="canonical wide")
        return {"wide_files": sorted(f for f in os.listdir(wdir) if f.endswith(".scorefile.gz"))}
    stage("wide", s_wide)

    for st in STAGES:
        if not m["stages"].get(st, {}).get("done"):
            die(f"stage {st} not complete")
    m["complete"] = True
    m["built_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    save_manifest(cdir, m)
    print(cdir)



# ----------------------------------------------------------------------------- score
def tree_rss_gb(root):
    """GB of resident memory of a process and all its descendants (reads /proc)."""
    todo, seen, kb = [root], set(), 0
    while todo:
        p = todo.pop()
        if p in seen:
            continue
        seen.add(p)
        try:
            for line in open(f"/proc/{p}/status"):
                if line.startswith("VmRSS:"):
                    kb += int(line.split()[1])
                    break
            for t in os.listdir(f"/proc/{p}/task"):
                todo += [int(x) for x in open(f"/proc/{p}/task/{t}/children").read().split()]
        except OSError:
            pass
    return kb * 1024 / 1e9


class Step:
    """Wall-time + peak container memory for one step (docker stats polling)."""
    log = []


def run_step(name, image, script, workdir, mem_gb=None, env=None):
    import threading
    os.makedirs(workdir, exist_ok=True)
    cname = f"fp_{os.getpid()}_{name}".replace("/", "_")
    cmd = docker_base(workdir, cname)
    if mem_gb:
        cmd += ["--memory", f"{mem_gb}g", "--memory-swap", f"{mem_gb}g"]
    for k, v in (env or {}).items():
        cmd += ["-e", f"{k}={v}"]
    cmd += [image, "bash", "-euo", "pipefail", "-c", script]
    peak = [0.0]
    stop = threading.Event()

    def poll():
        # peak resident set of the step: sum of VmRSS over the container's process tree (found via the
        # container init pid), max over 2 s samples.  (docker stats/top are not per-container here.)
        pid = None
        while not stop.wait(2.0):
            try:
                if pid is None:
                    o = subprocess.run(["docker", "inspect", "-f", "{{.State.Pid}}", cname], capture_output=True,
                                       text=True, timeout=10).stdout.strip()
                    pid = int(o) if o.isdigit() and int(o) > 0 else None
                    if pid is None:
                        continue
                peak[0] = max(peak[0], tree_rss_gb(pid))
            except Exception:
                pass
    th = threading.Thread(target=poll, daemon=True); th.start()
    t0 = time.time()
    with open(os.path.join(workdir, f"{name}.log"), "w") as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    stop.set(); th.join(timeout=15)
    wall = time.time() - t0
    Step.log.append({"step": name, "wall_s": round(wall, 1), "peak_rss_gb": round(peak[0], 2)})
    log(f"step {name}: {wall:.0f} s, peak {peak[0]:.1f} GB (rc={r.returncode})")
    if r.returncode != 0:
        die(f"step {name} failed rc={r.returncode}; see {workdir}/{name}.log", 3)


def load_cache(a):
    root = a.cache
    if os.path.isdir(os.path.join(root, "target")) is False and os.path.isdir(root):
        cands = sorted(glob.glob(os.path.join(root, "fp-*", "manifest.json")))
        mans = [json.load(open(c)) for c in cands]
        want = sorted(x for x in a.pgs.replace(",", " ").split() if x)
        hit = [c for c, mm in zip(cands, mans) if mm.get("complete") and mm.get("pgs") == want]
        if len(hit) != 1:
            die(f"no unique complete fastpath cache under {root} for this PGS list ({len(hit)} found); "
                f"run: bin/fastpath.py build ... (no fallback)")
        root = os.path.dirname(hit[0])
    mpath = os.path.join(root, "manifest.json")
    if not os.path.exists(mpath):
        die(f"fastpath cache {root} has no manifest.json")
    m = json.load(open(mpath))
    if not m.get("complete"):
        die(f"fastpath cache {root} is incomplete (build did not finish)")
    if os.path.realpath(a.ref) != m["reference_fasta"]:
        die(f"cache key mismatch: cache built with reference {m['reference_fasta']}, this run passes {os.path.realpath(a.ref)}")
    want = sorted(x for x in a.pgs.replace(",", " ").split() if x)
    if m["pgs"] != want:
        die(f"cache key mismatch: cache was built for {len(m['pgs'])} scores, this run asks for {len(want)} "
            f"(differ: {sorted(set(m['pgs']) ^ set(want))[:5]}...)")
    sfs = score_paths(a.scorefile_cache, want)
    cur = {p: sha256_file(f) for p, f in zip(want, sfs)}
    if cur != m["scorefile_sha256"]:
        bad = [p for p in want if cur[p] != m["scorefile_sha256"].get(p)]
        die(f"cache key mismatch: scoring file content changed for {bad[:5]}")
    if sha256_file(os.path.join(root, "target", "targets.tsv")) != m["target_sites_sha256"]:
        die("cache key mismatch: target site set file altered")
    pid = panel_identity(a.panel)
    if (pid["size"], pid["head_tail_sha256"]) != (m["panel"]["size"], m["panel"]["head_tail_sha256"]):
        die("cache key mismatch: reference panel differs from the one the cache was built with")
    tools = tool_versions()
    for k in ("utils_image", "plink2_image", "fraposa_image", "pgsc_calc_revision", "fastpath"):
        if tools[k] != m["tools"][k]:
            die(f"cache key mismatch: tool version {k}: built with {m['tools'][k]}, now {tools[k]}")
    if m["key"] != compute_key(want, cur, m["target_sites_sha256"], m["panel"], m["tools"]):
        die("cache key mismatch: manifest key does not recompute")
    for f in ("match/candidates.ipc.zst", "intersect/matched_variants.txt.gz", "intersect/filter_ids.txt.gz",
              "intersect/intersect_counts_ALL.txt", "formatted/log_scorefiles.json", "label/L_00.ipc", "label/summary_counts.ipc",
              "label/scoring_rows.ipc", "label/wide/idsplit.ipc", "label/wide/split_counts.ipc",
              "label/wide/passing.txt"):
        if not os.path.isfile(os.path.join(root, f)):
            die(f"fastpath cache {root} is missing {f}")
    return root, m


def _zstd_lines(path):
    pz = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True)
    for line in pz.stdout:
        yield line
    if pz.wait() != 0:
        die(f"zstd failed reading {path}")


def derive_intersect(canon_matched, canon_counts, canon_pvar, pvar, afreq, vmiss, out_matched, out_counts,
                     extra_dir, delta_runner=None, maf_filter=0.0, vmiss_filter=0.1, canon_dup_keys=0):
    """The sample's INTERSECT_VARIANTS output (matched_variants.txt.gz) from the cached canonical one.

    pgscatalog-intersect pairs reference and target variants by sorted allele key and takes AAF /
    F_MISS_DOSAGE from the target's own .afreq/.vmiss; the pairing of a target row is independent of
    the other target rows.  So the sample's table is
      (a) the canonical table restricted to target IDs present in the sample pvar, with AAF / F_MISS_DOSAGE
          / PCA_ELIGIBLE recomputed with the tool's own formula, plus
      (b) the real tool run on the sample-only rows (IDs absent from the canonical site set, e.g. indel
          representations bcftools emitted) against just the reference rows at their positions,
    merged in the tool's sort order.  Validated byte-for-byte against the real tool
    (bin/test_fastpath_unit.py, docs/fastpath.md).  Returns stats incl. the sample-only pvar path."""
    import csv, heapq
    canon_ids = set()
    for line in _zstd_lines(canon_pvar):
        if not line.startswith("#"):
            canon_ids.add(line.split("\t", 3)[2])
    tab, extra = {}, []
    n_t = 0
    pv_header = None
    with gzip.open(afreq, "rt") as fa, gzip.open(vmiss, "rt") as fm:
        fa_h, fm_h = next(fa), next(fm)
        for pl_ in _zstd_lines(pvar):
            if pl_.startswith("#"):
                if not pl_.startswith("##"):
                    pv_header = pl_
                continue
            n_t += 1
            v = pl_.rstrip("\n").split("\t")
            fa_l, fm_l = next(fa), next(fm)
            fa_ = fa_l.rstrip("\n").split("\t")
            fm_ = fm_l.rstrip("\n").split("\t")
            if not (v[2] == fa_[1] == fm_[0]):
                die("target pvar/afreq/vmiss rows are not aligned")
            if len(v[4].split(",")) != 1:
                die("sample pvar contains a multiallelic row; fastpath expects --max-alleles 2 input")
            if v[2] in canon_ids:
                tab[v[2]] = (str(float(fa_[4].split(",")[0])), fm_[1])
            elif v[2] != ".":
                extra.append((pl_, fa_l, fm_l, v))
    extra_pvar = None
    os.makedirs(extra_dir, exist_ok=True)
    with open(os.path.join(extra_dir, "absent_ids.txt"), "w") as fh:
        fh.write("".join(i + "\n" for i in sorted(canon_ids - tab.keys())))
    delta_rows, delta_keys = [], set()
    if extra:
        os.makedirs(extra_dir, exist_ok=True)
        pre = os.path.join(extra_dir, f"{TARGET_BUILD}_extra_ALL")
        with open(pre + ".pvar", "w") as fp, gzip.open(pre + ".afreq.gz", "wt") as fa, gzip.open(pre + ".vmiss.gz", "wt") as fm:
            fp.write(pv_header); fa.write(fa_h); fm.write(fm_h)
            for pl_, fa_l, fm_l, v in extra:
                fp.write(pl_); fa.write(fa_l); fm.write(fm_l)
        extra_pvar = pre + ".pvar"
        if delta_runner is not None:
            delta_rows = delta_runner(pre, [(v[0], int(v[1])) for _, _, _, v in extra])
    delta_keys = {r[0] for r in delta_rows}
    n_match = 0
    tmp_out = out_matched + ".canon.tmp"
    with gzip.open(canon_matched, "rt", newline="") as fi, gzip.open(tmp_out, "wt", newline="") as fo:
        rd = csv.reader(fi, delimiter="\t")
        wr = csv.writer(fo, delimiter="\t")
        hdr = next(rd)
        ix = {c: i for i, c in enumerate(hdr)}
        wr.writerow(hdr)
        i_t, i_aaf, i_fm = ix["ID_TARGET"], ix["AAF"], ix["F_MISS_DOSAGE"]
        for r in rd:
            t = tab.get(r[i_t])
            if t is None:
                continue
            r[i_aaf], r[i_fm] = t
            r[ix["IS_MA_TARGET"]] = "False"
            r[ix["SAME_REF"]] = "True"
            aaf = float(r[i_aaf]); maf = 1 - aaf if aaf > 0.5 else aaf
            # pgscatalog-intersect compares a python bool to the string "True" for indels, so an indel is
            # never PCA-eligible (reproduced on purpose; verified byte-for-byte against the tool)
            elig = (r[ix["IS_MA_REF"]] == "False") and (
                r[ix["IS_INDEL"]] == "False" and r[ix["STRANDAMB"]] == "False")
            elig = elig and (maf >= maf_filter) and (float(r[i_fm]) <= vmiss_filter)
            r[ix["PCA_ELIGIBLE"]] = str(elig)
            if r[0] in delta_keys:
                die(f"allele-key collision between a sample-only variant and a cached one ({r[0]}); "
                    f"the positional pairing of pgscatalog-intersect cannot be reproduced -- refusing")
            wr.writerow(r)
            n_match += 1
    if canon_dup_keys:
        die("cached canonical intersect has duplicate allele keys; fastpath cannot reproduce the pairing")
    # merge canonical-derived rows and delta rows in the tool's order (key, ID_REF, REF_REF)
    sk = lambda r: (r[0], r[1], r[2])
    with gzip.open(tmp_out, "rt", newline="") as fi, gzip.open(out_matched, "wt", newline="") as fo:
        rd = csv.reader(fi, delimiter="\t")
        wr = csv.writer(fo, delimiter="\t")
        hdr = next(rd)
        wr.writerow(hdr)
        for r in heapq.merge(rd, sorted(delta_rows, key=sk), key=sk):
            wr.writerow(r)
    os.remove(tmp_out)
    n_match += len(delta_rows)
    n_ref = open(canon_counts).read().split("\n")[1]
    open(out_counts, "w").write("\n".join(map(str, [n_t, n_ref, n_match])))
    return {"target_variants": n_t, "matched": n_match, "sample_only_rows": len(extra), "absent_ids": os.path.join(extra_dir, "absent_ids.txt"),
            "sample_only_matched": len(delta_rows), "extra_pvar": extra_pvar}


def relabel_awk(matched, scorefile, out):
    """== pgscatalog-relabel --col_from ID_TARGET --col_to ID_REF --target_col ID --combined (ID is column 1)."""
    prog = ("NR==FNR{ if (FNR>1) { sub(/\\r$/,\"\"); m[$7]=$2 } ; next } "
            "FNR==1{ printf \"%s\\r\\n\", $0; next } "
            "{ i=index($0,\"\\t\"); id=substr($0,1,i-1); if (!(id in m)) { print \"unmapped ID \" id > \"/dev/stderr\"; exit 3 } "
            "printf \"%s%s\\r\\n\", m[id], substr($0,i) }")
    subprocess.run(["/usr/bin/time", "-f", "%M", "-o", out + ".rss", "bash", "-c",
                    "set -euo pipefail; awk -F'\\t' \"$1\" <(zcat \"$2\") <(zcat \"$3\") | gzip -c > \"$4\"",
                    "_", prog, matched, scorefile, out], check=True)
    rss = int(open(out + ".rss").read().split()[-1]) * 1024 / 1e9
    os.remove(out + ".rss")
    return rss


def run_ancestry_chain(d, matched, tpf, pre, S, ref, refpre, king, ld, thr, memmb):
    """ANCESTRY_PROJECT of pgsc_calc (modules/local/ancestry/*.nf) for one sample: reference QC, LD-thin,
    intersect, relabel, bed export, orient, FRAPOSA PCA + OADP projection.  Independent of the scoring files."""
    os.makedirs(d, exist_ok=True)
    shutil.copy(ld, os.path.join(d, os.path.basename(ld)))      # /home is not mounted in the containers
    ld = os.path.join(d, os.path.basename(ld))
    f1 = os.path.join(d, "filter"); os.makedirs(f1)
    run_step("filter_variants", IMG_PLINK2,
             f"awk '($13 ~ /^True/) {{print $2}}' <(zcat {matched}) | gzip -c > shared.txt.gz\n"
             f"plink2 --threads {thr} --memory {memmb} --pfile {ref}/{refpre} vzs --remove {king} "
             f"--extract shared.txt.gz --max-alleles 2 --snps-only just-acgt --rm-dup exclude-all "
             f"--geno 0.1 --mind 0.1 --maf 0.05 --hwe 0.0001 --make-pgen vzs --allow-extra-chr "
             f"--autosome --freq zs --out {TARGET_BUILD}_reference\n"
             f"plink2 --threads {thr} --memory {memmb} --pfile vzs {TARGET_BUILD}_reference "
             f"--indep-pairwise 1000 50 0.05 --exclude range {ld} --out {refpre}_thinned\n"
             f"gzip -f *.prune.in *.prune.out", f1)
    f2 = os.path.join(d, "makebed_ref"); os.makedirs(f2)
    run_step("makebed_ref", IMG_PLINK2,
             f"plink2 --threads {thr} --memory {memmb} --seed 31 --pgen {f1}/{TARGET_BUILD}_reference.pgen "
             f"--psam {f1}/{TARGET_BUILD}_reference.psam --pvar {f1}/{TARGET_BUILD}_reference.pvar.zst "
             f"--make-bed --extract {f1}/{refpre}_thinned.prune.in.gz --out {TARGET_BUILD}_reference_extracted", f2)
    f3 = os.path.join(d, "intersect_thinned"); os.makedirs(f3)
    for e in (".pgen", ".psam", ".pvar.zst"):
        os.symlink(tpf + e, os.path.join(f3, pre + e))
    os.symlink(matched, os.path.join(f3, f"{S}_ALL_matched.txt.gz"))
    os.symlink(f"{f1}/{refpre}_thinned.prune.in.gz", os.path.join(f3, f"{refpre}_thinned.prune.in.gz"))
    # INTERSECT_THINNED (oadp/intersect_thinned.nf), single combined pfile
    run_step("intersect_thinned", IMG_PLINK2,
             f"awk 'BEGIN {{ i=0 }} NR!=FNR && FNR == 1 && i==0 {{ print; i++ }} NR==FNR {{ F1[$0]; next }} "
             f"$2 in F1 {{ print }}' <(gunzip -c {refpre}_thinned.prune.in.gz) <(gunzip -c {S}_ALL_matched.txt.gz) "
             f"> {S}_ALL_matched_thinned.txt\n"
             f"cut -f 7 -d ' ' {S}_ALL_matched_thinned.txt > {S}_shared.txt\n"
             f"echo {pre}.pgen {pre}.psam {pre}.pvar.zst | tr -s ' ' '\\n' | cut -f 1 -d '.' | uniq > ids.txt\n"
             f"cp ids.txt autosomes.txt\nmkdir -p extracted\n"
             f"while read f; do plink2 --threads {thr} --memory {memmb} --seed 31 --pfile $f vzs "
             f"--extract {S}_shared.txt --chr 1-22 --allow-extra-chr --make-pgen --sort-vars "
             f"--out extracted/${{f}}_extracted; done < autosomes.txt\n"
             f"mv extracted/*.p* .\nrm -r extracted autosomes.txt ids.txt\ngzip *.txt *.pvar", f3)
    f4 = os.path.join(d, "relabel_ids"); os.makedirs(f4)
    run_step("relabel_ids", IMG_UTILS,
             f"pgscatalog-relabel --maps {f3}/{S}_ALL_matched_thinned.txt.gz --col_from ID_TARGET --col_to ID_REF "
             f"--target_file {f3}/{pre}_extracted.pvar.gz --target_col ID --dataset {S}.pvar --verbose "
             f"--split --combined --outdir $PWD", f4)
    f5 = os.path.join(d, "makebed_target"); os.makedirs(f5)
    run_step("makebed_target", IMG_PLINK2,
             f"plink2 --threads {thr} --memory {memmb} --seed 31 --pgen {f3}/{pre}_extracted.pgen "
             f"--psam {f3}/{pre}_extracted.psam --pvar {f4}/{S}.pvar_ALL_relabelled.gz --make-bed "
             f"--out {TARGET_BUILD}_{S}\nsplit -l 50000 {TARGET_BUILD}_{S}.fam {S}_splitfam", f5)
    f6 = os.path.join(d, "orient"); os.makedirs(f6)
    run_step("orient", IMG_PLINK2,
             f"plink2 --threads {thr} --memory {memmb} --seed 31 --bed {f5}/{TARGET_BUILD}_{S}.bed "
             f"--fam {f5}/{TARGET_BUILD}_{S}.fam --bim {f5}/{TARGET_BUILD}_{S}.bim "
             f"--alt1-allele {f2}/{TARGET_BUILD}_reference_extracted.bim 5 2 --make-bed "
             f"--out {TARGET_BUILD}_{S}_oriented", f6)
    f7 = os.path.join(d, "fraposa_pca"); os.makedirs(f7)
    for e in (".bed", ".bim", ".fam"):
        os.symlink(f"{f2}/{TARGET_BUILD}_reference_extracted{e}", os.path.join(f7, f"{TARGET_BUILD}_reference_extracted{e}"))
    run_step("fraposa_pca", IMG_FRAPOSA,
             f"fraposa {TARGET_BUILD}_reference_extracted --method oadp --dim_ref 10", f7)
    f8 = os.path.join(d, "fraposa_project"); os.makedirs(f8)
    for src in glob.glob(f"{f7}/*"):
        os.symlink(src, os.path.join(f8, os.path.basename(src)))
    for e in (".bed", ".bim", ".fam"):
        os.symlink(f"{f6}/{TARGET_BUILD}_{S}_oriented{e}", os.path.join(f8, f"{TARGET_BUILD}_{S}_oriented{e}"))
    os.symlink(f"{f5}/{S}_splitfamaa", os.path.join(f8, f"{S}_splitfamaa"))
    run_step("fraposa_project", IMG_FRAPOSA,
             f"fraposa {TARGET_BUILD}_reference_extracted --method oadp --dim_ref 10 "
             f"--stu_filepref {TARGET_BUILD}_{S}_oriented --stu_filt_iid {S}_splitfamaa "
             f"--out {TARGET_BUILD}_{S}_oriented_{S}_splitfamaa", f8)
    return {"filter": f1, "ref_pca": f7, "project": f8}


def cmd_score(a):
    t_all = time.time()
    cache, man = load_cache(a)
    pgs = man["pgs"]
    S = a.sample
    out = os.path.abspath(a.outdir)
    prep = os.path.abspath(a.prepped)
    vcf = prep if prep.endswith(".vcf.gz") else prep + ".vcf.gz"
    if not os.path.isfile(vcf):
        die(f"prepped VCF not found: {vcf}")
    if os.environ.get("FASTPATH_WORK_ROOT"):
        W = os.path.join(os.environ["FASTPATH_WORK_ROOT"], f"fastpath_work_{S}")
    else:
        W = os.path.join(out, "fastpath_work")
    shutil.rmtree(W, ignore_errors=True); os.makedirs(W)
    thr, mem = str(a.threads), a.mem_gb
    memmb = min(int(a.mem_gb * 1000), 32000)   # plink2 workspace only; never changes results
    ref = os.path.join(cache, "ref")
    refpre = f"{TARGET_BUILD}_HGDP+1kGP_ALL"
    for ext in (".pgen", ".psam", ".pvar.zst"):
        if not os.path.isfile(os.path.join(ref, refpre + ext)):
            die(f"cache reference panel file missing: {refpre}{ext}")
    king = os.path.join(ref, f"{TARGET_BUILD}_HGDP+1kGP.king.cutoff.out.id")
    ld = os.path.join(PGSC_ASSETS, "assets", "ancestry", "high-LD-regions-hg38-GRCh38.txt")
    formatted = sorted(glob.glob(os.path.join(cache, "formatted", "normalised_*")))
    D = lambda n: os.path.join(W, n)
    tm = {}

    def timed(label, fn):
        t0 = time.time(); r = fn(); tm[label] = round(time.time() - t0, 1); return r

    # 1. PLINK2_VCF
    vdir = D("vcf"); os.makedirs(vdir)
    pre = f"{TARGET_BUILD}_{S}_ALL"
    run_step("plink2_vcf", IMG_PLINK2, plink2_vcf_cmd(vcf, pre, thr, mem), vdir, mem_gb=None)
    tpf = os.path.join(vdir, pre)

    from concurrent.futures import ThreadPoolExecutor
    pool = ThreadPoolExecutor(max_workers=6)
    idir = D("intersect"); os.makedirs(idir)
    matched = os.path.join(idir, f"{S}_ALL_matched.txt.gz")
    counts = os.path.join(idir, "intersect_counts_ALL.txt")

    def delta_runner(pre, positions):
        """real pgscatalog-intersect on the sample-only rows vs the reference rows at those positions"""
        d = D("intersect_delta"); os.makedirs(d)
        with open(os.path.join(d, "regions.tsv"), "w") as f:
            for c, p_ in sorted(set(positions)):
                f.write(f"{c}\t{p_}\t{p_}\n")
        body = subprocess.run(["tabix", "-R", os.path.join(d, "regions.tsv"),
                               os.path.join(ref, f"{refpre}.pvar.gz")], capture_output=True, text=True, check=True).stdout
        with open(os.path.join(d, "mini_ref.pvar"), "w") as f:
            f.write(open(os.path.join(ref, "pvar.header")).read()); f.write(body)
        try:
            run_step("intersect_delta", IMG_UTILS,
                     f"pgscatalog-intersect --ref {d}/mini_ref.pvar --target {pre}.pvar --chrom ALL --maf_target 0 "
                     f"--geno_miss 0.1 --outdir {d} -v", d)
        except SystemExit:
            if "ZeroDivisionError" in open(os.path.join(d, "intersect_delta.log")).read():
                return []                   # the tool divides by zero when nothing at all matches
            raise
        import csv
        with gzip.open(os.path.join(d, "matched_variants.txt.gz"), "rt", newline="") as f:
            rd = csv.reader(f, delimiter="\t"); next(rd)
            return [r for r in rd]
    t0 = time.time()
    dstat = derive_intersect(os.path.join(cache, "intersect", "matched_variants.txt.gz"),
                             os.path.join(cache, "intersect", "intersect_counts_ALL.txt"),
                             os.path.join(cache, "target", f"{TARGET_BUILD}_CANON_ALL.pvar.zst"),
                             tpf + ".pvar.zst", tpf + ".afreq.gz", tpf + ".vmiss.gz", matched, counts,
                             D("extra"), delta_runner,
                             canon_dup_keys=man["stages"]["intersect_stats"]["canonical_duplicate_keys"])
    Step.log.append({"step": "derive_intersect", "wall_s": round(time.time() - t0, 1),
                     **{k: v for k, v in dstat.items() if k != "extra_pvar"}})
    log(f"derive_intersect: { {k: v for k, v in dstat.items() if k != 'extra_pvar'} } in {time.time()-t0:.0f} s")

    mdir = D("match")
    fm_env = {"FASTPATH_FULL_LOG": "1" if a.full_log else "0"}

    def do_match():
        # delta match/merge against the cached labelled canonical candidates (bin/fastpath_delta.py)
        pre_cmd = ""
        extra_arg = "-"
        if dstat["extra_pvar"]:
            ed = os.path.join(mdir, "extra"); extra_arg = f"{ed}/extra.ipc"
            pre_cmd = (f"mkdir -p {ed}/w; pgscatalog-match --dataset {S} --scorefile " + " ".join(formatted) +
                       f" --target {dstat['extra_pvar']} --only_match --outdir {ed}/w -v; "
                       f"mv {ed}/w/*.ipc.zst {extra_arg}\n")
        run_step("match", IMG_UTILS,
                 f"export POLARS_MAX_THREADS={thr}\nmkdir -p {mdir}\n{pre_cmd}python {os.path.join(HERE, 'fastpath_delta.py')} "
                 f"{os.path.join(cache, 'label')} {dstat['absent_ids']} {extra_arg} "
                 f"{os.path.join(cache, 'intersect', 'filter_ids.txt.gz')} {mdir} {S} {MIN_OVERLAP}",
                 W, mem_gb=a.mem_gb)
    fut_match = pool.submit(do_match)

    fut_anc = pool.submit(run_ancestry_chain, D("anc"), matched, tpf, pre, S, ref, refpre, king, ld, thr, memmb)

    fut_match.result()
    anc = fut_anc.result()
    sc = D("score"); os.makedirs(sc)
    mfiles = sorted(glob.glob(os.path.join(mdir, f"{S}_ALL_*_*.scorefile.gz")))
    if not mfiles:
        die("no matched scoring files were written (all scores below min_overlap?)")

    # 3. relabel reference scoring files + afreq, in parallel
    rl = D("relabel"); os.makedirs(rl)
    afr = D("relabel_afreq"); os.makedirs(afr)

    def relabel_one(f):
        """RELABEL_SCOREFILES as a streaming awk map (checked byte for byte against pgscatalog-relabel in
        bin/test_fastpath_unit.py): first column ID_TARGET -> ID_REF through the sample's matched table,
        rows written CRLF like csv.DictWriter; an unmapped ID is a hard error (as the tool's KeyError)."""
        b = os.path.basename(f)
        w = os.path.join(rl, b.replace(".scorefile.gz", "")); os.makedirs(w)
        outn = "reference" + b[b.index("_"):]
        t0_ = time.time()
        rss = relabel_awk(matched, f, os.path.join(w, outn))
        Step.log.append({"step": "relabel_" + b.replace(".scorefile.gz", ""), "wall_s": round(time.time() - t0_, 1),
                         "peak_rss_gb": round(rss, 2)})
        log(f"step relabel_{b}: {time.time()-t0_:.0f} s")
        return os.path.join(w, outn)
    futs = [pool.submit(relabel_one, f) for f in mfiles]
    fut_afreq = None
    # RELABEL_AFREQ needs FILTER_VARIANTS' afreq (ancestry chain finished above)
    run_step("relabel_afreq", IMG_UTILS,
             f"pgscatalog-relabel --maps {matched} --col_from ID_REF --col_to ID_TARGET "
             f"--target_file {anc['filter']}/{TARGET_BUILD}_reference.afreq.zst --target_col ID "
             f"--dataset {S}.afreq --verbose --combined --outdir $PWD", afr)
    afreq_rl = os.path.join(afr, f"{S}.afreq_ALL_relabelled.gz")
    ref_sf = [f.result() for f in futs]

    # 4. PLINK2_SCORE: sample and reference panel
    def score_one(kind, sfile):
        b = os.path.basename(sfile)
        stem = b.replace(".scorefile.gz", "")          # <id>_ALL_additive_<n>
        w = os.path.join(sc, stem if kind == "sample" else stem); os.makedirs(w, exist_ok=True)
        n_scores = len(gzip.open(sfile, "rt").readline().rstrip("\n").split("\t")) - 2
        effect = stem.split("_")[2]
        cols = "header-read cols=+scoresums,+denom,+fid"
        rec = " recessive " if effect == "recessive" else ""
        dom = " dominant " if effect == "dominant" else ""
        args2 = " ".join(["zs", cols, "list-variants", "", rec, dom, ""])
        pf = tpf if kind == "sample" else f"{ref}/{refpre}"
        pvar = ".pvar.zst"
        if kind == "sample":
            gpath = f"{tpf}"
        else:
            gpath = f"{ref}/{refpre}"
        extra = f"--score-col-nums 3-{n_scores + 2}" if n_scores > 1 else ""
        run_step(f"score_{kind}_{stem}", IMG_PLINK2,
                 f"plink2 --threads {a.score_threads} --memory {memmb} --seed 31 --extract {sfile} "
                 f"--read-freq {afreq_rl} --allow-extra-chr --score {sfile} {args2} {extra} "
                 f"--pfile vzs {gpath} --out {stem}", w)
        return stem, w
    sfut = [pool.submit(score_one, "sample", f) for f in mfiles] + [pool.submit(score_one, "reference", f) for f in ref_sf]
    scored = [f.result() for f in sfut]

    # 5. SCORE_AGGREGATE
    ag = D("aggregate"); os.makedirs(ag)
    for (stem, w), kind in zip(scored, ["sample"] * len(mfiles) + ["reference"] * len(ref_sf)):
        outstem = stem if kind == "sample" else "reference" + stem[stem.index("_"):]
        for e in (".sscore.zst", ".sscore.vars"):
            os.symlink(os.path.join(w, stem + e), os.path.join(ag, outstem + e))
    for f in mfiles:
        os.symlink(f, os.path.join(ag, os.path.basename(f)))
    for f in ref_sf:
        os.symlink(f, os.path.join(ag, os.path.basename(f)))
    run_step("score_aggregate", IMG_UTILS,
             "pgscatalog-aggregate -s " + " ".join(sorted(os.path.basename(x) for x in glob.glob(ag + "/*.sscore.zst")))
             + " -o . -v --no-split --verify_variants", ag)

    # 6. ANCESTRY_ANALYSIS
    an = D("ancestry_analysis"); os.makedirs(os.path.join(an, "target_pcs")); os.makedirs(os.path.join(an, "ref_pcs"))
    tp = glob.glob(os.path.join(anc["project"], "*.pcs"))
    tp = [x for x in tp if os.path.basename(x).startswith(f"{TARGET_BUILD}_{S}_oriented_{S}_splitfamaa")]
    rp = glob.glob(os.path.join(anc["ref_pca"], "*.pcs"))
    if len(tp) != 1 or len(rp) != 1:
        die(f"expected one target and one reference .pcs, found {tp} / {rp}")
    os.symlink(tp[0], os.path.join(an, "target_pcs", "001.pcs"))
    os.symlink(rp[0], os.path.join(an, "ref_pcs", "1.pcs"))
    run_step("ancestry_analysis", IMG_UTILS,
             f"pgscatalog-ancestry-adjust -d {S} -r reference --psam {ref}/{refpre}.psam --ref_pcs ref_pcs/1.pcs "
             f"--target_pcs target_pcs/*.pcs -x {king} -p SuperPop -s {ag}/aggregated_scores.txt.gz -a RandomForest "
             f"--n_popcomp 5 -n empirical mean mean+var --n_normalization 4 --outdir . -v", an)

    # 7. SCORE_REPORT (produces report.html and pop_summary.csv, which grade_pgs reads)
    rp_dir = D("report"); os.makedirs(rp_dir)
    for f in glob.glob(os.path.join(an, f"{S}_*")):
        shutil.copy(f, rp_dir)
    shutil.copy(os.path.join(mdir, f"{S}_summary.csv"), rp_dir)
    shutil.copy(os.path.join(cache, "formatted", "log_scorefiles.json"), rp_dir)
    shutil.copy(counts, rp_dir)
    for f in ("report.qmd", "logo.css", "PGS_Logo.png", "pgs_header_background.png"):
        shutil.copy(os.path.join(PGSC_ASSETS, "assets", "report", f), rp_dir)
    run_step("score_report", "ghcr.io/pgscatalog/report:2-beta",
             "printf 'project:\\n  type: default\\n' > _quarto.yml\nexport TMPDIR=$PWD\n"
             "echo fastpath > command.txt\n"
             "printf 'keep_multiallelic: false\\nkeep_ambiguous   : false\\nmin_overlap      : 0.1\\n' > params.txt\n"
             f"quarto render report.qmd -M 'self-contained:true' -P score_path:{S}_pgs.txt.gz -P sampleset:{S} "
             f"-P run_ancestry:true -P reference_panel_name:HGDP+1kGP -P version:2.3.0 -o report.html", rp_dir)

    # 8. publish in pgsc_calc's layout
    so = os.path.join(out, "score"); shutil.rmtree(so, ignore_errors=True)
    os.makedirs(os.path.join(so, S, "match")); os.makedirs(os.path.join(so, S, "score")); os.makedirs(os.path.join(so, "pipeline_info"))
    for f in glob.glob(os.path.join(mdir, f"{S}_*")):
        shutil.copy(f, os.path.join(so, S, "match"))
    for f in glob.glob(os.path.join(rp_dir, "*")):
        b = os.path.basename(f)
        if b.startswith(f"{S}_") and b.endswith((".txt.gz", ".json.gz")) or b in ("pop_summary.csv", "report.html"):
            shutil.copy(f, os.path.join(so, S, "score"))
    total = time.time() - t_all
    json.dump({"fastpath_version": FASTPATH_VERSION, "cache": cache, "cache_key": man["key"], "sample": S,
               "total_wall_s": round(total, 1), "steps": Step.log},
              open(os.path.join(so, "pipeline_info", "fastpath_timings.json"), "w"), indent=1)
    if not a.keep_work:
        shutil.rmtree(W, ignore_errors=True)
    log(f"score done in {total:.0f} s -> {so}")


def plink2_vcf_cmd(vcf, out, threads, mem_gb):
    """PLINK2_VCF (modules/local/plink2_vcf.nf; ext.args '--new-id-max-allele-len 100 missing')."""
    mem_mb = min(int(mem_gb * 1000), 32000) if mem_gb else 16000
    return (f"plink2 --threads {threads} --memory {mem_mb} --set-all-var-ids '@:#:$r:$a' --max-alleles 2 "
            f"--freq --missing vcols=fmissdosage,fmiss --new-id-max-allele-len 100 missing "
            f"--vcf {vcf} --allow-extra-chr --chr 1-22, X, Y, XY "
            f"--make-pgen vzs pvar-cols=\"-xheader,-maybequal,-maybefilter,-maybeinfo,-maybecm\" --out {out}; "
            f"gzip -f {out}.vmiss; gzip -f {out}.afreq")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    b = sp.add_parser("build")
    b.add_argument("--pgs", required=True)
    b.add_argument("--scorefile-cache", required=True)
    b.add_argument("--panel", required=True)
    b.add_argument("--ref", required=True, help="reference FASTA (part of the cache id)")
    b.add_argument("--cache-root", default=os.path.join(ROOT, "cache", "fastpath"))
    b.add_argument("--threads", type=int, default=8)
    b.add_argument("--mem-gb", type=int, default=100)
    c = sp.add_parser("score")
    c.add_argument("--sample", required=True)
    c.add_argument("--prepped", required=True, help="prepped VCF (or its prefix) of the sample")
    c.add_argument("--outdir", required=True)
    c.add_argument("--pgs", required=True)
    c.add_argument("--scorefile-cache", required=True)
    c.add_argument("--panel", required=True)
    c.add_argument("--ref", required=True, help="reference FASTA the cache was built with (checked)")
    c.add_argument("--cache", default=os.path.join(ROOT, "cache", "fastpath"))
    c.add_argument("--threads", type=int, default=8)
    c.add_argument("--score-threads", type=int, default=8)
    c.add_argument("--mem-gb", type=int, default=100)
    c.add_argument("--full-log", action="store_true", help="also write the multi-GB per-variant <S>_log.csv.gz")
    c.add_argument("--keep-work", action="store_true")
    a = ap.parse_args()
    if a.cmd == "build":
        cmd_build(a)
    else:
        cmd_score(a)


if __name__ == "__main__":
    main()
