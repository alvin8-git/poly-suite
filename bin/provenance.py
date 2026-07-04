#!/usr/bin/env python3
"""Assemble provenance.json for a poly-suite run — the reproducibility record
(spec §2). Captures pipeline version, timestamp, reference, scores, ancestry
panel, tool versions (from pgsc_calc's versions.yml), and a content hash of the
graded contract so a report is reproducible to the release.

Importable (build/write) and runnable: python3 bin/provenance.py [results_dir]
"""
import os, sys, json, glob, csv, hashlib, subprocess
from datetime import datetime, timezone

POLY_SUITE_VERSION = "0.1"


def _sha256(path):
    if not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _tool_versions(results_dir):
    hits = sorted(glob.glob(f"{results_dir}/**/versions.yml", recursive=True))
    if not hits:
        return None
    txt = open(hits[0]).read()
    try:
        import yaml
        return yaml.safe_load(txt)
    except Exception:
        return {"_raw": txt}          # yaml not installed -> keep raw text


def _first_line(args):
    try:
        out = subprocess.run(args, capture_output=True, text=True, timeout=15).stdout
        for ln in out.splitlines():
            if ln.strip():
                return ln.strip()
    except Exception:
        return None
    return None


def _contract_scores(contract):
    scores, calibrated = [], False
    if os.path.exists(contract):
        with open(contract) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                scores.append(r["pgs_id"])
                if r.get("percentile") not in (None, "", "NA"):
                    calibrated = True
    return sorted(set(scores)), calibrated


def build(results_dir, ref_fasta=None, panel=None, sample=None, bam=None, params=None):
    contract = os.path.join(results_dir, "pgs_scores.tsv")
    scores, calibrated = _contract_scores(contract)
    return {
        "poly_suite_version": POLY_SUITE_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "target_build": "GRCh38",
        "sample": sample,
        "input_bam": bam,
        "reference_fasta": ref_fasta,
        "ancestry_panel": (panel if calibrated
                           else (panel or "none (uncalibrated — no --run_ancestry)")),
        "calibrated": calibrated,
        "pgs_scores": scores,
        "pgs_catalog_harmonization": "hmPOS GRCh38",
        "genotype_prep": "force-genotype at scoring loci (bcftools mpileup | call -A -C alleles)",
        "tool_versions": _tool_versions(results_dir),
        "host_tools": {
            "nextflow": _first_line(["/home/alvin/bin/nextflow", "-v"]),
            "samtools": _first_line(["samtools", "--version"]),
            "python": _first_line(["python3", "--version"]),
        },
        "contract_sha256": _sha256(contract),
        "parameters": params or {"min_overlap": 0.1, "target_build": "GRCh38"},
    }


def write(results_dir, **kw):
    prov = build(results_dir, **kw)
    out = os.path.join(results_dir, "provenance.json")
    with open(out, "w") as fh:
        json.dump(prov, fh, indent=2)
    return out, prov


if __name__ == "__main__":
    rd = sys.argv[1] if len(sys.argv) > 1 else "results"
    out, prov = write(rd)
    print(f"provenance -> {out}  (calibrated={prov['calibrated']}, "
          f"{len(prov['pgs_scores'])} scores, sha={str(prov['contract_sha256'])[:12]})")
