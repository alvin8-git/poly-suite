---
name: Bug report
about: A run failed or produced a result you believe is wrong
labels: bug
---

**What happened**
<!-- The command you ran (replace local paths with placeholders) and what went wrong. -->

**What you expected**

**Log excerpt**
<!-- The last lines of <outdir>/score.log or run.log. Remove hostnames, usernames and internal paths. -->

**Environment**
- poly-suite commit (`git rev-parse --short HEAD`):
- pgsc_calc revision (`nextflow info pgscatalog/pgsc_calc`):
- Nextflow version:
- Docker version, rootless or not:
- Input route: `--bam` / `--vcf` / `--reuse-prep`; `--panel` yes or no; `--fastpath` yes or no

**Does `tests/selftest.sh` pass? Does `tests/smoke/run.sh` pass?**

<!-- Please use only public samples (GIAB HG001-HG007) in examples. Do not attach private data. -->
