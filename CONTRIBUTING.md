# Contributing to poly-suite

Thank you for helping. poly-suite is a research and educational tool, so changes are judged first on whether
they keep scores correct and reproducible.

## Before you open a pull request

1. Run the unit tests: `tests/selftest.sh` (no network or Docker needed). It must end with `ALL SELFTESTS PASSED`.
2. If your change touches scoring, input handling or `bin/run.sh`, also run the smoke test:
   `tests/smoke/run.sh` (needs Docker and Nextflow; under a minute once the containers are pulled).
3. If a change alters any score, say so in the pull request and show the before and after values. Changes that
   only make things faster must leave every output identical. `bin/compare_fastpath.py` compares two score
   directories.
4. Add a dated entry to `docs/CHANGES.md` describing what changed and how it was checked.

## Data rules

- Use only public samples (GIAB HG001 to HG007) in tests, examples, screenshots and issue reports.
- Do not commit sample data, BAMs, VCFs or results. `data/`, `results/`, `work/`, `cache/` and `scratch/` are
  gitignored for this reason. The only committed genotype data is the small GIAB slice in `tests/smoke/`.
- Do not include hostnames, usernames or internal file paths in code, docs or logs you paste.

## Style

- Keep the output contract (`pgs_scores.tsv`, checked by `tests/validate_contract.py`) stable. A new column
  needs a matching change in OmniGen, so discuss it in an issue first.
- Plain, non-clinical wording in reports and docs.
- Diagrams in `docs/img/` are Mermaid sources (`*.mmd`) rendered with `docs/img/render.sh`.
