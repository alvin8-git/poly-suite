# Candidate traits — what else can be added, and the tradeoffs

A data-grounded survey of traits **not** in the current 61-trait launch set
([Documentation §6](Documentation.md#6-the-tiered-launch-set)), with the real
PGS Catalog picture and the tradeoff for each group. Grades below are what
poly-suite's `select_pgs.grade()` would assign the best available score (GWAS
size + ancestry-eval thresholds), computed live against the Catalog.

## The headline: availability is not the binding constraint

The instinct is "which traits have good enough PGS?" The data says: **most do.**
Biomarkers, most common cancers, several common diseases, and behavioral traits
all reach grade **A** on the best available score. So the decision is not gated
by evidence — it's gated by five other things:

1. **Absolute-risk convertibility.** poly-suite emits a natural-frequency "1 in N"
   only for traits with a verified per-SD effect **and** a baseline incidence in
   `resources/`. Today that's **9 traits** (CAD, breast, prostate, T2D, RA, IBD,
   AF, colorectal, Alzheimer). Every new *disease* added without that pair renders
   percentile-only — informative, but not the flagship figure.
2. **Actionability.** A raised genetic likelihood is worth showing only if a reader
   can do something (screening, lifestyle, "mention to your doctor"). A PGS for
   varicose veins or hearing loss is real but inert.
3. **Consent sensitivity.** Some high-grade traits (psychiatric, behavioral,
   substance use, educational attainment) carry disclosure and misuse risk that a
   general report should not surface without explicit opt-in — the same reason the
   current set gates 6 "sensitive" traits.
4. **Redundancy / compute.** Prep force-genotypes the *union* of scoring loci. Each
   added trait grows that union (prep + scoring time). Correlated traits (LDL is in
   the set; adding apoB, non-HDL, total cholesterol adds loci for little new signal).
5. **Interpretation burden.** 60 traits already. Every addition is another row a
   lay reader must triage; the marginal trait has to earn its scroll.

## The tranches

### A. Quantitative biomarkers — grade A, cheap, but a different report mode
`platelet count` A (N 634k), `red blood cell count` A (728k), `alanine
aminotransferase` A (437k), `glomerular filtration rate` A (**1.16M**, 91 scores),
`urate` A (458k), `systolic blood pressure` A (**1.36M**), `heart rate` A (239k),
`FEV1` A (321k), `bone mineral density` A (395k).

- **Why add:** enormous GWAS → trivially grade A, and several are genuinely
  actionable proxies — BMD (fracture/osteoporosis risk, distinct from the osteoporosis
  *diagnosis* score already in the set), SBP (hypertension/CAD), eGFR/urate (CKD/gout),
  ALT (fatty-liver). They also anchor the report's credibility ("here are things with
  million-sample evidence").
- **Tradeoffs:** these are **continuous measures, not diseases** — "absolute risk /
  1 in N" doesn't apply; they'd report as "genetic predisposition to higher/lower X."
  That needs a second report mode (a distribution + direction, no dot-array), or they
  clutter the disease-framed card. Many are also **downstream of traits already shown**
  (SBP↔hypertension, urate↔gout, BMD↔osteoporosis) — risk of double-counting the same
  biology. **Verdict: add BMD, eGFR, SBP as a small "biomarker" section with its own
  framing; skip the blood counts (grade A but near-zero lay actionability).**

### B. Common diseases, well-powered — the strongest disease additions
`peripheral arterial disease` A, `endometriosis` A (209k), `diverticular disease`
A (455k), `obstructive sleep apnea` A (359k), `insomnia` A (**1.33M**),
`age-related hearing impairment` A (254k).

- **Why add:** grade A, common, and several are actionable (PAD → cardiovascular
  workup; endometriosis → shortens a notoriously long diagnostic odyssey; OSA →
  treatable). Fills real gaps (no vascular-peripheral, no gynaecological-benign,
  no sleep trait today).
- **Tradeoffs:** no effect+baseline in `resources/` yet → percentile-only until
  curated (per-trait literature work, same pattern as the RA/IBD additions). Some are
  **sex-specific** (endometriosis) and need the sex plumbing (already present).
  Hearing loss and insomnia are high-grade but low-actionability. **Verdict: add PAD,
  endometriosis, OSA; curate their absolute-risk inputs so they get a "1 in N."**

### C. Common diseases, thin evidence — the grade gate does its job
`epilepsy` C (45k), `uterine fibroid` D, `kidney stone` D, `gastroesophageal
reflux disease` D, `varicose veins` D.

- **Why add:** completeness; common conditions people ask about.
- **Tradeoffs:** they land **C/D** and render with a loud caveat or number-less — i.e.
  they add rows without adding trustworthy findings. This is the system working as
  designed, but each is a row that mostly says "not enough evidence." **Verdict: hold.
  Re-score automatically when the Catalog gains a larger GWAS (the grade will lift itself).**

### D. Additional cancers — mostly grade A, high perceived value
`esophageal cancer` A (394k), `stomach cancer` A (208k), `cervical cancer` A (327k),
`glioma` A (104k), `multiple myeloma` A (**1.05M**); `renal cell carcinoma` C (42k).

- **Why add:** the set already leans into cancer (12 sites); these are the obvious
  missing common/again-screenable ones. Cancer PGS have SEER lifetime baselines readily
  available, so **absolute-risk curation is easiest here** (SEER gives the baseline; the
  score's own paper gives the per-SD OR).
- **Tradeoffs:** cancer findings carry the highest anxiety-per-row and the strongest
  "false precision" risk — a grade-A PGS still has modest AUROC (the report's new
  reported-perf column matters most here). Cervical cancer is HPV-driven; a PGS is a
  weak sideshow to the real (viral) cause — misleading without that context.
  **Verdict: add esophageal, gastric, cervical (with an HPV caveat); curate SEER
  baselines so they get a "1 in N." Hold renal cell (C).**

### E. Rarer autoimmune — thin, and the set is already autoimmune-heavy
`vitiligo` C (40k), `systemic sclerosis` C (27k), `Sjogren syndrome` D.

- **Tradeoffs:** small GWAS → C/D, and the launch set already has the well-powered
  autoimmune cluster (T1D, RA, IBD, coeliac, psoriasis, MS, SLE, ankylosing spondylitis,
  thyroid). These share HLA, so they're **correlated** — more autoimmune rows add little
  independent signal and dilute the strong ones. **Verdict: skip until a large GWAS lands.**

### F. Psychiatric / neuro additions — high grade, high sensitivity
`anxiety` A (333k), `post-traumatic stress disorder` A (196k). (OCD, anorexia, ALS
exist under other ontology labels; ALS and anorexia GWAS are mid-sized.)

- **Tradeoffs:** grade A is achievable, but psychiatric PGS are the **most consent- and
  stigma-sensitive** class, have poor positive predictive value at the individual level,
  and interact with the existing sensitive-tier gating (schizophrenia, MDD, bipolar,
  ADHD, autism are already opt-in). **Verdict: only inside the explicitly-gated sensitive
  tier, never in the default suite.**

### G. Behavioral / sensitive — grade A, but out of scope for a health report
`educational attainment` A (334k), `smoking initiation` A (**1.23M**, 44 scores),
`alcohol dependence` A (492k), `age at menopause` B (69k).

- **Tradeoffs:** technically pristine (huge GWAS), but educational-attainment / risk /
  substance PGS are the textbook **misuse-risk** traits — not health findings, high
  potential for discrimination, and a reputational liability for a clinical-adjacent tool.
  `age at menopause`/`menarche` are the exception: benign, useful reproductive-planning
  context, sex-specific. **Verdict: exclude the behavioral ones entirely; consider
  age-at-menopause/menarche in a reproductive sub-section.**

### H. Pharmacogenomics — a different tool, not this one
Warfarin dose, statin-induced myopathy, clopidogrel response, DPYD/TPMT, etc.

- **Tradeoffs:** these are **not polygenic-score-shaped** — they're few-variant,
  star-allele / guideline-driven (CPIC), and belong to the sibling **pgx-suite**, not a
  PGS report. Adding them here would blur the producer/consumer boundary. **Verdict:
  out of scope by design; defer to pgx-suite.**

## Cross-cutting cost of any addition

- **Prep union growth.** Each trait adds its scoring loci to the force-genotyped union
  (currently 12.6M). Biomarkers with million-sample GWAS can carry large variant sets;
  budget prep + scoring wall-clock accordingly.
- **Absolute-risk curation is the real work.** The score is free (Catalog); the honest
  "1 in N" needs a per-SD effect + a baseline with sources. That's ~an hour of literature
  work per trait and is the actual bottleneck to a *useful* (not just present) addition.
- **Redundancy audit.** Before adding, check it isn't a proxy for an existing trait
  (SBP/hypertension, urate/gout, BMD/osteoporosis, apoB/LDL).

## Recommended next tranche (highest value / lowest risk)

1. **Diseases with easy baselines + real actionability:** peripheral arterial disease,
   endometriosis, obstructive sleep apnea — curate effect+baseline so they get a "1 in N."
2. **Cancers with SEER baselines:** esophageal, gastric, cervical (HPV caveat).
3. **A small biomarker section (new report mode):** bone mineral density, eGFR, systolic
   blood pressure — reported as genetic predisposition (direction + percentile), not risk.

Hold the thin-evidence C/D diseases (auto-lift when the Catalog improves), keep
psychiatric/behavioral traits behind the sensitive-tier gate or out entirely, and leave
pharmacogenomics to pgx-suite.
