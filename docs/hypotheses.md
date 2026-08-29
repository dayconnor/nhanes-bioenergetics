# Project Bioenergetics — Hypotheses

Central question: Does Ray Peat's conversion-centric model of thyroid function hold up against national health data?

This document is the pre-registered methods record: hypotheses and exclusion
criteria as decided before/during analysis. Every item below was an open
question at some point during the project; all are now resolved, with a
pointer to where each decision is actually implemented. Results and
interpretation live in [`README.md`](../README.md), not here.

---

## Study population (applies to all hypotheses unless noted)
Non-pregnant adults (age ≥ 20), not taking thyroid medication, no self-reported
thyroid disease.
Cycles: NHANES 2007–2008, 2009–2010, 2011–2012 (the thyroid-panel cycles).

Exclusion variables (confirmed against real data, not assumed):
- Age: `RIDAGEYR >= 20`
- Pregnancy: `RIDEXPRG != 1` (DEMO; populated only for women aged 20-44, `NaN`
  for everyone else — `!=` correctly treats `NaN` as "not pregnant" since
  `NaN != 1` evaluates `True` in pandas)
- Self-reported thyroid disease: `MCQ160M != 1` (MCQ, "ever told you had a
  thyroid problem")
- Thyroid medication: any drug classified `THYROID HORMONES` or
  `ANTITHYROID AGENTS` in `RXQ_DRUG`'s `RXDDCN1B`, joined via `RXQ_RX`

Implemented identically in `pipeline/02_build_h1_sample.py`,
`03_build_h2_sample.py`, `04_build_h3_sample.py`, `05_build_h4_sample.py`.

---

## H1 — FLAGSHIP: conversion → triglycerides
Among the study population, a higher free-T3-to-free-T4 ratio (LBXT3F / LBXT4F)
is associated with LOWER fasting triglycerides (LBXTR), independent of age and sex
— and the body-size-independent portion of that association is tested by
additionally adjusting for BMI.

- Two models, reported as different questions, NOT one:
    (a) adjusted for age + sex → TOTAL association (includes path through body fat)
    (b) also adjusted for BMI -> body-size-INDEPENDENT association
- Sample/weight: triglycerides = FASTING subsample → `WTSAF2YR`, pooled by
  dividing by 3 (three cycles combined).
- Falsified if: the ratio's coefficient on triglycerides is null or positive
  after adjustment.
- **Result: falsified.** Coefficient positive and significant in both models
  (barely moved between them, +14.9 → +14.2), survives multiple-comparison
  correction. See `pipeline/02_build_h1_sample.py`, `notebooks/02_h1_weighting.ipynb`.
- Waist-circumference robustness swap (BMXWAIST in place of BMXBMI) was
  scoped as optional and was **not completed** — deprioritized given the
  project timeline, not attempted and abandoned. Would be the natural next
  step if this project is picked up again.

## H2 — COMPETING MARKERS: free T3 vs TSH
Among the study population, free T3 (LBXT3F) is a stronger predictor of
metabolic-syndrome status than TSH (LBXTSH1) is, both adjusted for age, sex, BMI.

- "Stronger" defined by three metrics: AUC (headline), AIC (corroboration),
  standardized coefficient (mechanism — both predictors z-scored before
  fitting so coefficients are directly comparable).
- Both models run on the identical n (verified explicitly, not assumed —
  both `LBXT3F` and `LBXTSH1` required non-null by the sample-building step).
- Weighted throughout (not left as unweighted future work — `svy` was
  already in place by the time H2 was built).
- Falsified if: TSH's AUC is equal to or higher than free T3's.
- **Result: ambiguous, not falsified by the literal AUC test** (0.8323 vs
  0.8304 — free T3 higher, but the margin is thin enough that no formal
  test of AUC difference was computed). AIC modestly favors free T3. The
  standardized-coefficient comparison cuts the other way: free T3's point
  estimate is larger but not significant (p=0.067); TSH's smaller effect is
  significant (p=0.028). See `pipeline/03_build_h2_sample.py`,
  `notebooks/03_h2_weighting.ipynb`.

## H3: free T3 -> inflammation
Among the study population, LOWER free T3 (LBXT3F) is associated with HIGHER
C-reactive protein (LBXCRP), independent of age, sex, and BMI.

- DATA CONSTRAINT: CRP (LBXCRP) exists ONLY in 2007–2008 and 2009–2010. H3
  runs on a 2-cycle sample, pooled weight `WTMEC2YR / 2` (not /3), full MEC
  exam sample (no fasting-subsample filter — CRP is non-fasting).
- CRP is heavily right-skewed (mean 0.42, median 0.19, max 20 mg/dL); primary
  model uses `log(LBXCRP)` (min=0.01, safely positive, no zero-handling
  needed). Raw-CRP model kept only as a less-trustworthy sensitivity check,
  not co-equal with the log model.
- Interpretation caveat (CENTRAL to this one): reverse causation via
  non-thyroidal illness / "euthyroid sick" syndrome would predict a spurious
  NEGATIVE correlation — it does not straightforwardly explain the POSITIVE
  finding actually observed, so no causal story is claimed in either
  direction.
- Falsified if: the fT3–CRP association is null or positive after adjustment.
- **Result: significant and positive (wrong direction) on the primary
  log-CRP model (p=0.027 raw) — but does not survive FDR correction
  (adjusted p=0.054) once tested alongside H1/H2/H4.** Report as "suggestive
  but not statistically robust," not "falsified with significance." See
  `pipeline/04_build_h3_sample.py`, `notebooks/04_h3_weighting.ipynb`.

## H4 — THE CHOLESTEROL CLAIM
Among the study population AND excluding lipid-lowering medication users, higher
free T3 (LBXT3F) is associated with LOWER total cholesterol (LBXTC), independent
of age, sex, and BMI.

- Sample/weight: total cholesterol is NON-fasting → full MEC sample,
  `WTMEC2YR`, pooled /3, larger n than H1/H2.
- **Lipid-lowering medication, resolved: the broad definition — everyone on
  any drug in `RXQ_DRUG`'s `ANTIHYPERLIPIDEMIC AGENTS` category** (statins,
  fibrates, ezetimibe, bile acid sequestrants, PCSK9 inhibitors), not statins
  alone. Chosen for consistency with the thyroid-med exclusion's scope and to
  avoid hand-picking a subset of drug classes by their own presumed potency —
  a pre-specifiable category-level rule is more defensible than a curated list.
  Accepted trade-off: narrows generalizability further than a statins-only
  exclusion would (older, sicker, more medicated population excluded),
  stated explicitly.
- PRIMARY model excludes these users (n=6,311). SENSITIVITY check keeps them,
  adjusts for an `on_lipid_lowering` covariate (n=7,697).
- Medication merge (multi-row-per-SEQN `RXQ_RX` → one per-person flag):
  implemented once, reused across the thyroid-med, MetS treatment-clause, and
  this lipid-lowering exclusion — same pattern each time.
- Falsified if: the fT3 coefficient on total cholesterol is null or positive
  after adjustment.
- **Result: null under both specifications** (primary p=0.432, sensitivity
  p=0.235, both positive-direction) — a robust null, not a fragile one. See
  `pipeline/05_build_h4_sample.py`, `notebooks/05_h4_weighting.ipynb`.

---

## Metabolic syndrome flag (needed for H2)

Definition: 2009 Harmonized criteria (Alberti et al., *Circulation* 2009,
verified directly against the primary source, not assumed), using the AHA/NHLBI
(ATP III) waist cutpoints. Positive = meets ≥ 3 of 5 criteria.

  1. Waist: ≥102 cm (men) / ≥88 cm (women)          [BMXWAIST, sex-branch on RIAGENDR]
  2. Triglycerides: ≥150 mg/dL                       [LBXTR — fasting]
  3. HDL: <40 (men) / <50 (women) mg/dL              [LBDHDD, sex-branch]
  4. Blood pressure: systolic ≥130 OR diastolic ≥85  [BPX — mean of up to 4 readings; a
     diastolic reading of exactly 0 is NHANES's "no audible sound" code, excluded
     from the average, not treated as a real value]
  5. Fasting glucose: ≥100 mg/dL                     [LBXGLU — fasting]

**Treatment clauses (resolved — confirmed directly from the Alberti et al.
Table 1):** drug treatment counts as an alternate positive
indicator for triglycerides/HDL (fibrates + niacin), blood pressure (any antihypertensive class), and
glucose (any antidiabetic agent, including insulin). No treatment clause
exists for waist circumference. Drug categories implemented in
`pipeline/03_build_h2_sample.py` (`THYROID_MED_CATEGORIES`,
`LIPID_LOWERING_DRUGS`, `ANTIHYPERTENSIVE_CATEGORIES`, `ANTIDIABETIC_CATEGORIES`).

**Missing-component logic (resolved):** rather than excluding anyone with any
missing component, each person is classified as definite-positive (≥3 known
positives already), definite-negative (even a best-case assumption on missing
components can't reach 3), or ambiguous (could go either way given
what's missing). Only the ambiguous group is excluded. In practice this is
56 of 3,459 people (1.6%). Implemented as `flag_component()` /
`build_metabolic_syndrome_flag()` in `pipeline/03_build_h2_sample.py`.

---

## H5 — Exploratory: "practical" thyroid markers (post-hoc, not pre-registered)

Added after H1-H4 were already analyzed, in response to the Peat-sphere objection that serum free T3/T4/TSH are poor proxies for
cellular thyroid status, and that resting pulse rate, basal body temperature,
total cholesterol, and prolactin are better. Body temperature and prolactin
are not available in these NHANES cycles; pulse rate (`BPXPLS`) and total
cholesterol (`LBXTC`) are, so those two are tested against the same three
outcomes as H1-H3 (triglycerides, metabolic syndrome, log-CRP), reusing each
hypothesis's exact population/weights for direct comparability.

This section is explicitly post-hoc and carries a real methodological risk
flagged before it was built: declaring the pre-registered markers "wrong"
specifically because they gave unfavorable results, then reaching for markers
the framework's own proponents prefer, is a textbook unfalsifiable move. This
section exists to test those markers honestly — including the possibility
they'd also come back null or wrong-direction — not to rescue H1-H4.
Corrected separately from H1-H4's multiplicity family for the same reason a
pre-registered and a post-hoc test shouldn't be pooled into one correction.

Pulse-rate models exclude beta-blocker users specifically (`RXDDCN1B ==
'BETA-ADRENERGIC BLOCKING AGENTS'`), which are narrower than H2's broad
antihypertensive category, since beta-blockers uniquely suppress heart rate
pharmacologically, independent of thyroid status.

**Result: pulse rate reinforces H1/H3's disconfirming direction rather than
rescuing the hypothesis** (all three pulse models significant, same "wrong"
direction; a descriptive calibration check confirms pulse correlates
positively with free T3, as expected physiologically, so the two markers are
telling a consistent story). **Cholesterol's result is directionally
Peat-consistent but uninformative**, where higher cholesterol predicting worse
MetS/CRP outcomes is already fully explained by conventional cardiometabolic
physiology with no thyroid mechanism required; this data cannot distinguish
that explanation from the Peat framework's. See
`pipeline/06_build_h5_markers.py`, `notebooks/07_h5_exploratory_markers.ipynb`.

---

## Variables used (all verified against real data during analysis, not assumed)
`LBXT3F`, `LBXT4F`, `LBXTSH1`, `LBXTR`, `LBXTC`, `LBDHDD`, `LBXGLU`, `LBXCRP`,
`BMXWAIST`, `BMXBMI`, `BPXPLS`, `RIDAGEYR`, `RIAGENDR`, `RIDEXPRG`,
`MCQ160M`, `WTMEC2YR`, `WTSAF2YR`, `SDMVSTRA`, `SDMVPSU`, `RXDDRGID`,
`RXDDCN1B` — plus `BPXSY1-4`/`BPXDI1-4` (blood pressure readings) and
`RXQ_RX`/`RXQ_DRUG` medication/class codes.
