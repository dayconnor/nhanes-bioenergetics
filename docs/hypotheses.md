# Project Bioenergetics hypotheses

Central question. Does Ray Peat's conversion-centric model of thyroid function
hold up against national health data?

This document is the methods record for the project. It sets out the four
pre-specified hypotheses, their falsification criteria, and the exclusion
rules, together with a pointer to where each decision is implemented. The
hypotheses and criteria were set before the models were fit. This document
was committed after the analyses were complete, however, so it is not a
timestamped pre-registration. Results are summarized here for reference, and
the full interpretation is in [`README.md`](../README.md).

## Study population

These rules apply to all hypotheses unless noted. The population is
non-pregnant adults aged 20 and older who are not taking thyroid medication
and do not report thyroid disease, drawn from NHANES 2007-2008, 2009-2010, and
2011-2012 (the cycles with a thyroid panel).

The exclusion variables were confirmed against the data.

- Age is `RIDAGEYR >= 20`.
- Pregnancy is `RIDEXPRG != 1`. The DEMO variable is populated only for women
  aged 20 to 44 and is missing for everyone else, and `NaN != 1` evaluates as
  true in pandas, so missing values are correctly treated as not pregnant.
- Self-reported thyroid disease is `MCQ160M != 1` ("ever told you had a
  thyroid problem").
- Thyroid medication is any drug in the `THYROID HORMONES` or
  `ANTITHYROID AGENTS` class, from `RXQ_RX` joined to `RXQ_DRUG`.
- One record (SEQN 53952) is excluded as a laboratory or recording error. The
  rule is free T3 above 8.4 pg/mL (twice the upper reference limit) with
  total T3 at or below 200 ng/dL, which cannot both be true physiologically.
  The pipeline flags the record (`ft3_implausible`) and the notebooks drop it,
  so that `notebooks/09_robustness_checks.ipynb` can report results with it
  retained.

Medication classes are matched against all four category slots in `RXQ_DRUG`
(`RXDDCN1B` to `RXDDCN4B`), and against the third-level subclass where a class
is only recorded there, so that combination products are classified
correctly. The shared rules are in `pipeline/nhanes_helpers.py` and are used
by `pipeline/02_build_h1_sample.py` through `pipeline/06_build_h5_markers.py`.

## Survey weights

The thyroid panel covered the full examination sample in 2007-2008 and a
random one-third subsample from 2009-2010 onward, which CDC weights with
`WTSA2YR`. Each analytic sample carries a thyroid-adjusted two-year weight,
`wt_thyroid_2yr`, which is pooled by dividing by the number of cycles.

- Non-fasting outcomes (H3, H4) use `WTMEC2YR` in 2007-2008 and `WTSA2YR`
  thereafter.
- Fasting outcomes (H1, H2) use `WTSAF2YR` in 2007-2008 and
  `WTSAF2YR * WTSA2YR / WTMEC2YR` thereafter. NHANES publishes no weight for
  the intersection of the fasting and thyroid subsamples, so this
  approximation treats the two subsamples as independent.

## H1. Conversion and triglycerides (flagship)

Among the study population, a higher free-T3-to-free-T4 ratio
(`LBXT3F / LBXT4F`) is associated with lower fasting triglycerides (`LBXTR`),
independent of age and sex. The association independent of body size is
tested by additionally adjusting for BMI.

- Two models answer different questions. Model (a) adjusts for age and sex and
  estimates the total association, including any path through body fat.
  Model (b) adds BMI and estimates the association independent of body size.
- Triglycerides are measured in the fasting subsample, so the base weight is
  `WTSAF2YR`.
- Falsified if the ratio's coefficient on triglycerides is null or positive
  after adjustment.
- Result. Falsified. The coefficient is positive and significant in both
  models (+20.3 in model a and +19.7 in model b, both p < 0.001) and survives
  multiple-comparison correction. See `pipeline/02_build_h1_sample.py` and
  `notebooks/02_h1_weighting.ipynb`.
- A robustness check replacing BMI with waist circumference (`BMXWAIST`) was
  scoped as optional and was not completed because of the project timeline.
  It would be the natural next step if the project is resumed.
- A post-hoc sensitivity check replaced free T3 and free T4 with total T3 and
  total T4 (`LBXTT3`/`LBXTT4`). Peat's writing is skeptical of free hormone
  measurements and favors total hormone levels, so this tests the same claim
  with the measure the framework considers more valid. The result has the
  same disconfirming direction, with a coefficient of +2.97 in model (a) and
  +2.99 in model (b) (p = 0.024 in both). See `notebooks/02_h1_weighting.ipynb`.

## H2. Competing markers, free T3 and TSH

Among the study population, free T3 (`LBXT3F`) is a stronger predictor of
metabolic syndrome than TSH (`LBXTSH1`), with both adjusted for age, sex, and
BMI.

- "Stronger" is defined by three metrics. AUC is the headline metric, AIC
  provides corroboration, and the standardized coefficient addresses
  mechanism. Both predictors are z-scored before fitting so that their
  coefficients are directly comparable.
- Both models run on the same participants. The sample-building step requires
  both `LBXT3F` and `LBXTSH1` to be present.
- All models are survey weighted.
- Falsified if TSH's AUC is equal to or higher than free T3's.
- Result. Supported as stated, with a small margin. Free T3 leads on all
  three metrics. Its AUC is 0.829 against 0.824 for TSH, its AIC is lower
  (3,275 against 3,307), and its standardized coefficient is larger and
  significant (0.38, p < 0.001) where TSH's is not (0.16, p = 0.058). Neither
  marker adds much to age, sex, and BMI alone (AUC 0.823). The free T3
  coefficient is positive, so higher free T3 is associated with higher odds
  of metabolic syndrome, which runs against the framework's premise. No formal
  test of the AUC difference was computed. Because the AUC criterion has no
  p-value, the free T3 coefficient represents H2 in the multiple-comparison
  correction. See `pipeline/03_build_h2_sample.py` and
  `notebooks/03_h2_weighting.ipynb`.

## H3. Free T3 and inflammation

Among the study population, lower free T3 (`LBXT3F`) is associated with higher
C-reactive protein (`LBXCRP`), independent of age, sex, and BMI.

- CRP exists only in 2007-2008 and 2009-2010, so H3 runs on a two-cycle sample
  pooled by dividing by 2. CRP is a non-fasting measure, so no fasting
  restriction applies.
- CRP is heavily right-skewed (mean 0.42, median 0.19, maximum 20 mg/dL). The
  primary model uses `log(LBXCRP)`, which needs no handling of zeros because
  the minimum is 0.01. A model on untransformed CRP is kept only as a less
  reliable sensitivity check.
- The central interpretive caveat is reverse causation through non-thyroidal
  illness ("euthyroid sick" syndrome). That mechanism would produce a
  spurious negative correlation and does not explain a positive one, so no
  causal interpretation is offered in either direction.
- Falsified if the free T3 association with CRP is null or positive after
  adjustment.
- Result. Falsified. The association is positive and significant in the
  primary log-CRP model (0.23 per pg/mL, p = 0.003) and survives correction
  (FDR-adjusted p = 0.004, Bonferroni-adjusted p = 0.012). The untransformed
  model is null (p = 0.63). With the excluded record retained, the primary
  association is not significant (p = 0.072). See
  `pipeline/04_build_h3_sample.py` and `notebooks/04_h3_weighting.ipynb`.

## H4. The cholesterol claim

Among the study population, excluding users of lipid-lowering medication,
higher free T3 (`LBXT3F`) is associated with lower total cholesterol
(`LBXTC`), independent of age, sex, and BMI.

- Total cholesterol is a non-fasting measure taken in the full examination
  sample, so H4 has a larger sample than H1 and H2.
- Lipid-lowering medication is defined broadly as any drug in the
  `ANTIHYPERLIPIDEMIC AGENTS` class (statins, fibrates, ezetimibe, bile acid
  sequestrants, niacin, and PCSK9 inhibitors) rather than statins alone. A
  category-level rule that can be stated in advance is more defensible than a
  hand-picked list of drug classes, and it matches the scope of the thyroid
  medication exclusion. The trade-off is reduced generalizability, because an
  older, sicker, and more medicated group is excluded.
- The primary model excludes these users (n = 6,283). The sensitivity model
  keeps them and adjusts for an `on_lipid_lowering` covariate (n = 7,696).
- Falsified if the free T3 coefficient on total cholesterol is null or
  positive after adjustment.
- Result. Falsified. The primary model is positive and not significant
  (4.5 mg/dL per pg/mL, p = 0.094). The sensitivity model is positive and
  significant (6.0 mg/dL per pg/mL, p = 0.021). The positive association is
  significant among euthyroid participants but absent in 2007-2008 alone, so
  it is not stable. See `pipeline/05_build_h4_sample.py` and
  `notebooks/05_h4_weighting.ipynb`.

## Metabolic syndrome flag (used in H2)

The definition follows the 2009 Harmonized criteria (Alberti et al.,
*Circulation* 2009, checked against the primary source) with the AHA/NHLBI
(ATP III) waist cutpoints. A participant is positive if they meet at least
three of five criteria.

1. Waist of at least 102 cm in men or 88 cm in women (`BMXWAIST`, by
   `RIAGENDR`).
2. Fasting triglycerides of at least 150 mg/dL (`LBXTR`).
3. HDL below 40 mg/dL in men or 50 mg/dL in women (`LBDHDD`, by `RIAGENDR`).
4. Systolic blood pressure of at least 130 or diastolic of at least 85 mmHg
   (`BPX`, the mean of up to four readings). A diastolic reading of exactly 0
   is NHANES's code for no audible sound and is excluded from the mean.
5. Fasting glucose of at least 100 mg/dL (`LBXGLU`).

Drug treatment counts as an alternate positive indicator, as specified in
Table 1 of Alberti et al. Fibrates and nicotinic acid count for the
triglyceride and HDL criteria, matched by ingredient so that combination
products are included. Any antihypertensive class counts for blood pressure,
including ACE inhibitors, angiotensin II inhibitors, beta-blockers, calcium
channel blockers, diuretics, aldosterone receptor antagonists, renin
inhibitors, vasodilators, centrally and peripherally acting antiadrenergic
agents, and antihypertensive combinations. Any antidiabetic agent, including
insulin, counts for glucose. Waist circumference has no treatment clause.

Participants with missing components are not excluded outright. Each is
classified as definite positive (at least three known positives), definite
negative (unable to reach three even if every missing component were
positive), or ambiguous. Only the ambiguous group is excluded, which in
practice is 57 of 3,458 participants (1.6%). This is implemented in
`flag_component()` and `build_metabolic_syndrome_flag()` in
`pipeline/03_build_h2_sample.py`.

## H5. Exploratory analysis of practical thyroid markers (post-hoc)

H5 was added after H1-H4 were analyzed and was not pre-specified. It responds
to the argument, common in Peat-oriented online discussion, that serum free
T3, free T4, and TSH are poor proxies for cellular thyroid status and that
resting pulse rate, basal body temperature, total cholesterol, and prolactin
are better. Body temperature and prolactin are not available in these NHANES
cycles. Pulse rate (`BPXPLS`) and total cholesterol (`LBXTC`) are, so both are
tested against the outcomes of H1-H3 (triglycerides, metabolic syndrome, and
log-CRP), reusing each hypothesis's population and weights for direct
comparability.

The methodological risk was identified before the section was built.
Declaring the pre-specified markers invalid because they gave unfavorable
results, and then turning to the markers the framework's proponents prefer,
would make the framework unfalsifiable. The section therefore tests these
markers on their own terms, with the possibility that they also come back
null or in the wrong direction, and is corrected as a separate family from
H1-H4.

Pulse-rate models exclude users of beta-blockers (`BETA-ADRENERGIC BLOCKING
AGENTS`, plus the `BETA BLOCKERS WITH THIAZIDES` combination subclass). This
is narrower than H2's antihypertensive category because beta-blockers
specifically suppress heart rate, independent of thyroid status. Cholesterol
models exclude users of any `ANTIHYPERLIPIDEMIC AGENTS` drug, since their
cholesterol reflects the medication.

Result. Pulse rate reinforces the disconfirming direction of H1 and H3. All
three pulse models are significant and point the same way, and a descriptive
calibration check confirms that pulse rises with free T3, as expected
physiologically, so the two markers tell a consistent story. Cholesterol's
results are consistent with the framework in direction but uninformative,
because conventional cardiometabolic physiology already explains why higher
cholesterol accompanies metabolic syndrome and higher CRP, with no thyroid
mechanism required. See `pipeline/06_build_h5_markers.py` and
`notebooks/07_h5_exploratory_markers.ipynb`.

## Variables used

All variables were checked against the data during analysis.

`LBXT3F`, `LBXT4F`, `LBXTT3`, `LBXTT4`, `LBXTSH1`, `LBXTPO`, `LBXTR`, `LBXTC`,
`LBDHDD`, `LBXGLU`, `LBXCRP`, `BMXWAIST`, `BMXBMI`, `BPXPLS`, `RIDAGEYR`,
`RIAGENDR`, `RIDEXPRG`, `MCQ160M`, `WTMEC2YR`, `WTSAF2YR`, `WTSA2YR`,
`SDMVSTRA`, `SDMVPSU`, `RXDDRGID`, `RXDDRUG`, `RXDDCN1B` to `RXDDCN4B`, and
`RXDDCN1C` to `RXDDCN4C`, plus the blood pressure readings `BPXSY1` to
`BPXSY4` and `BPXDI1` to `BPXDI4`, and the `RXQ_RX` and `RXQ_DRUG` medication
files.
