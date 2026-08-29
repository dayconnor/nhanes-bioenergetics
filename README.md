# Project Bioenergetics

An end-to-end analysis of the thyroid–metabolic axis in U.S. adults, using NHANES
public health data to test Ray Peat's bioenergetic hypotheses about thyroid
function.

Peat's framework argues that efficient thyroid hormone conversion (T4 -> T3) and
adequate free T3 drive efficient, pro-metabolic physiology, protecting against
elevated triglycerides, inflammation, and other markers of metabolic
dysfunction. This project asks: does that framework hold up against national
health survey data? The honest answer, detailed below, is
**mostly no** — three of four pre-registered hypotheses were not supported, and
the one significant result runs in the opposite direction from what the
framework predicts. A post-hoc exploratory section (H5) tests two "practical"
markers the framework's own proponents prefer over serum thyroid tests; it
doesn't change this picture.

## Data

[NHANES](https://www.cdc.gov/nchs/nhanes/) (National Health and Nutrition
Examination Survey), cycles **2007-2008, 2009-2010, 2011-2012** — the only
three cycles with a full thyroid panel (free T3, free T4, TSH, total T3/T4).
NHANES is a stratified, multistage probability sample of the U.S. civilian
population; all analyses use design-based survey weights, strata, and primary
sampling units (PSUs), not naive unweighted statistics.

Study population (all hypotheses): non-pregnant adults age ≥ 20, excluding
self-reported thyroid disease and thyroid/antithyroid medication users.
Hypothesis-specific exclusions (fasting subsample, lipid-lowering medication
use) are documented per hypothesis below and in
[`docs/hypotheses.md`](docs/hypotheses.md).

## Methodology

- **Survey-weighted regression** via [`svy`](https://svylab.com), using each cycle's appropriate weight
  (`WTSAF2YR` for fasting-subsample outcomes, `WTMEC2YR` for non-fasting),
  pooled across cycles by dividing by the number of cycles combined, which is the
  standard NHANES approach for multi-cycle analysis.
- **Confounder adjustment**: age, sex, and BMI at minimum; hypothesis-specific
  models add further adjustment (e.g. lipid-lowering medication use for H4).
- **Metabolic syndrome** (used in H2) is constructed from the 2009 Harmonized
  criteria ([Alberti et al., *Circulation* 2009](https://doi.org/10.1161/CIRCULATIONAHA.109.192644)),
  including the criteria's drug-treatment alternate indicators, with an
  explicit definite-positive / definite-negative / ambiguous framework for
  missing components rather than blanket-excluding anyone with any missing
  lab value.
- **Multiple-comparison correction** (FDR and Bonferroni) applied across all
  four hypotheses' primary tests, since evaluating several hypotheses at
  p < 0.05 without correction overstates how much evidence you actually have.
- **Reproducible pipeline**: raw NHANES `.XPT` files (git-ignored, not
  redistributed) -> cleaned/merged analytic samples (`pipeline/`, numbered by
  run order) -> exploratory/modeling notebooks (`notebooks/`) -> processed
  output (`data/processed/`, committed since raw data is not).

## Findings

| # | Hypothesis | Result | Direction |
|---|---|---|---|
| **H1** | Higher FT3/FT4 ratio -> lower triglycerides | **Falsified** — significant, opposite direction (survives multiple-comparison correction) | Wrong |
| **H2** | Free T3 outperforms TSH at predicting metabolic syndrome | **Ambiguous** — narrow edge to free T3 on AUC/AIC, but TSH has the statistically significant individual coefficient | Mixed |
| **H3** | Lower free T3 -> higher CRP (inflammation) | **Not supported after correction** — significant and opposite-direction on the primary model alone, but does not survive FDR/Bonferroni correction | Wrong, weak |
| **H4** | Higher free T3 -> lower total cholesterol | **Null** — no significant association, held under both a primary model (excluding lipid-lowering medication users) and a sensitivity check (including them, adjusted) | None |

**Only H1 remains statistically significant once corrected for testing four
hypotheses.** It is also the one result that runs clearly counter to the
Peat framework: a higher free-T3/free-T4 ratio is associated with *higher*,
not lower, triglycerides (pooled n=3,459; coefficient +14.9 mg/dL per unit
ratio, p<0.001, robust to BMI adjustment).

![Primary coefficient across all four hypotheses](docs/figures/summary_all_hypotheses.png)

![Raw vs. corrected p-values across all four hypotheses](docs/figures/multiplicity_correction.png)

### H1 — FT3/FT4 ratio and triglycerides
Two models were fit: age+sex-adjusted (total association) and
+BMI-adjusted (body-size-independent association), since BMI is a plausible
mediator on the causal pathway between thyroid function and triglycerides.
The coefficient barely moved between the two (14.9 to 14.2), indicating the
association is **not substantially explained by body size**.

### H2 — Free T3 vs. TSH as predictors of metabolic syndrome
Compared via three pre-registered metrics: AUC (0.832 vs. 0.830, visually 
indistinguishable on the ROC curve), AIC (a modest edge to
free T3), and standardized coefficient (free T3's point estimate is larger
but not significant at p=0.067; TSH's smaller effect is significant at
p=0.028). No clean winner. Free T3 edges ahead on ranking ability, TSH is
the one with the more statistically certain individual effect.

### H3 — Free T3 and CRP (inflammation)
CRP exists only in the 2007-2008 and 2009-2010 cycles, so this hypothesis
runs on a 2-cycle pooled sample (n=6,271) with its own pooled weight. 
CRP is heavily right skewed (mean 0.42, median 0.19, max 20
mg/dL). The primary model uses log-transformed CRP, which is standard
practice for this outcome and the statistically appropriate choice given the
skew. That model found free T3 *positively* associated with inflammation
(p=0.027 raw), the opposite of the hypothesized direction — though this
result does not survive correction for multiple comparisons (FDR-adjusted
p=0.054). A raw-CRP sensitivity check (less trustworthy given the skew) was
null in both directions. Note: this hypothesis's central caveat — that
illness can suppress T3 while independently raising CRP ("non-thyroidal
illness syndrome"), producing a spurious *negative* correlation — does not
straightforwardly explain a *positive* finding, so no causal story is claimed
either way. Cross-sectional data cannot establish causal direction.

### H4 — Free T3 and total cholesterol
Total cholesterol is measured on the full (non-fasting) MEC exam sample,
giving the largest pooled n of the four hypotheses (7,803). Lipid-lowering
medication users (the full `ANTIHYPERLIPIDEMIC AGENTS` drug class — statins,
fibrates, ezetimibe, bile acid sequestrants, PCSK9 inhibitors) were excluded
in the primary model (n=6,311) and included with an adjustment covariate in
a sensitivity check (n=7,697). Both gave a non-significant, positive-direction
coefficient, which translates to a null result.

![H4 primary vs sensitivity model](docs/figures/h4_primary_vs_sensitivity.png)

## H5 — Exploratory: "practical" thyroid markers (post-hoc)

Bioenergetic/Peat-sphere commentary on X/Twitter (e.g. Danny Roddy) argues that circulating
free T3, free T4, and TSH are poor proxies for tissue-level thyroid status,
preferring resting pulse rate, basal body temperature, total cholesterol, and
prolactin instead. Body temperature and prolactin are not available in these
NHANES cycles; **resting pulse rate (`BPXPLS`) and total cholesterol
(`LBXTC`, already used in H4) are**, so this section re-tests the same three
outcomes (triglycerides, metabolic syndrome, log-CRP) against these two
alternative markers, reusing H1/H2/H3's exact populations, weights, and
covariate adjustment for direct comparability.

**This section is explicitly post-hoc and exploratory**, added after seeing
H1-H4's results, not pre-registered. That ordering carries a real risk: within
the Peat framework, declaring the serum markers "wrong" specifically because
they gave unfavorable results — and then reaching for markers the framework's
own proponents prefer — is a textbook unfalsifiable move, not a data-based
correction. This section exists to test those alternative markers honestly,
not to rescue H1-H4; the results below are reported the same way regardless
of which direction they point. Corrected separately from H1-H4's multiplicity
family, since mixing a pre-registered confirmatory set with a post-hoc
exploratory one would muddle what each corrected p-value actually means.

**Pulse rate models exclude beta-blocker users** — a narrower flag than H2's
broad antihypertensive category, since beta-blockers specifically suppress
heart rate pharmacologically, independent of thyroid status, and would
confound any pulse-based result if left in.

| Test | Coefficient | p (raw) | p (FDR) | p (Bonferroni) | Direction |
|---|---|---|---|---|---|
| Calibration: pulse ~ free T3 (descriptive) | +2.02 | 0.033 | — | — | Higher T3 -> higher pulse (expected; thyroid hormone raises heart rate) |
| Pulse -> triglycerides | +1.26 | <0.001 | 0.0001 | 0.0005 | Wrong direction |
| Pulse -> metabolic syndrome | +0.026 | <0.001 | 0.0001 | 0.0005 | Wrong direction |
| Pulse -> log-CRP | +0.016 | <0.001 | 0.0001 | 0.0005 | Wrong direction |
| Cholesterol -> metabolic syndrome | +0.0035 | 0.0066 | 0.0066 | 0.0330 | Directionally Peat-consistent — see caveat |
| Cholesterol -> log-CRP | +0.0024 | <0.001 | 0.0001 | 0.0005 | Directionally Peat-consistent — see caveat |

All five formal tests survive correction, but the two marker types support different conclusions:

![H5 exploratory markers](docs/figures/h5_exploratory_markers.png)

**Pulse rate reinforces the disconfirmation rather than rescuing the
hypothesis.** All three pulse models point the same "wrong" direction as the
original free-T3 findings. The calibration check shows this is coherent, not
contradictory: free T3 correlates positively with pulse (which is expected, since thyroid
hormone raises heart rate, standard endocrinology, not Peat-specific), and
both markers independently predict *worse*, not better, metabolic outcomes.
Two physiologically-linked markers agreeing is stronger evidence against the
hypothesis than one marker alone.

**Cholesterol's result requires a hard caveat, not a headline.** Higher
cholesterol predicting worse MetS/CRP outcomes is fully explained by ordinary
cardiometabolic physiology — dyslipidemia clusters with insulin resistance
and low-grade inflammation as a matter of routine metabolic biology, with no
reference to thyroid function required. This data cannot distinguish
"cholesterol reflects thyroid status, which drives the outcome" from
"cholesterol is an independent conventional risk marker unrelated to
thyroid" — both predict the identical observed pattern. Reporting this
result as support for the Peat framework without this caveat would
misrepresent what the data can actually show, and is the specific failure
mode this section was built to test rather than commit.

## Limitations

- **Cross-sectional data cannot establish causal direction.** Every
  association reported here could run in either direction, or both could be
  driven by an unmeasured third factor. This is most acute for H3, where
  reverse causation (illness suppressing T3) is a live alternative
  explanation even though it doesn't cleanly fit the direction actually
  observed.
- **Different hypotheses use different sub-populations and weights** (fasting
  subsample vs. full exam sample), by necessity of which lab values require
  fasting — documented per hypothesis above, not an inconsistency.
- **The thyroid panel's subsample design changed after 2007-2008** — CDC
  moved from near-full-sample coverage to a random ~1/3 subsample starting in
  2009-2010, which affects the relative weight of cycles in every pooled
  analysis in this project.
- Multiple-comparison correction was applied across one primary p-value per
  hypothesis; H1 and H4 each generated a second (sensitivity/alternate-model)
  p-value not included in that correction, to avoid double-counting a single
  underlying test.

## Repo structure

```
data/raw/          NHANES .XPT files, git-ignored, regenerated by downloading
                    the files listed in each pipeline script
data/processed/     Cleaned, merged, pooled analytic samples (one CSV per
                    hypothesis, committed)
pipeline/           Reproducible sample-building scripts, numbered by run
                    order (02-06; sample selection only, not modeling)
notebooks/          Exploration, survey-weighted modeling, and visualization
                    (numbered 01-08)
docs/hypotheses.md  The four pre-registered, falsifiable hypotheses and study
                    population definitions, written before analysis began
docs/figures/       Exported presentation-quality charts (PNG), sourced from
                    notebooks/08_summary_visualizations.ipynb
```

## Stack

Python, pandas, Polars, [`svy`](https://svylab.com) (survey-weighted
regression), pyreadstat, statsmodels, seaborn/matplotlib, scikit-learn (AUC/
ROC), SQLite (planned), Streamlit (dashboard, planned).

## Status

All four hypotheses, plus the H5 exploratory extension, are analyzed,
visualized, and multiple-comparison corrected, with presentation-quality
summary charts (`docs/figures/`, built with a validated colorblind-safe
palette). Remaining work: an interactive dashboard, and a public write-up
site (lower priority than the analysis itself, per this project's own
priorities).
