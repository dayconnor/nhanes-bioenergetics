import pandas as pd

# RXQ_DRUG lists up to four therapeutic categories per drug. Second-level class
# names live in RXDDCN1B-RXDDCN4B and third-level subclass names in
# RXDDCN1C-RXDDCN4C. Combination products are often filed under a combination
# class first (e.g. amlodipine/atorvastatin is ANTIHYPERTENSIVE COMBINATIONS in
# slot 1 and ANTIHYPERLIPIDEMIC AGENTS in slot 2), so every slot is checked.
class_cols = ['RXDDCN1B', 'RXDDCN2B', 'RXDDCN3B', 'RXDDCN4B']
subclass_cols = ['RXDDCN1C', 'RXDDCN2C', 'RXDDCN3C', 'RXDDCN4C']

thyroid_med_classes = ['THYROID HORMONES', 'ANTITHYROID AGENTS']
antihyperlipidemic_classes = ['ANTIHYPERLIPIDEMIC AGENTS']
antidiabetic_classes = ['ANTIDIABETIC AGENTS']
antihypertensive_classes = [
    'ANGIOTENSIN CONVERTING ENZYME (ACE) INHIBITORS', 'ANGIOTENSIN II INHIBITORS',
    'ANGIOTENSIN RECEPTOR BLOCKERS AND NEPRILYSIN INHIBITORS', 'ANTIHYPERTENSIVE COMBINATIONS',
    'BETA-ADRENERGIC BLOCKING AGENTS', 'CALCIUM CHANNEL BLOCKING AGENTS', 'DIURETICS',
    'ALDOSTERONE RECEPTOR ANTAGONISTS', 'RENIN INHIBITORS', 'VASODILATORS',
    'ANTIADRENERGIC AGENTS, CENTRALLY ACTING', 'ANTIADRENERGIC AGENTS, PERIPHERALLY ACTING',
]
beta_blocker_classes = ['BETA-ADRENERGIC BLOCKING AGENTS']
# Beta-blocker/thiazide combinations (e.g. atenolol; chlorthalidone) carry only
# ANTIHYPERTENSIVE COMBINATIONS at the class level, so they are caught here.
beta_blocker_subclasses = ['BETA BLOCKERS WITH THIAZIDES']

# Metabolic syndrome TG/HDL treatment clause (Alberti et al. 2009): fibrates and
# nicotinic acid. Matched per ingredient so that combination products such as
# niacin; simvastatin are included.
mets_lipid_ingredients = ['CLOFIBRATE', 'GEMFIBROZIL', 'FENOFIBRATE', 'FENOFIBRIC ACID', 'NIACIN']


def merge_drug_classes(rxq_rx_df, rxq_drug_df):
    drug_cols = ['RXDDRGID', 'RXDDRUG'] + class_cols + subclass_cols
    return pd.merge(
        left=rxq_rx_df[['SEQN', 'RXDDRGID']], right=rxq_drug_df[drug_cols],
        on='RXDDRGID', how='left',
    )


def seqns_on(rx_merged, classes=(), subclasses=(), ingredients=()):
    in_class = rx_merged[class_cols].isin(classes).any(axis=1)
    in_subclass = rx_merged[subclass_cols].isin(subclasses).any(axis=1)
    drug_ingredients = rx_merged['RXDDRUG'].fillna('').str.split('; ')
    has_ingredient = drug_ingredients.apply(lambda parts: any(p in ingredients for p in parts))
    return set(rx_merged.loc[in_class | in_subclass | has_ingredient, 'SEQN'])


def thyroid_subsample_weight(df, cycle_label, base_weight):
    # The thyroid panel covered the full MEC exam sample in 2007-2008 but only a
    # random one-third subsample from 2009-2010 on, which CDC weights with
    # WTSA2YR. Using WTMEC2YR or WTSAF2YR alone for the later cycles ignores that
    # subsampling and shrinks those cycles to about a third of their correct
    # share of a pooled estimate.
    #
    # base_weight='WTMEC2YR' (non-fasting outcomes): WTSA2YR is the exact weight.
    # base_weight='WTSAF2YR' (fasting outcomes): NHANES publishes no weight for
    # the fasting x thyroid intersection, so the fasting weight is inflated by
    # each person's inverse thyroid-subsampling factor (WTSA2YR / WTMEC2YR).
    # This treats the two subsamples as independent and brings each cycle's
    # weighted total back in line with 2007-2008's.
    if cycle_label == '2007-2008':
        return df[base_weight]
    if base_weight == 'WTMEC2YR':
        return df['WTSA2YR']
    return df['WTSAF2YR'] * df['WTSA2YR'] / df['WTMEC2YR']


def implausible_free_t3(df):
    # Data-quality exclusion: free T3 above twice the upper reference limit
    # (4.2 pg/mL) while total T3 is within its reference range (<= 200 ng/dL)
    # cannot both be true physiologically, so the record is treated as a
    # laboratory or recording error. Concordantly high values (high free and
    # total T3 with suppressed TSH) are kept as undiagnosed hyperthyroidism.
    return (df['LBXT3F'] > 8.4) & (df['LBXTT3'] <= 200)
