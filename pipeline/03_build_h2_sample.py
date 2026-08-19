import pathlib

import numpy as np
import pandas as pd
import pyreadstat

repo_root = pathlib.Path(__file__).resolve().parents[1]
raw_dir = repo_root / 'data' / 'raw'
processed_dir = repo_root / 'data' / 'processed'

thyroid_med_categories = ['THYROID HORMONES', 'ANTITHYROID AGENTS']
lipid_lowering_drugs = ['CLOFIBRATE', 'GEMFIBROZIL', 'FENOFIBRATE', 'FENOFIBRIC ACID', 'NIACIN']
antihypertensive_categories = [
    'ANGIOTENSIN CONVERTING ENZYME (ACE) INHIBITORS', 'ANGIOTENSIN II INHIBITORS',
    'ANGIOTENSIN RECEPTOR BLOCKERS AND NEPRILYSIN INHIBITORS', 'ANTIHYPERTENSIVE COMBINATIONS',
    'BETA-ADRENERGIC BLOCKING AGENTS', 'CALCIUM CHANNEL BLOCKING AGENTS', 'DIURETICS',
]
antidiabetic_categories = ['ANTIDIABETIC AGENTS']

cycles = {
    '2007-2008': {
        'thyroid_path': raw_dir / '2007-2008' / 'THYROD_E.xpt',
        'trigly_path': raw_dir / '2007-2008' / 'TRIGLY_E.xpt',
        'demo_path': raw_dir / '2007-2008' / 'DEMO_E.xpt',
        'mcq_path': raw_dir / '2007-2008' / 'MCQ_E.xpt',
        'rxq_rx_path': raw_dir / '2007-2008' / 'RXQ_RX_E.xpt',
        'bmx_path': raw_dir / '2007-2008' / 'BMX_E.xpt',
        'hdl_path': raw_dir / '2007-2008' / 'HDL_E.xpt',
        'glu_path': raw_dir / '2007-2008' / 'GLU_E.xpt',
        'bpx_path': raw_dir / '2007-2008' / 'BPX_E.xpt',
    },
    '2009-2010': {
        'thyroid_path': raw_dir / '2009-2010' / 'THYROD_F.xpt',
        'trigly_path': raw_dir / '2009-2010' / 'TRIGLY_F.xpt',
        'demo_path': raw_dir / '2009-2010' / 'DEMO_F.xpt',
        'mcq_path': raw_dir / '2009-2010' / 'MCQ_F.xpt',
        'rxq_rx_path': raw_dir / '2009-2010' / 'RXQ_RX_F.xpt',
        'bmx_path': raw_dir / '2009-2010' / 'BMX_F.xpt',
        'hdl_path': raw_dir / '2009-2010' / 'HDL_F.xpt',
        'glu_path': raw_dir / '2009-2010' / 'GLU_F.xpt',
        'bpx_path': raw_dir / '2009-2010' / 'BPX_F.xpt',
    },
    '2011-2012': {
        'thyroid_path': raw_dir / '2011-2012' / 'THYROD_G.xpt',
        'trigly_path': raw_dir / '2011-2012' / 'TRIGLY_G.xpt',
        'demo_path': raw_dir / '2011-2012' / 'DEMO_G.xpt',
        'mcq_path': raw_dir / '2011-2012' / 'MCQ_G.xpt',
        'rxq_rx_path': raw_dir / '2011-2012' / 'RXQ_RX_G.xpt',
        'bmx_path': raw_dir / '2011-2012' / 'BMX_G.xpt',
        'hdl_path': raw_dir / '2011-2012' / 'HDL_G.xpt',
        'glu_path': raw_dir / '2011-2012' / 'GLU_G.xpt',
        'bpx_path': raw_dir / '2011-2012' / 'BPX_G.xpt',
    },
}


def average_blood_pressure(bpx_df):
    sy_cols = ['BPXSY1', 'BPXSY2', 'BPXSY3', 'BPXSY4']
    di_cols = ['BPXDI1', 'BPXDI2', 'BPXDI3', 'BPXDI4']

    sy = bpx_df[sy_cols]
    # A diastolic reading of exactly 0 means no audible sound was detected.
    # It's NHANES's missing-data code, not a real value, so it must be excluded
    # before averaging or it would drag the average down artificially.
    di = bpx_df[di_cols].replace(0, np.nan)

    out = bpx_df[['SEQN']].copy()
    out['bp_systolic'] = sy.mean(axis=1)
    out['bp_diastolic'] = di.mean(axis=1)
    return out


def flag_component(lab_value, meets_threshold, on_treatment):
    # 1 = positive, 0 = negative, NaN = unknown (lab missing, not on treatment)
    result = pd.Series(np.nan, index=lab_value.index, dtype='float')
    result[meets_threshold & lab_value.notna()] = 1
    result[~meets_threshold & lab_value.notna()] = 0
    result[on_treatment] = 1
    return result


def build_metabolic_syndrome_flag(df):
    waist_meets = ((df['RIAGENDR'] == 1) & (df['BMXWAIST'] >= 102)) | ((df['RIAGENDR'] == 2) & (df['BMXWAIST'] >= 88))
    tg_meets = df['LBXTR'] >= 150
    hdl_meets = ((df['RIAGENDR'] == 1) & (df['LBDHDD'] < 40)) | ((df['RIAGENDR'] == 2) & (df['LBDHDD'] < 50))
    bp_meets = (df['bp_systolic'] >= 130) | (df['bp_diastolic'] >= 85)
    glucose_meets = df['LBXGLU'] >= 100

    components = pd.DataFrame({
        'waist': flag_component(df['BMXWAIST'], waist_meets, pd.Series(False, index=df.index)),
        'tg': flag_component(df['LBXTR'], tg_meets, df['on_lipid_lowering']),
        'hdl': flag_component(df['LBDHDD'], hdl_meets, df['on_lipid_lowering']),
        'bp': flag_component(df['bp_systolic'], bp_meets, df['on_antihypertensive']),
        'glucose': flag_component(df['LBXGLU'], glucose_meets, df['on_antidiabetic']),
    })

    known_positive = (components == 1).sum(axis=1)
    missing = components.isna().sum(axis=1)

    definite_positive = known_positive >= 3
    definite_negative = (known_positive + missing) < 3

    mets = pd.Series(np.nan, index=df.index)
    mets[definite_positive] = 1
    mets[definite_negative & ~definite_positive] = 0
    return mets, components


def medication_seqns(rxq_rx_df, rxq_drug_df):
    merged = pd.merge(left=rxq_rx_df, right=rxq_drug_df, on='RXDDRGID', how='left', suffixes=('_rx', '_drug'))

    thyroid_med = set(merged.loc[merged['RXDDCN1B'].isin(thyroid_med_categories), 'SEQN'])
    antihypertensive = set(merged.loc[merged['RXDDCN1B'].isin(antihypertensive_categories), 'SEQN'])
    antidiabetic = set(merged.loc[merged['RXDDCN1B'].isin(antidiabetic_categories), 'SEQN'])
    lipid_lowering = set(merged.loc[merged['RXDDRUG_drug'].isin(lipid_lowering_drugs), 'SEQN'])

    return thyroid_med, antihypertensive, antidiabetic, lipid_lowering


def load_and_filter_cycle(paths, rxq_drug_df):
    thyroid_df, _ = pyreadstat.read_xport(str(paths['thyroid_path']))
    trigly_df, _ = pyreadstat.read_xport(str(paths['trigly_path']))
    demo_df, _ = pyreadstat.read_xport(str(paths['demo_path']))
    mcq_df, _ = pyreadstat.read_xport(str(paths['mcq_path']))
    rxq_rx_df, _ = pyreadstat.read_xport(str(paths['rxq_rx_path']))
    bmx_df, _ = pyreadstat.read_xport(str(paths['bmx_path']))
    hdl_df, _ = pyreadstat.read_xport(str(paths['hdl_path']))
    glu_df, _ = pyreadstat.read_xport(str(paths['glu_path']))
    bpx_df, _ = pyreadstat.read_xport(str(paths['bpx_path']))

    merged_df = pd.merge(left=thyroid_df, right=trigly_df, on='SEQN', how='inner')
    merged_df = pd.merge(left=merged_df, right=demo_df, on='SEQN', how='left')

    # Study population: thyroid panel + fasting subsample + age >= 20 + not pregnant,
    # same base exclusions as H1 (hypotheses.md: applies to all hypotheses unless noted).
    analytic_df = merged_df.dropna(subset=['LBXT3F', 'LBXTSH1', 'LBXTR'])
    analytic_df = analytic_df[analytic_df['WTSAF2YR'] > 0]
    analytic_df = analytic_df[analytic_df['RIDAGEYR'] >= 20]
    analytic_df = analytic_df[analytic_df['RIDEXPRG'] != 1]

    analytic_df = pd.merge(left=analytic_df, right=mcq_df, on='SEQN', how='left')
    analytic_df = analytic_df[analytic_df['MCQ160M'] != 1]

    thyroid_med_seqns, antihypertensive_seqns, antidiabetic_seqns, lipid_lowering_seqns = medication_seqns(rxq_rx_df, rxq_drug_df)
    analytic_df = analytic_df[~analytic_df['SEQN'].isin(thyroid_med_seqns)]

    analytic_df = pd.merge(left=analytic_df, right=bmx_df[['SEQN', 'BMXWAIST', 'BMXBMI']], on='SEQN', how='left')
    analytic_df = pd.merge(left=analytic_df, right=hdl_df[['SEQN', 'LBDHDD']], on='SEQN', how='left')
    analytic_df = pd.merge(left=analytic_df, right=glu_df[['SEQN', 'LBXGLU']], on='SEQN', how='left')

    bp_avg = average_blood_pressure(bpx_df)
    analytic_df = pd.merge(left=analytic_df, right=bp_avg, on='SEQN', how='left')

    analytic_df['on_lipid_lowering'] = analytic_df['SEQN'].isin(lipid_lowering_seqns)
    analytic_df['on_antihypertensive'] = analytic_df['SEQN'].isin(antihypertensive_seqns)
    analytic_df['on_antidiabetic'] = analytic_df['SEQN'].isin(antidiabetic_seqns)

    analytic_df['metabolic_syndrome'], components = build_metabolic_syndrome_flag(analytic_df)
    analytic_df = pd.concat([analytic_df, components.add_prefix('mets_')], axis=1)

    return analytic_df


def build_h2_sample():
    rxq_drug_df, _ = pyreadstat.read_xport(str(raw_dir / 'RXQ_DRUG.xpt'))

    pooled_frames = []
    for cycle_label, paths in cycles.items():
        analytic_df = load_and_filter_cycle(paths, rxq_drug_df)
        ambiguous = analytic_df['metabolic_syndrome'].isna().sum()
        print(f"{cycle_label}: n = {len(analytic_df)}, metabolic_syndrome ambiguous/excluded = {ambiguous}")
        pooled_frames.append(analytic_df.assign(cycle=cycle_label))

    pooled = pd.concat(pooled_frames, ignore_index=True)
    n_before_mets_drop = len(pooled)
    pooled_with_mets = pooled.dropna(subset=['metabolic_syndrome'])
    print(f"pooled n (before MetS-ambiguous drop) = {n_before_mets_drop}")
    print(f"pooled n (MetS classifiable)           = {len(pooled_with_mets)}")
    print(f"MetS prevalence in classifiable sample  = {pooled_with_mets['metabolic_syndrome'].mean():.3f}")
    return pooled_with_mets


if __name__ == '__main__':
    pooled = build_h2_sample()
    processed_dir.mkdir(parents=True, exist_ok=True)
    out_path = processed_dir / 'h2_analytic_sample.csv'
    pooled.to_csv(out_path, index=False)
    print(f"saved to {out_path}")
