import pathlib

import pandas as pd
import pyreadstat

from nhanes_helpers import (
    antihyperlipidemic_classes, implausible_free_t3, merge_drug_classes, seqns_on,
    thyroid_med_classes, thyroid_subsample_weight,
)

repo_root = pathlib.Path(__file__).resolve().parents[1]
raw_dir = repo_root / 'data' / 'raw'
processed_dir = repo_root / 'data' / 'processed'

cycles = {
    '2007-2008': {
        'thyroid_path': raw_dir / '2007-2008' / 'THYROD_E.xpt',
        'tchol_path': raw_dir / '2007-2008' / 'TCHOL_E.xpt',
        'demo_path': raw_dir / '2007-2008' / 'DEMO_E.xpt',
        'mcq_path': raw_dir / '2007-2008' / 'MCQ_E.xpt',
        'rxq_rx_path': raw_dir / '2007-2008' / 'RXQ_RX_E.xpt',
        'bmx_path': raw_dir / '2007-2008' / 'BMX_E.xpt',
    },
    '2009-2010': {
        'thyroid_path': raw_dir / '2009-2010' / 'THYROD_F.xpt',
        'tchol_path': raw_dir / '2009-2010' / 'TCHOL_F.xpt',
        'demo_path': raw_dir / '2009-2010' / 'DEMO_F.xpt',
        'mcq_path': raw_dir / '2009-2010' / 'MCQ_F.xpt',
        'rxq_rx_path': raw_dir / '2009-2010' / 'RXQ_RX_F.xpt',
        'bmx_path': raw_dir / '2009-2010' / 'BMX_F.xpt',
    },
    '2011-2012': {
        'thyroid_path': raw_dir / '2011-2012' / 'THYROD_G.xpt',
        'tchol_path': raw_dir / '2011-2012' / 'TCHOL_G.xpt',
        'demo_path': raw_dir / '2011-2012' / 'DEMO_G.xpt',
        'mcq_path': raw_dir / '2011-2012' / 'MCQ_G.xpt',
        'rxq_rx_path': raw_dir / '2011-2012' / 'RXQ_RX_G.xpt',
        'bmx_path': raw_dir / '2011-2012' / 'BMX_G.xpt',
    },
}


def medication_seqns(rxq_rx_df, rxq_drug_df):
    rx_merged = merge_drug_classes(rxq_rx_df, rxq_drug_df)
    thyroid_med = seqns_on(rx_merged, classes=thyroid_med_classes)
    lipid_lowering = seqns_on(rx_merged, classes=antihyperlipidemic_classes)
    return thyroid_med, lipid_lowering


def load_and_filter_cycle(cycle_label, paths, rxq_drug_df):
    thyroid_df, _ = pyreadstat.read_xport(str(paths['thyroid_path']))
    tchol_df, _ = pyreadstat.read_xport(str(paths['tchol_path']))
    demo_df, _ = pyreadstat.read_xport(str(paths['demo_path']))
    mcq_df, _ = pyreadstat.read_xport(str(paths['mcq_path']))
    rxq_rx_df, _ = pyreadstat.read_xport(str(paths['rxq_rx_path']))
    bmx_df, _ = pyreadstat.read_xport(str(paths['bmx_path']))

    merged_df = pd.merge(left=thyroid_df, right=tchol_df, on='SEQN', how='inner')
    merged_df = pd.merge(left=merged_df, right=demo_df, on='SEQN', how='left')

    # H4 study population: thyroid panel + total cholesterol measured + age >= 20
    # + not pregnant. Total cholesterol is non-fasting (same as H3), no
    # fasting-subsample filter here, and WTMEC2YR (general exam weight) is the
    # right weight, not WTSAF2YR.
    analytic_df = merged_df.dropna(subset=['LBXT3F', 'LBXTC'])
    analytic_df = analytic_df[analytic_df['RIDAGEYR'] >= 20]
    analytic_df = analytic_df[analytic_df['RIDEXPRG'] != 1]
    # Flagged rather than dropped so notebooks/09_robustness_checks.ipynb can
    # show results with the record retained; every other notebook drops it.
    analytic_df = analytic_df.assign(ft3_implausible=implausible_free_t3(analytic_df))

    analytic_df = pd.merge(left=analytic_df, right=mcq_df, on='SEQN', how='left')
    analytic_df = analytic_df[analytic_df['MCQ160M'] != 1]

    thyroid_med_seqns, lipid_lowering_seqns = medication_seqns(rxq_rx_df, rxq_drug_df)
    analytic_df = analytic_df[~analytic_df['SEQN'].isin(thyroid_med_seqns)]

    analytic_df = pd.merge(left=analytic_df, right=bmx_df[['SEQN', 'BMXBMI']], on='SEQN', how='left')

    # Lipid-lowering-med users are NOT dropped here, since H4 needs both a primary
    # model (excluding them) and a sensitivity check (keeping them, adjusted),
    # so both groups must survive into the saved sample. The flag lets the
    # modeling step decide which rows to use for which model.
    analytic_df['on_lipid_lowering'] = analytic_df['SEQN'].isin(lipid_lowering_seqns)

    # WTMEC2YR in 2007-2008, WTSA2YR (thyroid subsample) in 2009-2012.
    analytic_df['wt_thyroid_2yr'] = thyroid_subsample_weight(analytic_df, cycle_label, 'WTMEC2YR')

    return analytic_df


def build_h4_sample():
    rxq_drug_df, _ = pyreadstat.read_xport(str(raw_dir / 'RXQ_DRUG.xpt'))

    pooled_frames = []
    for cycle_label, paths in cycles.items():
        analytic_df = load_and_filter_cycle(cycle_label, paths, rxq_drug_df)
        n_lipid_lowering = analytic_df['on_lipid_lowering'].sum()
        print(f"{cycle_label}: n = {len(analytic_df)}, on lipid-lowering meds = {n_lipid_lowering}")
        pooled_frames.append(analytic_df.assign(cycle=cycle_label))

    pooled = pd.concat(pooled_frames, ignore_index=True)
    print(f"pooled n = {len(pooled)}")
    print(f"pooled n, primary model (excludes lipid-lowering meds) = {(~pooled['on_lipid_lowering']).sum()}")
    return pooled


if __name__ == '__main__':
    pooled = build_h4_sample()
    processed_dir.mkdir(parents=True, exist_ok=True)
    out_path = processed_dir / 'h4_analytic_sample.csv'
    pooled.to_csv(out_path, index=False)
    print(f"saved to {out_path}")
