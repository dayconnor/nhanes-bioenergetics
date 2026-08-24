import pathlib

import pandas as pd
import pyreadstat

repo_root = pathlib.Path(__file__).resolve().parents[1]
raw_dir = repo_root / 'data' / 'raw'
processed_dir = repo_root / 'data' / 'processed'

thyroid_med_categories = ['THYROID HORMONES', 'ANTITHYROID AGENTS']

# H3 only has CRP in 2007-2008 and 2009-2010.
cycles = {
    '2007-2008': {
        'thyroid_path': raw_dir / '2007-2008' / 'THYROD_E.xpt',
        'crp_path': raw_dir / '2007-2008' / 'CRP_E.xpt',
        'demo_path': raw_dir / '2007-2008' / 'DEMO_E.xpt',
        'mcq_path': raw_dir / '2007-2008' / 'MCQ_E.xpt',
        'rxq_rx_path': raw_dir / '2007-2008' / 'RXQ_RX_E.xpt',
        'bmx_path': raw_dir / '2007-2008' / 'BMX_E.xpt',
    },
    '2009-2010': {
        'thyroid_path': raw_dir / '2009-2010' / 'THYROD_F.xpt',
        'crp_path': raw_dir / '2009-2010' / 'CRP_F.xpt',
        'demo_path': raw_dir / '2009-2010' / 'DEMO_F.xpt',
        'mcq_path': raw_dir / '2009-2010' / 'MCQ_F.xpt',
        'rxq_rx_path': raw_dir / '2009-2010' / 'RXQ_RX_F.xpt',
        'bmx_path': raw_dir / '2009-2010' / 'BMX_F.xpt',
    },
}


def load_and_filter_cycle(paths, rxq_drug_df):
    thyroid_df, _ = pyreadstat.read_xport(str(paths['thyroid_path']))
    crp_df, _ = pyreadstat.read_xport(str(paths['crp_path']))
    demo_df, _ = pyreadstat.read_xport(str(paths['demo_path']))
    mcq_df, _ = pyreadstat.read_xport(str(paths['mcq_path']))
    rxq_rx_df, _ = pyreadstat.read_xport(str(paths['rxq_rx_path']))
    bmx_df, _ = pyreadstat.read_xport(str(paths['bmx_path']))

    merged_df = pd.merge(left=thyroid_df, right=crp_df, on='SEQN', how='inner')
    merged_df = pd.merge(left=merged_df, right=demo_df, on='SEQN', how='left')

    # H3 study population: thyroid panel + CRP measured + age >= 20 + not pregnant.
    # No fasting-subsample filter, CRP is measured on the full MEC exam
    # sample, not restricted to the fasting subsample like H1/H2's triglycerides
    # and glucose were, so WTMEC2YR (general exam weight) is the right weight,
    # not WTSAF2YR.
    analytic_df = merged_df.dropna(subset=['LBXT3F', 'LBXCRP'])
    analytic_df = analytic_df[analytic_df['RIDAGEYR'] >= 20]
    analytic_df = analytic_df[analytic_df['RIDEXPRG'] != 1]

    analytic_df = pd.merge(left=analytic_df, right=mcq_df, on='SEQN', how='left')
    analytic_df = analytic_df[analytic_df['MCQ160M'] != 1]

    rxq_merge = pd.merge(left=rxq_rx_df, right=rxq_drug_df, on='RXDDRGID', how='left')
    thyroid_med_seqns = rxq_merge[rxq_merge['RXDDCN1B'].isin(thyroid_med_categories)]['SEQN'].unique()
    analytic_df = analytic_df[~analytic_df['SEQN'].isin(thyroid_med_seqns)]

    analytic_df = pd.merge(left=analytic_df, right=bmx_df[['SEQN', 'BMXBMI']], on='SEQN', how='left')

    return analytic_df


def build_h3_sample():
    rxq_drug_df, _ = pyreadstat.read_xport(str(raw_dir / 'RXQ_DRUG.xpt'))

    pooled_frames = []
    for cycle_label, paths in cycles.items():
        analytic_df = load_and_filter_cycle(paths, rxq_drug_df)
        print(f"{cycle_label}: n = {len(analytic_df)}")
        pooled_frames.append(analytic_df.assign(cycle=cycle_label))

    pooled = pd.concat(pooled_frames, ignore_index=True)
    print(f"pooled n = {len(pooled)}")
    return pooled


if __name__ == '__main__':
    pooled = build_h3_sample()
    processed_dir.mkdir(parents=True, exist_ok=True)
    out_path = processed_dir / 'h3_analytic_sample.csv'
    pooled.to_csv(out_path, index=False)
    print(f"saved to {out_path}")
