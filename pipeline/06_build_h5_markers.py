import pathlib

import pandas as pd
import pyreadstat

from nhanes_helpers import antihyperlipidemic_classes, beta_blocker_classes, beta_blocker_subclasses, merge_drug_classes, seqns_on

repo_root = pathlib.Path(__file__).resolve().parents[1]
raw_dir = repo_root / 'data' / 'raw'
processed_dir = repo_root / 'data' / 'processed'

# H5 (exploratory, post-hoc): pulse rate and total cholesterol as
# community-preferred "functional" thyroid markers, tested against the
# same outcomes/populations/weights as H1-H3. Not a replacement for the
# primary hypotheses -- an addition, kept in its own multiplicity family.
cycles = {
    '2007-2008': {
        'bpx_path': raw_dir / '2007-2008' / 'BPX_E.xpt',
        'tchol_path': raw_dir / '2007-2008' / 'TCHOL_E.xpt',
        'rxq_rx_path': raw_dir / '2007-2008' / 'RXQ_RX_E.xpt',
    },
    '2009-2010': {
        'bpx_path': raw_dir / '2009-2010' / 'BPX_F.xpt',
        'tchol_path': raw_dir / '2009-2010' / 'TCHOL_F.xpt',
        'rxq_rx_path': raw_dir / '2009-2010' / 'RXQ_RX_F.xpt',
    },
    '2011-2012': {
        'bpx_path': raw_dir / '2011-2012' / 'BPX_G.xpt',
        'tchol_path': raw_dir / '2011-2012' / 'TCHOL_G.xpt',
        'rxq_rx_path': raw_dir / '2011-2012' / 'RXQ_RX_G.xpt',
    },
}


def build_cycle_markers(paths, rxq_drug_df):
    bpx_df, _ = pyreadstat.read_xport(str(paths['bpx_path']))
    tchol_df, _ = pyreadstat.read_xport(str(paths['tchol_path']))
    rxq_rx_df, _ = pyreadstat.read_xport(str(paths['rxq_rx_path']))

    markers_df = pd.merge(left=bpx_df[['SEQN', 'BPXPLS']], right=tchol_df[['SEQN', 'LBXTC']], on='SEQN', how='outer')

    rx_merged = merge_drug_classes(rxq_rx_df, rxq_drug_df)
    beta_blocker_users = seqns_on(rx_merged, classes=beta_blocker_classes, subclasses=beta_blocker_subclasses)
    markers_df['on_beta_blocker'] = markers_df['SEQN'].isin(beta_blocker_users)

    # Cholesterol is only meaningful as a thyroid marker when it is not being
    # lowered pharmacologically, so cholesterol models exclude anyone on any
    # antihyperlipidemic agent (the same broad class H4 uses).
    antihyperlipidemic_users = seqns_on(rx_merged, classes=antihyperlipidemic_classes)
    markers_df['on_antihyperlipidemic'] = markers_df['SEQN'].isin(antihyperlipidemic_users)

    return markers_df


def build_h5_markers():
    rxq_drug_df, _ = pyreadstat.read_xport(str(raw_dir / 'RXQ_DRUG.xpt'))

    pooled_frames = []
    for cycle_label, paths in cycles.items():
        markers_df = build_cycle_markers(paths, rxq_drug_df)
        print(f"{cycle_label}: n = {len(markers_df)}, on beta-blocker = {markers_df['on_beta_blocker'].sum()}")
        pooled_frames.append(markers_df.assign(cycle=cycle_label))

    pooled = pd.concat(pooled_frames, ignore_index=True)
    print(f"pooled n = {len(pooled)}")
    return pooled


if __name__ == '__main__':
    pooled = build_h5_markers()
    processed_dir.mkdir(parents=True, exist_ok=True)
    out_path = processed_dir / 'h5_markers.csv'
    pooled.to_csv(out_path, index=False)
    print(f"saved to {out_path}")
