import pathlib

import requests

repo_root = pathlib.Path(__file__).resolve().parents[1]
raw_dir = repo_root / 'data' / 'raw'

base_url = 'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public'

# Cycle label -> (first year of cycle used in CDC's URL, file-name suffix)
cycles = {
    '2007-2008': ('2007', 'E'),
    '2009-2010': ('2009', 'F'),
    '2011-2012': ('2011', 'G'),
}

# Every component file read by pipeline/02-06. CRP was only collected through
# 2009-2010, so it is skipped for 2011-2012.
components = ['DEMO', 'MCQ', 'RXQ_RX', 'BMX', 'BPX', 'THYROD', 'TRIGLY', 'HDL', 'GLU', 'TCHOL', 'CRP']
crp_cycles = ['2007-2008', '2009-2010']


def download(url, out_path):
    if out_path.exists():
        print(f"exists, skipping: {out_path.relative_to(repo_root)}")
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    # CDC serves an HTML error page with status 200 for some missing files,
    # so check for the SAS transport header before saving.
    if not response.content.startswith(b'HEADER RECORD'):
        raise ValueError(f"not an XPT file: {url}")
    out_path.write_bytes(response.content)
    print(f"downloaded: {out_path.relative_to(repo_root)}")


def download_all():
    for cycle_label, (year, suffix) in cycles.items():
        for component in components:
            if component == 'CRP' and cycle_label not in crp_cycles:
                continue
            file_name = f'{component}_{suffix}.xpt'
            download(f'{base_url}/{year}/DataFiles/{file_name}', raw_dir / cycle_label / file_name)

    # The prescription drug lookup table is a single cross-cycle file that CDC
    # files under its 1988 directory.
    download(f'{base_url}/1988/DataFiles/RXQ_DRUG.xpt', raw_dir / 'RXQ_DRUG.xpt')


if __name__ == '__main__':
    download_all()
