import typer

import pandas as pd
from pathlib import Path
from pandarallel import pandarallel


def main(r1_paf: Path, r2_paf: Path, probe_fa: Path, count_file: Path, threads: int = 4) -> None:
    if count_file.is_file():
        probe_df = pd.read_csv(count_file, sep="\t")
    else:
        probe_list = [each.strip().replace('>', '') for each in open(probe_fa) if each.startswith('>')]
        probe_df = pd.DataFrame(probe_list, columns=['probe'])
        pandarallel.initialize(progress_bar=True, nb_workers=threads)
        r1_df = pd.read_csv(r1_paf, usecols=[0, 5], header=None, names=['read_name', 'probe'], sep="\t")
        r2_df = pd.read_csv(r2_paf, usecols=[0, 5], header=None, names=['read_name', 'probe'], sep='\t')
        r1_df['read_prefix'] = r1_df['read_name'].parallel_map(lambda x: x.split('/')[0])
        r2_df['read_prefix'] = r2_df['read_name'].parallel_map(lambda x: x.split('/')[0])
        merged_df = pd.concat([r1_df, r2_df])
        merged_df.drop_duplicates(subset=['probe', 'read_prefix'])
        probe_count = pd.DataFrame(merged_df['probe'].value_counts()).reset_index()
        probe_df = probe_df.merge(probe_count, how='left')
        probe_df.fillna(0, inplace=True)
        probe_df['count'] = probe_df['count'].astype('int')
        probe_df.to_csv(count_file, sep="\t", index=False)
    plot = probe_df.plot.hist(column = ["count"], bins=50)
    fig = plot.get_figure()
    plot_path = count_file.with_suffix('.hist.pdf')
    fig.savefig(str(plot_path))

if __name__ == "__main__":
    typer.run(main)