import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

fig, axes = plt.subplots(2, 2, figsize=(15, 12))

df = pd.concat([
    pd.read_csv(f)
    for f in snakemake.input['csvs']
])
df['sample'] = df.apply(lambda row: '_'.join(row[['antibody', 'condition']]), axis=1)


def _bed_lengths(f, meta):
    bed = pd.read_csv(f, header=None, sep='\t', comment='t')
    header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'] # , 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
    bed.columns = header[:len(df.columns)]
    return pd.DataFrame({
        'width': bed['chromEnd'] - bed['chromStart'],
        'sample': '_'.join(meta[['antibody', 'condition']]),
        'type': meta['type']
    })

peak_widths = pd.concat([_bed_lengths(f, df.iloc[0]) for i, f in enumerate(snakemake.input['beds'])])

# number of peaks
axes[0, 0].set_title('#(Peaks)')
sns.swarmplot(data=df, x='sample', hue='type', y='total_peaks', ax=axes[0, 0])
axes[0, 0].set_xticklabels(axes[0,0].get_xticklabels(), rotation=30, ha='right')

# width of peaks
axes[0, 1].set_title('Width of peaks')
sns.violinplot(data=peak_widths, x='sample', y='width', hue='type')

# reproduced peaks
axes[1, 0].set_title('Reproducibility of peaks across replicates')
sns.swarmplot(data=df, x='sample', hue='type', y='repl_percent', ax=axes[1, 0])

# frips
axes[1, 1].set_title('Percentage of fragments in peaks')
sns.swarmplot(data=df, x='sample', hue='type', y='frips', ax=axes[1, 1])
axes[1, 1].set_xticklabels(axes[1, 1].get_xticklabels(), rotation=30)

plt.tight_layout()

fig.savefig(snakemake.output['png'])
fig.savefig(snakemake.output['svg'])
