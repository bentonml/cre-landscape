###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2021.11.14
#   updated   | 2023.03.22
#
###

import pandas as pd
import numpy  as np

from scipy import stats


### // constants and paths // ###
DATA_PATH = '../dat'           # relative to slurm dir
GEN_DATA_PATH = '../../../../data'
gene_to_ensembl = pd.read_table(f'{GEN_DATA_PATH}/dna/hg19/hgnc_to_ensembl.txt')
### \\


### // functions \\ ###
def sim(row):
    return float(len(row.peakachu_targets.intersection(row.abc_targets)) / len(row.peakachu_targets.union(row.abc_targets)))

def rel_sim(row):
    p = len(row.peakachu_targets)
    h = len(row.abc_targets)
    if p < h:
        relmax = p / len(row.peakachu_targets.union(row.abc_targets))
    else:
        relmax = h / len(row.peakachu_targets.union(row.abc_targets))
    if relmax == 0:
        return 0
    return float(len(row.peakachu_targets.intersection(row.abc_targets)) / len(row.peakachu_targets.union(row.abc_targets))) / relmax

def peakachu_len(row):
    return len(row.peakachu_targets)

def abc_len(row):
    return len(row.abc_targets)
### \\


### // compare peakachu loops to ABC pairs \\ ###
all_tis = []

for tis in ['liver', 'heart', 'muscle', 'ovary', 'pancreas', 'spleen']:
    print(f'{tis}')

    col_names=['chrom', 'start', 'end', 'gene_loop', 'gtex_tpm', 'hk', 'lof_intol', 'essential', 'rel_entropy', 'enh_rel_entropy', 'enh_length',
               'phastcons_overlap', 'tisspec_enh', 'tissue', 'chrom2', 'start2', 'end2', 'name', 'class', 'TargetGene', 'TargetGeneTSS',
               'TargetGeneIsExpressed', 'distance', 'ABC.Score', 'CellType', 'overlap']

    df = pd.read_table(f'{DATA_PATH}/{tis}_peakachuloop_abc_intersect_wo.bed', sep='\t', header=None, names=col_names)
    df = df.filter(['chrom', 'start', 'end', 'gene_loop', 'name', 'class', 'TargetGene', 'TargetGeneIsExpressed', 'ABC.Score'])
    df = df.merge(gene_to_ensembl, left_on='TargetGene', right_on='Approved symbol', how='left')
    df = df.drop(columns=['Approved symbol', 'HGNC ID', 'Approved name', 'Status'])
    df = df.rename(columns={'Ensembl ID(supplied by Ensembl)':'gene_abc'})

    print(f'CREs intersecting ABC enhancer-gene pair: {df.shape[0]}')
    print(f'Number of rows remove b/c no Ensembl mapping: {df.shape[0] - df.dropna().shape[0]}')

    # drop genes that cannot be mapped to ensembl IDs
    df = df.dropna()

    loop = df.groupby(['chrom', 'start', 'end']).gene_loop.agg(lambda x: set(x)).reset_index().rename(columns={'gene_loop':'peakachu_targets'})
    abc  = df.groupby(['chrom', 'start', 'end']).gene_abc.agg(lambda x: set(x)).reset_index().rename(columns={'gene_abc':'abc_targets'})

    loop_abc_merge_df = loop.merge(abc)
    loop_abc_merge_df['peakachu_len'] = loop_abc_merge_df.apply(peakachu_len, axis=1)
    loop_abc_merge_df['abc_len'] = loop_abc_merge_df.apply(abc_len, axis=1)
    loop_abc_merge_df['similarity'] = loop_abc_merge_df.apply(sim, axis=1)
    loop_abc_merge_df['rel_similarity'] = loop_abc_merge_df.apply(rel_sim, axis=1)

    print(loop_abc_merge_df.filter(['peakachu_len', 'abc_len', 'similarity', 'rel_similarity']).describe())
    all_tis.append(loop_abc_merge_df.filter(['peakachu_len', 'abc_len', 'similarity', 'rel_similarity']).assign(tissue=tis))

all_df = pd.concat(all_tis)
all_df.to_csv(f'../res/all_tis_stats.tsv', sep='\t', header=True, index=False)

