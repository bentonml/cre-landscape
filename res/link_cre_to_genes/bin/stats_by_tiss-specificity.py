###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#
###

import os
import pandas as pd
import numpy  as np

from pybedtools import BedTool
from datetime import date

import seaborn as sns
import matplotlib.pyplot as plt

### // constants and paths // ###
EXP_DATA_PATH = f'../dat/2022-01-07'
RES_PATH = '..'
FIG_PATH = f'{RES_PATH}/fig/{str(date.today())}'

# create a date stamped dir for files
if not os.path.isdir(FIG_PATH):
    os.makedirs(FIG_PATH)

thresh='Q05'          # for hic q-values
t = 0.05              # for hic q-values

# plotting options
fmt='pdf'
rc = {'font.size':14, 'axes.titlesize':16,'axes.labelsize':14, 'legend.fontsize': 14,
      'xtick.labelsize': 14, 'ytick.labelsize': 14}

tis_order = ['spleen','liver', 'heart_left_ventricle', 'brain_hippocampus', 'lung', 'pancreas',
             'brain_prefrontal_cortex', 'psoas_muscle', 'small_intestine', 'ovary']

def facet_titles(axes):
    ax = axes.flatten()
    ax[0].set_title('Spleen')
    ax[1].set_title('Liver')
    ax[2].set_title('Heart')
    ax[3].set_title('Hippocampus')
    ax[4].set_title('Lung')
    ax[5].set_title('Pancreas')
    ax[6].set_title('Prefrontal cortex')
    ax[7].set_title('Muscle')
    ax[8].set_title('Small intestine')
    ax[9].set_title('Ovary')

def facet_x_axis(axes, title):
    for ax in axes.flatten():
        ax.set_xlabel(title)

def read_data(landscape_def):
    df_lst = []
    for tis in tis_order:
        if landscape_def == 'contact':
            infile = f'{EXP_DATA_PATH}/{tis}_gene_with_enh-count_frac-enh-spec_hicQ05-linked.tsv'
        elif landscape_def == 'loop':
            infile = f'{EXP_DATA_PATH}/{tis}_gene_with_enh-count_frac-enh-spec_peakachuloop-linked.tsv'
        enh_num_by_gene = pd.read_table(infile)
        enh_num_by_gene = enh_num_by_gene.rename(columns={tis:'gtex_exp'})
        enh_num_by_gene = enh_num_by_gene.assign(tissue=tis)
        df_lst.append(enh_num_by_gene)
    all_tis = pd.concat(df_lst)
    return all_tis

for landscape_def in ['loop', 'contact']:
    # create dataframe of enhancer number by gene with all tissues
    all_tis = read_data(landscape_def)
    all_tis_exp = all_tis[all_tis['exp']==1].copy(deep=True)

    # plot rel_entropy vs. housekeeping gene
    g = sns.catplot(x='hk', y='rel_entropy', col='tissue', data=all_tis_exp,
                order=[0, 1], palette=['tab:gray', 'tab:green'], col_order=tis_order,
                kind='box', col_wrap=5, showfliers=True, fliersize=1)
    for ax in g.axes.flatten():
        ax.axhline(y=0.6, c='.3', linestyle='--')
        ax.axhline(y=0.4, c='.3', linestyle='--')
        ax.axhline(y=0.2, c='.3', linestyle='--')
    facet_titles(g.axes)
    facet_x_axis(g.axes, '')
    plt.savefig(f'{FIG_PATH}/all_{landscape_def}_rel-entropy_bytissue_v_hkstatus_boxplot.{fmt}', format=fmt, dpi=400)
    plt.close()

    # calculate different thresholds
    all_tis_exp = all_tis_exp.assign(ts_06=lambda x: np.where(x.rel_entropy > 0.6, 1, 0),
                                     ts_05=lambda x: np.where(x.rel_entropy > 0.5, 1, 0),
                                     ts_04=lambda x: np.where(x.rel_entropy > 0.4, 1, 0),
                                     ts_03=lambda x: np.where(x.rel_entropy > 0.3, 1, 0),
                                     ts_02=lambda x: np.where(x.rel_entropy > 0.2, 1, 0))

    df06 = (all_tis_exp.groupby(['tissue', 'ts_06'], as_index=False).target_gene.count()
                .pivot(index='tissue', columns='ts_06', values='target_gene')
                .rename(columns={0:'ts06_broad', 1:'ts06_spec'})
                .assign(prop_tis06=lambda x: x.ts06_spec / x.ts06_broad))

    df05 = (all_tis_exp.groupby(['tissue', 'ts_05'], as_index=False).target_gene.count()
                .pivot(index='tissue', columns='ts_05', values='target_gene')
                .rename(columns={0:'ts05_broad', 1:'ts05_spec'})
                .assign(prop_tis05=lambda x: x.ts05_spec / x.ts05_broad))

    df04 = (all_tis_exp.groupby(['tissue', 'ts_04'], as_index=False).target_gene.count()
                .pivot(index='tissue', columns='ts_04', values='target_gene')
                .rename(columns={0:'ts04_broad', 1:'ts04_spec'})
                .assign(prop_tis04=lambda x: x.ts04_spec / x.ts04_broad))

    df03 = (all_tis_exp.groupby(['tissue', 'ts_03'], as_index=False).target_gene.count()
                .pivot(index='tissue', columns='ts_03', values='target_gene')
                .rename(columns={0:'ts03_broad', 1:'ts03_spec'})
                .assign(prop_tis03=lambda x: x.ts03_spec / x.ts03_broad))

    df02 = (all_tis_exp.groupby(['tissue', 'ts_02'], as_index=False).target_gene.count()
                .pivot(index='tissue', columns='ts_02', values='target_gene')
                .rename(columns={0:'ts02_broad', 1:'ts02_spec'})
                .assign(prop_tis02=lambda x: x.ts02_spec / x.ts02_broad))

    full_df = pd.concat([df02, df03, df04, df05, df06], axis=1)
    full_df.to_csv(f'{RES_PATH}/tissue_specific_stats_{landscape_def}_{str(date.today())}.txt',
                    sep='\t', index=True, header=True)
