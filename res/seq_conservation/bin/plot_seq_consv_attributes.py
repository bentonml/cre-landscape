###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2021.04.28
#
#   assumes pwd is here with hydrogen:
#       /dors/capra_lab/users/bentonml/cross_species_gain_loss/bin/slurm/
###

import os
import pandas as pd
import numpy  as np
from datetime import date
from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm


### // constants, paths, functions \\ ###
CRE_PATH = '../../link_cre_to_genes/dat/2022-05-10'
LIN_PATH = '../dat/linsight'
DAT_PATH = '../dat/'
FIG_PATH = f'../fig/{str(date.today())}'

# create a date stamped dir for files
if not os.path.isdir(FIG_PATH):
        os.makedirs(FIG_PATH)

fmt = 'pdf'
rc = {'font.size':14, 'axes.titlesize':16,'axes.labelsize':14, 'legend.fontsize': 14,
      'xtick.labelsize': 14, 'ytick.labelsize': 14, 'font.sans-serif':'Arial', 'font.family':'sans-serif'}

tis_order = ['spleen','liver', 'heart_left_ventricle', 'brain_hippocampus', 'lung', 'pancreas',
             'brain_prefrontal_cortex', 'psoas_muscle', 'small_intestine', 'ovary']
tis_name = {'ovary':'Ovary', 'psoas_muscle': 'Muscle', 'heart_left_ventricle':'Heart',
            'lung':'Lung', 'spleen':'Spleen', 'small_intestine':'Small intestine',
            'liver':'Liver', 'pancreas':'Pancreas', 'brain_prefrontal_cortex':'Prefrontal cortex',
            'brain_hippocampus':'Hippocampus'}

linsight_cols = ['chrom', 'start', 'end', 'landscape_cre_num', 'landscape_tisspec_cre_num',
                 'landscape_cre_quartile', 'landscape_cre_quartile_id',
                 'chr_lin', 'start_lin', 'end_lin', 'linsight_score', 'overlap']
### \\


for lndscp_type in ['peakachuloop', 'hicQ05']:
    print(lndscp_type)
    spearman_lin, pearson_lin = [], []
    spearman_pc, pearson_pc = [], []

    for tis in tis_order:
        # read in landscape attr per enh
        lndscp_size = (pd.read_table(f'{CRE_PATH}/{tis}_enh_to_gene_{lndscp_type}-linked.tsv')
                         .filter(['enh_chrom', 'enh_start', 'enh_end', 'target_gene'])
                         .drop_duplicates()
                         .groupby('target_gene', as_index=False)
                         .agg({'enh_chrom':'count'})
                         .rename(columns={'enh_chrom':'num_cre'}))

        cre_df = (pd.read_table(f'{CRE_PATH}/{tis}_enh_to_gene_{lndscp_type}-linked.tsv')
                    .assign(phastcons_prop=lambda x: x.phastcons_overlap / x.enh_length)
                    .filter(['enh_chrom', 'enh_start', 'enh_end', 'target_gene', 'gtex_tpm', 'enh_rel_entropy', 'rel_entropy', 'phastcons_prop'])
                    .drop_duplicates()
                    .merge(lndscp_size, validate='m:1')
                 )

        lin = (pd.read_table(f'{LIN_PATH}/{tis}_{lndscp_type}_cre_x_linsight.bed', names=linsight_cols, na_values='.')
                 .dropna()
                 .groupby(['chrom', 'start', 'end'], as_index=False)
                 .agg({'linsight_score':'mean', 'overlap':'sum'})
                 .merge(cre_df, left_on=['chrom', 'start', 'end'], right_on=['enh_chrom', 'enh_start', 'enh_end'], validate='m:m')
              )

        lin_grp = lin.groupby(['target_gene'], as_index=False).agg({'linsight_score':['mean', 'median', 'max', 'std'], 'num_cre':'mean', 'overlap':'sum', 'phastcons_prop':['mean', 'median']})
        lin_grp.columns = ['_'.join(col).strip('_') for col in lin_grp.columns.values]

        edges = pd.qcut(lin_grp.num_cre_mean, q=4, retbins=True)[1]
        labels = [f'{edges[0]:.0f}-{edges[1]:.0f}', f'{edges[1]+1:.0f}-{edges[2]:.0f}', f'{edges[2]+1:.0f}-{edges[3]:.0f}', f'{edges[3]+1:.0f}+']
        lin_grp['landscape_cre_quartile'] = pd.qcut(lin_grp.num_cre_mean, q=4, labels=labels)

        with sns.plotting_context("paper", rc=rc):
            fig, axes = plt.subplots(1, 2, figsize=(7.5, 3), sharey=False)
            sns.despine(fig=fig)
            sns.violinplot(x='landscape_cre_quartile', y='linsight_score_mean', data=lin_grp, palette='Blues', ax=axes[0], cut=0)
            axes[0].set(xlabel='# CREs', ylabel='Mean LINSIGHT Score')
            sns.violinplot(x='landscape_cre_quartile', y='phastcons_prop_mean', data=lin_grp, palette='Blues', ax=axes[1], cut=0)
            axes[1].set(xlabel='# CREs', ylabel='Mean % PhastCons')
            plt.tight_layout()
            plt.savefig(f'{FIG_PATH}/{str(date.today())}_{tis}_{lndscp_type}_linsight-mean_prop-phastcons_v_cre-quartile_violinplot.{fmt}', format=fmt, dpi=400)
            plt.close()

        print(tis)
        print(f"linsight\t{lin_grp.groupby('landscape_cre_quartile').linsight_score_mean.median()}")
        print(f"phastcons\t{lin_grp.groupby('landscape_cre_quartile').phastcons_prop_mean.median()}")

        corr, p = stats.spearmanr(lin_grp.linsight_score_mean, lin_grp.num_cre_mean)
        spearman_lin.append(f'{tis}\tmean\t{corr}\t{p}')
        corr, p = stats.spearmanr(lin_grp.phastcons_prop_mean, lin_grp.num_cre_mean)
        spearman_pc.append(f'{tis}\tmean\t{corr}\t{p}')

        corr, p = stats.pearsonr(lin_grp.linsight_score_mean, lin_grp.num_cre_mean)
        pearson_lin.append(f'{tis}\tmean\t{corr}\t{p}')
        corr, p = stats.pearsonr(lin_grp.phastcons_prop_mean, lin_grp.num_cre_mean)
        pearson_pc.append(f'{tis}\tmean\t{corr}\t{p}')

    print('spearman linsight')
    for line in spearman_lin:
        print(line)
    print('pearson linsight')
    for line in pearson_lin:
        print(line)

    print('spearman phastcons')
    for line in spearman_pc:
        print(line)

    print('spearman phastcons')
    for line in pearson_pc:
        print(line)
