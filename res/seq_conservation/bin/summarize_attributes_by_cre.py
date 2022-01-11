###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2021.02.23
#
#   assumes pwd is here with hydrogen:
#       /dors/capra_lab/users/bentonml/cross_species_gain_loss/bin/slurm/
###


import pandas as pd
import numpy  as np
from datetime import date
from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm


### // constants, paths, functions \\ ###
CRE_PATH = '../../link_cre_to_genes/dat/2022-01-07'
DAT_PATH = '../dat'
FIG_PATH = '../fig'

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
    print(f'tissue\ttest\trho\tp')

    for tis in tis_order:
        ldf = (pd.read_table(f'{DAT_PATH}/linsight/{tis}_{lndscp_type}_cre_x_linsight.bed', names=linsight_cols, na_values='.')
                .drop(columns=['chr_lin', 'start_lin', 'end_lin'])
                .fillna(0)
                .groupby(['chrom', 'start', 'end', 'landscape_cre_num'], as_index=False)
                .agg({'landscape_cre_quartile_id':'mean', 'linsight_score':['mean', 'median', 'max'], 'overlap':'sum'}))

        ldf.columns = ['_'.join(col).strip('_') for col in ldf.columns.values]

        # kept the enh associated with the smallest landscape
        ldf_min = ldf.sort_values('landscape_cre_num', ascending=True).drop_duplicates(subset=['chrom', 'start', 'end'], keep='first')
        ldf_min['landscape_cre_decile'] = pd.qcut(ldf_min.landscape_cre_num, q=10)

        # read in landscape attr per enh
        cre_df = pd.read_table(f'{CRE_PATH}/{tis}_enh_to_gene_peakachuloop-linked.tsv')
        landscape_count = (cre_df.groupby(['target_gene'], as_index=False).enh_start.count()).rename(columns={'enh_start':'cre_num'})
        cre_df = cre_df.merge(landscape_count, how='left', validate='m:1')

        # generate final dataframe
        df = cre_df.filter(['enh_chrom', 'enh_start', 'enh_end', 'cre_num', 'gtex_tpm', 'phastcons_overlap', 'enh_length', 'enh_rel_entropy']).assign(phastcons_prop=lambda x: x.phastcons_overlap / x.enh_length)
        ddf = df.merge(ldf_min.filter(['chrom', 'start', 'end', 'landscape_cre_num', 'linsight_score_mean']), left_on=['enh_chrom', 'enh_start', 'enh_end'], right_on=['chrom', 'start', 'end']).drop_duplicates()
        ddf['cre_num_quartile'] = pd.qcut(ddf.cre_num, q=4, labels=[1,2,3,4])

        # kruskal wallis with linsight mean
        ts, p = stats.kruskal(*[group['linsight_score_mean'].values for name, group in ldf.groupby('landscape_cre_quartile_id_mean')])
        means = ldf.groupby('landscape_cre_quartile_id_mean').linsight_score_mean.mean()

        with sns.plotting_context("paper", rc=rc):
            fig, axes = plt.subplots(1, 2, figsize=(8, 4), sharey=False)
            fig.suptitle(f'{tis_name[tis]}')
            sns.despine(fig=fig)
            sns.violinplot(x='cre_num_quartile', y='linsight_score_mean', data=ddf, palette='Blues', ax=axes[0], cut=0)
            axes[0].set(xlabel='CRE landscape quartile', ylabel='Mean LINSIGHT Score')
            sns.violinplot(x='cre_num_quartile', y='phastcons_prop', data=ddf, palette='Blues', ax=axes[1], cut=0)
            axes[1].set(xlabel='CRE landscape quartile', ylabel='% conserved bp')
            plt.tight_layout()
            plt.savefig(f'{FIG_PATH}/{str(date.today())}_{tis}_{lndscp_type}_linsight-mean_prop-phastcons_v_cre-quartile_violinplot.{fmt}', format=fmt, dpi=400)
            plt.close()

        rho, p = stats.spearmanr(ddf.cre_num, ddf.linsight_score_mean)
        print(f'{tis}\t# CRE v. Mean LINSIGHT\t{rho}\t{p}')
        rho, p = stats.spearmanr(ddf.cre_num, ddf.phastcons_prop)
        print(f'{tis}\t# CRE v. % PhastCons\t{rho}\t{p}')
        print(f'{tis_name[tis]}\t{means[1]}\t{means[2]}\t{means[3]}\t{means[4]}\t{ts}\t{p}\t{p*20}\n')

