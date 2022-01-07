###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2021.11.12
#
#   plot gene sets matched on covariates (expression level)
###

import os
import pandas as pd
import numpy  as np
import seaborn as sns
import matplotlib.pyplot as plt

from datetime import date
from scipy import stats

plt.switch_backend('agg')  # add to save plots non-interactively

### // constants and paths \\ ###
DORS = '/dors/capra_lab/users/bentonml/cre_landscape'
#TODO: update paths
EXP_DAT_PATH = f'{DORS}/res/link_enh_to_genes/dat/2021-07-30'
CV_FILE_PATH = f'{DORS}/res/exp_variation/dat'
RES_PATH = f'{DORS}/res/exp_variation/fig/{str(date.today())}'


tis_order = ['spleen','liver', 'heart_left_ventricle', 'brain_hippocampus', 'lung', 'pancreas',
             'brain_prefrontal_cortex', 'psoas_muscle', 'small_intestine', 'ovary']
tis_names = ['Spleen','Liver', 'Heart', 'Hippocampus', 'Lung', 'Pancreas',
             'Prefrontal cortex', 'Muscle', 'Small intestine', 'Ovary']
tis_dict  = {'ovary':'Ovary', 'psoas_muscle':'Muscle', 'heart_left_ventricle':'Heart',
             'lung':'Lung', 'spleen':'Spleen', 'small_intestine':'Small intestine',
             'pancreas':'Pancreas', 'liver':'Liver', 'brain_prefrontal_cortex':'Prefrontal cortex',
             'brain_hippocampus':'Hippocampus'}

# plotting options
fmt='pdf'
rc = {'font.size':14, 'axes.titlesize':16,'axes.labelsize':14, 'legend.fontsize': 14,
      'xtick.labelsize': 14, 'ytick.labelsize': 14}

def read_data(landscape_def):
    df_lst = []
    for tis in tis_order:
        if landscape_def == 'contact':
            infile = f'{EXP_DAT_PATH}/{tis}_gene_with_enh-count_frac-enh-spec_hicQ05-linked.tsv'
        elif landscape_def == 'loop':
            infile = f'{EXP_DAT_PATH}/{tis}_gene_with_enh-count_frac-enh-spec_peakachuloop-linked.tsv'
        enh_num_by_gene = pd.read_table(infile)
        enh_num_by_gene = enh_num_by_gene.assign(rel_entropy_bins=lambda x: pd.cut(x['rel_entropy'], bins=10))
        enh_num_by_gene = enh_num_by_gene.assign(enh_rel_entropy_bins=lambda x: pd.cut(x['enh_rel_entropy'], bins=5))
        enh_num_by_gene = enh_num_by_gene.assign(exp_nocat=lambda x: np.where((x.exp == 1) & (x.hk == 0) & (x.lof_intol == 0) & (x.essential == 0), 1, 0))
        enh_num_by_gene = enh_num_by_gene.assign(exp_broad=lambda x: np.where((x.exp == 1) & (x.tis_spec == 0), 1, 0))
        enh_num_by_gene = enh_num_by_gene.assign(log2_exp=lambda x: np.log2(x[tis] + 1))
        enh_num_by_gene = enh_num_by_gene.rename(columns={tis:'gtex_exp'})
        enh_num_by_gene = enh_num_by_gene.assign(tissue=tis)
        df_lst.append(enh_num_by_gene)

    all_tis = pd.concat(df_lst)
    tisspec_bins = ['(0,0.1]', '(0.1,0.2]', '(0.2,0.3]', '(0.3,0.4]', '(0.4,0.5]',
                    '(0.5,0.6]', '(0.6,0.7]', '(0.7,0.8]', '(0.8,0.9]', '(0.9,1.0]']
    all_tis['frac_tisspec_enh_bins'] = pd.cut(all_tis.frac_tisspec_enh, bins=10,
                                              labels=tisspec_bins)
    return all_tis

def read_cv_data(tis, thresh='0.8'):
    return pd.read_table(f'{CV_FILE_PATH}/{tis}_loess_cv_resid.csv', sep=',')

#TODO: update to add hic contact
for landscape_def in ['loop']:
    print(landscape_def)

    for tis in tis_order:
        df = read_cv_data(tis)

        # create dataframe of enhancer number by gene with all tissues
        all_tis = read_data(landscape_def)
        all_tis_exp = all_tis[all_tis['exp']==1]

        df_merge = df.merge(all_tis_exp, left_on='Name', right_on='target_gene', how='inner')
        df_merge = df_merge.assign(enh_num_quartile=lambda x: pd.qcut(x['enh_num'], q=4, labels=[1,2,3,4]))
        df_merge = df_merge.assign(log2_tpm_quartile=lambda x: pd.qcut(x['median_log2_tpm'], q=4, labels=[1,2,3,4]))

        with sns.plotting_context("paper", rc=rc):
            g = sns.lmplot(x='median_log2_tpm', y='log2cv', data=df_merge, lowess=True, scatter_kws={'alpha':0.3, 'color':'.3'}, line_kws={'color':'tab:red'})
            g.set_xlabels('Median log2 Expression')
            g.set_ylabels('CV')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_def}_log2exp_v_cv_lmplot.{fmt}', format=fmt, dpi=400)
            plt.close()

            g = sns.lmplot(x='median_log2_tpm', y='loess_resid', data=df_merge, lowess=True, scatter_kws={'alpha':0.3, 'color':'.3'}, line_kws={'color':'tab:red'})
            g.set_xlabels('Median log2 Expression')
            g.set_ylabels('Expression Variation')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_def}_log2exp_v_expVar_lmplot.{fmt}', format=fmt, dpi=400)
            plt.close()

        rho, p = stats.spearmanr(df_merge.enh_num, df_merge.loess_resid)
        print(f'{tis}\t# CRE\t{rho}\t{p}')
        rho, p = stats.spearmanr(df_merge.rel_entropy, df_merge.loess_resid)
        print(f'{tis}\tRelative entropy\t{rho}\t{p}')
        rho, p = stats.spearmanr(df_merge.frac_tisspec_enh, df_merge.loess_resid)
        print(f'{tis}\t% tissue-specific CRE\t{rho}\t{p}')

        gt0 = df_merge.filter(['frac_phastcons', 'loess_resid']).dropna()
        rho, p = stats.spearmanr(gt0.frac_phastcons, gt0.loess_resid)
        print(f'{tis}\t% PhastCons\t{rho}\t{p}')

        with sns.plotting_context("paper", rc=rc):
            df_merge = df_merge.rename(columns={'loess_resid':'Expression variation', 'rel_entropy':'Tissue-specificity (gene)',
                                        'enh_num':'# CREs', 'frac_tisspec_enh':'% tissue-specific CREs', 'frac_phastcons':'% PhastCons'})
            corr = df_merge.filter(['Expression variation', 'Tissue-specificity (gene)', '# CREs', '% tissue-specific CREs', '% PhastCons']).corr('spearman')
            mask = np.triu(np.ones_like(corr, dtype=bool))
            f, ax = plt.subplots(figsize=(11, 9))
            cmap = sns.diverging_palette(240, 10, as_cmap=True)
            g = sns.heatmap(corr, mask=mask, cmap=cmap, vmax=1, center=0, vmin=-1,
                            square=True, linewidths=.5, cbar_kws={"shrink": .5}, annot=True)
            g.set_xticklabels(g.get_xticklabels(), rotation = 30, horizontalalignment='right')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_def}_expVar_heatmap.{fmt}', format=fmt, dpi=400)
            plt.close()
