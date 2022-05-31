###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2021.11.12
#
#   plot gene sets matched on covariates (expression level)
#   tis_spec threshold = 0.3
###

import os
import pandas as pd
import numpy  as np
import seaborn as sns
import matplotlib.pyplot as plt

from datetime import date
from scipy.stats import spearmanr
import statsmodels.api as sm


plt.switch_backend('agg')  # add to save plots non-interactively

### // constants and paths \\ ###
EXP_DAT_PATH = f'../../link_cre_to_genes/dat/2022-05-10'
CV_FILE_PATH = f'../dat'
RES_PATH = f'../fig/{str(date.today())}'

# create a date stamped dir for files
if not os.path.isdir(RES_PATH):
        os.makedirs(RES_PATH)

tis_order = ['spleen','liver', 'heart_left_ventricle', 'brain_hippocampus', 'lung', 'pancreas',
             'brain_prefrontal_cortex', 'psoas_muscle', 'small_intestine', 'ovary']
tis_names = ['Spleen','Liver', 'Heart', 'Hippocampus', 'Lung', 'Pancreas',
             'Prefrontal cortex', 'Muscle', 'Small intestine', 'Ovary']
tis_dict  = {'ovary':'Ovary', 'psoas_muscle':'Muscle', 'heart_left_ventricle':'Heart',
             'lung':'Lung', 'spleen':'Spleen', 'small_intestine':'Small intestine',
             'pancreas':'Pancreas', 'liver':'Liver', 'brain_prefrontal_cortex':'Prefrontal cortex',
             'brain_hippocampus':'Hippocampus'}

tisspec_thresh = 0.3

# plotting options
fmt='pdf'
rc = {'font.size':14, 'axes.titlesize':16,'axes.labelsize':14, 'legend.fontsize': 14,
      'xtick.labelsize': 14, 'ytick.labelsize': 14}

def spearman_pval(x,y):
    return spearmanr(x,y)[1]

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
        enh_num_by_gene = enh_num_by_gene.assign(tis_spec=lambda x: np.where(x.rel_entropy > tisspec_thresh, 1, 0))
        enh_num_by_gene = enh_num_by_gene.assign(exp_nocat=lambda x: np.where((x.exp == 1) & (x.hk == 0) & (x.lof_intol == 0) & (x.essential == 0), 1, 0))
        enh_num_by_gene = enh_num_by_gene.assign(exp_broad=lambda x: np.where((x.exp == 1) & (x.tis_spec == 0), 1, 0))
        enh_num_by_gene = enh_num_by_gene.assign(log2_exp=lambda x: np.log2(x[tis] + 1))
        enh_num_by_gene = enh_num_by_gene.rename(columns={tis:'gtex_exp'})
        enh_num_by_gene = enh_num_by_gene.assign(tissue=tis)
        df_lst.append(enh_num_by_gene)

    all_tis = pd.concat(df_lst)
    return all_tis

def read_cv_data(tis, thresh='0.8'):
    return pd.read_table(f'{CV_FILE_PATH}/{tis}_loess_cv_resid.csv', sep=',')

for landscape_def in ['loop', 'contact']:
    print(landscape_def)

    # create dataframe of enhancer number by gene with all tissues
    all_tis = read_data(landscape_def)

    for tis in tis_order:
        df = read_cv_data(tis)
        tis_exp = all_tis[all_tis['exp']==1].query(f'tissue=="{tis}"')
        df_merge = df.merge(tis_exp, left_on='Name', right_on='target_gene', how='inner', validate='1:1')

        with sns.plotting_context("paper", rc=rc):
            plt.figure(figsize=(4,4))
            g = sns.regplot(x='median_log2_tpm', y='log2cv', data=df_merge,
                           lowess=True, scatter_kws={'alpha':0.5, 'color':'.2', 's':1},
                           line_kws={'color':'tab:red'})
            plt.xlabel('log2(Expression)')
            plt.ylabel('log2(CV)')
            sns.despine()
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_def}_log2exp_v_cv_lmplot.{fmt}', format=fmt, dpi=400)
            plt.close()

        with sns.plotting_context("paper", rc=rc):
            plt.figure(figsize=(4,4))
            g = sns.regplot(x='median_log2_tpm', y='loess_resid', data=df_merge,
                           lowess=False, scatter_kws={'alpha':0.5, 'color':'.2', 's':1},
                           line_kws={'color':'tab:red'}, ci=False)
            plt.xlabel('log2(Expression)')
            plt.ylabel('Expression Variation')
            sns.despine()
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_def}_log2exp_v_expVar_lmplot.{fmt}', format=fmt, dpi=400)
            plt.close()

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


            gt0 = df_merge[df_merge['# CREs'] > 0]
            corr = gt0.filter(['Expression variation', 'Tissue-specificity (gene)', '# CREs', '% tissue-specific CREs', '% PhastCons']).corr('spearman')
            mask = np.triu(np.ones_like(corr, dtype=bool))
            f, ax = plt.subplots(figsize=(11, 9))
            cmap = sns.diverging_palette(240, 10, as_cmap=True)
            g = sns.heatmap(corr, mask=mask, cmap=cmap, vmax=1, center=0, vmin=-1,
                            square=True, linewidths=.5, cbar_kws={"shrink": .5}, annot=True)
            g.set_xticklabels(g.get_xticklabels(), rotation = 30, horizontalalignment='right')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_def}_expVar_heatmap_gt0.{fmt}', format=fmt, dpi=400)
            plt.close()

            print(f'{tis}, CRE > 0, spearman')
            with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 150):
                print(corr)

            print(f'{tis}, CRE > 0, p-values')
            with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 150):
                corr_p = gt0.filter(['Expression variation', 'Tissue-specificity (gene)', '# CREs', '% tissue-specific CREs', '% PhastCons']).corr(method=spearman_pval)
                print(corr_p)


for landscape_def in ['loop', 'contact']:
    dfs = []
    # create dataframe of enhancer number by gene with all tissues
    all_tis = read_data(landscape_def)

    # create dataframe of all tissues
    for tis in tis_order:
        df = read_cv_data(tis)
        tis_exp = all_tis[all_tis['exp']==1].query(f'tissue=="{tis}"')
        df_merge = df.merge(tis_exp, left_on='Name', right_on='target_gene', how='inner', validate='1:1')
        df_merge = df_merge.rename(columns={'loess_resid':'Expression variation', 'rel_entropy':'Tissue-specificity (gene)',
                                    'enh_num':'# CREs', 'frac_tisspec_enh':'% tissue-specific CREs', 'frac_phastcons':'% PhastCons'})
        corr = df_merge.filter(['Expression variation', 'Tissue-specificity (gene)', '# CREs', '% tissue-specific CREs', '% PhastCons']).corr('spearman')
        d = pd.DataFrame(corr['Expression variation']).T.assign(Tissue=tis)
        dfs.append(d)

    # combine all tissues
    df_corr = pd.concat(dfs).set_index('Tissue').drop(columns='Expression variation')


    # plot heatmap of all tissues
    with sns.plotting_context("paper", rc=rc):
        f, ax = plt.subplots(figsize=(12, 9))
        cmap = sns.diverging_palette(240, 10, as_cmap=True)
        g = sns.heatmap(df_corr, cmap=cmap, vmax=1, center=0, vmin=-1,
                        square=True, linewidths=.5, cbar_kws={"shrink": .5}, annot=True)
        g.set_xticklabels(g.get_xticklabels(), rotation = 30, horizontalalignment='right')
        g.set_yticklabels(tis_names)
        g.set_ylabel('')
        plt.tight_layout()
        plt.savefig(f'{RES_PATH}/{landscape_def}_expVar_heatmap_alltis.{fmt}', format=fmt, dpi=400)
        plt.close()
