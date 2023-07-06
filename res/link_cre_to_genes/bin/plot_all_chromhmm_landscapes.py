###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2020.07.22
#
#   manuscript plots for link_enh_to_genes output, both landscape types
###

import os
import logging
import pandas as pd
import numpy  as np
import seaborn as sns
import matplotlib.pyplot as plt
import scikit_posthocs as sp

from datetime import date
from scipy import stats

plt.switch_backend('agg')  # add to save plots non-interactively

### // constants and paths \\ ###
EXP_DATA_PATH = f'../dat/2022-12-20'
RES_PATH = f'../fig/{str(date.today())}'
#RES_PATH = f'../fig/2022-05-10'

# create a date stamped dir for files
if not os.path.isdir(RES_PATH):
    os.makedirs(RES_PATH)

tis_order = ['spleen','liver', 'heart_left_ventricle', 'brain_hippocampus', 'lung', 'pancreas',
             'brain_prefrontal_cortex', 'psoas_muscle', 'small_intestine', 'ovary']

tis_names = ['Spleen','Liver', 'Heart', 'Hippocampus', 'Lung', 'Pancreas',
             'Prefrontal cortex', 'Muscle', 'Small intestine', 'Ovary']

# plotting options
fmt='pdf'
rc = {'font.size':14, 'axes.titlesize':16,'axes.labelsize':14, 'legend.fontsize': 14,
      'xtick.labelsize': 14, 'ytick.labelsize': 14}

# set up log file to record data about run
logging.basicConfig(filename=f'{RES_PATH}/README', level=logging.INFO)

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
            infile = f'{EXP_DATA_PATH}/{tis}_gene_with_chromhmm-enh-count_frac-enh-spec_hicQ05-linked.tsv'
        elif landscape_def == 'loop':
            infile = f'{EXP_DATA_PATH}/{tis}_gene_with_chromhmm-enh-count_frac-enh-spec_peakachuloop-linked.tsv'
        enh_num_by_gene = pd.read_table(infile)
        enh_num_by_gene = enh_num_by_gene.assign(rel_entropy_bins=lambda x: pd.cut(x['rel_entropy'], bins=10))
        enh_num_by_gene = enh_num_by_gene.assign(enh_rel_entropy_bins=lambda x: pd.cut(x['enh_rel_entropy'], bins=5))
        enh_num_by_gene = enh_num_by_gene.assign(exp_nocat=lambda x: np.where((x.exp == 1) & (x.hk == 0) & (x.lof_intol == 0) & (x.essential == 0), 1, 0))
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
### \\


logging.info(f'Creating enh_num_by_gene dataframe for all tissues from {EXP_DATA_PATH}')

for landscape_def in ['loop', 'contact']:
    logging.info(f'Using {landscape_def}-based definition')
    ### // create dataframe of enhancer number by gene with all tissues \\ ###
    all_tis = read_data(landscape_def)
    all_tis_exp = all_tis[all_tis['exp']==1].copy(deep=True)
    all_tis_exp.head()

    # df for plots by gene type using hue
    enh_by_gene_anno = (all_tis.filter(['target_gene', 'tissue', 'exp', 'hk', 'lof_intol', 'essential', 'exp_nocat', 'enh_num'])
                                        .set_index(['target_gene', 'tissue', 'enh_num', 'exp'])
                                        .stack()
                                        .reset_index()
                                        .rename(columns={'level_4':'anno', 0:'val'})
                                        .query('val == 1')
                                        .merge(all_tis.filter(['target_gene', 'tissue', 'cds_length', 'gene_length', 'enh_rel_entropy', 'frac_tisspec_enh_bins', 'frac_phastcons']), how='inner', validate='m:1'))
    ### \\


    ### // plots of distribution of number of enhancers \\ ###
    logging.info('Plot distribution of the number of enhancers per gene')

    logging.info('Number of genes per tissue (exp/not):')
    logging.info(all_tis.groupby(['tissue', 'exp']).target_gene.count())

    with sns.plotting_context("paper", rc=rc):
        ### fig: distribution of enhancer number per gene, exp/not, no shared axis
        g = sns.FacetGrid(all_tis, col='tissue', col_wrap=5, sharex=False, sharey=False,
                          col_order=tis_order, height=4)
        g = g.map(sns.histplot, 'enh_num', color='#363737', bins=25, kde=False)
        g.fig.subplots_adjust(hspace=0.4, wspace=.15)
        facet_titles(g.axes)
        facet_x_axis(g.axes, 'Number of CREs')
        sns.despine()
        plt.savefig(f'{RES_PATH}/all_{landscape_def}_chromhmm-enh-num_bytissue_distplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot of number of enhancers by gene type (exp/not), no fliers
        g = sns.catplot(x='exp', y='enh_num', col='tissue', data=enh_by_gene_anno,
                        palette=['tab:gray', 'tab:blue'], col_order=tis_order,
                        kind="box", col_wrap=5, sharex=False, sharey=True,
                        showfliers=False, fliersize=1)
        g.fig.subplots_adjust(hspace=0.25, wspace=.15)
        facet_titles(g.axes)
        facet_x_axis(g.axes, '')
        g.set_ylabels('Number of CREs')
        g.set_xticklabels(['Not Expressed', 'Expressed'])
        sns.despine()
        plt.savefig(f'{RES_PATH}/all_{landscape_def}_expvchromhmm-enh-num_bytissue_boxplot_nofliers.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot of enh num by expression status for all tissues
        fig = plt.figure(figsize=(8,10))
        g = sns.boxplot(data=enh_by_gene_anno,
                        x='enh_num', y='tissue', hue='exp', order=tis_order,
                        palette=['tab:gray', 'tab:blue'], showfliers=False)
        g.set_yticklabels(tis_names)
        g.set_xlabel('Number of CREs')
        g.set_ylabel('')
        g.legend(handles=g.legend_.legendHandles, bbox_to_anchor=(1.01, 1), frameon=False,
                 labels=['Not expressed', 'Expressed'])
        sns.despine()
        plt.tight_layout()
        plt.savefig(f'{RES_PATH}/all_{landscape_def}_expvchromhmm-enh-num_bytissue_boxplot_hue.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: H_boxplot of enh num by expression status for all tissues
        fig = plt.figure(figsize=(10,6))
        g = sns.boxplot(data=enh_by_gene_anno,
                        y='enh_num', x='tissue', hue='exp', order=tis_order,
                        palette=['tab:gray', 'tab:blue'], showfliers=False)
        g.set_xticklabels(tis_names, rotation=30, horizontalalignment='right')
        g.set_ylabel('Number of CREs')
        g.set_xlabel('')
        g.legend(handles=g.legend_.legendHandles, bbox_to_anchor=(1.01, 1), frameon=False,
                 labels=['Not expressed', 'Expressed'])
        sns.despine()
        plt.tight_layout()
        plt.savefig(f'{RES_PATH}/all_{landscape_def}_expvchromhmm-enh-num_bytissue_H_boxplot_hue.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: H_boxplot of enh num by expression status for all tissues, fliers
        fig = plt.figure(figsize=(10,6))
        g = sns.boxplot(data=enh_by_gene_anno, order=tis_order,
                        y='enh_num', x='tissue', hue='exp', flierprops=dict(marker='o', markersize=1),
                        palette=['tab:gray', 'tab:blue'], showfliers=True)
        g.set_xticklabels(tis_names, rotation=30, horizontalalignment='right')
        g.set_ylabel('Number of CREs')
        g.set_xlabel('')
        g.legend(handles=g.legend_.legendHandles, bbox_to_anchor=(1.01, 1), frameon=False,
                 labels=['Not expressed', 'Expressed'])
        sns.despine()
        plt.tight_layout()
        plt.savefig(f'{RES_PATH}/all_{landscape_def}_expvchromhmm-enh-num_bytissue_H_boxplot_hue_fliers.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: H_boxplot of enh num by expression status for all tissues, # CRE > 0
        fig = plt.figure(figsize=(10,6))
        g = sns.boxplot(data=enh_by_gene_anno.query('enh_num > 0'),
                        y='enh_num', x='tissue', hue='exp', order=tis_order,
                        palette=['tab:gray', 'tab:blue'], showfliers=False)
        g.set_xticklabels(tis_names, rotation=30, horizontalalignment='right')
        g.set_ylabel('Number of CREs')
        g.set_xlabel('')
        g.legend(handles=g.legend_.legendHandles, bbox_to_anchor=(1.01, 1), frameon=False,
                 labels=['Not expressed', 'Expressed'])
        sns.despine()
        plt.tight_layout()
        plt.savefig(f'{RES_PATH}/all_{landscape_def}_expvchromhmm-enh-num_bytissue_H_boxplot_hue_gt0.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: H split violinplot of enh num by expression status for all tissues
        fig = plt.figure(figsize=(10,6))
        g = sns.violinplot(data=enh_by_gene_anno, order=tis_order,
                        y='enh_num', x='tissue', hue='exp', inner="quart", linewidth=1, cut=0,
                        palette=['tab:gray', 'tab:blue'], split=True)
        g.set_xticklabels(tis_names, rotation=30, horizontalalignment='right')
        g.set_ylabel('Number of CREs')
        g.set_xlabel('')
        g.legend(handles=g.legend_.legendHandles, bbox_to_anchor=(1.01, 1), frameon=False,
                 labels=['Not expressed', 'Expressed'])
        sns.despine()
        plt.tight_layout()
        plt.savefig(f'{RES_PATH}/all_{landscape_def}_expvchromhmm-enh-num_bytissue_H_violinplot_hue.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

    ### fig/table: MWU test of number of enhancers exp vs. not
    logging.info('# CRE by expressed v. not expressed genes')
    for tis_name in tis_order:
        ts, p = stats.mannwhitneyu(all_tis[all_tis['exp']==0].query(f'tissue=="{tis_name}"').enh_num,
                                   all_tis[all_tis['exp']==1].query(f'tissue=="{tis_name}"').enh_num, alternative='less')
        logging.info(f'{tis_name}, p = {p:.3}, corrected = {p < (0.05/10)}')

    logging.info('Not expressed')
    logging.info(all_tis[all_tis['exp']==0].groupby('tissue').enh_num.median())
    logging.info('Expressed')
    logging.info(all_tis[all_tis['exp']==1].groupby('tissue').enh_num.median())

    logging.info('Not expressed, CRE >0')
    logging.info(all_tis[(all_tis['exp']==0) & (all_tis.enh_num>0)].groupby('tissue').enh_num.median())
    logging.info('Expressed, CRE > 0')
    logging.info(all_tis[(all_tis['exp']==1) & (all_tis.enh_num>0)].groupby('tissue').enh_num.median())
    ### end_fig
    ### \\


    ### // plots of distribution of fraction tissue-specific enhancers \\ ###
    logging.info('Generate dataframe of the gene annotations mapped to enhancer attributes (Osterwalder-style)')

    # create expressed-only version
    enh_by_expgene_anno = (all_tis_exp
                               .filter(['target_gene', 'tissue', 'exp', 'hk', 'lof_intol', 'essential', 'exp_nocat', 'enh_num'])
                               .assign(expressed=lambda x: x.exp)
                               .set_index(['target_gene', 'tissue', 'enh_num', 'exp'])
                               .stack()
                               .reset_index()
                               .rename(columns={'level_4':'anno', 0:'val'})
                               .query('val == 1')
                               .merge(all_tis_exp.filter(['target_gene', 'gtex_exp', 'cds_length', 'gene_length', 'tissue', 'enh_rel_entropy', 'frac_tisspec_enh', 'frac_phastcons']), how='inner', validate='m:1'))

    enh_by_expgene_anno = enh_by_expgene_anno.assign(num_tisspec_enh=lambda x: x.enh_num * x.frac_tisspec_enh)
    enh_by_expgene_anno = enh_by_expgene_anno.assign(gtex_exp_log2=lambda x: np.log2(x.gtex_exp))

    logging.info('Plot enhancer attributes by gene annotation (Osterwalder-style)')

    with sns.plotting_context("paper", rc=rc):
        ### fig: boxplot of conserved enhancer bp by annotation, exp only, fliers
        fig = plt.figure(figsize=(8,10))
        g = sns.boxplot(x='frac_phastcons', hue='anno', y='tissue', data=enh_by_expgene_anno.query('anno=="hk" | anno=="lof_intol" | anno=="exp_nocat"'),
                        palette=['tab:green', 'tab:orange', 'tab:blue'],
                        hue_order=['hk', 'lof_intol', 'exp_nocat'], order=tis_order,
                        showfliers=True, fliersize=1)
        g.set_yticklabels(tis_names)
        g.set_xlabel('Fraction of conserved CRE bp')
        g.set_ylabel('')
        g.legend(handles=g.legend_.legendHandles, frameon=False, loc='upper right',
                 labels=['Housekeeping', 'LoF intolerant', 'Expressed'], title='Gene Category')
        sns.despine()
        plt.savefig(f'{RES_PATH}/all_{landscape_def}_genetype_noess_v_fracconsv_bytissue_exponly_fliers_boxplot_hue_chromhmm.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot of enhancer number by annotation, exp only, no fliers
        fig = plt.figure(figsize=(8,10))
        g = sns.boxplot(x='enh_num', hue='anno', y='tissue', data=enh_by_expgene_anno.query('anno=="hk" | anno=="lof_intol" | anno=="exp_nocat"'),
                        palette=['tab:green', 'tab:orange', 'tab:blue'],
                        hue_order=['hk', 'lof_intol', 'exp_nocat'], order=tis_order,
                        showfliers=False, fliersize=1)
        g.set_yticklabels(tis_names)
        g.set_xlabel('Number of CREs')
        g.set_ylabel('')
        g.legend(handles=g.legend_.legendHandles, frameon=False, loc='upper right',
                 labels=['Housekeeping', 'LoF intolerant', 'Expressed'], title='Gene Category')
        sns.despine()
        plt.savefig(f'{RES_PATH}/all_{landscape_def}_genetype_noess_v_chromhmm-enhnum_bytissue_exponly_nofliers_boxplot_hue.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot of enhancer number by annotation, exp only, no fliers
        fig = plt.figure(figsize=(8,10))
        g = sns.boxplot(x='enh_num', hue='anno', y='tissue', data=enh_by_expgene_anno.query('anno=="hk" | anno=="essential" | anno=="lof_intol" | anno=="exp_nocat"'),
                        palette=['tab:green', 'tab:red', 'tab:orange', 'tab:blue'],
                        hue_order=['hk', 'essential', 'lof_intol', 'exp_nocat'], order=tis_order,
                        showfliers=False, fliersize=1)
        g.set_yticklabels(tis_names)
        g.set_xlabel('Number of CREs')
        g.set_ylabel('')
        g.legend(handles=g.legend_.legendHandles, frameon=False, loc='upper right',
                 labels=['Housekeeping', 'Essential', 'LoF intolerant', 'Expressed'], title='Gene Category')
        sns.despine()
        plt.savefig(f'{RES_PATH}/all_{landscape_def}_genetype_v_chromhmm-enhnum_bytissue_exponly_nofliers_boxplot_hue.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig


    ### fig: boxplot of enh num by log2 expression for each tissue
    for idx, tis in enumerate(tis_order):
        with sns.plotting_context("paper", rc=rc):
            nonzero = all_tis_exp.query(f'enh_num>0 & tissue=="{tis}"').assign(exp_bins=lambda x: pd.qcut(x['log2_exp'], q=4, labels=[1,2,3,4]))
            g = sns.catplot(x='exp_bins', y='enh_num', data=nonzero,
                            kind='box', notch=True, order=[1,2,3,4],
                            palette='Greys', dodge=False,
                            showfliers=False, width=.6, height=6, aspect=.7)
            g.set_xlabels('Expression Quartile')
            g.set_ylabels('Number of CREs')
            g.fig.suptitle(f'{tis_names[idx]}', x=0.6, y=0.92)
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_def}_chromhmm-enh-numXlog2exp_bytissue_boxplot.{fmt}', format=fmt, dpi=400)
            plt.close()
    ### end_fig

    ### fig/table: log exp v. # CRE stats to file
    logging.info('Correlation between log2 expression and CRE decile, test 1st v. 4th')
    for tis_name in tis_order:
        nonzero = all_tis_exp.query(f'enh_num>0 & tissue=="{tis_name}"').assign(exp_bins=lambda x: pd.qcut(x['log2_exp'], q=4, labels=[1,2,3,4]))
        rho, corrp = stats.spearmanr(nonzero.query(f'tissue=="{tis_name}"').enh_num, nonzero.query(f'tissue=="{tis_name}"').log2_exp)
        ts, p = stats.mannwhitneyu(nonzero.query(f'tissue=="{tis_name}" & exp_bins == 1').enh_num,
                                   nonzero.query(f'tissue=="{tis_name}" & exp_bins == 4').enh_num)
        logging.info(f'{tis_name}, 1stv4th p = {p:.3}, spearmanR = {rho:.3}, p = {corrp:.3}')

        logging.info('Median # CRE per CRE quartile')
        logging.info(nonzero.groupby(['tissue', 'exp_bins']).enh_num.median().to_string())
    ### end_fig

    ### fig/table: log stats to file, kw # CRE
    logging.info('KW between CRE # and gene type (hk, lof, expressed)')
    for tis in tis_order:
        logging.info(f'Running enhancer number Kruskal Wallis on {tis}')
        ts, p = stats.kruskal(*[group['enh_num'].values for name, group in enh_by_expgene_anno.query(f'tissue=="{tis}" & anno!="expressed" & anno!="essential"').groupby('anno')])
        logging.info(enh_by_expgene_anno.query(f'tissue=="{tis}" & anno!="expressed" & anno!="essential"').groupby('anno').mean())
        logging.info(f'{tis} test_stat = {ts}, p = {p}')

        if p < 0.05:
            logging.info(f'Running post-hoc Dunn test on {tis}')
            p = sp.posthoc_dunn(enh_by_expgene_anno.query(f'tissue=="{tis}" & anno!="expressed" & anno!="essential"'), val_col='enh_num', group_col='anno', sort=True, p_adjust='fdr_bh')
            logging.info(p)
            logging.info(sp.sign_table(p))
    ### end_fig
    ### \\
