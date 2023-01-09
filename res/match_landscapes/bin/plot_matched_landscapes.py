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
from scipy import stats
import statsmodels.api as sm

plt.switch_backend('agg')  # add to save plots non-interactively

### // constants and paths \\ ###
CHROMHMM = True
if CHROMHMM:
    EXP_DAT_PATH = f'../../link_cre_to_genes/dat/2022-12-20'
    MAT_DAT_PATH = f'../dat/2023-01-09'
else:
    EXP_DAT_PATH = f'../../link_cre_to_genes/dat/2022-05-10'
    MAT_DAT_PATH = f'../dat/2022-05-20'

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

if CHROMHMM:
    landscape_label = f'chromhmm_{landscape_def}'
else:
    landscape_label = landscape_def

def read_data(landscape_def, chmm):
    df_lst = []
    for tis in tis_order:
        if chmm:
            if landscape_def == 'contact':
                infile = f'{EXP_DATA_PATH}/{tis}_gene_with_chromhmm-enh-count_frac-enh-spec_hicQ05-linked.tsv'
            elif landscape_def == 'loop':
                infile = f'{EXP_DATA_PATH}/{tis}_gene_with_chromhmm-enh-count_frac-enh-spec_peakachuloop-linked.tsv'
        else:
            if landscape_def == 'contact':
                infile = f'{EXP_DATA_PATH}/{tis}_gene_with_enh-count_frac-enh-spec_hicQ05-linked.tsv'
            elif landscape_def == 'loop':
                infile = f'{EXP_DATA_PATH}/{tis}_gene_with_enh-count_frac-enh-spec_peakachuloop-linked.tsv'
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
    tisspec_bins = ['(0,0.1]', '(0.1,0.2]', '(0.2,0.3]', '(0.3,0.4]', '(0.4,0.5]',
                    '(0.5,0.6]', '(0.6,0.7]', '(0.7,0.8]', '(0.8,0.9]', '(0.9,1.0]']
    all_tis['frac_tisspec_enh_bins'] = pd.cut(all_tis.frac_tisspec_enh, bins=10,
                                              labels=tisspec_bins)
    return all_tis

def read_matched_data(tis, landscape_type, gene_type, chmm):
    if chmm:
        return pd.read_table(f'{MAT_DAT_PATH}/matched_{gene_type}_{tis}_{landscape_type}_chromhmm.tsv', sep='\t')
    return pd.read_table(f'{MAT_DAT_PATH}/matched_{gene_type}_{tis}_{landscape_type}.tsv', sep='\t')
### \\

for landscape_def in ['loop', 'contact']:
    print(landscape_def)
    # create dataframe of enhancer number by gene with all tissues
    all_tis = read_data(landscape_def, CHROMHMM)
    all_tis_exp = all_tis[all_tis['exp']==1]


    ### // full hk and lof sets \\ ###
    for tis in tis_order:
        df = read_matched_data(tis=tis, gene_type='hk', landscape_type=landscape_def, chmm=CHROMHMM)
        merged = df.merge(all_tis_exp.query(f'tissue=="{tis}"'), how='left')

        with sns.plotting_context("paper", rc=rc):
            ### fig: plot boxplot of hk v. expressed by cre_num, matched, fliers
            g = sns.catplot(x='anno', y='enh_num', data=merged, kind='box',
                            palette=['tab:green', 'tab:blue'], order=['Housekeeping', 'Expressed'], showfliers=True)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').enh_num, merged.query(f'anno=="Expressed"').enh_num)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_ylabels('Number of CREs')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXenh-num_hk_matched_boxplot_fliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: ecdf of hk v. expressed by cre_num, matched
            g = sns.displot(hue='anno', y='enh_num', data=merged, kind='ecdf',
                            palette=['tab:green', 'tab:blue'], hue_order=['Housekeeping', 'Expressed'],
                            aspect=1, legend=False)
            g.set_xlabels('')
            g.set_ylabels('Number of CREs')
            plt.legend(loc='upper left', frameon=False, labels=['Housekeeping', 'Expressed'])
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXenh-num_hk_matched_ecdf.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: plot boxplot of hk v. expressed by cre_num, matched, notch, nofliers
            g = sns.catplot(x='anno', y='enh_num', data=merged, kind='box',
                            notch=True, order=['Housekeeping', 'Expressed'],
                            hue='anno', hue_order=['Housekeeping', 'Expressed'],
                            palette=['tab:green', 'tab:blue'], dodge=False,
                            showfliers=False, width=.65, height=4.5, aspect=.8)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').enh_num, merged.query(f'anno=="Expressed"').enh_num)
            if p < 0.01:
              plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.1e}', horizontalalignment='center')
            else:
              plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_xticklabels(['HK', 'Expressed'])
            g.set_ylabels('Number of CREs')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXenh-num_hk_matched_boxplot_nofliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: plot boxplot of hk v. expressed by % phastcons, matched, fliers
            g = sns.catplot(x='anno', y='frac_phastcons', data=merged, kind='box',
                            palette=['tab:green', 'tab:blue'], order=['Housekeeping', 'Expressed'], showfliers=True)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').frac_phastcons, merged.query(f'anno=="Expressed"').frac_phastcons)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_ylabels('% PhastCons')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXphastcons_hk_matched_boxplot_fliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: plot boxplot of hk v. expressed by % phastcons, matched, nofliers
            g = sns.catplot(x='anno', y='frac_phastcons', data=merged, kind='box',
                            notch=True, order=['Housekeeping', 'Expressed'],
                            hue='anno', hue_order=['Housekeeping', 'Expressed'],
                            palette=['tab:green', 'tab:blue'], dodge=False,
                            showfliers=False, width=.65, height=4.5, aspect=.8)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').frac_phastcons, merged.query(f'anno=="Expressed"').frac_phastcons)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_xticklabels(['HK', 'Expressed'])
            g.set_ylabels('% PhastCons')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXphastcons_hk_matched_boxplot_nofliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig


        df = read_matched_data(tis=tis, gene_type='lof', landscape_type=landscape_def, chmm=CHROMHMM)
        merged = df.merge(all_tis_exp.query(f'tissue=="{tis}"'), how='left')

        with sns.plotting_context("paper", rc=rc):
            ### fig: plot boxplot of lof v. expressed by cre_num, matched, fliers
            g = sns.catplot(x='anno', y='enh_num', data=merged, kind='box',
                            palette=['tab:orange', 'tab:blue'], order=['LoF Intolerant', 'Expressed'], showfliers=True)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').enh_num, merged.query(f'anno=="Expressed"').enh_num)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_ylabels('Number of CREs')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXenh-num_lof_matched_boxplot_fliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: ecdf of lof v. expressed by cre_num, matched
            g = sns.displot(hue='anno', y='enh_num', data=merged, kind='ecdf',
                            palette=['tab:orange', 'tab:blue'], hue_order=['LoF Intolerant', 'Expressed'],
                            aspect=1, legend=False)
            g.set_xlabels('')
            g.set_ylabels('Number of CREs')
            plt.legend(loc='upper left', frameon=False, labels=['Expressed', 'LoF Intolerant'])
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXenh-num_lof_matched_ecdf.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: plot boxplot of lof v. expressed by cre_num, matched, nofliers
            g = sns.catplot(x='anno', y='enh_num', data=merged, kind='box',
                            palette=['tab:orange', 'tab:blue'], order=['LoF Intolerant', 'Expressed'], showfliers=False)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').enh_num, merged.query(f'anno=="Expressed"').enh_num)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_ylabels('Number of CREs')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXenh-num_lof_matched_boxplot_nofliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: plot boxplot of lof v. expressed by cre_num, matched, notch, nofliers
            g = sns.catplot(x='anno', y='enh_num', data=merged, kind='box', notch=True,
                        hue='anno', hue_order=['LoF Intolerant', 'Expressed'], palette=['tab:orange', 'tab:blue'],
                        dodge=False, showfliers=False, width=.65, height=4.5, aspect=.8)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').enh_num, merged.query(f'anno=="Expressed"').enh_num)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_ylabels('Number of CREs')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXenh-num_lof_matched_boxplot_nofliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: plot boxplot of hk v. expressed by % phastcons, matched, fliers
            g = sns.catplot(x='anno', y='frac_phastcons', data=merged, kind='box',
                            palette=['tab:orange', 'tab:blue'], order=['LoF Intolerant', 'Expressed'], showfliers=True)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').frac_phastcons, merged.query(f'anno=="Expressed"').frac_phastcons)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_ylabels('% PhastCons')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXphastcons_lof_matched_boxplot_fliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: plot boxplot of lof v. expressed by % phastcons, matched, notch, nofliers
            g = sns.catplot(x='anno', y='frac_phastcons', data=merged, kind='box',
                            notch=True, order=['LoF Intolerant', 'Expressed'],
                            hue='anno', hue_order=['LoF Intolerant', 'Expressed'],
                            palette=['tab:orange', 'tab:blue'], dodge=False,
                            showfliers=False, width=.65, height=4.5, aspect=.8)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').frac_phastcons, merged.query(f'anno=="Expressed"').frac_phastcons)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_ylabels('% PhastCons')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXphastcons_lof_matched_boxplot_nofliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig
    ### \\


    ### // CRE > 0 hk and lof sets \\ ###
    for tis in tis_order:
        df = read_matched_data(tis=tis, gene_type='hk_gt0', landscape_type=landscape_def, chmm=CHROMHMM)
        merged = df.merge(all_tis_exp.query(f'tissue=="{tis}"'), how='left')

        with sns.plotting_context("paper", rc=rc):
            ### fig: plot boxplot of hk v. expressed by cre_num, matched, fliers
            g = sns.catplot(x='anno', y='enh_num', data=merged, kind='box',
                            palette=['tab:green', 'tab:blue'], order=['Housekeeping', 'Expressed'], showfliers=True)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').enh_num, merged.query(f'anno=="Expressed"').enh_num)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_ylabels('Number of CREs')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXenh-num_hk_gt0_matched_boxplot_fliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: plot boxplot of hk v. expressed by cre_num, matched, notch, nofliers
            g = sns.catplot(x='anno', y='enh_num', data=merged, kind='box',
                            notch=True, order=['Housekeeping', 'Expressed'],
                            hue='anno', hue_order=['Housekeeping', 'Expressed'],
                            palette=['tab:green', 'tab:blue'], dodge=False,
                            showfliers=False, width=.65, height=4.5, aspect=.8)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').enh_num, merged.query(f'anno=="Expressed"').enh_num)
            if p < 0.01:
              plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.1e}', horizontalalignment='center')
            else:
              plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_xticklabels(['HK', 'Expressed'])
            g.set_ylabels('Number of CREs')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXenh-num_hk_gt0_matched_boxplot_nofliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: ecdf of hk v. expressed by cre_num, matched
            g = sns.displot(hue='anno', y='enh_num', data=merged, kind='ecdf',
                            palette=['tab:green', 'tab:blue'], hue_order=['Housekeeping', 'Expressed'],
                            aspect=1, legend=False)
            g.set_xlabels('')
            g.set_ylabels('Number of CREs')
            plt.legend(loc='upper left', frameon=False, labels=['Housekeeping', 'Expressed'])
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXenh-num_hk_gt0_matched_ecdf.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: plot boxplot of hk v. expressed by % phastcons, matched, nofliers
            g = sns.catplot(x='anno', y='frac_phastcons', data=merged, kind='box',
                            notch=True, order=['Housekeeping', 'Expressed'],
                            hue='anno', hue_order=['Housekeeping', 'Expressed'],
                            palette=['tab:green', 'tab:blue'], dodge=False,
                            showfliers=False, width=.65, height=4.5, aspect=.8)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').frac_phastcons, merged.query(f'anno=="Expressed"').frac_phastcons)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_xticklabels(['HK', 'Expressed'])
            g.set_ylabels('% PhastCons')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXphastcons_hk_gt0_matched_boxplot_nofliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: violinplot of hk v. expressed by % tisspec CRE, matched
            g = sns.violinplot(x='anno', y='frac_tisspec_enh', data=merged, cut=0, inner='quart',
                               palette=['tab:green', 'tab:blue'], order=['Housekeeping', 'Expressed'])
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').frac_tisspec_enh, merged.query(f'anno=="Expressed"').frac_tisspec_enh)
            if p < 0.01:
                plt.text(0.5, g.get_ylim()[1], f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.get_ylim()[1], f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabel('')
            g.set_ylabel('Fraction tissue-specific CRES')
            sns.despine()
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXfrac-tispec-cre_hk_gt0_matched_violinplot.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: ecdf of hk v. expressed by %tis_spec CRE, matched
            g = sns.displot(hue='anno', y='frac_tisspec_enh', data=merged, kind='ecdf',
                            palette=['tab:green', 'tab:blue'], hue_order=['Housekeeping', 'Expressed'],
                            aspect=1, legend=False)
            g.set_xlabels('')
            g.set_ylabels('Fraction of tissue-specifc CREs')
            plt.legend(loc='upper left', frameon=False, labels=['Expressed', 'Housekeeping'])
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXfrac-tispec-cre_hk_gt0_matched_ecdf.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

        df = read_matched_data(tis=tis, gene_type='lof_gt0', landscape_type=landscape_def, chmm=CHROMHMM)
        merged = df.merge(all_tis_exp.query(f'tissue=="{tis}"'), how='left')

        with sns.plotting_context("paper", rc=rc):
            ### fig: plot boxplot of lof v. expressed by cre_num, matched, notch, nofliers
            g = sns.catplot(x='anno', y='enh_num', data=merged, kind='box', notch=True,
                        hue='anno', hue_order=['LoF Intolerant', 'Expressed'], palette=['tab:orange', 'tab:blue'],
                        dodge=False, showfliers=False, width=.65, height=4.5, aspect=.8)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').enh_num, merged.query(f'anno=="Expressed"').enh_num)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1]+1, f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_ylabels('Number of CREs')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXenh-num_lof_gt0_matched_boxplot_nofliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: ecdf of lof v. expressed by cre_num, matched
            g = sns.displot(hue='anno', y='enh_num', data=merged, kind='ecdf',
                            palette=['tab:orange', 'tab:blue'], hue_order=['LoF Intolerant', 'Expressed'],
                            aspect=1, legend=False)
            g.set_xlabels('')
            g.set_ylabels('Number of CREs')
            plt.legend(loc='upper left', frameon=False, labels=['Expressed', 'LoF Intolerant'])
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXenh-num_lof_gt0_matched_ecdf.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: plot boxplot of lof v. expressed by % phastcons, matched, fliers
            g = sns.catplot(x='anno', y='frac_phastcons', data=merged, kind='box',
                            palette=['tab:orange', 'tab:blue'], order=['LoF Intolerant', 'Expressed'], showfliers=True)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').frac_phastcons, merged.query(f'anno=="Expressed"').frac_phastcons)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_ylabels('% PhastCons')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXphastcons_lof_gt0_matched_boxplot_fliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: plot boxplot of lof v. expressed by % phastcons, matched, notch, nofliers
            g = sns.catplot(x='anno', y='frac_phastcons', data=merged, kind='box',
                            notch=True, order=['LoF Intolerant', 'Expressed'],
                            hue='anno', hue_order=['LoF Intolerant', 'Expressed'],
                            palette=['tab:orange', 'tab:blue'], dodge=False,
                            showfliers=False, width=.65, height=4.5, aspect=.8)
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').frac_phastcons, merged.query(f'anno=="Expressed"').frac_phastcons)
            if p < 0.01:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.ax.get_ylim()[1], f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabels('')
            g.set_ylabels('% PhastCons')
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXphastcons_lof_gt0_matched_boxplot_nofliers.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: violinplot of hk v. expressed by % tisspec CRE, matched
            g = sns.violinplot(x='anno', y='frac_tisspec_enh', data=merged, cut=0, inner='quart',
                               palette=['tab:orange', 'tab:blue'], order=['LoF Intolerant', 'Expressed'])
            ts, p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').frac_tisspec_enh, merged.query(f'anno=="Expressed"').frac_tisspec_enh)
            if p < 0.01:
                plt.text(0.5, g.get_ylim()[1], f'p = {p:.1e}', horizontalalignment='center')
            else:
                plt.text(0.5, g.get_ylim()[1], f'p = {p:.2f}', horizontalalignment='center')
            g.set_xlabel('')
            g.set_ylabel('Fraction tissue-specific CRES')
            sns.despine()
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXfrac-tispec-cre_lof_gt0_matched_violinplot.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: ecdf of lof v. expressed by %tis_spec CRE, matched
            g = sns.displot(hue='anno', y='frac_tisspec_enh', data=merged, kind='ecdf',
                            palette=['tab:orange', 'tab:blue'], hue_order=['LoF Intolerant', 'Expressed'],
                            aspect=1, legend=False)
            g.set_xlabels('')
            g.set_ylabels('Fraction of tissue-specifc CREs')
            plt.legend(loc='upper left', frameon=False, labels=['Expressed', 'LoF Intolerant'])
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/{tis}_{landscape_label}_categoryXfrac-tispec-cre_lof_gt0_matched_ecdf.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig
    ### \\

    ### // write p value table \\ ###
    hrows = []
    lrows = []

    print('')
    print(f'tisue\tmed_cre_num_expressed\tmed_cre_num_housekeeping\tmed_frac_phastcons_expressed\tmed_frac_phastcons_housekeeping')
    for tis in tis_order:
        df = read_matched_data(tis=tis, gene_type='hk', landscape_type=landscape_def, chmm=CHROMHMM)
        merged = df.merge(all_tis_exp.query(f'tissue=="{tis}"'), how='left')

        # print medians
        cre_med = merged.groupby('anno').enh_num.median()
        phast_med = merged.groupby('anno').frac_phastcons.median()
        print(f'{tis}\t{cre_med["Expressed"]}\t{cre_med["Housekeeping"]}\t{phast_med["Expressed"]}\t{phast_med["Housekeeping"]}')

        # print MWU stats
        _, cre_p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').enh_num, merged.query(f'anno=="Expressed"').enh_num)
        _, phast_p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').frac_phastcons, merged.query(f'anno=="Expressed"').frac_phastcons)
        _, frac_p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').frac_tisspec_enh, merged.query(f'anno=="Expressed"').frac_tisspec_enh)
        hrows.append([tis, 'housekeeping', cre_p, phast_p, frac_p])

    print('')
    print(f'tisue\tmed_cre_num_expressed\tmed_cre_num_lof\tmed_frac_phastcons_expressed\tmed_frac_phastcons_lif')
    for tis in tis_order:
        df = read_matched_data(tis=tis, gene_type='lof', landscape_type=landscape_def, chmm=CHROMHMM)
        merged = df.merge(all_tis_exp.query(f'tissue=="{tis}"'), how='left')

        # print medians
        cre_med = merged.groupby('anno').enh_num.median()
        phast_med = merged.groupby('anno').frac_phastcons.median()
        print(f'{tis}\t{cre_med["Expressed"]}\t{cre_med["LoF Intolerant"]}\t{phast_med["Expressed"]}\t{phast_med["LoF Intolerant"]}')

        # print MWU stats
        _, cre_p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').enh_num, merged.query(f'anno=="Expressed"').enh_num)
        _, phast_p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').frac_phastcons, merged.query(f'anno=="Expressed"').frac_phastcons)
        _, frac_p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').frac_tisspec_enh, merged.query(f'anno=="Expressed"').frac_tisspec_enh)
        lrows.append([tis, 'lof_intolerant', cre_p, phast_p, frac_p])

    hdf = pd.DataFrame(hrows, columns=['tissue', 'genetype', 'cre_num_p', 'phastcons_p', 'frac_tisspec_p'])
    hdf = hdf.assign(cre_num_fdr=lambda x: sm.stats.multipletests(x.cre_num_p, alpha=0.05, method='fdr_bh')[1],
                     phastcons_fdr=lambda x: sm.stats.multipletests(x.phastcons_p, alpha=0.05, method='fdr_bh')[1],
                     frac_tisspec_fdr=lambda x: sm.stats.multipletests(x.frac_tisspec_p, alpha=0.05, method='fdr_bh')[1])
    hdf.to_csv(f'../res/housekeeping_{landscape_label}_df.out', sep='\t', index=False, header=True)

    ldf = pd.DataFrame(lrows, columns=['tissue', 'genetype', 'cre_num_p', 'phastcons_p', 'frac_tisspec_p'])
    ldf = ldf.assign(cre_num_fdr=lambda x: sm.stats.multipletests(x.cre_num_p, alpha=0.05, method='fdr_bh')[1],
                     phastcons_fdr=lambda x: sm.stats.multipletests(x.phastcons_p, alpha=0.05, method='fdr_bh')[1],
                     frac_tisspec_fdr=lambda x: sm.stats.multipletests(x.frac_tisspec_p, alpha=0.05, method='fdr_bh')[1])
    ldf.to_csv(f'../res/lof-intol_{landscape_label}_df.out', sep='\t', index=False, header=True)
    ### \\

    ### // write p value table, CRE > 0 \\ ###
    hrows = []
    lrows = []

    print('')
    print(f'tisue\tmed_cre_num_expressed\tmed_cre_num_housekeeping\tmed_frac_phastcons_expressed\tmed_frac_phastcons_housekeeping')
    for tis in tis_order:
        df = read_matched_data(tis=tis, gene_type='hk_gt0', landscape_type=landscape_def, chmm=CHROMHMM)
        merged = df.merge(all_tis_exp.query(f'tissue=="{tis}"'), how='left')

        # print medians
        cre_med = merged.groupby('anno').enh_num.median()
        phast_med = merged.groupby('anno').frac_phastcons.median()
        print(f'{tis}\t{cre_med["Expressed"]}\t{cre_med["Housekeeping"]}\t{phast_med["Expressed"]}\t{phast_med["Housekeeping"]}')

        # print MWU stats
        _, cre_p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').enh_num, merged.query(f'anno=="Expressed"').enh_num)
        _, phast_p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').frac_phastcons, merged.query(f'anno=="Expressed"').frac_phastcons)
        _, frac_p = stats.mannwhitneyu(merged.query(f'anno=="Housekeeping"').frac_tisspec_enh, merged.query(f'anno=="Expressed"').frac_tisspec_enh)
        hrows.append([tis, 'housekeeping', cre_p, phast_p, frac_p])

    print('')
    print(f'tisue\tmed_cre_num_expressed\tmed_cre_num_lof\tmed_frac_phastcons_expressed\tmed_frac_phastcons_lif')
    for tis in tis_order:
        df = read_matched_data(tis=tis, gene_type='lof_gt0', landscape_type=landscape_def, chmm=CHROMHMM)
        merged = df.merge(all_tis_exp.query(f'tissue=="{tis}"'), how='left')

        # print medians
        cre_med = merged.groupby('anno').enh_num.median()
        phast_med = merged.groupby('anno').frac_phastcons.median()
        print(f'{tis}\t{cre_med["Expressed"]}\t{cre_med["LoF Intolerant"]}\t{phast_med["Expressed"]}\t{phast_med["LoF Intolerant"]}')

        # print MWU stats
        _, cre_p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').enh_num, merged.query(f'anno=="Expressed"').enh_num)
        _, phast_p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').frac_phastcons, merged.query(f'anno=="Expressed"').frac_phastcons)
        _, frac_p = stats.mannwhitneyu(merged.query(f'anno=="LoF Intolerant"').frac_tisspec_enh, merged.query(f'anno=="Expressed"').frac_tisspec_enh)
        lrows.append([tis, 'lof_intolerant', cre_p, phast_p, frac_p])

    hdf = pd.DataFrame(hrows, columns=['tissue', 'genetype', 'cre_num_p', 'phastcons_p', 'frac_tisspec_p'])
    hdf = hdf.assign(cre_num_fdr=lambda x: sm.stats.multipletests(x.cre_num_p, alpha=0.05, method='fdr_bh')[1],
                     phastcons_fdr=lambda x: sm.stats.multipletests(x.phastcons_p, alpha=0.05, method='fdr_bh')[1],
                     frac_tisspec_fdr=lambda x: sm.stats.multipletests(x.frac_tisspec_p, alpha=0.05, method='fdr_bh')[1])
    hdf.to_csv(f'../res/housekeeping_{landscape_label}_gt0_df.out', sep='\t', index=False, header=True)

    ldf = pd.DataFrame(lrows, columns=['tissue', 'genetype', 'cre_num_p', 'phastcons_p', 'frac_tisspec_p'])
    ldf = ldf.assign(cre_num_fdr=lambda x: sm.stats.multipletests(x.cre_num_p, alpha=0.05, method='fdr_bh')[1],
                     phastcons_fdr=lambda x: sm.stats.multipletests(x.phastcons_p, alpha=0.05, method='fdr_bh')[1],
                     frac_tisspec_fdr=lambda x: sm.stats.multipletests(x.frac_tisspec_p, alpha=0.05, method='fdr_bh')[1])
    ldf.to_csv(f'../res/lof-intol_{landscape_label}_gt0_df.out', sep='\t', index=False, header=True)
    ### \\
