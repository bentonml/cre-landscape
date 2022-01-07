###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2021.05.27
#
#   assumes pwd is here with hydrogen:
#       /dors/capra_lab/users/bentonml/cross_species_gain_loss/bin/slurm/
#           .. EXP_DATA_PATH = '../../results/link_enh_to_genes/dat/2021-01-26'
#           .. RES_PATH = '../../results/link_enh_to_genes/fig/' + str(date.today())
###

import os
import logging
import pandas as pd
import numpy  as np
import seaborn as sns
import matplotlib.pyplot as plt

from pybedtools import BedTool
from datetime import date
from scipy import stats

plt.switch_backend('agg')  # add to save plots non-interactively

### // constants and paths \\ ###
EXP_DATA_PATH = '../dat/2021-01-26'      # relative to current dir
RES_PATH = f'../fig/{str(date.today())}'

# create a date stamped dir for files
if not os.path.isdir(RES_PATH):
    os.makedirs(RES_PATH)

tis_order = ['spleen','liver', 'heart_left_ventricle', 'brain_hippocampus', 'lung', 'pancreas',
             'brain_prefrontal_cortex', 'psoas_muscle', 'small_intestine', 'ovary']

fmt='pdf'
rc = {'font.size':14, 'axes.titlesize':16,'axes.labelsize':14, 'legend.fontsize': 14,
      'xtick.labelsize': 14, 'ytick.labelsize': 14, 'font.sans-serif':'Arial', 'font.family':'sans-serif'}


# set up log file to record data about run
logging.basicConfig(filename=f'{RES_PATH}/README', level=logging.INFO)
logging.info('Running on %s with peakachu loops.', str(date.today()))
### \\


### // functions \\ ###
def plot_dist_by_genetype_facetline(x, y, dat, tis, detail='', res_path='', fmt='pdf', rc=rc):
    ''' Plot line plot of the distribution of [y] stratified by gene type. CI is
    95% confidence interval from bootstrapping.
    output name: f'{res_path}/{str(date.today())}_{tis}_peakachuloop_{x}v{y}_{detail}bygenetype_lineplot.{fmt}'
    '''
    with sns.plotting_context("paper", rc=rc):
        fig, axes = plt.subplots(2, 2, figsize=(10,15), sharey='row')
        sns.lineplot(x=x, y=y, hue='hk_first',        data=dat, palette=['tab:gray', 'tab:green'], ax=axes[0,0])
        sns.lineplot(x=x, y=y, hue='essential_third', data=dat, palette=['tab:gray', 'tab:red'], ax=axes[0,1])
        sns.lineplot(x=x, y=y, hue='lof_intol_second', data=dat, palette=['tab:gray', 'tab:orange'], ax=axes[1,0])
        sns.lineplot(x=x, y=y, hue='tis_spec', data=dat, palette=['tab:gray', 'tab:blue'], ax=axes[1,1])
        sns.despine()
        plt.savefig(f'{res_path}/{tis}_peakachuloop_{x}v{y}_{detail}bygenetype_lineplot.{fmt}', format=fmt, dpi=400)
        plt.close()

def plot_dist_by_genetype_boxplot(x, y, dat, tis, detail='', res_path='', fmt='pdf', rc=rc):
    ''' Plot box plot of the distribution of [y] stratified by gene type.
    output name: f'{res_path}/{str(date.today())}_{tis}_peakachuloop_{x}v{y}_{detail}bygenetype_boxplot.{fmt}'
    '''
    with sns.plotting_context("paper", rc=rc):
        fig, axes = plt.subplots(2, 2, figsize=(10,15), sharey='row')
        sns.boxplot(x=x, y=y, hue='hk_first',        data=dat, palette=['tab:gray', 'tab:green'], ax=axes[0,0])
        sns.boxplot(x=x, y=y, hue='essential_third', data=dat, palette=['tab:gray', 'tab:red'], ax=axes[0,1])
        sns.boxplot(x=x, y=y, hue='lof_intol_second', data=dat, palette=['tab:gray', 'tab:orange'], ax=axes[1,0])
        sns.boxplot(x=x, y=y, hue='tis_spec',  data=dat, palette=['tab:gray', 'tab:blue'], ax=axes[1,1])
        sns.despine()
        plt.savefig(f'{res_path}/{tis}_peakachuloop_{x}v{y}_{detail}bygenetype_boxplot.{fmt}', format=fmt, dpi=400)
        plt.close()

def plot_dist_by_genetype_boxenplot(y, dat, tis, detail='', res_path='', fmt='pdf', rc=rc):
    ''' Plot boxenplot of the distribution of [y] stratified by gene type.
    output name: f'{res_path}/{str(date.today())}_{tis}_peakachuloop_{detail}bygenetype_boxenplot.{fmt}'
    '''
    with sns.plotting_context("paper", rc=rc):
        fig, axes = plt.subplots(2, 2, figsize=(10,15), sharey='row')
        sns.boxenplot(x='hk_first',        y=y, data=dat, palette=['tab:gray', 'tab:green'], ax=axes[0,0])
        sns.boxenplot(x='essential_third', y=y, data=dat, palette=['tab:gray', 'tab:red'], ax=axes[0,1])
        sns.boxenplot(x='lof_intol_second', y=y, data=dat, palette=['tab:gray', 'tab:orange'], ax=axes[1,0])
        sns.boxenplot(x='tis_spec',  y=y, data=dat, palette=['tab:gray', 'tab:blue'], ax=axes[1,1])
        sns.despine()
        plt.savefig(f'{res_path}/{tis}_peakachuloop_{detail}bygenetype_boxenplot.{fmt}', format=fmt, dpi=400)
        plt.close()
### \\


for tis in tis_order:
    logging.info('Running %s', tis)

    # create a tissue subdir within date stamped dir
    RES_TIS_PATH = f'{RES_PATH}/{tis}'
    if not os.path.isdir(RES_TIS_PATH):
        os.makedirs(RES_TIS_PATH)

    logging.info('Setting results directory to %s', RES_TIS_PATH)

    ### // read in dataframe \\ ###
    EXP_DATA_PATH = '../dat/2021-01-26'
    infile = f'{EXP_DATA_PATH}/{tis}_gene_with_enh-count_frac-enh-spec_peakachuloop-linked.tsv'
    logging.info('Read %s', infile)

    enh_num_by_gene = pd.read_table(infile)
    enh_num_by_gene = enh_num_by_gene.assign(rel_entropy_bins=lambda x: pd.cut(x['rel_entropy'], bins=10))
    enh_num_by_gene = enh_num_by_gene.assign(enh_rel_entropy_bins=lambda x: pd.cut(x['enh_rel_entropy'], bins=5))
    enh_num_by_gene = enh_num_by_gene.assign(exp_nocat=lambda x: np.where((x.exp == 1) & (x.hk == 0) & (x.lof_intol == 0) & (x.essential == 0), 1, 0))
    enh_num_by_gene = enh_num_by_gene.assign(hk_first=lambda x: np.where((x.hk == 1), 1, 0))
    enh_num_by_gene = enh_num_by_gene.assign(lof_intol_second=lambda x: np.where((x.hk == 0) & (x.lof_intol == 1), 1, 0))
    enh_num_by_gene = enh_num_by_gene.assign(essential_third=lambda x: np.where((x.hk == 0) & (x.lof_intol == 0) & (x.essential == 1), 1, 0))
    ### \\

    print(enh_num_by_gene.head())


    ### // plot distribution of CRE number and entropy values \\ ###
    logging.info('Plot enh num distribution for %s', infile)

    with sns.plotting_context("paper", rc=rc):
        ### fig: distribution of CRE number per gene
        sns.displot(x='enh_num', data=enh_num_by_gene, kind='hist')
        sns.despine()
        plt.xlabel('Number of CREs')
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_enh-num_distplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot of CRE number per gene, hue exp/not, no fliers
        plt.figure(figsize=(8,8))
        sns.boxplot(x='exp', y='enh_num', data=enh_num_by_gene.query('enh_num>0'), palette=['tab:gray', 'tab:blue'], showfliers=False)
        sns.despine()

        # difference of medians
        ts, p = stats.mannwhitneyu(enh_num_by_gene.query('enh_num>0 & exp==0').enh_num, enh_num_by_gene.query('enh_num>0 & exp==1').enh_num)
        yloc = enh_num_by_gene.query('enh_num>0').enh_num.describe()['75%'] * 2.1
        plt.text(-0.4, yloc, f'p = {p:.3e} (MWU)', horizontalalignment='left', size='medium', color='black')

        plt.xlabel('Expressed')
        plt.ylabel('Number of CREs')
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_enh-num_gt0_boxplot_nofliers.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

    logging.info('Plot entropy distributions for %s', infile)

    with sns.plotting_context("paper", rc=rc):
        ### fig: distribution of gene entropy values
        sns.displot(x='rel_entropy', data=enh_num_by_gene, kind='hist')
        sns.despine()
        plt.xlabel('Relative entropy')
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_gene-entropy_distplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        enh_num_by_gene[enh_num_by_gene.rel_entropy.notna()]
        ### fig: distribution of gene entropy per gene, hue exp/not
        sns.displot(x='rel_entropy', data=enh_num_by_gene[enh_num_by_gene.rel_entropy.notna()], hue='exp',
                    hue_order=[1, 0], palette=['tab:blue', 'tab:gray'], kind='kde')
        sns.despine()
        plt.xlabel('Relative entropy')
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_gene-entropy_byexp_kdeplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: distribution of average CRE entropy per gene
        sns.displot(x='enh_rel_entropy', data=enh_num_by_gene[enh_num_by_gene.rel_entropy.notna()], kind='hist')
        sns.despine()
        plt.xlabel('Relative entropy')
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_enh-entropy_distplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: distribution of average CRE entropy per gene, hue exp/not
        sns.displot(x='enh_rel_entropy', data=enh_num_by_gene[enh_num_by_gene.enh_rel_entropy.notna()], hue='exp',
                    hue_order=[1, 0], palette=['tab:blue', 'tab:gray'], kind='kde')
        sns.despine()
        plt.xlabel('Relative entropy')
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_enh-entropy_byexp_kdeplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig
    ### \\


    ### // plot distribution of CRE number by tss bin \\ ###
    logging.info('Plot enh num distribution by tss bin for %s', infile)
    with sns.plotting_context("paper", rc=rc):
        ### fig: box plot of CRE number by tss bin, hue exp/not, no fliers
        plt.figure(figsize=(8,8))
        sns.boxplot(x='tss_bins', y='enh_num', hue='exp', data=enh_num_by_gene, palette=['tab:gray', 'tab:blue'], showfliers=False)
        sns.despine()
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_tssvenh_binned_boxplot_nofliers.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig
    ### \\


    ### // plot gene entropy by CRE number \\ ###
    logging.info('Plot gene entropy by CRE number for %s', infile)

    with sns.plotting_context("paper", rc=rc):
        ### fig: scatterplot of CRE num by relative entropy bin
        plt.figure(figsize=(12,8))
        sns.scatterplot(x='rel_entropy', y='enh_num', data=enh_num_by_gene, alpha=0.1)
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_rel-entropyvenh-num_scatter.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot of CRE num by relative entropy bin (10), exp and not
        plt.figure(figsize=(12,8))
        ax = sns.boxplot(x='rel_entropy_bins', y='enh_num', data=enh_num_by_gene.query('enh_num>0'), width=.75, showfliers=False, color='#866f85')
        ax = sns.swarmplot(x='rel_entropy_bins', y='enh_num', data=enh_num_by_gene.query('enh_num>0'), color='.3')
        _ = ax.set_xticklabels(['(0,0.1]', '(0.1,0.2]', '(0.2,0.3]', '(0.3,0.4]', '(0.4,0.5]',
                                '(0.5,0.6]', '(0.6,0.7]', '(0.7,0.8]', '(0.8,0.9]', '(0.9,1.0]'])
        plt.xlabel('Relative entropy')
        plt.ylabel('Number of CREs')
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_rel-entropyvenh-numenh_gt0_expandnot_boxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot of CRE num by relative entropy bin (10), exp only
        plt.figure(figsize=(12,8))
        ax = sns.violinplot(x='rel_entropy_bins', y='enh_num', data=enh_num_by_gene.query('enh_num>0 & exp==1'), color='tab:blue', scale='area')
        _ = ax.set_xticklabels(['(0,0.1]', '(0.1,0.2]', '(0.2,0.3]', '(0.3,0.4]', '(0.4,0.5]',
                                '(0.5,0.6]', '(0.6,0.7]', '(0.7,0.8]', '(0.8,0.9]', '(0.9,1.0]'])
        plt.xlabel('Relative entropy')
        plt.ylabel('Number of CREs')
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_rel-entropyvenh-num_gt0_exponly_boxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: violinplot of CRE num by relative entropy bin (10), exp only
        plt.figure(figsize=(12,8))
        ax = sns.violinplot(x='rel_entropy_bins', y='enh_num', data=enh_num_by_gene.query('enh_num>0 & exp==1'), color='tab:blue', scale='area')
        _ = ax.set_xticklabels(['(0,0.1]', '(0.1,0.2]', '(0.2,0.3]', '(0.3,0.4]', '(0.4,0.5]',
                                '(0.5,0.6]', '(0.6,0.7]', '(0.7,0.8]', '(0.8,0.9]', '(0.9,1.0]'])
        plt.xlabel('Relative entropy')
        plt.ylabel('Number of CREs')
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_rel-entropyvenh-num_gt0_exponly_violinplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: scatter lmplot of gene entropy vs CRE number, hue exp/not
        sns.lmplot(x='rel_entropy', y='enh_num', hue='exp', col='exp', data=enh_num_by_gene,
                   scatter_kws={'alpha': 0.1}, palette=['tab:gray', 'tab:blue'], truncate=True)
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_tis-spec_lmplot_order1_scatter.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: binned lmplot of gene entropy vs CRE number, hue exp/not
        sns.lmplot(x='rel_entropy', y='enh_num', hue='exp', col='exp', data=enh_num_by_gene,
                   palette=['tab:gray', 'tab:blue'], x_bins=10, truncate=True)
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_tis-spec_lmplot_order1_binned_equalobs.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: binned lmplot of gene entropy vs CRE number, hue exp/not, eq spaced
        sns.lmplot(x='rel_entropy', y='enh_num', hue='exp', col='exp', data=enh_num_by_gene,
                   palette=['tab:gray', 'tab:blue'], x_bins=[.1,.3,.5,.7,.9], truncate=True)
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_tis-spec_lmplot_order1_binned_equaldist.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig
    ### \\


    ### // plot CRE number distributions by gene type \\ ###
    logging.info('Plot CRE number by gene type for %s', infile)

    with sns.plotting_context("paper", rc=rc):
        ### fig: boxplot of CRE number by exp/not, hue tis_spec/not for gene
        plt.figure(figsize=(8,8))
        sns.boxplot(x='exp', y='enh_num', hue='tis_spec', data=enh_num_by_gene, palette=['tab:gray','#aec7e8'])
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_enh-numvtisspec-gene_expandnot_boxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

    #### figs: lineplots of CRE number by gene type, exp only and exp/not, by tss
    #plot_dist_by_genetype_facetline(x='tss', y='enh_num', dat=enh_num_by_gene, tis=tis, res_path=RES_TIS_PATH, fmt=fmt)
    #plot_dist_by_genetype_facetline(x='tss_bins', y='enh_num', dat=enh_num_by_gene, tis=tis, detail='binned_', res_path=RES_TIS_PATH, fmt=fmt)
    #plot_dist_by_genetype_facetline(x='tss', y='enh_num', dat=enh_num_by_gene.query('exp > 0'), tis=tis, detail='exp_only_', res_path=RES_TIS_PATH, fmt=fmt)
    #plot_dist_by_genetype_facetline(x='tss_bins', y='enh_num', dat=enh_num_by_gene.query('exp > 0'), tis=tis, detail='exp_only_binned_', res_path=RES_TIS_PATH, fmt=fmt)
    #### end_fig

    #### figs: boxplots of CRE number by gene type, exp only and exp/not, by tss
    #plot_dist_by_genetype_boxplot(x='tss', y='enh_num', dat=enh_num_by_gene, tis=tis, res_path=RES_TIS_PATH, fmt=fmt)
    #plot_dist_by_genetype_boxplot(x='tss_bins', y='enh_num', dat=enh_num_by_gene, tis=tis, detail='binned_', res_path=RES_TIS_PATH, fmt=fmt)
    #plot_dist_by_genetype_boxplot(x='tss', y='enh_num', dat=enh_num_by_gene.query('exp > 0'), tis=tis, detail='exp_only_', res_path=RES_TIS_PATH, fmt=fmt)
    #plot_dist_by_genetype_boxplot(x='tss_bins', y='enh_num', dat=enh_num_by_gene.query('exp > 0'), tis=tis, detail='exp_only_binned_', res_path=RES_TIS_PATH, fmt=fmt)
    #### end_fig
    ### \\


    ### // plot osterwalder18 style boxplots ### \\
    logging.info('Creating enh_num_by_gene_anno dataframe for %s', infile)
    enh_by_gene_anno = (enh_num_by_gene.filter(['target_gene', 'exp', 'hk', 'lof_intol', 'essential', 'tis_spec', 'exp_nocat', 'hk_first', 'lof_intol_second', 'essential_third', 'enh_num'])
                                        .assign(expressed=lambda x: x.exp)
                                        .set_index(['target_gene', 'enh_num', 'exp'])
                                        .stack()
                                        .reset_index()
                                        .rename(columns={'level_3':'anno', 0:'val'})
                                        .query('val == 1')
                                        .merge(enh_num_by_gene.filter(['target_gene', 'enh_rel_entropy', 'frac_tisspec_enh', 'tss']), how='inner', validate='m:1'))

    logging.info('Plot Osterwalder18 style CRE number boxplots for %s', infile)

    with sns.plotting_context("paper", rc=rc):
        ### fig: boxplot for CRE num by gene type, both exp/not
        plt.figure(figsize=(10,7))
        sns.boxplot(x='enh_num', y='anno', data=enh_by_gene_anno,
                    palette=['tab:green', 'tab:red', 'tab:orange', 'tab:blue'],
                    order=['hk_first', 'essential_third', 'lof_intol_second', 'exp_nocat'], orient='h', fliersize=1)
        sns.despine()
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_enh-numbygenetype_expandnon_boxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot for CRE num by gene type, exp only
        plt.figure(figsize=(10,7))
        sns.boxplot(x='enh_num', y='anno', data=enh_by_gene_anno.query('exp > 0'),
                    palette=['tab:green', 'tab:red', 'tab:orange', 'tab:blue'],
                    order=['hk_first', 'essential_third', 'lof_intol_second', 'exp_nocat'], orient='h', fliersize=1)
        sns.despine()
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_enh-numbygenetype_exponly_boxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot for CRE num by gene type, tis-spec, exp only
        plt.figure(figsize=(10,7))
        sns.boxplot(x='enh_num', y='anno', data=enh_by_gene_anno.query('exp > 0'),
                    palette=['#aec7e8', 'tab:blue'],
                    order=['tis_spec', 'expressed'], orient='h', fliersize=1)
        sns.despine()
        plt.xlabel('Number of CREs')
        plt.ylabel('')
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_enh-numbygenetype_exponly_justtis-spec_boxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot for CRE num by gene type, exp only, no fliers
        plt.figure(figsize=(10,7))
        sns.boxplot(x='enh_num', y='anno', data=enh_by_gene_anno.query('exp > 0 & enh_num > 0'),
                    palette=['tab:green', 'tab:red', 'tab:orange', 'tab:blue'],
                    order=['hk_first', 'essential_third', 'lof_intol_second', 'exp_nocat'],
                    orient='h', fliersize=1, showfliers=False)
        sns.despine()
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_enh-numbygenetype_exponly_boxplot_enhgt0_nofliers.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

    logging.info('Plot Osterwalder18 style CRE entropy boxplots for %s', infile)

    with sns.plotting_context("paper", rc=rc):
        ### fig: boxplot for CRE entropy by gene type, exp only
        plt.figure(figsize=(10,7))
        sns.boxplot(x='enh_rel_entropy', y='anno', data=enh_by_gene_anno.query('exp > 0'),
                    palette=['tab:green', 'tab:red', 'tab:orange', 'tab:blue'],
                    order=['hk_first', 'essential_third', 'lof_intol_second', 'exp_nocat'], orient='h', fliersize=2)
        sns.despine()
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_enh-entropybygenetype_exponly_boxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot for CRE entropy by gene type, exp only, with swarm
        plt.figure(figsize=(10,7))
        sns.boxplot(x='enh_rel_entropy', y='anno', data=enh_by_gene_anno.query('exp > 0'),
                    palette=['tab:green', 'tab:red', 'tab:orange', 'tab:blue'],
                    order=['hk_first', 'essential_third', 'lof_intol_second', 'exp_nocat'], orient='h', showfliers=False)
        sns.swarmplot(x='enh_rel_entropy', y='anno', data=enh_by_gene_anno.query('exp > 0'), size=3, color='.3',
                      linewidth=0, order=['hk_first', 'essential_third', 'lof_intol_second', 'tis_spec', 'exp_nocat'])
        sns.despine()
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_enh-entropybygenetype_exponly_swarmboxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

    logging.info('Plot Osterwalder18 style fraction tissue-specific CRE boxplots for %s', infile)

    with sns.plotting_context("paper", rc=rc):
        ### fig: boxplot for frac tissue specific CREs by gene type, exp only
        plt.figure(figsize=(10,7))
        sns.boxplot(x='frac_tisspec_enh', y='anno', data=enh_by_gene_anno.query('exp > 0'),
                    palette=['tab:green', 'tab:red', 'tab:orange', 'tab:blue'],
                    order=['hk_first', 'essential_third', 'lof_intol_second', 'exp_nocat'], orient='h', fliersize=2)
        sns.despine()
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_frac-spec-enhbygenetype_exponly_boxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot for frac tissue specific CREs by gene type, exp only, with swarm
        plt.figure(figsize=(10,7))
        sns.boxplot(x='frac_tisspec_enh', y='anno', data=enh_by_gene_anno.query('exp > 0'),
                    palette=['tab:green', 'tab:red', 'tab:orange', 'tab:blue'],
                    order=['hk_first', 'essential_third','lof_intol_second', 'exp_nocat'], orient='h', showfliers=False)
        sns.swarmplot(x='frac_tisspec_enh', y='anno', data=enh_by_gene_anno.query('exp > 0'), size=3, color='.3',
                      linewidth=0, order=['hk_first', 'essential_third', 'lof_intol_second', 'tis_spec', 'exp_nocat'])
        sns.despine()
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_frac-spec-enhbygenetype_exponly_swarmboxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot for CRE num by gene type, tis-spec, exp only
        plt.figure(figsize=(10,7))
        sns.boxplot(x='frac_tisspec_enh', y='anno', data=enh_by_gene_anno.query('exp > 0'),
                    palette=['#aec7e8', 'tab:blue'],
                    order=['tis_spec', 'expressed'], orient='h', fliersize=1)
        sns.despine()
        plt.xlabel('Fraction of tissue-specific CREs')
        plt.ylabel('')
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_frac-spec-enhbygenetype_exponly_justtis-spec_boxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot for frac tissue specific CREs by gene type, exp only, 1 tss
        plt.figure(figsize=(10,7))
        sns.boxplot(x='frac_tisspec_enh', y='anno', data=enh_by_gene_anno.query('exp > 0 & tss == 1'),
                    palette=['tab:green', 'tab:red', 'tab:orange', 'tab:blue'],
                    order=['hk_first', 'essential_third', 'lof_intol_second', 'exp_nocat'], orient='h', showfliers=False)
        sns.despine()
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_frac-spec-enhbygenetype_exponly_1tss_boxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig

        ### fig: boxplot for frac tissue specific CREs by gene type, exp only, 1 tss, with swarm
        plt.figure(figsize=(10,7))
        sns.boxplot(x='frac_tisspec_enh', y='anno', data=enh_by_gene_anno.query('exp > 0 & tss == 1'),
                    palette=['tab:green', 'tab:red', 'tab:orange', 'tab:blue'],
                    order=['hk_first', 'essential_third', 'lof_intol_second', 'exp_nocat'], orient='h', showfliers=False)
        sns.swarmplot(x='frac_tisspec_enh', y='anno', data=enh_by_gene_anno.query('exp > 0 & tss == 1'), size=3, color='.3',
                      linewidth=0, order=['hk_first', 'essential_third', 'lof_intol_second', 'exp_nocat'])
        sns.despine()
        plt.savefig(f'{RES_TIS_PATH}/{tis}_peakachuloop_frac-spec-enhbygenetype_exponly_1tss_swarmboxplot.{fmt}', format=fmt, dpi=400)
        plt.close()
        ### end_fig
    ### \\
