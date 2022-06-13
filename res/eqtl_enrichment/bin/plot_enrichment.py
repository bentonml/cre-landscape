###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2020.07.21
#   updated   | 2020.07.22
#
#   assumes pwd is here with hydrogen:
#       /dors/capra_lab/users/bentonml/cross_species_gain_loss/bin/slurm/
###

import os
import pandas as pd
import numpy  as np
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import date
from scipy.stats import ttest_ind, mannwhitneyu


### // constants, paths, functions \\ ###
RES_PATH = '../res'
FIG_PATH = '../fig'
DATE = '2022-06-13'

# create a date stamped dir for figures
if not os.path.isdir(FIG_PATH):
    os.makedirs(FIG_PATH)

tissues = ['spleen','liver', 'heart_left_ventricle', 'brain_hippocampus', 'lung', 'pancreas',
           'brain_prefrontal_cortex', 'psoas_muscle', 'small_intestine', 'ovary']
tis_names = ['Spleen','Liver', 'Heart', 'Hippocampus', 'Lung', 'Pancreas',
             'Prefrontal cortex', 'Muscle', 'Small intestine', 'Ovary']

fmt='pdf'
rc = {'font.size':14, 'axes.titlesize':16,'axes.labelsize':14, 'legend.fontsize': 14,
      'xtick.labelsize': 14, 'ytick.labelsize': 14}
### \\


### // functions \\ ###
def read_enrichment_output_quartile(tissue, landscape_type, date, path):
    with open(f'{path}/{date}_{tissue}_eqtl_{landscape_type}_enrichment.out', 'r') as res_file:
        res = []

        for q in ['q1', 'q2', 'q3', 'q4']:
            next(res_file)  # skip command line
            next(res_file)  # skip header line
            #next(res_file)  # skip iterations line
            line = res_file.readline().strip('\n').split('\t')
            res.append([tissue, q] + line)
            # next(res_file)  # skip the timing lines
            # next(res_file)
            # next(res_file)
            # next(res_file)
    return res


def create_enrichment_df_quartile(landscape_type, date=DATE, path=RES_PATH, tissues=tissues):
    res = []
    for tis in tissues:
        d = read_enrichment_output_quartile(tissue=tis, landscape_type=landscape_type, date=date, path=path)
        res += d
    return pd.DataFrame(res, columns=['tissue', 'landscape_quartile',
                                      'obs', 'exp', 'std_dev', 'fold_change', 'p'])


def create_enrichment_df_full(tissues, landscape_type, date, path):
    res = []
    for tis in tissues:
        with open(f'{path}/{date}_{tis}_eqtl_{landscape_type}_enrichment.out', 'r') as res_file:
            next(res_file)  # skip command line
            next(res_file)  # skip header line
            line = res_file.readline().strip('\n').split('\t')
            res.append([tis] + line)
    return pd.DataFrame(res, columns=['tissue', 'obs', 'exp', 'std_dev', 'fold_change', 'p'])


def format_enrichment_df(dd):
    dd['obs'] = pd.to_numeric(dd.obs)
    dd['exp'] = pd.to_numeric(dd.exp)
    dd['std_dev'] = pd.to_numeric(dd.std_dev)
    dd['fold_change'] = pd.to_numeric(dd.fold_change)
    dd['p'] = pd.to_numeric(dd.p)
    num_tests = dd.shape[0]

    dd = dd.assign(bonf_signif=lambda x: np.where(x.p < (0.05/num_tests), '*', 'NS'),
                   signif=lambda x: np.where(x.p < 0.05, '*', 'NS'),
                   log2FoldChange=lambda x: np.log2(x.fold_change))
    return dd


def read_counts(filename):
    with open(filename, 'r') as infile:
        obs = int(infile.readline().strip())

        # expects data in counts format: [ obs ] [ exp \t exp \t ... ]
        exp_list = infile.readline().strip().split('\t')

    return obs, np.array(exp_list, dtype=int)


def calculate_fc_distribution(obs, exp_list):
    return np.array([(obs + 1.0) / (exp + 1.0) for exp in exp_list])


def calculate_empirical_p(obs_fc_a, fc_list_a, obs_fc_b, fc_list_b):
    delta_final = obs_fc_a - obs_fc_b
    delta_dist  = [fc_list_a[i] - fc_list_b[i] for i in range(len(fc_list_a))]

    # sum number of differences >= to observed difference
    p_sum = sum(1 for delta in delta_dist if abs(delta) >= abs(delta_final))

    return (p_sum + 1.0) / (len(delta_dist) + 1.0)


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


def add_hline(axes):
    for ax in axes.flatten():
        ax.axhline(0, c='k')

def add_vline(axes):
    for ax in axes.flatten():
        ax.axvline(0, c='k')
### \\

for landscape in ['peakachuloop', 'hicQ05']:
    # plot for CRE enrichment without quartile, all tissues on same plot
    eqtl_full = create_enrichment_df_full(tissues, landscape_type=landscape, date=DATE, path=RES_PATH)
    eqtl_full = format_enrichment_df(eqtl_full)

    # barplot of eQTL enrichment with CREs, all tissues
    with sns.plotting_context("paper", rc=rc):
        g = sns.catplot(x='log2FoldChange', y='tissue', data=eqtl_full, kind='bar',
                        hue='bonf_signif', hue_order=['*', 'NS'], dodge=False,
                        palette=['#2171b5', '#bdd7e7'])
        g.set(xlim=(-0.1, None), ylabel='', xlabel=r'log$_2$(Fold Change)')
        g.set_yticklabels(tis_names)
        g._legend.set_title('Significance')
        add_vline(g.axes)
        plt.savefig(f'{FIG_PATH}/{str(date.today())}_{landscape}_eqtl_full_barplot.{fmt}', format=fmt, dpi=400)
        plt.close()

    expected_counts = []
    for tis in tissues:
        obs, exp_list = read_counts(f'../dat/{tis}_{landscape}_cre_full.counts')
        c = pd.DataFrame(exp_list, columns=['expected_count'])
        c['tissue'] = tis
        expected_counts.append(c[c.expected_count >= 0])
    count_df = pd.concat(expected_counts)

    # boxplot of eQTL enrichment with CREs, all tissues (showing expected counts)
    with sns.plotting_context("paper", rc=rc):
        g = sns.boxplot(x='expected_count', y='tissue', data=count_df, showfliers=False, color='.5')
        sns.scatterplot(x='obs', y='tissue', data=eqtl_full, marker='*', s=100, color='tab:red')
        g.set(ylabel='', xlabel=r'Overlap (bp)')
        g.set_yticklabels(tis_names)
        sns.despine()
        plt.savefig(f'{FIG_PATH}/{str(date.today())}_{landscape}_eqtl_full_boxplot_with_counts.{fmt}', format=fmt, dpi=400)
        plt.close()

    # read in data and generate table
    eqtl_df = create_enrichment_df_quartile(landscape, date='2022-06-01', path=RES_PATH)
    eqtl_df = format_enrichment_df(eqtl_df)
    eqtl_df['quartile'] = eqtl_df.landscape_quartile.map({'q1':1, 'q2':2, 'q3':3, 'q4':4})

    # bar plot of enrichment by quartile, hue for significance
    with sns.plotting_context("paper", rc=rc):
        g = sns.catplot(x='landscape_quartile', y='log2FoldChange', data=eqtl_df, kind='bar',
                        col='tissue', col_wrap=5, col_order=tissues, hue='bonf_signif',
                        hue_order=['*', 'NS'], dodge=False, palette=['#2171b5', '#bdd7e7'],
                        sharey=False, sharex=False)
        g.fig.subplots_adjust(hspace=0.5, wspace=.25)
        facet_titles(g.axes)
        facet_x_axis(g.axes, 'CRE landscape quartile')
        g.set(ylabel=r'log$_2$(Fold Change)', ylim=(-0.1, None))
        g.set_xticklabels(['1', '2', '3', '4'])
        g._legend.set_title('')
        add_hline(g.axes)
        plt.savefig(f'{FIG_PATH}/{str(date.today())}_{landscape}_eqtl_barplot.{fmt}', format=fmt, dpi=400)
        plt.close()

    expected_counts = []
    for tis in tissues:
        for q in range(1,5):
            obs, exp_list = read_counts(f'../dat/{tis}_{landscape}_cre_quart{q}.counts')
            c = pd.DataFrame(exp_list, columns=['expected_count'])
            c['observed'] = obs
            c['tissue'] = tis
            c['quartile'] = q
            expected_counts.append(c[c.expected_count >= 0])
    count_df = pd.concat(expected_counts)
    count_df['log2FC'] = np.log2(count_df.observed / count_df.expected_count)

    # boxplot for enrichment showing all fold changes
    with sns.plotting_context("paper", rc=rc):
        g = sns.catplot(x='quartile', y='log2FC', col='tissue', col_wrap=5, data=count_df,
                    palette='Greys', kind='box', sharey=False, sharex=False, showfliers=False)
        for i, ax in enumerate(g.axes.flatten()):
            min_fc = min(count_df.query(f'tissue=="{tissues[i]}"').log2FC)
            ax.set_ylim((min(0, min_fc-0.1), None))
        g.fig.subplots_adjust(hspace=0.25, wspace=.25)
        g.set(ylabel=r'log$_2$(Fold Change)', xlabel='CRE landscape quartile')
        facet_titles(g.axes)
        plt.savefig(f'{FIG_PATH}/{str(date.today())}_{landscape}_eqtl_quartile_boxplot.{fmt}', format=fmt, dpi=400)
        plt.close()

    # violin plot for enrichment with CREs, all tissues (showing expected counts)
    with sns.plotting_context("paper", rc=rc):
        g = sns.catplot(x='quartile', y='expected_count', col='tissue', col_wrap=5, data=count_df, color='.5', kind='violin', sharey=False)
        for i, ax in enumerate(g.axes.flatten()):
            sns.scatterplot(x=eqtl_df.query(f'tissue=="{tissues[i]}"').quartile-1, y='obs',
                            data=eqtl_df.query(f'tissue=="{tissues[i]}"'), marker='*', s=100, color='tab:red', ax=ax)
        g.set(ylabel='Overlap (bp)', xlabel='CRE landscape quartile')
        facet_titles(g.axes)
        plt.savefig(f'{FIG_PATH}/{str(date.today())}_{landscape}_eqtl_quartile_violinplot_with_counts.{fmt}', format=fmt, dpi=400)
        plt.close()

    for tis in tissues:
        obs_a, exp_list_a = read_counts(f'../dat/{tis}_{landscape}_cre_quart1.counts')
        obs_b, exp_list_b = read_counts(f'../dat/{tis}_{landscape}_cre_quart4.counts')

        fc_list_a = calculate_fc_distribution(obs_a, exp_list_a)
        fc_list_b = calculate_fc_distribution(obs_b, exp_list_b)

        # ts, p = mannwhitneyu(fc_list_a, fc_list_b, equal_var=False)
        ts, p = mannwhitneyu(fc_list_a, fc_list_b, alternative="greater")
        print(f'{tis}: {np.median(fc_list_a):.2f} v. {np.median(fc_list_b):.2f}  {ts} (p = {p:.3e})')
