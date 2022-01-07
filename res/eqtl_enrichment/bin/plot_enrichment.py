###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2020.07.21
#   updated   | 2020.07.22
#
#   assumes pwd is here with hydrogen:
#       /dors/capra_lab/users/bentonml/cross_species_gain_loss/bin/slurm/
###

import pandas as pd
import numpy  as np
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import date


### // constants, paths, functions \\ ###
RES_PATH = '../res'
FIG_PATH = '../fig'

# create a date stamped dir for figures
if not os.path.isdir(FIG_PATH):
    os.makedirs(FIG_PATH)

tissues = ['spleen','liver', 'heart_left_ventricle', 'brain_hippocampus', 'lung', 'pancreas',
           'brain_prefrontal_cortex', 'psoas_muscle', 'small_intestine', 'ovary']

fmt='pdf'
rc = {'font.size':14, 'axes.titlesize':16,'axes.labelsize':14, 'legend.fontsize': 14,
      'xtick.labelsize': 14, 'ytick.labelsize': 14}
### \\


### // functions \\ ###
def read_enrichment_output(tissue, landscape_type, path=RES_PATH):
    #TODO: update path date
    with open(f'{path}/2021-12-17_{tissue}_eqtl_{landscape_type}_enrichment.out', 'r') as res_file:
        res = []

        for q in ['q1', 'q2', 'q3', 'q4']:
            next(res_file)  # skip command line
            next(res_file)  # skip header line
            next(res_file)  # skip iterations line
            line = res_file.readline().strip('\n').split('\t')
            res.append([tissue, q] + line)
            next(res_file)  # skip the timing lines
            next(res_file)
            next(res_file)
            next(res_file)
    return res

def create_enrichment_df(landscape_type, path=RES_PATH):
    res = []
    for tis in tissues:
        d = read_enrichment_output(tis, landscape_type, path)
        res += d
    return pd.DataFrame(res, columns=['tissue', 'landscape_quartile',
                                      'obs', 'exp', 'std_dev', 'fold_change', 'p'])

def format_enrichment_df(dd):
    dd['obs'] = pd.to_numeric(dd.obs)
    dd['exp'] = pd.to_numeric(dd.exp)
    dd['std_dev'] = pd.to_numeric(dd.std_dev)
    dd['fold_change'] = pd.to_numeric(dd.fold_change)
    dd['p'] = pd.to_numeric(dd.p)

    dd = dd.assign(bonf_signif=lambda x: np.where(x.p < (0.05/40), '*', 'NS'),
                   signif=lambda x: np.where(x.p < 0.05, '*', 'NS'),
                   log2FoldChange=lambda x: np.log2(x.fold_change))
    return dd

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
### \\

for landscape in ['peakachuloop', 'hicQ05']:
    # read in data and generate table
    eqtl_df = create_enrichment_df(landscape, RES_PATH)
    eqtl_df = format_enrichment_df(eqtl_df)

    with sns.plotting_context("paper", rc=rc):
        g = sns.catplot(x='landscape_quartile', y='log2FoldChange', data=eqtl_df, kind='bar',
                        col='tissue', col_wrap=5, col_order=tissues, hue='bonf_signif',
                        hue_order=['*', 'NS'], dodge=False, palette=['#2171b5', '#bdd7e7'],
                        sharex=False)
        g.fig.subplots_adjust(hspace=0.5, wspace=.15)
        facet_titles(g.axes)
        facet_x_axis(g.axes, 'CRE Landscape Quartile')
        g.set_xticklabels(['1', '2', '3', '4'])
        g._legend.set_title('')
        add_hline(g.axes)
        plt.savefig(f'{FIG_PATH}/{str(date.today())}_{landscape}_eqtl_barplot.{fmt}', format=fmt, dpi=400)
        plt.close()

