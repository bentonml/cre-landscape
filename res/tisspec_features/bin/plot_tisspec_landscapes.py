###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2022.02.18
#
#   manuscript plots for tis_spec matching output, both landscape types
###

import os
import logging
import pandas as pd
import numpy  as np
import seaborn as sns
import matplotlib.pyplot as plt

from datetime import date
from scipy import stats

plt.switch_backend('agg')  # add to save plots non-interactively

EXP_DAT_PATH = f'../../link_cre_to_genes/dat/2022-03-08'
MAT_DAT_PATH = f'../dat/2022-03-08'
RES_PATH = f'../fig/{str(date.today())}'

# create a date stamped dir for files
if not os.path.isdir(RES_PATH):
    os.makedirs(RES_PATH)

# set up log file to record data about run
logging.basicConfig(filename=f'{RES_PATH}/README', level=logging.INFO)

tisspec_thresh = 0.3
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
        enh_num_by_gene = enh_num_by_gene.assign(tis_spec=lambda x: np.where(x.rel_entropy > tisspec_thresh, 1, 0))
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

def read_matched_data(tis, landscape_type, dtype):
    return pd.read_table(f'{MAT_DAT_PATH}/matched_{dtype}{tis}_{landscape_type}.tsv', sep='\t')
### \\

for cretype in ['', 'gt0_']:
    for landscape_def in ['loop', 'contact']:
        # create dataframe of enhancer number by gene with all tissues
        all_tis = read_data(landscape_def)
        all_tis_exp = all_tis[all_tis['exp']==1]

        ### // full matched sets \\ ###
        df_lst = []
        for tis in tis_order:
            dd = (read_matched_data(tis=tis, landscape_type=landscape_def, dtype=cretype)
                    .merge(all_tis_exp.query(f'tissue=="{tis}"'), how='left')
                    .filter(['tissue', 'anno', 'gtex_exp_log2', 'target_gene', 'enh_num', 'frac_tisspec_enh', 'frac_phastcons'])
                 )
            df_lst.append(dd)
        all_tis_merged = pd.concat(df_lst)

        logging.info(all_tis_merged.groupby('tissue').anno.value_counts())
        ### \\

        with sns.plotting_context("paper", rc=rc):
            ### fig: H boxplot of all tissues by frac-tissue-spec enhancers, hue = exp type
            fig = plt.figure(figsize=(10,6))
            g = sns.boxplot(data=all_tis_merged,
                            y='frac_tisspec_enh', x='tissue', hue='anno', order=tis_order,
                            palette={'Broad':'#fee5d9', 'Tissue-specific':'#de2d26'},
                            showfliers=False)
            g = sns.stripplot(data=all_tis_merged, x='tissue', y='frac_tisspec_enh',
                              hue='anno', dodge=True, order=tis_order, palette=['.3', '.3'])
            g.set_xticklabels(tis_names, rotation=30, horizontalalignment='right')
            g.set_ylabel('Fraction tissue-specific CREs')
            g.set_xlabel('')
            g.legend(handles=g.legend_.legendHandles[0:2], frameon=False, bbox_to_anchor=(1, 1.15))
            sns.despine()
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/all_{landscape_def}_{cretype}frac-tisspec-enhXtisspec_bytissue_boxplot_withpoints.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: H boxplot of all tissues by frac-tissue-spec enhancers, hue = exp type
            fig = plt.figure(figsize=(10,6))
            g = sns.boxplot(data=all_tis_merged,
                            y='frac_tisspec_enh', x='tissue', hue='anno', order=tis_order,
                            palette={'Broad':'#fee5d9', 'Tissue-specific':'#de2d26'},
                            showfliers=False)
            g.set_xticklabels(tis_names, rotation=30, horizontalalignment='right')
            g.set_ylabel('Fraction tissue-specific CREs')
            g.set_xlabel('')
            g.legend(handles=g.legend_.legendHandles[0:2], frameon=False, bbox_to_anchor=(1, 1.15))
            sns.despine()
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/all_{landscape_def}_{cretype}frac-tisspec-enhXtisspec_bytissue_boxplot.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: H boxplot of all tissues by frac-tissue-spec enhancers, hue = exp type
            fig = plt.figure(figsize=(10,6))
            g = sns.boxplot(data=all_tis_merged,
                            y='enh_num', x='tissue', hue='anno', order=tis_order,
                            palette={'Broad':'#fee5d9', 'Tissue-specific':'#de2d26'},
                            showfliers=False)
            g.set_xticklabels(tis_names, rotation=30, horizontalalignment='right')
            g.set_ylabel('# of CREs')
            g.set_xlabel('')
            g.legend(handles=g.legend_.legendHandles, frameon=False, bbox_to_anchor=(1, 1.15))
            sns.despine()
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/all_{landscape_def}_{cretype}num-enhXtisspec_bytissue_boxplot.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig

            ### fig: H boxplot of all tissues by frac-tissue-spec enhancers, hue = exp type
            fig = plt.figure(figsize=(10,6))
            g = sns.boxplot(data=all_tis_merged,
                            y='frac_phastcons', x='tissue', hue='anno', order=tis_order,
                            palette={'Broad':'#fee5d9', 'Tissue-specific':'#de2d26'},
                            showfliers=False)
            g.set_xticklabels(tis_names, rotation=30, horizontalalignment='right')
            g.set_ylabel('Fraction of conserved CRE bp')
            g.set_xlabel('')
            g.legend(handles=g.legend_.legendHandles, frameon=False, bbox_to_anchor=(1, 1.15))
            sns.despine()
            plt.tight_layout()
            plt.savefig(f'{RES_PATH}/all_{landscape_def}_{cretype}frac-phastconsXtisspec_bytissue_H_boxplot_hue.{fmt}', format=fmt, dpi=400)
            plt.close()
            ### end_fig


        ### fig/table: MWU test of fraction of tissue-specific enhancers
        logging.info(f'MWU test of fraction of tissue-specific enhancers, {cretype}{landscape_def}')
        for tis_name in tis_order:
            ts, p = stats.mannwhitneyu(all_tis_merged.query('anno=="Tissue-specific"').query(f'tissue=="{tis_name}"').frac_tisspec_enh,
                                       all_tis_merged.query('anno=="Broad"').query(f'tissue=="{tis_name}"').frac_tisspec_enh)
            logging.info(f'{tis_name}, p = {p:.3}')

        logging.info(all_tis_merged.groupby(['tissue', 'anno']).frac_tisspec_enh.median())
        ### end_fig

        ### fig/table: MWU test of number of enhancers
        logging.info(f'MWU test of number of enhancers, {cretype}{landscape_def}')
        for tis_name in tis_order:
            ts, p = stats.mannwhitneyu(all_tis_merged.query('anno=="Tissue-specific"').query(f'tissue=="{tis_name}"').enh_num,
                                       all_tis_merged.query('anno=="Broad"').query(f'tissue=="{tis_name}"').enh_num)
            logging.info(f'{tis_name}, p = {p:.3}')

        logging.info(all_tis_merged.groupby(['tissue', 'anno']).enh_num.median())
        ### end_fig

        ### fig/table: MWU test of fraction of phastcons bp (in cre)
        logging.info(f'MWU test of fraction of phastcons bp, {cretype}{landscape_def}')
        for tis_name in tis_order:
            ts, p = stats.mannwhitneyu(all_tis_merged.query('anno=="Tissue-specific"').query(f'tissue=="{tis_name}"').frac_phastcons,
                                       all_tis_merged.query('anno=="Broad"').query(f'tissue=="{tis_name}"').frac_phastcons)
            logging.info(f'{tis_name}, p = {p:.3}')

        logging.info(all_tis_merged.groupby(['tissue', 'anno']).frac_phastcons.median())
        ### end_fig
