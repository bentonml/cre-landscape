###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2022.02.18
#
#  tis_spec threshold = 0.3
###

import os
import pandas as pd
import numpy  as np
from datetime import date


### // constants and paths \\ ###
EXP_DAT_PATH = f'../../link_cre_to_genes/dat/2022-03-08'
OUT_DAT_PATH = f'../dat/{str(date.today())}'

# create a date stamped dir for files
if not os.path.isdir(OUT_DAT_PATH):
    os.makedirs(OUT_DAT_PATH)

landscape_def = ['loop', 'contact']
tisspec_thresh = 0.3

tis_order = ['ovary', 'psoas_muscle', 'heart_left_ventricle', 'lung', 'spleen',
             'small_intestine', 'pancreas', 'liver', 'brain_prefrontal_cortex', 'brain_hippocampus']

def read_data(landscape_def, EXP_DATA_PATH):
    df_lst = []
    for tis in tis_order:
        if landscape_def == 'contact':
            infile = f'{EXP_DATA_PATH}/{tis}_gene_with_enh-count_frac-enh-spec_hicQ05-linked.tsv'
        elif landscape_def == 'loop':
            infile = f'{EXP_DATA_PATH}/{tis}_gene_with_enh-count_frac-enh-spec_peakachuloop-linked.tsv'
        enh_num_by_gene = pd.read_table(infile)
        enh_num_by_gene = enh_num_by_gene.assign(exp_nocat=lambda x: np.where((x.exp == 1) & (x.hk == 0) & (x.lof_intol == 0) & (x.essential == 0), 1, 0))
        enh_num_by_gene = enh_num_by_gene.assign(tis_spec=lambda x: np.where(x.rel_entropy > tisspec_thresh, 1, 0))
        enh_num_by_gene = enh_num_by_gene.assign(exp_broad=lambda x: np.where((x.exp == 1) & (x.tis_spec == 0), 1, 0))
        enh_num_by_gene = enh_num_by_gene.assign(log2_exp=lambda x: np.log2(x[tis] + 1))
        enh_num_by_gene = enh_num_by_gene.rename(columns={tis:'gtex_exp'})
        enh_num_by_gene = enh_num_by_gene.assign(tissue=tis)
        df_lst.append(enh_num_by_gene)

    all_tis = pd.concat(df_lst)
    return all_tis
### \\


for landscape_def in ['loop', 'contact']:
    ### // create dataframe of enhancer number by gene with all tissues \\ ###
    all_tis = read_data(landscape_def, EXP_DAT_PATH)
    all_tis_exp = all_tis[all_tis['exp']==1].copy(deep=True)

    # create expressed-only version
    enh_by_expgene_anno = (all_tis_exp
                               .filter(['target_gene', 'tissue', 'exp', 'hk', 'lof_intol', 'essential', 'tis_spec', 'exp_nocat', 'exp_broad_nohk', 'exp_broad', 'enh_num'])
                               .assign(expressed=lambda x: x.exp)
                               .set_index(['target_gene', 'tissue', 'enh_num', 'exp'])
                               .stack()
                               .reset_index()
                               .rename(columns={'level_4':'anno', 0:'val'})
                               .query('val == 1')
                               .merge(all_tis_exp.filter(['target_gene', 'gtex_exp', 'cds_length', 'gene_length', 'tissue', 'enh_rel_entropy', 'frac_tisspec_enh', 'frac_phastcons', 'tss']), how='inner', validate='m:1'))

    enh_by_expgene_anno = enh_by_expgene_anno.assign(num_tisspec_enh=lambda x: x.enh_num * x.frac_tisspec_enh)
    enh_by_expgene_anno = enh_by_expgene_anno.assign(gtex_exp_log2=lambda x: np.log2(x.gtex_exp))

    for tis in tis_order:
        # calculate dataframes for each tissue
        ts = (enh_by_expgene_anno
                .query(f'tissue=="{tis}" & (anno=="tis_spec" | anno=="exp_broad")')
                .filter(['anno', 'gtex_exp_log2', 'target_gene', 'enh_num']))
        ts.to_csv(f'{OUT_DAT_PATH}/{tis}_{landscape_def}.tsv', sep='\t', index=False, header=True)

        # calculate dataframe for each tissue, CRE > 0
        ts = (enh_by_expgene_anno
                .query(f'tissue=="{tis}" & enh_num>0 & (anno=="tis_spec" | anno=="exp_broad")')
                .filter(['anno', 'gtex_exp_log2', 'target_gene', 'enh_num']))
        ts.to_csv(f'{OUT_DAT_PATH}/gt0_{tis}_{landscape_def}.tsv', sep='\t', index=False, header=True)

