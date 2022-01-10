###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2021.11.12
#
#   compare gene sets matched on covariates (expression level)
###

import os
import pandas as pd
import numpy  as np
from datetime import date

### // constants and paths \\ ###
EXP_DAT_PATH = f'../../link_cre_to_genes/dat/2022-01-07'
OUT_DAT_PATH = f'../dat/match_landscapes/dat/{str(date.today())}'

# create a date stamped dir for files
if not os.path.isdir(OUT_DAT_PATH):
    os.makedirs(OUT_DAT_PATH)

landscape_def = ['loop', 'contact']

tis_order = ['ovary', 'psoas_muscle', 'heart_left_ventricle', 'lung', 'spleen',
             'small_intestine', 'pancreas', 'liver', 'brain_prefrontal_cortex', 'brain_hippocampus']

def read_data(landscape_def):
    df_lst = []
    for tis in tis_order:
        if landscape_def == 'contact':
            infile = f'{EXP_DAT_PATH}/{tis}_gene_with_enh-count_frac-enh-spec_hicQ05-linked.tsv'
        elif landscape_def == 'loop':
            infile = f'{EXP_DAT_PATH}/{tis}_gene_with_enh-count_frac-enh-spec_peakachuloop-linked.tsv'
        enh_num_by_gene = pd.read_table(infile)
        enh_num_by_gene = enh_num_by_gene.assign(exp_nocat=lambda x: np.where((x.exp == 1) & (x.hk == 0) & (x.lof_intol == 0) & (x.essential == 0), 1, 0))
        enh_num_by_gene = enh_num_by_gene.assign(exp_broad=lambda x: np.where((x.exp == 1) & (x.tis_spec == 0), 1, 0))
        enh_num_by_gene = enh_num_by_gene.rename(columns={tis:'gtex_exp'})
        enh_num_by_gene = enh_num_by_gene.assign(tissue=tis)
        df_lst.append(enh_num_by_gene)
    all_tis = pd.concat(df_lst)
    return all_tis
### \\


for landscape in landscape_def:
    # read all data for landscape type
    all_tis = read_data(landscape)
    all_tis_exp = all_tis[all_tis['exp']==1].copy(deep=True)

    # create expressed gene only dataframe
    enh_by_expgene_anno = (all_tis_exp
                               .filter(['target_gene', 'tissue', 'exp', 'hk', 'rel_entropy', 'lof_intol', 'essential', 'tis_spec', 'exp_nocat', 'exp_broad_nohk', 'exp_broad', 'enh_num'])
                               .assign(expressed=lambda x: x.exp)
                               .set_index(['target_gene', 'tissue', 'enh_num', 'exp', 'rel_entropy'])
                               .stack()
                               .reset_index()
                               .rename(columns={'level_5':'anno', 0:'val'})
                               .query('val == 1')
                               .merge(all_tis_exp.filter(['target_gene', 'gtex_exp', 'cds_length', 'gene_length', 'tissue', 'enh_rel_entropy', 'frac_tisspec_enh', 'frac_phastcons', 'tss']), how='inner', validate='m:1')
                               .assign(gtex_exp_log2=lambda x: np.log2(x.gtex_exp)))

    for tis in tis_order:
        # calculate hk and lof dataframes for each tissue
        hk = (enh_by_expgene_anno
                .query(f'tissue=="{tis}" & (anno=="hk" | anno=="exp_nocat")')
                .filter(['anno', 'gtex_exp_log2', 'target_gene', 'enh_num']))
        hk.to_csv(f'{OUT_DAT_PATH}/hk_{tis}_{landscape}.tsv', sep='\t', index=False, header=True)

        lo = (enh_by_expgene_anno
                .query(f'tissue=="{tis}" & (anno=="lof_intol" | anno=="exp_nocat")')
                .filter(['anno', 'gtex_exp_log2', 'target_gene', 'enh_num']))
        lo.to_csv(f'{OUT_DAT_PATH}/lof_{tis}_{landscape}.tsv', sep='\t', index=False, header=True)

        # calculate hk and lof dataframes for each tissue, CRE > 0
        hk = (enh_by_expgene_anno
                .query(f'tissue=="{tis}" & enh_num>0 & (anno=="hk" | anno=="exp_nocat")')
                .filter(['anno', 'gtex_exp_log2', 'target_gene', 'enh_num']))
        hk.to_csv(f'{OUT_DAT_PATH}/hk_gt0_{tis}_{landscape}.tsv', sep='\t', index=False, header=True)

        lo = (enh_by_expgene_anno
                .query(f'tissue=="{tis}" & enh_num>0 & (anno=="lof_intol" | anno=="exp_nocat")')
                .filter(['anno', 'gtex_exp_log2', 'target_gene', 'enh_num']))
        lo.to_csv(f'{OUT_DAT_PATH}/lof_gt0_{tis}_{landscape}.tsv', sep='\t', index=False, header=True)
