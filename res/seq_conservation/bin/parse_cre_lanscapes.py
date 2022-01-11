###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#
#   create CRE quartiles for each landscape definition
###

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# relative to bin/ dir
DATA_PATH = '../../link_cre_to_genes/dat/2022-01-07'
RES_PATH = '../dat'

final_cols = ['enh_chrom', 'enh_start', 'enh_end', 'landscape_cre_num',
              'landscape_tisspec_cre_num', 'landscape_cre_quartile', 'landscape_cre_quartile_id']

def label_cre_quartiles(landscape_type):
    for tis in ['ovary', 'psoas_muscle', 'heart_left_ventricle', 'lung', 'spleen',
                'small_intestine', 'pancreas', 'liver', 'brain_prefrontal_cortex', 'brain_hippocampus']:

        # read in CRE and landscape data
        cre = pd.read_table(f'{DATA_PATH}/{tis}_enh_to_gene_{landscape_type}-linked.tsv')
        lndscp = (cre.groupby('target_gene', as_index=False)
                     .agg({'enh_chrom':'count', 'enh_rel_entropy':'mean', 'tisspec_enh':'sum'})
                     .rename(columns={'enh_rel_entropy':'landscape_enh_rel_entropy',
                                      'enh_chrom':'landscape_cre_num',
                                      'tisspec_enh':'landscape_tisspec_cre_num'}))

        # separate landscapes into quartiles
        lndscp['landscape_cre_quartile'] = pd.qcut(lndscp.landscape_cre_num, 4)
        lndscp['landscape_cre_quartile_id'] = pd.qcut(lndscp.landscape_cre_num, 4, [1, 2, 3, 4])

        # merge with CRE locations
        lndscp_merge = cre.merge(lndscp, how='left', validate='m:1')
        lndscp_merge = lndscp_merge.filter(final_cols).drop_duplicates()

        lndscp_merge.to_csv(f'{RES_PATH}/{tis}_{landscape_type}_cre.bed', sep='\t', index=False, header=False)

label_cre_quartiles(landscape_type='peakachuloop')
label_cre_quartiles(landscape_type='hicQ05')

