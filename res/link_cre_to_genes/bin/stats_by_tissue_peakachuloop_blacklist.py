###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2021.07.23
#
#   assumes pwd is located here with hydrogen:
#       /dors/capra_lab/users/bentonml/cross_species_gain_loss/bin/slurm/
#
#  tis_spec threshold = 0.6
###

import os
import logging
import pandas as pd
import numpy  as np

from pybedtools import BedTool
from datetime import date
from scipy import stats

### // constants and paths // ###
DATA_PATH = '../../../dat'           # relative to current dir
GEN_DATA_PATH = '../../../../data'
DORS_DATA_PATH = '../../../../../../data'
RES_PATH = '../'

tisspec_thresh = 0.6

tis_to_loop = {'ovary':'Schmitt_2016.Ovary.hg19.peakachu-merged.loops',
               'psoas_muscle':'Schmitt_2016.Psoas.hg19.peakachu-merged.loops',
               'heart_left_ventricle':'Leung_2015.VentricleLeft.hg19.peakachu-merged.loops',
               'lung':'Schmitt_2016.Lung.hg19.peakachu-merged.loops',
               'spleen':'Schmitt_2016.Spleen.hg19.peakachu-merged.loops',
               'small_intestine':'Schmitt_2016.Bowel_Small.hg19.peakachu-merged.loops',
               'pancreas':'Schmitt_2016.Pancreas.hg19.peakachu-merged.loops',
               'liver':'Leung_2015.Liver.hg19.peakachu-merged.loops',
               'brain_prefrontal_cortex':'Schmitt_2016.Cortex_DLPFC.hg19.peakachu-merged.loops',
               'brain_hippocampus':'Schmitt_2016.Hippocampus.hg19.peakachu-merged.loops'}

# set up log file to record data about run
logging.basicConfig(filename=f'{RES_PATH}/README_{str(date.today())}', level=logging.INFO)
logging.info('STATS ON BLACKLIST LOOPS')
logging.info('Running on %s with peakachu loops.', str(date.today()))
### \\


### // functions // ###
def intersect_loop_with_enh(loop_df, enh_df):
    ''' Intersect all windows with a Hi-C loop with the enhancer annotations for
    enhancers in a given tissue.
    '''
    n = ['chrom', f'start', f'end', 'enh_chrom', 'enh_start', 'enh_end', 'enh_rel_entropy']
    return loop_df.intersect(BedTool.from_dataframe(enh_df), wo=True).cut([0,1,2,3,4,5,9]).to_dataframe(names=n)


def intersect_loop_with_gene(loop_df, tss, tissue):
    ''' Intersect all windows with a Hi-C loop with the TSS annotations for
    genes in a given tissue.
    '''
    n = ['chrom', 'start', 'end', 'target_gene', tissue, 'hk', 'lof_intol', 'essential', 'rel_entropy']
    return loop_df.intersect(tss, wo=True).cut([0,1,2,6,7,8,9,10,11]).to_dataframe(names=n)


def count_per_loop(df, gene=True):
    ''' Count the number of annotations per window. For enhancers calculate the
    mean enhancer entropy for all enhancers in that window as well (NaN if none).
    '''
    res = df.groupby(['chrom', 'start', 'end'], as_index=False)
    if gene:
        res = (res.count()
                  .filter([f'chrom', f'start', f'end', 'target_gene'])
                  .rename(columns={'target_gene':f'g_count'}))
    else:
        res = (res.agg({'enh_chrom':'count', 'enh_rel_entropy':'mean'})
                  .filter([f'chrom', f'start', f'end', 'enh_chrom', 'enh_rel_entropy'])
                  .rename(columns={'enh_chrom':f'e_count', 'enh_rel_entropy':f'e_entropy'}))
    return res


def create_df_anno(loop_df, ec, gc):
    ''' Create full dataframe with all overlapping annotations. Fills in NaN counts
    with 0 since these are the locations with no overlapping annotations.
    '''
    res = (loop_df.merge(ec, how='left', validate='m:1')
                  .merge(gc, how='left', validate='m:1'))
    res['e_count'] = res.e_count.fillna(0)
    res['g_count'] = res.g_count.fillna(0)
    return res


def create_enh_num_by_gene_df(enh_to_gene_df, gene_anno_df):
    ''' Creates a dataframe of single Hi-C windows with all gene and enhancer
    annotations. Calculates the fraction of tissue specific enhancers in a landscape,
    specified as enhancers with entropy > 0.6. Genes are counted as tissue specific
    if they have entropy > 0.6 as well.
    '''
    res = (enh_to_gene_df.groupby('target_gene', as_index=False)
                         .agg({'enh_chrom':'count', 'tss':'mean', 'enh_rel_entropy':'mean', 'tisspec_enh':'sum'}))

    # add in genes that don't have any assigned enhancers
    res = (gene_anno_df.merge(res, how='left', left_on='name', right_on='target_gene')
                       .drop(columns=['target_gene'])
                       .rename(columns={'enh_chrom':'enh_num', 'name':'target_gene'}))

    # fill with 0 for genes with no enhancers or other tss in window
    res['enh_num'] = res.enh_num.fillna(0)
    res['tisspec_enh'] = res.tisspec_enh.fillna(0)
    res['tss'] = res.tss.fillna(0)

    # because all windows with tss count themselves, counteract blanket nan to zero
    res = res.assign(tss=lambda x: np.where(x.tss == 0, 1, x.tss),
                     tss_bins=lambda x: pd.cut(x['tss'], bins=[0,1,2,3,4,8], labels=['1', '2', '3', '4','5+']),
                     tis_spec=lambda x: np.where(x.rel_entropy > tisspec_thresh, 1, 0),
                     frac_tisspec_enh=lambda x: np.where(x.enh_num > 0, x.tisspec_enh/x.enh_num, np.NaN))
    return res
### \\


### // read general datasets // ###
# read blacklist and set of hg19 ensembl genes
blacklist = BedTool(f'{GEN_DATA_PATH}/dna/hg19/hg19_blacklist_gap.bed').sort().merge()
gene = BedTool(f'{GEN_DATA_PATH}/dna/hg19/ensembl/Hs_GRCh37-75_filter.bed')

# save locations from hg19(grch37.75) ensembl minus blacklist
gene_loc = (gene.intersect(blacklist, v=True)
                .to_dataframe(names=['chrom', 'start', 'end', 'name', 'gene_type', 'strand'])
                .assign(tss_start=lambda x: x.start - 1, tss_end=lambda x: x.start))

# read gtex rna-seq expression values for all tissue, save schmitt et al 2016 subset
hpa = (pd.read_table(f'{GEN_DATA_PATH}/human_protein_atlas/rna_tissue_gtex.tsv')
         .pivot(index='Gene', columns='Tissue', values='TPM')
         .dropna(axis=1)
         .filter(['ovary', 'skeletal muscle', 'heart muscle', 'lung', 'spleen',
                  'small intestine', 'pancreas', 'liver', 'cerebral cortex', 'hippocampal formation'])
         .rename(columns={'skeletal muscle':'psoas_muscle', 'heart muscle':'heart_left_ventricle',
                          'small intestine':'small_intestine','cerebral cortex':'brain_prefrontal_cortex',
                          'hippocampal formation':'brain_hippocampus'}))

# read/parse list of housekeeping genes from eisenberg et al 2013
hk = (pd.read_table(f'{GEN_DATA_PATH}/dna/hg19/eisenberg13_housekeeping_ensemblhg19.bed',
                    header=None, names=['hkchrom', 'hkstart', 'hkend', 'gene'])
            .assign(hk=1)
            .filter(['gene', 'hk']))

# read/parse list of loss of function intolerant genes from gnomad
lof = (pd.read_table(f'{GEN_DATA_PATH}/dna/hg19/gnomad.v2.1.1.lof_metrics.by_gene.bed')
            .assign(lof_intol=lambda x: np.where(x.oe_lof_upper < 0.35, 1, 0))
            .rename(columns={'gene_id':'gene'})
            .filter(['gene', 'lof_intol']))

# read/parse list of essential genes from hart et al 2015
essen = (pd.read_table(f'{GEN_DATA_PATH}/dna/hg19/hart2015_hgnc_core_essential_genes.txt',
                            header=None, names=['symbol', 'hgnc'])
               .drop(columns='symbol')
               .merge(pd.read_table(f'{GEN_DATA_PATH}/dna/hg19/hgnc_to_ensembl.txt')
                         .rename(columns={'HGNC ID':'hgnc',
                                          'Approved symbol':'symbol',
                                          'Approved name':'name',
                                          'Status':'status',
                                          'Ensembl ID(supplied by Ensembl)':'ensembl_id'})
                         .drop(columns=['symbol', 'name', 'status']))
               .rename(columns={'ensembl_id':'gene'})
               .drop(columns='hgnc')
               .assign(essential=1))

# read gene entropy values from file
gene_spec = (pd.read_table(f'{DATA_PATH}/gene-specificity_GTEx_TPM_Schmitt16tissues.tsv')
               .filter(['Gene', 'rel_entropy'])
               .rename(columns={'Gene':'gene'}))

# read enhancer entropy values from file
enhs_filter_entropy = pd.read_table(f'{DATA_PATH}/2020-07-17_enhs_filter_entropy_220bp.tsv')
enhs_filter_entropy = (BedTool.from_dataframe(enhs_filter_entropy)
                              .intersect(blacklist, v=True)
                              .to_dataframe(names=['chrom', 'start', 'end', 'tis_code', 'tis_name', 'old_len', 'entropy']))
### \\

for tis in ['ovary', 'psoas_muscle', 'heart_left_ventricle', 'lung', 'spleen', 'small_intestine', 'pancreas', 'liver', 'brain_prefrontal_cortex', 'brain_hippocampus']:
    hpa_tis = (hpa.reset_index()
                  .filter(['Gene', tis])
                  .rename(columns={'Gene':'gene'})
                  .assign(exp=lambda x: np.where(x[tis] > 1, 1, 0)))

    tis_enhs = enhs_filter_entropy[enhs_filter_entropy.tis_name == tis]

    gene_anno = (gene_loc.filter(['chrom', 'start', 'end', 'name', 'gene_type', 'tss_start', 'tss_end'])
                         .query("gene_type.str.contains('protein_coding')")
                         .merge(hpa_tis, left_on='name', right_on='gene', how='inner')
                         .merge(gene_spec, left_on='name', right_on='gene', how='inner')
                         .merge(hk, left_on='name', right_on='gene', how='left')
                         .merge(lof, left_on='name', right_on='gene', how='left')
                         .merge(essen, left_on='name', right_on='gene', how='left')
                         .drop(columns=['gene_x', 'gene_y', 'gene_type', 'gene'])
                         .drop_duplicates()
                         .replace(np.NaN, 0))

    hic_pairs = pd.read_table(f'{DATA_PATH}/peakachu/{tis_to_loop[tis]}',
                    names=['rchrom', 'rstart', 'rend', 'cchrom', 'cstart', 'cend'])

    # filtering out negative lengths for now, could update to flip coords
    hic_regions = (hic_pairs.filter(['rchrom', 'rstart', 'cend'])
                            .rename(columns={'rchrom':'chrom', 'rstart':'start', 'cend':'end'})
                            .assign(length=lambda x: x.end - x.start)
                            .query('length > 0'))

    tss = BedTool.from_dataframe(gene_anno.filter(['chrom', 'tss_start', 'tss_end', 'name', tis, 'hk', 'lof_intol', 'essential', 'rel_entropy']))
    blacklistloops = BedTool.from_dataframe(hic_regions).sort().cut([0,1,2]).intersect(blacklist, u=True)

    # intersection of row/col loc and enhancer locations
    loop_enh = intersect_loop_with_enh(blacklistloops, tis_enhs)

    # intersection of row/col loc and gene id/type
    loop_gene = intersect_loop_with_gene(blacklistloops, tss, tissue=tis)

    # count the number of enhancers or genes per location
    enh_count = count_per_loop(loop_enh, gene=False)
    gene_count = count_per_loop(loop_gene, gene=True)

    # create full dataframe with all overlapping annotations
    hic_anno = create_df_anno(blacklistloops.to_dataframe(names=['chrom', 'start', 'end']), enh_count, gene_count)

    logging.info(f'{tis}')
    logging.info(f'CRE > 0:  {hic_anno.query("e_count>0").shape[0]}')
    logging.info(f'Gene > 0: {hic_anno.query("g_count>0").shape[0]}')
    logging.info(f'Link > 0: {hic_anno.query("e_count>0 & g_count>0").shape[0]}, {sum(hic_anno.query("e_count>0 & g_count>0").g_count)} genes, {sum(hic_anno.query("e_count>0 & g_count>0").e_count)} CREs')
    logging.info(hic_anno.filter(['e_count', 'g_count']).describe())
