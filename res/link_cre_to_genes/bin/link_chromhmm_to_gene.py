###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2019.12.10
#
#  tis_spec threshold = 0.6 (CRE)
###

import os
import logging
import pandas as pd
import numpy  as np

from pybedtools import BedTool
from datetime import date
from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt

### // constants and paths // ###
DATA_PATH = '../../../dat'           # relative to current dir
GEN_DATA_PATH = '../../../../data'
DORS_DATA_PATH = '../../../../../../data'
RES_PATH = '../'
FIG_PATH = f'{RES_PATH}/fig/{str(date.today())}'

# create a date stamped dir for files
if not os.path.isdir(FIG_PATH):
    os.makedirs(FIG_PATH)

tisspec_thresh = 0.6  # for tissue-specificity metric (CRE)
thresh='Q05'          # for hic q-values
t = 0.05              # for hic q-values

# figure options
fmt='pdf'
rc = {'font.size':14, 'axes.titlesize':16,'axes.labelsize':14, 'legend.fontsize': 14,
      'xtick.labelsize': 14, 'ytick.labelsize': 14}

tis_to_tad = {'ovary':'Ovary_Schmitt2016-raw_TADs_hg19from38.bed',
              'psoas_muscle':'Muscle_Psoas_Donor-PO1-raw_TADs.txt',
              'heart_left_ventricle':'VentricleLeft_STL003_Leung2015-raw_TADs.txt',
              'lung':'Lung_Donor-LG1-raw_TADs.txt',
              'spleen':'Spleen_Donor-PX1-raw_TADs.txt',
              'small_intestine':'Bowel_Small_Donor-SB2-raw_TADs.txt',
              'pancreas':'Pancreas_Donor-PA2-raw_TADs.txt',
              'liver':'Liver_STL011_Leung2015-raw_TADs_hg19from38.bed',
              'brain_prefrontal_cortex':'Cortex_DLPFC_Donor-CO-raw_TADs.txt',
              'brain_hippocampus':'Hippocampus_Schmitt2016-raw_TADs_hg19from38.bed'}

tis_to_hic = {'ovary':'OV', 'psoas_muscle':'PO', 'heart_left_ventricle':'LV',
              'lung':'LG', 'spleen':'SX', 'small_intestine':'SB', 'pancreas':'PA',
              'liver':'LI', 'brain_prefrontal_cortex':'CO', 'brain_hippocampus':'HC'}

constants = ['chr1$', 'chr2$', 'chr3$', 'chr4$', 'chr5$', 'chr6$', 'chr7$', 'chr8$',
             'chr9$', 'chr10$', 'chr11$', 'chr12$', 'chr13$', 'chr14$', 'chr15$',
             'chr16$', 'chr17$', 'chr18$', 'chr19$', 'chr20$', 'chr21$', 'chr22$', 'chrX$']

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
logging.info(f'Link CREs to genes ({str(date.today())})')
### \\


### // functions // ###
def create_window_df(tad):
    ''' Creates 1Mb and 40kb non-overlapping windows tiled across the genome (hg19).
    Joins these with provided TAD annotations in bed format.
    '''
    window01mb = BedTool().window_maker(genome='hg19', w=1000000)
    window40kb = BedTool().window_maker(genome='hg19', w=40000)

    # create dataframe of 40kb and 1Mb windows
    joint_windows = (window40kb.intersect(window01mb, wao=True)
                               .cut([0,1,2,3,4,5])
                               .to_dataframe()
                               .assign(window_id=lambda x: x.groupby(['name', 'score', 'strand'], as_index=False).ngroup())
                               .drop(columns=['name', 'score', 'strand'])
                               .assign(window40_id=lambda x: x.groupby(['chrom']).cumcount()))

    # create dataframe of windows with TAD intersections
    by_tad = (window40kb.intersect(tad, wao=True)
                       .cut([0,1,2,3,4,5])
                       .to_dataframe()
                       .assign(tad_id=lambda x: x.groupby(['name', 'score', 'strand'], as_index=False).ngroup() - 1)
                       .drop(columns=['name', 'score', 'strand']))

    # return merged dataframe of windows and TAD annotations
    return joint_windows.merge(by_tad, how='left')


def create_hic_df(filename, pattern, bed, threshold=0.05):
    ''' Creates a dataframe of Hi-C interactions less than the provided Q-value
    threshold.
    '''
    return (pd.read_table(filename)
              .assign(RowID=lambda x: x.RowID - 1, ColumnID=lambda x: x.ColumnID- 1)
              .query('QValue < ' + str(threshold))
              .drop(columns=['ExpectedCount'])
              .merge(bed[bed.chrom == pattern[0:-1]], left_on='RowID', right_on='window40_id', how='left')
              .rename(columns={'chrom':'rchrom',
                              'start':'rstart',
                              'end':'rend',
                              'window_id':'rwindow',
                              'tad_id':'rtad'})
              .drop(columns=['window40_id'])
              .merge(bed[bed.chrom == pattern[0:-1]], left_on='ColumnID', right_on='window40_id', how='left')
              .rename(columns={'chrom':'cchrom',
                              'start':'cstart',
                              'end':'cend',
                              'window_id':'cwindow',
                              'tad_id':'ctad',
                              'ObservedCount':'obs_count',
                              'QValue':'q'})
              .filter(['rchrom', 'rstart', 'rend', 'rwindow', 'rtad',
                      'cchrom', 'cstart', 'cend', 'cwindow', 'ctad', 'obs_count', 'q'])
                .assign(distance=lambda x: x.cstart - x.rend))


def create_hic_pairs_df(tad_bed, hic_patterns, hic_code, hic_path, thresh):
    ''' Creates dataframe of paired locations connected by significant Hi-C contacts
    for a given tissues and threshold.

    tad_bed : bedfile object of tad locations
    hic_patterns : regex chromosome patterns for Hi-C matrix files from FitHiC
    hic_code : tissue code for Hi-C matrix files from FitHiC
    hic_path : location of Hi-C matrix files from FitHiC
    thresh : Qvalue threshold to determine significant contacts
    '''
    hic_dfs = []

    # run through each chromosome and calculate number of contacts
    for pattern in hic_patterns:
        # create dataframe of windows (40kb and 1Mb)
        bed = create_window_df(tad_bed)

        # read from FitHiC output files to get significant contacts
        hic_pattern = 'chr23' if pattern[0:-1] == 'chrX' else pattern[0:-1]
        filename = f'{GEN_DATA_PATH}/chromatin_structure/schmitt16/GSE87112/FitHiC_primary_cohort/FitHiC_output.{hic_code}_{hic_pattern}.sparse.matrix.gz'
        hic = create_hic_df(filename=filename, pattern=pattern, bed=bed, threshold=thresh)

        hic_dfs.append(hic.assign(sig_connect=lambda x: np.where(x.obs_count > 0, 1, 0),
                         under_1mb=lambda x: np.where(x.distance <= 1000 * 1000, 1, 0),
                         in_same_win=lambda x: np.where(x.rwindow == x.cwindow, 1, 0),
                         in_same_tad=lambda x: np.where((x.rtad != -1) & (x.ctad != -1) & (x.rtad == x.ctad), 1, 0)))

    return pd.concat(hic_dfs)


def intersect_hic_with_enh(hic_pair, enh_df, row=True):
    ''' Intersect all windows with a Hi-C contact with the enhancer annotations for
    enhancers in a given tissue.
    '''
    pfx = 'r' if row else 'c'
    n = [f'{pfx}chrom', f'{pfx}start', f'{pfx}end', 'enh_chrom', 'enh_start', 'enh_end', 'enh_rel_entropy', 'phastcons_overlap']
    return hic_pair.intersect(BedTool.from_dataframe(enh_df), wo=True).cut([0,1,2,3,4,5,9,10]).to_dataframe(names=n)


def intersect_hic_with_gene(hic_pair, tss, tissue, row=True):
    ''' Intersect all windows with a Hi-C contact with the TSS annotations for
    genes in a given tissue.
    '''
    pfx = 'r' if row else 'c'
    n = [f'{pfx}chrom', f'{pfx}start', f'{pfx}end', 'target_gene', tissue, 'hk', 'lof_intol', 'essential', 'rel_entropy']
    return hic_pair.intersect(tss, wo=True).cut([0,1,2,6,7,8,9,10,11]).to_dataframe(names=n)


def count_per_hic(df, row=True, gene=True):
    ''' Count the number of annotations per window, either in the row or column
    of the Hi-C matrix. For enhancers calculate the mean enhancer entropy for all
    enhancers in that window as well (NaN if no enhancers).
    '''
    pfx = 'r' if row else 'c'
    res = df.groupby([f'{pfx}chrom', f'{pfx}start', f'{pfx}end'], as_index=False)
    if gene:
        res = (res.count()
                  .filter([f'{pfx}chrom', f'{pfx}start', f'{pfx}end', 'target_gene'])
                  .rename(columns={'target_gene':f'{pfx}g_count'}))
    else:
        res = (res.agg({'enh_chrom':'count', 'enh_rel_entropy':'mean'})
                  .filter([f'{pfx}chrom', f'{pfx}start', f'{pfx}end', 'enh_chrom', 'enh_rel_entropy'])
                  .rename(columns={'enh_chrom':f'{pfx}e_count', 'enh_rel_entropy':f'{pfx}e_entropy'}))
    return res


def create_hic_df_anno(hic_pairs_df, rec, rgc, cec, cgc):
    ''' Create full dataframe with all overlapping annotations. Fills in NaN counts
    with 0 since these are the locations with no overlapping annotations.
    '''
    res = (hic_pairs_df.drop(columns=['obs_count', 'q', 'distance', 'rwindow', 'cwindow'])
                       .merge(rec, how='left', validate='m:1')
                       .merge(rgc, how='left', validate='m:1')
                       .merge(cec, how='left', validate='m:1')
                       .merge(cgc, how='left', validate='m:1'))
    res['re_count'] = res.re_count.fillna(0)
    res['ce_count'] = res.ce_count.fillna(0)
    res['rg_count'] = res.rg_count.fillna(0)
    res['cg_count'] = res.cg_count.fillna(0)
    return res


def intersect_loop_with_enh(loop_df, enh_df):
    ''' Intersect all windows with a Hi-C loop with the enhancer annotations for
    enhancers in a given tissue.
    '''
    n = ['chrom', f'start', f'end', 'enh_chrom', 'enh_start', 'enh_end', 'enh_rel_entropy', 'phastcons_overlap']
    return loop_df.intersect(BedTool.from_dataframe(enh_df), wo=True).cut([0,1,2,3,4,5,9,10]).to_dataframe(names=n)


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
    return res.drop_duplicates()


def create_enh_num_by_gene_df(enh_to_gene_df, gene_anno_df):
    ''' Creates a dataframe of single Hi-C windows with all gene and enhancer
    annotations. Calculates the fraction of tissue specific enhancers in a landscape,
    specified as enhancers with entropy > 0.6.
    '''
    res = (enh_to_gene_df.groupby('target_gene', as_index=False)
                         .agg({'enh_chrom':'count', 'enh_rel_entropy':'mean', 'tisspec_enh':'sum', 'enh_length':'sum', 'phastcons_overlap':'sum'}))

    # add in genes that don't have any assigned enhancers
    res = (gene_anno_df.merge(res, how='left', left_on='name', right_on='target_gene', validate='1:1')
                       .drop(columns=['target_gene'])
                       .rename(columns={'enh_chrom':'enh_num', 'name':'target_gene'}))

    # fill with 0 for genes with no enhancers or other tss in window
    res['enh_num'] = res.enh_num.fillna(0)
    res['tisspec_enh'] = res.tisspec_enh.fillna(0)

    # because all windows with tss count themselves, counteract blanket nan to zero
    res = res.assign(frac_tisspec_enh=lambda x: np.where(x.enh_num > 0, x.tisspec_enh/x.enh_num, np.NaN),
                     frac_phastcons=lambda x: np.where(x.enh_num > 0, x.phastcons_overlap/x.enh_length, np.NaN))
    return res


def to_csv_new_subdir(df, subdir, outname):
    ''' Saves provided dataframe to the subdirectory and filename; creates the
    tree for the subdirectory if it does not already exist.
    '''
    os.makedirs(subdir, exist_ok=True)
    full_name = os.path.join(subdir, outname)
    df.to_csv(full_name, sep='\t', index=False, header=True)
### \\


### // read general datasets // ###
# read blacklist and set of hg19 ensembl genes
blacklist = BedTool(f'{GEN_DATA_PATH}/dna/hg19/hg19_blacklist_gap.bed').sort().merge()
gene = BedTool(f'{GEN_DATA_PATH}/dna/hg19/ensembl/Hs_GRCh37-75_filter.bed')

# save locations from hg19(grch37.75) ensembl minus blacklist
gene_loc = (gene.intersect(blacklist, v=True)
                .to_dataframe(names=['chrom', 'start', 'end', 'name', 'gene_type', 'strand'])
                .assign(tss_start=lambda x: x.start - 1, tss_end=lambda x: x.start,
                        gene_length=lambda x: x.end - x.start))

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
         .rename(columns={'gene_id':'gene', 'oe_lof_upper':'loeuf'})
         .filter(['gene', 'lof_intol', 'loeuf']))

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

# read phastcons conserved elements
phastcons = BedTool(f'{GEN_DATA_PATH}/phastcons/phastConsElements46way_All_merged.bed')

# read enhancer entropy values from file
enhs_filter_entropy = pd.read_table(f'{DATA_PATH}/chromhmm/2022-12-20_chromhmm_entropy_220bp.tsv')
enhs_filter_entropy = (BedTool.from_dataframe(enhs_filter_entropy)
                              .intersect(blacklist, v=True)
                              .intersect(phastcons, wao=True)
                              .cut([0,1,2,3,4,5,6,10])
                              .to_dataframe(names=['chrom', 'start', 'end', 'tis_code', 'tis_name', 'elength', 'entropy', 'phastcons_overlap'])
                              .groupby(['chrom', 'start', 'end', 'tis_code', 'tis_name', 'elength', 'entropy'], as_index=False)
                              .sum())

# filter enhancers to remove ones in the top 1% (> ~5kb)
upper = enhs_filter_entropy.elength.quantile(q=.99)
len_fltr_enh = enhs_filter_entropy.query(f'elength < {upper}')
### \\

for tis in ['ovary', 'psoas_muscle', 'heart_left_ventricle', 'lung', 'spleen', 'small_intestine', 'pancreas', 'liver', 'brain_prefrontal_cortex', 'brain_hippocampus']:
    logging.info(f'Running {tis}')

    ### // read tissue specific datasets \\ ###
    # read tissue specific datasets
    hpa_tis = (hpa.reset_index()
                  .filter(['Gene', tis])
                  .rename(columns={'Gene':'gene'})
                  .assign(exp=lambda x: np.where(x[tis] > 1, 1, 0)))

    tads = BedTool(f'{GEN_DATA_PATH}/chromatin_structure/3d_genome_browser/hg19/{tis_to_tad[tis]}')

    tis_enhs = len_fltr_enh[len_fltr_enh.tis_name == tis]
    ### \\

    ### // combine relevant columns to single dataframe \\ ###
    # save all genes to single dataframe, lose genes without data in GTEx, n = 17313
    gene_anno = (gene_loc.filter(['chrom', 'start', 'end', 'name', 'gene_type', 'tss_start', 'tss_end', 'gene_length'])
                         .query("gene_type.str.contains('protein_coding')")
                         .merge(hpa_tis, left_on='name', right_on='gene', how='inner')
                         .merge(gene_spec, left_on='name', right_on='gene', how='inner')
                         .merge(hk, left_on='name', right_on='gene', how='left')
                         .merge(lof, left_on='name', right_on='gene', how='left')
                         .merge(essen, left_on='name', right_on='gene', how='left')
                         .drop(columns=['gene_x', 'gene_y', 'gene_type', 'gene'])
                         .drop_duplicates()
                         .replace(np.NaN, 0))
    ### \\

    ### // CONTACT-BASED LANDSCAPE \\ ###
    logging.info(f'Generate contact-based landscape for {tis}')
    ### create hic contact dataframe (hic_pairs)
    hic_pairs = create_hic_pairs_df(tad_bed=tads, hic_patterns=constants,
                                    hic_code=tis_to_hic[tis], hic_path=GEN_DATA_PATH,
                                    thresh=t)

    ### intersections with 40kb windows
    row = BedTool.from_dataframe(hic_pairs.filter(['rchrom', 'rstart', 'rend']).drop_duplicates())
    col = BedTool.from_dataframe(hic_pairs.filter(['cchrom', 'cstart', 'cend']).drop_duplicates())
    tss = BedTool.from_dataframe(gene_anno.filter(['chrom', 'tss_start', 'tss_end', 'name', tis, 'hk', 'lof_intol', 'essential', 'rel_entropy']))

    # intersection of row/col loc and enhancer locations
    row_enh = intersect_hic_with_enh(row, tis_enhs, row=True)
    col_enh = intersect_hic_with_enh(col, tis_enhs, row=False)

    # intersection of row/col loc and gene id/type
    row_gene = intersect_hic_with_gene(row, tss, tis, row=True)
    col_gene = intersect_hic_with_gene(col, tss, tis, row=False)

    # count the number of enhancers or genes per location
    row_enh_count = count_per_hic(row_enh, row=True, gene=False)
    col_enh_count = count_per_hic(col_enh, row=False, gene=False)
    row_gene_count = count_per_hic(row_gene, row=True, gene=True)
    col_gene_count = count_per_hic(col_gene, row=False, gene=True)

    # create full dataframe with all overlapping annotations
    hic_pairs_anno = create_hic_df_anno(hic_pairs, row_enh_count, row_gene_count, col_enh_count, col_gene_count)

    ### assign enhancers to targets
    enh_to_gene = (pd.concat([hic_pairs_anno.query('re_count > 0 & cg_count > 0').merge(col_gene).merge(row_enh).drop(columns=['rg_count', 'cg_count']),
                              hic_pairs_anno.query('ce_count > 0 & rg_count > 0').merge(row_gene).merge(col_enh).drop(columns=['cg_count', 'rg_count'])], sort=False)
                     .drop(columns=['rtad', 'ctad', 'in_same_win', 'under_1mb'])
                     .filter(['enh_chrom', 'enh_start', 'enh_end', 'sig_connect', 'in_same_tad', 'target_gene', tis, 'hk', 'core', 'lof_intol', 'essential', 'rel_entropy', 'enh_rel_entropy', 'phastcons_overlap'])
                     .assign(tisspec_enh=lambda x: np.where(x.enh_rel_entropy > tisspec_thresh, 1, 0),
                             enh_length=lambda x: x.enh_end - x.enh_start)
                     .drop_duplicates()
                     .query('in_same_tad > 0'))

    # save dataframe of number of enhancers linked to each gene
    enh_num_by_gene = create_enh_num_by_gene_df(enh_to_gene, gene_anno)
    to_csv_new_subdir(enh_num_by_gene, f'{RES_PATH}/dat/{str(date.today())}', f'{tis}_gene_with_chromhmm-enh-count_frac-enh-spec_hicQ05-linked.tsv')

    # save dataframe of enhancer to gene links
    enh_to_gene = enh_to_gene.rename(columns={tis:'gtex_tpm'}).assign(tissue=tis)
    to_csv_new_subdir(enh_to_gene, f'{RES_PATH}/dat/{str(date.today())}', f'{tis}_chromhmm-enh_to_gene_hicQ05-linked.tsv')
    ### \\

    ### // LOOP-BASED LANDSCAPE \\ ###
    logging.info(f'Generate loop-based landscape for {tis}')
    ### create hic contact dataframe (hic_pairs)
    hic_pairs = pd.read_table(f'{DATA_PATH}/peakachu/{tis_to_loop[tis]}',
                    names=['rchrom', 'rstart', 'rend', 'cchrom', 'cstart', 'cend'])

    # filtering out negative lengths for now, could update to flip coords
    hic_regions = (hic_pairs.filter(['rchrom', 'rstart', 'cend'])
                            .rename(columns={'rchrom':'chrom', 'rstart':'start', 'cend':'end'})
                            .assign(length=lambda x: x.end - x.start)
                            .query('length > 0'))

    ### intersections with 40kb windows
    loops = BedTool.from_dataframe(hic_regions).cut([0,1,2]).intersect(blacklist, v=True).sort()
    tss = BedTool.from_dataframe(gene_anno.filter(['chrom', 'tss_start', 'tss_end', 'name', tis, 'hk', 'lof_intol', 'essential', 'rel_entropy']))

    logging.info(f'{loops.count()} total loops')

    # intersection of row/col loc and enhancer locations
    loop_enh = intersect_loop_with_enh(loops, tis_enhs).assign(enh_length=lambda x: x.enh_end-x.enh_start)

    # intersection of row/col loc and gene id/type
    loop_gene = intersect_loop_with_gene(loops, tss, tissue=tis)
    genes_in_loops = loop_gene.target_gene.drop_duplicates()

    # count the number of enhancers or genes per location
    enh_count = count_per_loop(loop_enh, gene=False)
    gene_count = count_per_loop(loop_gene, gene=True)

    # create full dataframe with all overlapping annotations
    hic_anno = create_df_anno(loops.to_dataframe(names=['chrom', 'start', 'end']), enh_count, gene_count)

    ### fig: distplot for gene/enh per loop, length of loops in tissue
    estat = hic_anno.e_count.describe()
    ehi = (1.5 * (estat['75%'] - estat['25%'])) + estat['75%']
    gstat = hic_anno.g_count.describe()
    ghi = (1.5 * (gstat['75%'] - gstat['25%'])) + gstat['75%']

    with sns.plotting_context("paper", rc=rc):
        fig, axs = plt.subplots(nrows=3, figsize=(4,10))
        sns.histplot(hic_anno.query(f'e_count <= {ehi}').e_count, color='.3', kde=False, discrete=True, ax=axs[0])
        axs[0].set_xlabel('# CRE / Loop')
        sns.histplot(hic_anno.query(f'g_count <= {ghi}').g_count, color='.3', kde=False, discrete=True, ax=axs[1])
        axs[1].set_xlabel('# Gene / Loop')
        sns.histplot((hic_anno.end-hic_anno.start)/1000, kde=False, log_scale=True, color='.3', ax=axs[2])
        axs[2].set_xlabel('Length (kb)')
        sns.despine()
        plt.tight_layout()
        plt.savefig(f'{FIG_PATH}/{tis}_peakachuloop_anno-per-loop_distplot_chromhmm.{fmt}', format=fmt, dpi=400)
        plt.close()
    ### end_fig

    stat = (hic_anno.merge(loop_gene, how='left', validate='1:m')
                    .drop_duplicates()
                    .assign(exp=lambda x: np.where(x[tis]>1, 1, 0))
                    .groupby(['chrom', 'start', 'end'], as_index=False)
                    .agg({'g_count':'mean', 'e_count':'mean', 'exp':'sum'}))

    logging.info(f'Average loop length: {((hic_anno.end-hic_anno.start)/1000).describe()[1]:.2f} kb, Median loop length: {((hic_anno.end-hic_anno.start)/1000).describe()[5]:.2f} kb')
    logging.info('Loops with no genes')
    logging.info(f'{stat.query("g_count==0").shape[0]}, {stat.query("g_count==0").shape[0]/stat.shape[0]*100:.2f}%')
    logging.info('Loops with at least one gene and no enhancers')
    logging.info(f'{stat.query("g_count>0 & e_count==0").shape[0]}, {stat.query("g_count>0 & e_count==0").shape[0]/stat.shape[0]*100:.2f}%')
    logging.info('Loops with only not-expressed genes')
    logging.info(f'{stat.query("g_count>0 & exp==0").shape[0]}, {stat.query("g_count>0 & exp==0").shape[0]/stat.shape[0]*100:.2f}%')
    logging.info('Loops with only not-expressed genes and at least one enhancer')
    logging.info(f'{stat.query("g_count>0 & exp==0 & e_count>0").shape[0]}, {stat.query("g_count>0 & exp==0 & e_count>0").shape[0]/stat.shape[0]*100:.2f}%')
    logging.info(hic_anno.filter(['e_count', 'g_count']).describe())

    ### assign enhancers to targets
    enh_to_gene = (hic_anno.query('e_count > 0 & g_count > 0').merge(loop_gene).merge(loop_enh).rename(columns={'e_count':'enh_num'}).drop(columns=['e_entropy', 'g_count'])
                         .filter(['enh_chrom', 'enh_start', 'enh_end', 'target_gene', tis, 'hk', 'core', 'lof_intol', 'essential', 'rel_entropy', 'enh_rel_entropy', 'enh_length', 'phastcons_overlap'])
                         .assign(tisspec_enh=lambda x: np.where(x.enh_rel_entropy > tisspec_thresh, 1, 0))
                         .drop_duplicates())

    logging.info(f'Saving enh_num_by_gene dataframe for {tis}')
    # save dataframe of number of enhancers linked to each gene, only those in loops (possible to link)
    enh_num_by_gene = create_enh_num_by_gene_df(enh_to_gene, gene_anno[gene_anno.name.isin(genes_in_loops)]).drop_duplicates()
    to_csv_new_subdir(enh_num_by_gene, f'{RES_PATH}/dat/{str(date.today())}', f'{tis}_gene_with_chromhmm-enh-count_frac-enh-spec_peakachuloop-linked.tsv')

    logging.info(f'Saving enh_num_by_gene dataframe for {tis}')
    # save dataframe of enhancer to gene links
    enh_to_gene = enh_to_gene.rename(columns={tis:'gtex_tpm'}).assign(tissue=tis).drop_duplicates()
    to_csv_new_subdir(enh_to_gene, f'{RES_PATH}/dat/{str(date.today())}', f'{tis}_chromhmm-enh_to_gene_peakachuloop-linked.tsv')
    ### \\
