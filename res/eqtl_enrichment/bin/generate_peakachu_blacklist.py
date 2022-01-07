# marylaurenbenton, 2021
# generate blacklist file for peakachu loops for all tissues

import pandas as pd
from pybedtools import BedTool

#TODO: update paths
DATA_PATH = '../../../../data/'        # relative to current dir
GEN_DATA_PATH = '../../../../../data'


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

blacklist = BedTool(f'{GEN_DATA_PATH}/dna/hg19/hg19_blacklist_gap.bed').sort().merge()

for tis in tis_to_loop:
    hic_pairs = pd.read_table(f'{DATA_PATH}/peakachu/{tis_to_loop[tis]}', names=['rchrom', 'rstart', 'rend', 'cchrom', 'cstart', 'cend'])

    # filtering out negative lengths for now, could update to flip coords
    hic_regions = (hic_pairs.filter(['rchrom', 'rstart', 'cend'])
                            .rename(columns={'rchrom':'chrom', 'rstart':'start', 'cend':'end'})
                            .assign(length=lambda x: x.end - x.start)
                            .query('length > 0'))

    loops = BedTool.from_dataframe(hic_regions).cut([0,1,2]).intersect(blacklist, v=True).sort()
    loops.complement(genome='hg19').moveto(f'../dat/{tis}_blacklist.bed')

