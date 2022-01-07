###
#   name      | mary lauren benton
#   conda env | enh_gain-loss
#   created   | 2021.04.26
#
###


import pandas as pd
import numpy  as np
from datetime import date
from scipy import stats


### // constants, paths, functions \\ ###
GTEX_PATH = '/dors/capra_lab/users/bentonml/data/gtex/v8/expression'
DAT_PATH = '../dat'
RES_FILE = f'../res/{str(date.today())}_calc-cv.out'

exp_above_thresh = [0.8, 1.0]
exp_thresh = 1

tis_name = {'ovary':'Ovary', 'psoas_muscle': 'Muscle', 'heart_left_ventricle':'Heart',
            'lung':'Lung', 'spleen':'Spleen', 'small_intestine':'Small intestine',
            'liver':'Liver', 'pancreas':'Pancreas', 'brain_prefrontal_cortex':'Prefrontal cortex',
            'brain_hippocampus':'Hippocampus'}

tissue_to_smtsd = {'brain_prefrontal_cortex':'Brain - Frontal Cortex (BA9)',
                   'psoas_muscle':'Muscle - Skeletal',
                   'ovary':'Ovary',
                   'lung':'Lung',
                   'spleen':'Spleen',
                   'pancreas':'Pancreas',
                   'small_intestine':'Small Intestine - Terminal Ileum',
                   'heart_left_ventricle':'Heart - Left Ventricle',
                   'liver':'Liver',
                   'brain_hippocampus':'Brain - Hippocampus'}

def coef_var(x):
    ''' Calculates the coefficient of variation for a dataframe row.
    '''
    return np.std(x) / np.mean(x)

### \\


with open(RES_FILE, 'w') as outfile:
    # read sample ids
    id_map = pd.read_csv(f'{GTEX_PATH}/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t')

    for tis in tis_name.keys():
        outfile.write(f'{tis}\n')

        # merge sample id with tissue names
        tis_ids = list(id_map[id_map.SMTSD == tissue_to_smtsd[tis]]['SAMPID'])

        exp = (pd.read_csv(f'{GTEX_PATH}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct',
                           sep='\t', skiprows=2, usecols=lambda x: x in ['Name']+tis_ids)
                           .set_index('Name'))

        for required_prop in exp_above_thresh:
            # filter out lowly expressed genes— >exp_thresh TPM in at least 80%
            tot_individ = exp.shape[1]
            exp_filter = exp[exp.apply(lambda x: (sum(np.abs(x) > exp_thresh) / tot_individ) >= required_prop, axis=1)]
            outfile.write(f'{tis} (>{exp_thresh}): {tot_individ} samples, {exp_filter.shape[0]} genes\n')

            # transform the TPM (log2) to better account for heteroscedasticity
            exp_filter = np.log2(exp_filter+0.01)

            # calculate the coefficient of variation
            cv = list(exp_filter.apply(coef_var, axis=1))

            # calculate the median log2 expression (TPM)
            median_log2_tpm = list(exp_filter.median(axis=1))

            # create final df
            final = exp_filter.reset_index().filter(['Name'])
            final['Name'] = final.Name.str.split('.').str[0]
            final['median_log2_tpm'] = median_log2_tpm
            final['cv'] = cv

            # save intermediate files
            final.to_csv(f'{DAT_PATH}/{str(date.today())}_{tis}_{exp_thresh}_gt{required_prop}_gene_by_cv.tsv', sep='\t', index=False)
