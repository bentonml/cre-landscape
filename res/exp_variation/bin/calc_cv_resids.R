# marylaurenbenton, 2021
# calculated the loess residuals from cv ~ median_tpm to remove the confounding
# of expression level on variability
# process described previously in sigalova et al. 2020
# ml GCC/10.2.0  OpenMPI/4.0.5 R

library(stats)
library(ggplot2)

DAT_PATH <- '../dat/'
FIG_PATH <- '../fig/'

tis_names <- c('ovary', 'psoas_muscle', 'heart_left_ventricle', 'lung', 'spleen',
               'small_intestine', 'pancreas', 'liver', 'brain_prefrontal_cortex',
               'brain_hippocampus')

for (tissue in tis_names) {
    print(tissue)
    df  <- read.csv(paste(DAT_PATH, '2021-08-04_', tissue, '_1_gt0.8_gene_by_cv.tsv', sep=''), sep='\t')

    # log scale the cv -- b/c heteroscedasticity 
    df$log2cv <- log2(df$cv)

    # qc plot
    ggplot(df, aes(median_log2_tpm, log2cv)) + 
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    theme_classic()
    ggsave(filename=paste(FIG_PATH, tissue, '_loess_orig.pdf'), device='pdf')

    # run the loess model
    loess_mod <- loess(log2cv ~ median_log2_tpm, data=df, span=0.25, degree=1)

    # qc plot 2
    ggplot(df, aes(median_log2_tpm, loess_mod$residuals)) + 
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    theme_classic()
    ggsave(filename=paste(FIG_PATH, tissue, '_loess_resid.pdf'), device='pdf')

    # add the residuals as column in df
    df$loess_resid <- loess_mod$residuals

    # write residuals to a file for downstream analysis
    write.csv(df, paste(DAT_PATH, tissue, '_loess_cv_resid.csv', sep=''), row.names=FALSE, quote=FALSE)
}

