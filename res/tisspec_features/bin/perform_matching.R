###
#   name      | mary lauren benton
#   created   | 2022.02.18
#
#   use matchit in R to generate gene sets matched on covariates
#   requires: $ ml GCC/10.2.0  OpenMPI/4.0.5  R/4.0.5
###

library(MatchIt)
library(readr)
suppressMessages(library(dplyr))

# constants and paths
INDATE  <- '2022-03-08'
RESDATE <- '2022-03-08'

DORS <- '/dors/capra_lab/users/bentonml/cre_landscape'
DATA <- paste(DORS, '/res/tisspec_features/dat/', INDATE, sep='')
RESDATA <- paste(DORS, '/res/tisspec_features/dat/', RESDATE, sep='')

landscapes <- c('loop', 'contact')
tissues <- c('ovary', 'psoas_muscle', 'heart_left_ventricle', 'lung', 'spleen',
             'small_intestine', 'pancreas', 'liver', 'brain_prefrontal_cortex', 'brain_hippocampus')

# plots to assess match fit
plot_match_data <- function(mdat, t, dtype, ltype, plot) {
  pdf(paste('../fig/',RESDATE, '_', t, '_', dtype, ltype, '_match_', plot, '.pdf', sep=''))
  plot(mdat, type=plot)
  dev.off()
}

# generate fit plots and print summary stats
return_match_statistics <- function(mdat, tis, dtype, ltype) {
  print(paste(dtype, ltype, ':', tis))
  print(summary(matched, standardize = TRUE))
  plot_match_data(matched, t=tis, dtype=dtype, ltype=ltype, plot='ecdf')
  plot_match_data(matched, t=tis, dtype=dtype, ltype=ltype, plot='jitter')
}

for (landscape in landscapes) {
    for (tis in tissues) {
      # housekeeping v. expressed genes
      infile <- paste(DATA, '/', tis, '_', landscape, '.tsv', sep='')
      ts <- read.csv(infile, sep='\t')
      ts$anno <- factor(ts$anno, levels=c("exp_broad", "tis_spec"), labels=c('Broad','Tissue-specific'))
      matched <- matchit(anno ~ gtex_exp_log2, data=ts, estimand="ATT", replace=FALSE, ratio=10, caliper=0.1)

      return_match_statistics(matched, tis=tis, ltype=landscape, dtype='')
      mdata <- match.data(matched, data=ts, distance="prop.score")
      write_tsv(mdata, paste(RESDATA, '/matched_', tis, '_', landscape, '.tsv', sep=''))

      # housekeeping v. expressed genes, CRE > 0
      infile <- paste(DATA, '/gt0_', tis, '_', landscape, '.tsv', sep='')
      ts <- read.csv(infile, sep='\t')
      ts$anno <- factor(ts$anno, levels=c("exp_broad", "tis_spec"), labels=c('Broad','Tissue-specific'))
      matched <- matchit(anno ~ gtex_exp_log2, data=ts, estimand="ATT", replace=FALSE, ratio=10, caliper=0.1)

      return_match_statistics(matched, tis=tis, dtype='gt0_', ltype=landscape)
      mdata <- match.data(matched, data=ts, distance="prop.score")
      write_tsv(mdata, paste(RESDATA, '/matched_gt0_', tis, '_', landscape, '.tsv', sep=''))
    }
}
