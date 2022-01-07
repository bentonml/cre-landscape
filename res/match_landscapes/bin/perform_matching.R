###
#   name      | mary lauren benton
#   created   | 2021.11.12
#   updated   | 2021.12.15
#
#   use matchit in R to generate gene sets matched on covariates
#   requires: $ ml GCC/10.2.0  OpenMPI/4.0.5  R/4.0.5
###

library(MatchIt)
library(readr)
suppressMessages(library(dplyr))

# constants and paths
#TODO: update path dates
INDATE <- '2021-12-15'
RESDATE <- '2021-12-15'

DORS <- '/dors/capra_lab/users/bentonml/cre_landscape'
DATA <- paste(DORS, '/res/match_landscapes/dat/', INDATE, sep='')
RESDATA <- paste(DORS, '/res/match_landscapes/dat/', RESDATE, sep='')

tissues <- c('ovary', 'psoas_muscle', 'heart_left_ventricle', 'lung', 'spleen',
             'small_intestine', 'pancreas', 'liver', 'brain_prefrontal_cortex', 'brain_hippocampus')

# plots to assess match fit
plot_match_data <- function(mdat, t, dtype, ltype, plot) {
  pdf(paste('../fig/',RESDATE, '_', t,'_', dtype, '_', ltype, '_match_', plot, '.pdf', sep=''))
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
      infile <- paste(DATA, '/hk_', tis, '_', landscape, '.tsv', sep='')
      hk <- read.csv(infile, sep='\t')
      hk$anno <- factor(hk$anno, levels=c("exp_nocat", "hk"), labels=c('Expressed', 'Housekeeping'))
      matched <- matchit(anno ~ gtex_exp_log2, data=hk, estimand="ATT", replace=FALSE, ratio=2, caliper=0.1)

      return_match_statistics(matched, tis=tis, dtype='hk')
      mdata <- match.data(matched, data=hk, distance="prop.score")
      write_tsv(mdata, paste(RESDATA, '/matched_hk_', tis, '_', landscape, '.tsv', sep=''))

      # lof intolerant v. expressed genes
      infile <- paste(DATA, '/lof_', tis, '_', landscape, '.tsv', sep='')
      lo <- read.csv(infile, sep='\t')
      lo$anno <- factor(lo$anno, levels=c("exp_nocat", "lof_intol"), labels=c('Expressed', 'LoF Intolerant'))
      matched <- matchit(anno ~ gtex_exp_log2, data=lo, estimand="ATT", replace=FALSE, ratio=2, caliper=0.1)

      return_match_statistics(matched, tis=tis, dtype='lof')
      mdata <- match.data(matched, data=lo, distance="prop.score")
      write_tsv(mdata, paste(RESDATA, '/matched_lof_', tis, '_', landscape, '.tsv', sep=''))
      
      # housekeeping v. expressed genes, CRE > 0
      infile <- paste(DATA, '/hk_gt0_', tis, '_', landscape, '.tsv', sep='')
      hk <- read.csv(infile, sep='\t')
      hk$anno <- factor(hk$anno, levels=c("exp_nocat", "hk"), labels=c('Expressed', 'Housekeeping'))
      matched <- matchit(anno ~ gtex_exp_log2, data=hk, estimand="ATT", replace=FALSE, ratio=2, caliper=0.1)

      return_match_statistics(matched, tis=tis, dtype='hk_gt0')
      mdata <- match.data(matched, data=hk, distance="prop.score")
      write_tsv(mdata, paste(RESDATA, '/matched_hk_gt0_', tis, '_', landscape, '.tsv', sep=''))

      # lof intolerant v. expressed genes
      infile <- paste(DATA, '/lof_gt0_', tis, '_', landscape, '.tsv', sep='')
      lo <- read.csv(infile, sep='\t')
      lo$anno <- factor(lo$anno, levels=c("exp_nocat", "lof_intol"), labels=c('Expressed', 'LoF Intolerant'))
      matched <- matchit(anno ~ gtex_exp_log2, data=lo, estimand="ATT", replace=FALSE, ratio=2, caliper=0.1)

      return_match_statistics(matched, tis=tis, dtype='lof_gt0')
      mdata <- match.data(matched, data=lo, distance="prop.score")
      write_tsv(mdata, paste(RESDATA, '/matched_lof_gt0_', tis, '_', landscape, '.tsv', sep=''))
    }
}
