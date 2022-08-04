# Cis-regulatory landscape size, constraint, and tissue-specificity associate with gene function and expression

The code segments below will run the custom scripts from each analysis section.

### Generate CRE landscapes for genes
```
# pwd: res/link_cre_to_genes/bin
python link cre_to_gene.py
python plot_all_landscapes.py
```

### Summary statistics for CRE landscapes
```
# pwd: res/link_cre_to_genes/bin

python calc_cre_stats_by_tissue.py
python stats_by_tissue_peakachuloop_blacklist.py
python stats_by_tiss-specificity.py
```

### Tissue-specific gene analysis
```
# pwd: res/tisspec_features/bin
./run_tisspec_analysis
```

### Housekeeping and loss-of-function intolerant gene analysis
```
# pwd: res/match_landscapes/bin
./run_matching_analysis
```

### DNA sequence conservation analysis
```
# pwd: res/seq_conservation/bin
./run_seq_conserv_analysis
```

### Expression variability analysis
```
# pwd: res/exp_variation/bin
./run_variability_anlaysis
```

### eQTL enrichment analysis
```
# pwd: res/eqtl_enrichment/bin
./run_enrichment_local_loop
./run_enrichment_local_contact
```




