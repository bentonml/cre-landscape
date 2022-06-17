#!/bin/bash
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=4:30:00
#
#   requires:
#       ml Anaconda3
#       source activate enh_gain-loss
#
#   uses scripts:
#       ../resources/bin/calculate_enrichment.py 
###

# load modules
module load Anaconda3/4.4.0
source activate enh_gain-loss 

run_enrichment='/dors/capra_lab/users/bentonml/resources/bin/calculate_enrichment.py'
var_path='/dors/capra_lab/users/bentonml/cross_species_gain_loss/data/variants'
cre_path='../dat'
iters=1000


# gtex v8 eQTL 
time python $run_enrichment -i $iters -n 16 -b $cre_path/${tis}_blacklist.bed --print_counts_to $cre_path/${tis}_${lndscp}_cre_quart1.counts $cre_path/${tis}_${lndscp}_cre_quart1.bed $var_path/${tis}_gtex_v8_hg19_signif.bed
time python $run_enrichment -i $iters -n 16 -b $cre_path/${tis}_blacklist.bed --print_counts_to $cre_path/${tis}_${lndscp}_cre_quart2.counts $cre_path/${tis}_${lndscp}_cre_quart2.bed $var_path/${tis}_gtex_v8_hg19_signif.bed
time python $run_enrichment -i $iters -n 16 -b $cre_path/${tis}_blacklist.bed --print_counts_to $cre_path/${tis}_${lndscp}_cre_quart3.counts $cre_path/${tis}_${lndscp}_cre_quart3.bed $var_path/${tis}_gtex_v8_hg19_signif.bed
time python $run_enrichment -i $iters -n 16 -b $cre_path/${tis}_blacklist.bed --print_counts_to $cre_path/${tis}_${lndscp}_cre_quart4.counts $cre_path/${tis}_${lndscp}_cre_quart4.bed $var_path/${tis}_gtex_v8_hg19_signif.bed

