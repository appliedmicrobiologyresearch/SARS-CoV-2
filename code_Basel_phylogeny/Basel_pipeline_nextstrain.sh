#!/bin/bash

# accompanying code for "SARS-CoV-2 phylogeny during the early outbreak in the Basel area, Switzerland:
# import and spread dominated by a single B.1 lineage variant (C15324T)" 
# authored by Dr. Madlen Stange, stange.madlen@gmail.com

# the following steps assume that all genome consensus sequences (own and GISAID) are stored in sequences.fasta and all metadata as retrieved from GISAID is stored in metadata.tsv

# step 1: filter sequences by date 
# run rscript r_filter_by_date.r which outputs ID list that contains IDs for the set time window: filteredIDs.txt and metadata_acutePhase.tsv
Rscript r_filter_by_date.r
# use the ID list to filter for sequences from sequnces.fasta and output to filtered_time.fasta
awk 'FNR==NR{a[$0];next} /^>/{val=$0;sub(/^>/,"",val);flag=val in a?1:0} flag' filteredIDs.txt sequences.fasta > filtered_time.fasta

# step 2: filter by location using nextstrain software
# exclude samples from Romania have lots of letter "u" in their sequences what messes with subsequent analyses, as well as non-human host samples
# exclude: Switzerland/42* sequences , because those are the published USB Nanopore sequences on GISAID but we will use our polished Illumina ones
source activate nextstrain
augur filter \
--sequences data/filtered_time.fasta \
--metadata data/metadata_acutePhase.tsv \
--include config/ncov_include.txt \
--exclude config/ncov_exclude.txt \
--exclude-where country=Romania host!=human \
--output results/filtered_time_geo30_vJuly3.fasta \
--group-by country year month \
--sequences-per-group 30

# step 3: align sequences to Wuhan-reference using MAFFT v.7.467
mafft --thread 12 --reorder --keeplength --mapout --kimura 1 --addfragments filtered_time_geo30_vJuly3.fasta --auto /config/NC_045512_SARSCov2_ref.fasta > aligned_mafft_time_geo30_vJuly3.fasta

# step 4: mask homoplasic sites and trim ends, see Table S2 
augur mask --sequences results/aligned_mafft_time_geo30_vJuly3.fasta \
--mask config/maskHomoplasic_2020June19.txt \
--mask-from-beginning 50 \
--mask-from-end 99 \
--output results/masked_aligned_mafft_time_geo30_vJuly3.fasta

# step 5: construct phylogenetic tree
augur tree \
--alignment results/masked_aligned_mafft_time_geo30_vJuly3.fasta \
--nthreads 20 \
--substitution-model GTR+G \
--output results/masked_aligned_mafft_time_geo30_vJuly3_tree_raw.nwk

# step 6: time-calibrate tree; root to Wuhan/Hu-1/2019 and Wuhan/WH01/2019
augur refine \
--tree results/masked_aligned_mafft_time_geo30_vJuly3_tree_raw.nwk  \
--alignment results/masked_aligned_mafft_time_geo30_vJuly3.fasta \
--metadata data/metadata_acutePhase.tsv \
--output-tree results/masked_aligned_mafft_time_geo30_vJuly3_tree.nwk \
--output-node-data results/masked_aligned_mafft_time_geo30_vJuly3_branch_lengths.json \
--timetree \
--clock-rate 0.0008  \
--clock-std-dev 0.0004  \
--coalescent skyline \
--date-confidence \
--date-inference marginal \
--clock-filter-iqd 4  \
--root Wuhan/Hu-1/2019 Wuhan/WH01/2019

# step 7: ancestral state reconstruction for region, country, and country of exposure 
augur traits \
--tree results/masked_aligned_mafft_time_geo30_vJuly3_tree.nwk \
--metadata data/metadata_acutePhase.tsv \
--output results/masked_aligned_mafft_time_geo30_vJuly3_exposure_traits.json \
--sampling-bias-correction 2.5 \
--columns region_exposure country_exposure \

augur traits \
--tree results/masked_aligned_mafft_time_geo30_vJuly3_tree.nwk \
--metadata data/metadata_acutePhase.tsv \
--output results/masked_aligned_mafft_time_geo30_vJuly3_traits.json \
--sampling-bias-correction 2.5 \
--columns region country \

# step 8: infer nucleotide mutations per branch  
augur ancestral \
--tree results/masked_aligned_mafft_time_geo30_vJuly3_tree.nwk \
--alignment results/masked_aligned_mafft_time_geo30_vJuly3.fasta \
--output-node-data results/masked_aligned_mafft_time_geo30_vJuly3_nt_muts.json \
--inference joint

# step 9: translate nucleotide mutations to amino acid changes and annotate branches
augur translate \
--tree results/masked_aligned_mafft_time_geo30_vJuly3_tree.nwk  \
--ancestral-sequences results/masked_aligned_mafft_time_geo30_vJuly3_nt_muts.json \
--reference-sequence config/NC_045512_SARSCov2_ref.gb \
--output results/masked_aligned_mafft_time_geo30_vJuly3_aa_muts.json

# step 10: infer GISAID clades
augur clades \
--tree results/masked_aligned_mafft_time_geo30_vJuly3_tree.nwk \
--mutations results/masked_aligned_mafft_time_geo30_vJuly3_nt_muts.json  results/masked_aligned_mafft_time_geo30_vJuly3_aa_muts.json \
--clades config/clades.tsv \
--output-node-data results/masked_aligned_mafft_time_geo30_vJuly3_augur_clades.json

# step 11: export as json to be viewed with auspice
augur export v2 \
--tree results/masked_aligned_mafft_time_geo30_vJuly3_tree.nwk \
--metadata data/metadata_acutePhase.tsv \
--colors config/colours.tsv \
--color-by-metadata country divison location country_exposure region age pangolin_lineage GISAID_clade \
--auspice-config config/auspice_config.json \
--node-data results/masked_aligned_mafft_time_geo30_vJuly3_branch_lengths.json \
    results/masked_aligned_mafft_time_geo30_vJuly3_traits.json \
    results/masked_aligned_mafft_time_geo30_vJuly3_nt_muts.json \
    results/masked_aligned_mafft_time_geo30_vJuly3_aa_muts.json \
    results/masked_aligned_mafft_time_geo30_vJuly3_exposure_traits.json \
--lat-longs config/lat_longs.tsv \
--output auspice/ncov_augur_Basel_vJuly7_updatedColours.json  \
--title "Spread of SARS-CoV-2 in Basel Feb.26-Mar23"

source deactivate
