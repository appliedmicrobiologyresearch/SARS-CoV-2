# SARS-CoV-2

The project space accompanies the manuscript "SARS-CoV-2 phylogeny during the early outbreak in the Basel area, Switzerland: import and spread dominated by a single B.1 lineage variant (C15324T)" (preprint https://www.medrxiv.org/content/10.1101/2020.09.01.20186155v1). It is divided into two parts, the first contains the code for the COVID-19 Genome Analysis Pipieline (COVGAP). COVGAP is used to assemble consensus sequences from Illumina paired-end reads and to call genomic variants, which can be analysed in downstream applications. The second part contains the code that was used to produce the global phylogeny presented in the manuscript. 

## COVGAP pipeline:

A downloadable version of COVGAP as described in Stange et al., 2021 is available for linux and Mac at https://github.com/appliedmicrobiologyresearch/covgap.


## Dependencies for part 2: phylogenetic analysis of Basel sequences in global context
  - global genomes and metadata were downloaded from GISAID (https://www.gisaid.org/)
  - we filtered all genomes (GISAID and COVGAP produced consensus sequences) prior to the final analysis to contain less than 10% Ns (ambigious characters) using this biophyton script: https://biopython.org/wiki/Sequence_Cleaner
  - we time filter the reads using R and packages tidyr, dplyr, and readr. In theory, nextstrain offers a time filter by dd-mm-yyyy however, for some reason only mm-yyyy worked and was not sufficient for our purpose to filter by specific dates. We dropped all GISAID genomes that did not provide information to the day of sampling.
  - for installation of the nextstrain pipeline follow the instructions here https://github.com/nextstrain/ncov
  - nextstrain offers an alignment option within the pipeline, however we chose to use mafft v.7.467 (download here https://mafft.cbrc.jp/alignment/software/) and run it outside the nextstrain pipeline as this version can handle several thousands of genomes and is still very fast and accurate.
