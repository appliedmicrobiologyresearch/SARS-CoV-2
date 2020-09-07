#FILTER DATES BEFORE SUBSAMPLING. NEXTSTRAIN FILTER BY YEAR BUT NOT BY YYYY-MM-DD ALTHOUGH THE MANUAL SAYS SO
#in R
#***
library(tidyr)
library(dplyr)
library(readr)

master= read.delim("metadata.tsv",sep="\t", header = TRUE, fill=TRUE)
  master$date <- as.Date(master$date, "%Y-%m-%d")

# set date range: accute phase - 2019 Dec01 and 2020 Mar23
master.accutePhase=    master %>%
    select(strain,virus,gisaid_epi_isl,genbank_accession,date,region,country,division,location,region_exposure,country_exposure,division_exposure,segment,length,host,age,sex,pangolin_lineage,GISAID_clade,legacy_clade_membership,originating_lab,submitting_lab,authors,url,title,paper_url,date_submitted) %>%
    filter(date >= as.Date("2019-12-01") & date <= as.Date("2020-03-23"))

# save ID list to filter from original fasta file
write.table(master.accutePhase$strain, "filteredIDs.txt", sep="\t", col.names = F, row.names = F, quote = F)

# save new metadata file
write_tsv(as.data.frame(master.accutePhase), "metadata_acutePhase.tsv", na = "", append = FALSE, quote_escape = "double")
