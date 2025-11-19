# Supported genes
The complete list of PGx genes that StarPhase will call is somewhat dynamic.
The database generation component queries CPIC, PharmVar, and IMGT/HLA for genes and alleles.
While some genes are hard-coded for retrieval, many of the CPIC and PharmVar gene lists are dynamically pulled from API endpoints.
This means that if a new gene is added to the upstream database, StarPhase may automatically add it when the database is next constructed.
The following list contains the genes that were included with a particular database update.

Last updated: v2.0.0

# Genes by source
## CPIC and PharmVar
The database build can be configured to accept either CPIC or PharmVar as the upstream data source.
We recommend reviewing using `pbstarphase db-stat` to determine the upstream database for a given database.
See [pbstarphase_20251030.db_stat.txt](../data/v2.0.0/pbstarphase_20251030.db_stat.txt) for an example of the gene lists sourced from CPIC or PharmVar.

## PharmVar only
CYP2D6

## IMGT/HLA
HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1,
HLA-DQA1, HLA-DQB1, HLA-DRB1, HLA-DRB3, HLA-DRB4,
HLA-DRB5
