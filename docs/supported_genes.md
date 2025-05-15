# Supported genes
The complete list of PGx genes that StarPhase will call is somewhat dynamic.
The database generation component queries CPIC, PharmVar, and IMGT/HLA for genes and alleles.
While some genes are hard-coded for retrieval, many of the CPIC and PharmVar gene lists are dynamically pulled from API endpoints.
This means that if a new gene is added to the upstream database, StarPhase may automatically add it when the database is next constructed.
The following list contains the genes that were included with a particular database update.

Last updated: v1.4.0

# Genes by source
## CPIC
ABCG2, CACNA1S, CFTR,
CYP2B6, CYP2C19, CYP2C9, CYP3A5, CYP4F2,
DPYD, G6PD, IFNL3, MT-RNR1, NUDT15,
RYR1, SLCO1B1, TPMT, UTG1A1, VKORC1

## PharmVar
CYP1A2, CYP2A13, CYP2C8, CYP2D6, CYP3A4,
NAT2

## IMGT/HLA
HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1,
HLA-DQA1, HLA-DQB1, HLA-DRB1, HLA-DRB3, HLA-DRB4,
HLA-DRB5
