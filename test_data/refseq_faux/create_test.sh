# download the RefSeq file from:
# https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
zcat ./GRCh38_latest_genomic.gff.gz | \
    grep -e "#!" -e "genome=chromosome" -e "HLA-A" -e "HLA-B" \
        -e "HLA-DPA1" -e "HLA-DPB1" \
        -e "HLA-DQA1" -e "HLA-DQB1" \
        -e "HLA-DRB1" -e "HLA-DRB3" -e "HLA-DRB4" -e "HLA-DRB5" \
        -e "NAT2" -e "CACNA1S" | \
    gzip > refseq_small.gff.gz
