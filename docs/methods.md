# Methods
## CPIC genes
pb-StarPhase uses a relatively simple process to convert a VCF file into star-allele calls.
The high-level process is outlined below:

1. **Load the CPIC database** - This step loads the JSON database that was created from querying the CPIC API. Each haplotype is loaded into memory, and all variants in the haplotype are normalized. Any variants that fail to normalize are discarded from consideration.
2. **Search the VCF file for normalized variants** - Given a set of normalized variant coordinates, this step will search the corresponding VCF region for variant calls that match those from the database. This step also normalizes the genotypes and tracks any relevant phasing information, if available.
3. **Call diplotypes** - This step assigns all homozygous variants to both haplotypes. Heterozygous variants are then assigned combinatorially, and the resulting haplotype pairs are compared to the database. If both haplotypes match a known haplotype from the database, then the diplotype is recorded as a solution. Importantly, this method is phase-aware such that variants on the same phase set ID maintain the correct pairings. Additionally, variants on different phase set IDs are treated as effectively unphased relative to each other.
4. **Report findings** - Each diplotype solution as well as relevant variant information is saved into the output JSON.

## _HLA-A_ and _HLA-B_
The diplotyping of HLA genes is more complicated due to frequent mismappings of other HLA sequences to the wrong positions in the genome.
This can have the downstream effect of generating spurious or incorrect variant calling and/or phasing in the VCF file.
As a result, pb-StarPhase uses a different, alignment-based approach to diplotype the HLA genes:

1. **Load the HLA database** - This step loads the JSON database that was created from querying the IMGTHLA GitHub repo. Each haplotype sequence is loaded into memory for the following steps. 
2. **Load fully spanning reads** - This step will parse the alignment file (BAM) and extract all reads that **fully span** the gene region of interest. We note that the haplotypes stored in IMGTHLA do not have the same mapping coordinates, so we have defined the regions of interest as follows (subject to change as we iterate on performance):
    * _HLA-A_ - chr6:29942254-29945870
    * _HLA-B_ - chr6:31353362-31357442
3. **Score each HLA haplotype against each read** - This step will align each cDNA and DNA sequence against each read and save the edit distance of the mapping against that read. By default, any HLA haplotypes that are missing a DNA sequence are ignored. Any read with an edit distance that is too high will be ignored (this removes mismapped reads).
4. **Call best diplotype** - This step performs a ranked-choice vote where each read mapping is allowed one vote. The process iteratively excludes candidate haplotypes based on the number of votes as well as the edit distance for tie-breaking. If all candidates have a "high" vote fraction (default: >=10%), then a more complicated full comparison is used to exclude a candidate. Once only two candidates remain, the algorithm compares the votes to determine whether it should be reported as heterozygous for both candidates or homozygous for the dominant candidate.
5. **Report findings** - In contrast to the CPIC genes, only a single diplotype result is ever reported for the HLA genes. Additionally, relevant mapping information is saved into the output JSON.
