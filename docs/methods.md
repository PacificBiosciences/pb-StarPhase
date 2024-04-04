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

## _CYP2D6_
Diplotyping of _CYP2D6_ is further complicated due to the presence of homologous gene _CYP2D7_ as well as hybrids between the two genes.
This means that an individual may have 3 or more alleles present on a single haplotype, and 6 or more in a given diplotype.
Additionally, _CYP2D6_ is a much larger region than the HLA genes and requires the use of partially spanning reads to help resolve the diplotype.
To address this problem, pb-StarPhase uses an assembly-like approach to resolve _CYP2D6_:

1. **Load the CYP2D6 database** - This step loads the JSON database that was created from querying PharmVar. Each haplotype is loaded as a VCF file, distinguishing between variants tagged as having an impact ("VI" in VCF) and those that are not.
2. **Load spanning reads** - This step will parse the alignment file (BAM) and extract all reads that span regions of interest **by at least 25% (by default)**. The algorithm will search for reads within the following region of interest, and extract any signatures corresponding to _CYP2D6_, _CYP2D7_, hybrid alleles, or deletion alleles (*5).
    * Full extraction region: GRCh38, chr22:42122691-42145903
3. **Generate the detected alleles** - This step will look at all of the extracted signatures and generate multiple consensus sequences to explain the observations.
4. **Label the detected alleles** - Each allele is compared to the database and assigned a label corresponding to _CYP2D7_, a hybrid identifier, or a full _CYP2D6_ star allele. Incomplete alleles are labeled as "UNKNOWN" and discarded in the following steps.
5. **Generate allele chains** - The signatures from each read are then compared to each of the labeled alleles to form "chains" (i.e., partial haplotypes). An example chain might look like this, `0_CYP2D6*4 -> 1_CYP2D6*68 -> 2_CYP2D7`, indicating a CYP2D6 allele, followed by a hybrid allele (*68), followed by CYP2D7. Note that these chains are created in reference order (GRCh38).
6. **Identify the best full chain pair** - Finally, the chain observations are used to generate candidate full haplotypes. These haplotypes are then paired into diplotypes and scored based on how well they explain the observed chains and the number of observed copies of each allele. The best scoring pair is reported as the diplotype for CYP2D6.
