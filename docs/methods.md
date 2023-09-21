# Methods
PharmGOAT uses a relatively simple process to convert from a VCF file to star-allele calls.
The high-level process is outlined below:

1. Load the CPIC database - This step loads the JSON database that was created from querying the CPIC API. Each haplotype is loaded into memory, and all variants in the haplotype are normalized. Any that fail to normalize are discarded from consideration.
2. Search the VCF for normalized variants - Given a set of normalized variant coordinates, this step will search the corresponding VCF region for variant calls that match those from the database. This steps also normalizes the genotypes and tracks any relevant phasing information, if available.
3. Call diplotypes - This step assigns all homozygous variants to both haplotypes. Heterozygous variants are then assigned combinatorially, and the resulting haplotype pairs are compared to the database. If both haplotypes match a known haplotype from the database, then the diplotype is recorded as a solution. Importantly, this method is phase aware such that variants on the same phase set ID maintain the correct pairings. Additionally, variants on different phase set IDs are treated as effectively unphased relative to each other.
4. Report findings - Each diplotype solution as well as relevant variant information is saved into the output JSON.
