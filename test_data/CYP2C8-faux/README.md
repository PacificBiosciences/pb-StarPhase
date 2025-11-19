# CYP2C8 faux test data
This sub-folder contains a fake database and VCF files corresponding to cases we want to test that include sub-allele definitions.
These examples were hand-crafted to test edge cases.

## Files
* `database.json` - Contains fake allele definitions. Core alleles are *2 and *3. It also includes three sub-alleles of *2.
* `suballele_match.vcf.gz` - Contains the ideal scenario where both haplotypes exactly match a known sub-allele.
* `core_match.vcf.gz` - Slightly less ideal scenario where both haplotypes match a core allele. However, one of them does not match a defined sub-allele, so only core alleles should be reported in the main output. Inexact allele matches should be populated here.
* `inexact_match.vcf.gz` - Even less ideal scenario where a haplotype does not match a core allele. In this instance, only the inexact allele matches will have meaningful information for a user.
