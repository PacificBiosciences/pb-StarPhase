# v0.9.0
## Changes
- Added the ability to store _CYP2D6_ variant information from PharmVar in our database
- Released an updated database with CYP2D6 information
- Added the ability to call _CYP2D6_ using the updated database and an aligned BAM file for WGS datasets
- Added a CLI option to output realigned _CYP2D6_ segments to a small BAM file for IGV viewing and debugging
- The version of pb-StarPhase is now recorded as the `pbstarphase_version` in the output JSON

# v0.8.2
## Fixed
- Fixed a bug where failure to align an HLA allele to a read would cause an unhandled panic

# v0.8.1
## Fixed
- Fixed a bug where the "\*" was not pre-pended to HLA alleles in the output files

# v0.8.0
## Changes
- Added the ability to store HLA haplotype sequences from IMGTHLA in our database
- Added the ability to read and write from a gzip-compressed (.gz) database file
- Released an updated database with HLA haplotype sequences
- Added the ability to call _HLA-A_ and _HLA-B_ using an HLA-aware database file and an aligned BAM file
- Added several CLI parameters to support HLA calling

# v0.7.3
## Fixed
- Fixed a bug where the previous binary would segfault when trying to run `pbstarphase build`
- Pre-compiled static binary build system has changed to be Docker-based

# v0.7.2
Initial release.
