# v0.10.1
## Changes
- Updated the crate for _CYP2D6_ generation, improving alignment offset seeding and reducing over-splitting for shorter reads (e.g., targeted sequencing)
- Added a secondary penalty for using a non-unique chain pair during the _CYP2D6_ chaining step. This tends to reduce errors caused by ambiguous chains in WGS datasets.

## Fixed
- Fixed an issue with _CYP2D6_ where an "UNKNOWN" cluster could erroneously generate an extra allele when two defined alleles came from the same HPC-cluster
- Fixed an issue with _CYP2D6_ where alleles excluding from normalization could form their own unpenalized haplotype (e.g., a *68 allele by itself)
- Fixed a panic caused by typing an empty _CYP2D6_ sequence
- Fixed a panic caused by no usable reads in the _CYP2D6_ search region

# v0.10.0
## Changes
- Added support for calling _CYP2D6_ from targeted sequencing data
  - In general, accuracy for targeted datasets is less than that of WGS. This is largely due to difficulties with capture that lead to decreased coverage of hybrid or duplicated alleles.
  - We recommend using two additional parameters when using targeted sequencing data: `--infer-connections --normalize-d6-only`
- Added two new CLI options to support targeted sequencing datasets:
  - `--infer-connections` - If set, pb-StarPhase will infer allele connections that are not observed in the dataset but common in the population. For example, *4 and *68 are commonly found together, as are *10 and *36. This option is recommended when reads are too short to directly span from one allele to the next.
  - `--normalize-d6-only` - If set, pb-StarPhase will only normalize the copy numbers using the _CYP2D6_ alleles (i.e., excluding any discovered _CYP2D7_ alleles). This option is recommended when coverage of the _CYP2D7_ alleles is inconsistent relative to the _CYP2D6_ alleles.

## Fixed
- Fixed a reporting issue in the PharmCAT TSV where brackets were missing from combination alleles

# v0.9.1
## Changes
- The CLI settings log output has been updated for easier human readability
- Exposed three new CLI options that influence how the _CYP2D6_ algorithm works: `--min-consensus-fraction`, `--min-consensus-count`, and `--dual-max-ed-delta`. Most users should not need to modify the defaults.

## Fixed
- Improved the allele chaining algorithm to reduce errors with hybrid _CYP2D6_ allele copy-number counts
- Fixed the _CYP2D6_ star allele output to be in coding order (i.e. reverse relative to before)
- Fixed the PharmCAT TSV output to use simple _CYP2D6_ representation

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
