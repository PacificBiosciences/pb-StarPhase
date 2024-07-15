# v0.11.3
## Fixed
- Updated the `build` mode to account for IMGT-HLA's new database format that was released with `v3.57.0-alpha`

# v0.11.2
## Changes
- Added a penalty for CYP2D6 haplotypes composed of only hybrid alleles that are expected to be found paired with another allele
- Relaxed the requirements for CYP2D6 chain inferrence to allow for more inferred connections when observed connections are missing

## Fixed
- Fixed a bug where chains connecting from some region to REP6 were not being penalized (e.g., REP6 should always start a chain)
- Fixed a residual bug where a poorly mapping read in an HLA region would fail to match any alleles in the database and cause a panic

# v0.11.1
## Fixed
- Fixed a bug where empty haplotypes were treated as valid candidates and then reported as part of the diplotype
- Revised inferred chain logic to remove erroneous inferrences and add missing inferrences
- Fixed a bug where a poorly mapping read in an HLA region would fail to match any alleles in the database and cause a panic

# v0.11.0
## Changes
- The underlying methodology has been significantly altered to improve _CYP2D6_ diplotyping in targeted sequencing
  - The core _CYP2D6_ regions have been shrunk to contain just the regions containing variants that define the allele
  - Additional regions have been added solely for the purpose of linking _CYP2D6_ alleles: REP6, REP7, spacer, and "link_region" (region between CYP2D6 and CYP2D7 typically)
  - The chaining algorithm has been altered to account for the additional regions above. A "normal" haplotype chain is expected to have the following order of regions in GRCh38: REP6 -> CYP2D6 -> link_region -> REP7 -> spacer -> CYP2D7
  - The scoring of diplotype chain pairs is now based on a unified scoring scheme that accounts for: 1) edit distance of observations to the chain pair, 2) likelihood of the chain pair based on allele coverage and multinomial, 3) lasso penalty for duplications, and 4) penalty for unexpected chain events (see above "normal" chain)
  - The debug BAM will now output a haplotype block for each identified region (e.g., REP6 will have its own block)

# v0.10.2
## Changes
- The CLI has been modified such that the VCF file is now optional. If a VCF file is not provided, all variant-based diplotyping will be skipped and those genes will be absent from all output files. If no VCF or BAM files are provided, pb-StarPhase will generate a user error message.
- Two new options control the genes that are diplotyped: `--include-set` and `--exclude-set`. Only one of these options can be specified at a time. Both accept a plain text file with one gene name per line. If `--include-set` is specified, then only the genes in the given file will be diplotyped. If `--exclude-set` is specified, then all genes will be diplotyped _except_ the ones in the given file.
- Underlying CYP2D6 consensus algorithm was updated for greater compute efficiency

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
