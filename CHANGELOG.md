# v1.4.0
## Changes
- Adds in ability to query PharmVar API for additional genes to diplotype
- Database file has been adjusted slightly to add additional source information (i.e., CPIC v. PharmVar source)
- Release of new database file with the additional genes: `data/v1.4.0/pbstarphase_20250515.json.gz`
  - Includes new genes from PharmVar API changes: ["CYP1A2", "CYP2A13", "CYP2C8", "CYP3A4", "NAT2"]

# v1.3.2
## Changes
- Added statistics on HLA diplotype calling to the `hla_debug.json` file

## Fixed
- Instead of throwing an error when encountered a mapped read with no sequence in the HLA regions, the tool will now ignore the read and report a warning in the log

# v1.3.1
## Changes
- Adds additional structural variants to the default database construction: CYP2B6\*29, CYP4F2\*16, SLCO1B1\*48, and SLCO1B1\*49
- Release of new database file with the additional structural variants: `data/v1.3.1/pbstarphase_20250224.json.gz`

# v1.3.0
## Changes
- Adds support for externally called structural variants and corresponding PGx haplotypes
  - A new option `--sv-vcf` can be used to provide a structural variant VCF file from pbsv or sawfish
  - Updated the database build to include _CYP2C19_ known whole-gene (\*36) and partial-gene deletions (\*37) based on PharmVar definitions

# v1.2.0
## Changes
- Added a new debug output file: `cyp2d6_alleles.json`
  - This file contains additional details on the database variants that were identified in each reported _CYP2D6_ allele
  - Full details are specified in the [Debug Outputs](./docs/debug_outputs.md#cyp2d6-alleles)

# v1.1.0
## Changes
- HLA data configuration has been automated to support future haplotype additions
  - The `pbstarphase build` command now requires a reference genome FASTA file for hg38; this file is specified via `--reference`
  - `pbstarphase build` will now automatically pull the latest RefSeq file; coordinates for the MANE transcript are extracted from this resource
  - Configurations for HLA genes have been completely updated to a new structure; if using an older database file, the existing defaults from v1.0.0 will be loaded for HLA-A and HLA-B
- Additional HLA genes are now reported by StarPhase: HLA-C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1, -DRB3, -DRB4, and -DRB5
- HLA algorithm has been updated to accomodate additional HLA genes
  - Reads that only partially span an HLA locus are now supported (at least 50% overlap with a gene)
  - Reads are now pulled from all specified HLA regions and assigned to a single gene based on closest database match
  - Each batch of assigned reads is run through the HLA diplotyping process independently
  - cDNA consensus step has been removed and replaced with a HPC consensus step (similar to CYP2D6); haplotype label assignment still uses cDNA sequence
  - Additional logic has been added to support genes that are commonly absent or hemizygous (HLA-DRB3, -DRB4, and -DRB5) in individuals
- Debug folder updates:
  - `hla_debug.json` has been updated to include target and query unmapped base counts
  - `read_debug.json` has been added to the outputs. This file includes the best mapping of each read to an HLA allele.
  - `hla_igv_custom` has been added to the outputs. This is similar to the previous `cyp2d6_igv_custom`, but contains the assembled haplotypes for the HLA genes. Consensus sequence, database sequences, and user specified sequences are mapped to the custom assemblies. The IGV options on the session have been modified to reasonable defaults for comparing consensus to the reads.
- A new version of the database file has been uploaded to support the above changes (`v1.1.0/pbstarphase_20250110.json.gz`)

# v1.0.1
## Fixed
- Fixed an issue where a homozygous deletion in CYP2D6 (\*5/\*5) was not considered a valid diplotype when using `--normalize-d6-only`

# v1.0.0
## Changes
- Source code is included in the GitHub repository
- `LICENSE.md` has been updated to reflect terms and conditions of source code usage

# v0.14.2
## Fixed
- Fixed an issue where a CYP2D6 graph with no edges would lead to a panic during SVG graph visualization

# v0.14.1
## Changes
- Added a new output folder through the `--output-debug` option: `cyp2d6_igv_custom`. This folder contains an XML file describing an IGV session as well as the supporting data files to visualize full length alignments through the two constructed CYP2D6 haplotypes. For details, see the updated user docs.
- Released an updated database file: `data/v0.14.1/pbstarphase_20240826.json.gz`

## Fixed
- Fixed some off-by-one errors in the coordinates for miscellaneous extra regions in CYP2D6. These caused tiny overlaps in the output BAM file described above, but otherwise did not have an impact on diplotyping accuracy.

# v0.14.0
## Changes
- HLA allele labeling has been updated to improve 4th-field accuracy: When two potential definitions are compared, we now restrict the initial comparison to _only_ the shared regions of the two haplotype sequence definitions (this is often different, especially for DNA sequences). In the event of a tie, we revert to the full-length allele definitions.
- The HLA database configuration has been updated to include strand information for HLA genes. Defaults for _HLA-A_ and _HLA-B_ are set, so no database update is required. This modification will show in the next database release.
- HLA debug consensus outputs will now be output on the strand the gene is located to improve matching to IMGT/HLA sequences. For example, _HLA-A_ is already on the forward strand so no change will be made. In contrast, _HLA-B_ is on the reverse strand so the consensus sequences will be reverse complemented in the output FASTA file.
- **Breaking change**: _CYP2D6_ and the HLA genes now share a single debug BAM file through the `--output-debug` option: `debug_consensus.bam`
  - The previous debug file for _CYP2D6_, `cyp2d6_consensus.bam`, has been removed from the outputs. The mappings from this file have been moved into the new `debug_consensus.bam` file.
  - For both HLA genes, the BAM file contains alignments of the HLA consensus sequences and corresponding read sequences used to generate the consensus. Additionally, if the assigned haplotypes have DNA sequences in the database, those sequences are also aligned for comparison purposes.
  - Previously deprecated option `--debug-hla-target` has been repurposed to allow for specification of additional HLA haplotypes to get mapped in this debug BAM. As with the assigned haplotypes, these must have a DNA sequence in the database to get mapped.

## Fixed
- Fixed an issue where CDF filter was not filtering properly for HLA genes
- Fixed the CLI option syntax for `--hla-require-dna`
- Removed deprecated `--output-cyp2d6-bam` option from the list of CLI options, this is now part of the `--output-debug` files

# v0.13.3
## Fixed
- Replaced a panic with an error message when low coverage datasets fail to identify any CYP2D6 haplotypes to chain together. These will have a "NO_MATCH" diplotype in the results.
- Fixed a bug where duplicate consensus sequences in CYP2D6 could create a panic, duplicates are now flagged as FalseAlleles and ignored.

# v0.13.2
## Fixed
- Adjusted the alignment parameters for HLA mapping to reduce errors caused by soft-clipping of alignments near the end of a haplotype
- Replaced a panic with an error message when variants are found in unexpected states of zygosity and phase (e.g., phased homozygous)
- Debug messages for HLA calling have been adjusted to improve log reviewability

# v0.13.1
## Changes
- For HLA genes, StarPhase would previously ignore any HLA allele definitions that were missing a DNA sequence in the database. StarPhase now allows these partial HLA allele definitions by default.
- A new option was added to enable the previous behavior: `--hla-require-dna`. If this option is enabled, any HLA allele definition that is missing a DNA sequence will be ignored and never reported in StarPhase outputs.

## Fixed
- Fixed an issue where a _CYP2D6_ deletion allele (\*5) could be reported on the same haplotype as another allele. While this is biologically possible (e.g., deletion of one \*10 in a "\*10x2" haplotype), it is not considered a valid star-allele at this time. This combination will still show up in the debug log files, but it will get filtered in final reporting. For example: a "\*10+\*5" haplotype will now get reported as "\*10".

# v0.13.0
## Changes
- The algorithm for _HLA-A_ and _HLA-B_ has been modified to use a consensus-based approach to solve the alleles, a simpler version of the algorithm for _CYP2D6_.
  - CLI options related to consensus generation now control both HLA and _CYP2D6_ calling. These have been moved into a separate category on the CLI labeled "Consensus (HLA and CYP2D6)".
  - In internal tests, these changes slightly improved the accuracy of 4th-field entries in the HLA calls (2nd- and 3rd-field were unaffected). Additionally, the approach significantly reduced compute time requirements, averaging ~10% of CPU time required for v0.12.0.
  - With this change, the `--threads` option does not provide any benefit to the current algorithms. It has been deprecated, but may be added again if future optimizations allow it.
  - The `--max-error-rate` default has been adjusted for comparison to just the reference allele for each HLA gene, with a new default of 0.07 (previously 0.05).
  - Previous option `--min-allele-fraction` for HLA has been removed. The consensus option `--min-consensus-fraction` is used instead.
- Added a new option, `--output-debug`, that will create a debug folder with multiple additional files that are primarily for debugging the results from HLA and CYP2D6 calling, but may be useful for researchers. This folder is subject to change as the underlying methods develop. Some of the initial files included:
  - `consensus_{GENE}.fa` - Contains the full consensus sequences generated for a given `{GENE}`. Currently, this is only for HLA genes and _CYP2D6_.
  - `cyp2d6_consensus.bam` - Contains mapped substrings from the reads that were used to generate CYP2D6 consensus sequences. The phase set tag (PS) indicates which consensus the sequence was a part of. Useful for visualizing how the consensus ran and whether there are potential errors.
  - `cyp2d6_link_graph.svg` - A graphical representation of the connections present between CYP2D6 consensus segments.
  - `hla_debug.json` - Contains the summary mapping information of each database entry to the generated HLA consensus sequences.

## Fixed
- Fixed an issue with `build` where CPIC genes with no known chromosome would cause an error and exit. These entries are now ignored with a warning.
- Fixed an off-by-one error in the HLA gene region start coordinates. This has been corrected in the latest database release: `data/v0.13.0/pbstarphase_20240730.json.gz`

# v0.12.0
## Changes
- Hard-coded coordinates (GRCh38) for HLA-A, HLA-B, and CYP2D6 calling have been moved into the database file. Prior database versions without this new config information will automatically load the previously hard-coded values. This configuration is provided for transparency and experimentation in other reference coordinate systems, we do not recommend or support changing the provided default values.
- Released an updated database with the updated config format and name pattern change: `data/v0.12.0/pbstarphase_20240716.json.gz`
- Changed CYP2D6 consensus merging component to merge on sub-alleles instead of core alleles. 
    - To support the this change, consensuses with ambiguous CYP2D6 assignments (e.g. equal matches to "*4.001" and "*4.015") are labeled as unknown for the purpose of merging with similar consensus alleles prior to generating a final consensus set.
    - Internal tests showed this combination of changes led to increased sensitivity for sub-allele identification without injecting errors at the core allele or diplotype level.

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
