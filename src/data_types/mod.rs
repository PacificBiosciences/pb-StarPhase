
/// Contains serialization for generic alleles that only have a DNA definition
pub mod alleles;
/// Contains the coordinates functionality
pub mod coordinates;
/// Contains HGVS parser and formated output
pub mod hgvs;
/// Contains mapping stats for individual alignments
pub mod mapping;
/// Contains utility to normalize a variant to a standard definition that is non-ambiguous
pub mod normalized_variant;
/// Contains the data types for simple and complex PGx diplotypes and haplotypes
pub mod pgx_diplotype;
/// Contains region variant labeling and relationship types
pub mod region_variants;
/// Contains the final JSON output of the diplotype caller and related data types
pub mod starphase_json;
