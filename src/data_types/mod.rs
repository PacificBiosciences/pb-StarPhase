
/// Contains serialization for generic alleles that only have a DNA definition
pub mod alleles;
/// Contains the coordinates functionality
pub mod coordinates;
/// Contains serialization for CPIC API result types
pub mod cpic_api_results;
/// Contains definitions related to our underlying database of genes -> alleles -> variants
pub mod database;
/// Contains gene definition information parsed primarily from RefSeq; coordinates, strand, exons
pub mod gene_definition;
/// Contains mapping stats for individual alignments
pub mod mapping;
/// Contains utility to normalize a variant to a standard definition that is non-ambiguous
pub mod normalized_variant;
/// Contains definitions related to the representation of a final diplotype
pub mod pgx_diplotypes;
