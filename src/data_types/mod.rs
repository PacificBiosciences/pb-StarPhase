
/// Contains serialization for generic alleles that only have a DNA definition
pub mod alleles;
/// Contains the coordinates functionality
pub mod coordinates;
/// Contains serialization for CPIC API result types
pub mod cpic_api_results;
/// Contains definitions related to our underlying database of genes -> alleles -> variants
pub mod database;
/// Constants that are hard-coded and typically written to the database file
pub mod db_const;
/// Contains gene definition information parsed primarily from RefSeq; coordinates, strand, exons
pub mod gene_definition;
/// Contains HGVS parser and formated output
pub mod hgvs;
/// Contains mapping stats for individual alignments
pub mod mapping;
/// Contains utility to normalize a variant to a standard definition that is non-ambiguous
pub mod normalized_variant;
/// Contains serialization for PharmVar API result types
pub mod pharmvar_api_results;
/// Contains definitions related to the representation of a final diplotype
pub mod pgx_diplotypes;
/// Contains functionality for support CPIC structural variants
pub mod pgx_structural_variants;
