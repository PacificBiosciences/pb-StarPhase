
/// Contains serialization for CPIC API result types
pub mod cpic_api_results;
/// Contains builder functionality for database construction
pub mod db_config;
/// Constants that are hard-coded and typically written to the database file
pub mod db_const;
/// Contains gene definition information parsed primarily from RefSeq; coordinates, strand, exons
pub mod gene_definition;
/// Contains serialization for PharmVar API result types
pub mod pharmvar_api_results;
/// Contains definitions related to our underlying database of genes -> alleles -> variants
pub mod pgx_database;
/// Contains functionality for support CPIC structural variants
pub mod pgx_structural_variants;
