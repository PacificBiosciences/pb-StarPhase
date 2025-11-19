
use serde::Deserialize;

// CPIC API quickstart: https://github.com/cpicpgx/cpic-data/wiki
// CPI API full book: https://documenter.getpostman.com/view/1446428/Szt78VUJ?version=latest
// Useful postgrest reference: https://postgrest.org/en/v7.0.0/api.html#horizontal-filtering-rows

/// This captures a full CPIC definition, we only parse the elements we need though
#[derive(Debug, Deserialize)]
pub struct CpicAlleleDefinition {
    /// The gene for this definition, good for sanity checking mostly
    #[serde(alias = "genesymbol")]
    pub gene_symbol: String,
    /// The name of this allele
    #[serde(alias = "name")]
    pub allele_name: String,
    /// True if this is the reference allele
    #[serde(alias = "matchesreferencesequence")]
    pub is_reference: bool,
    /// True if this is a structural variant allele, which we will ignore
    #[serde(alias = "structuralvariation")]
    pub is_sv: bool,
    /// the list of variants for this allele
    #[serde(alias = "allele_location_value")]
    pub variants: Vec<CpicVariantDefinition>
}

/// This captures a CPIC variant allele definition
#[derive(Debug, Deserialize)]
pub struct CpicVariantDefinition {
    /// The allele for this variant
    #[serde(alias = "variantallele")]
    pub variant_allele: String,
    /// The sequence location definition for the variant
    pub sequence_location: CpicSequenceLocationDefinition
}

/// This captures a specific CPIC 
#[derive(Debug, Deserialize)]
pub struct CpicSequenceLocationDefinition {
    /// The unique ID for this variant
    pub id: u64,
    /// The name of this variant
    pub name: String,
    /// This is typically the g. for the variant
    #[serde(alias = "chromosomelocation")]
    pub gdot: String,
    /// The DBSNP ID (e.g. rs##)
    #[serde(alias = "dbsnpid")]
    pub dbsnp_id: Option<String>,
    /// The build 38 position of the variant, 1-based
    pub position: usize
}
