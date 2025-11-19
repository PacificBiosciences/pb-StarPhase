
use serde::Deserialize;

// PharmVar API quickstart: https://www.pharmvar.org/documentation
// Example query: https://www.pharmvar.org/api-service/genes/NAT2

/// This captures a full PharmVar gene definition, we only parse the elements we need though
#[derive(Clone, Debug, Deserialize)]
pub struct PharmvarGeneDefinition {
    /// The gene for this definition, good for sanity checking mostly
    #[serde(alias = "geneSymbol")]
    pub gene_symbol: String,
    /// the list of alleles
    pub alleles: Vec<PharmvarAlleleDefinition>
}

/// This captures a full PharmVar allele definition, we only parse the elements we need though
#[derive(Clone, Debug, Deserialize)]
pub struct PharmvarAlleleDefinition {
    /// The gene for this definition, good for sanity checking mostly
    #[serde(alias = "geneSymbol")]
    pub gene_symbol: String,
    /// The allele name (usually star-allele)
    #[serde(alias = "alleleName")]
    pub allele_name: String,
    /// If set, then this is a sub-allele and the value is the core allele
    #[serde(alias = "coreAllele")]
    pub core_allele: Option<String>,
    /// Seems to be either "Core" or "Sub"
    #[serde(alias = "alleleType")]
    pub allele_type: String,
    /// List of variants included
    pub variants: Vec<PharmvarVariantDefinition>
}

impl PharmvarAlleleDefinition {
    /// Returns the star allele associated with this allele.
    /// The PharmVar defs usually have the gene name on them, this just strips that after checking.
    pub fn star_allele(&self) -> &str {
        if self.allele_name.starts_with(&self.gene_symbol) {
            &self.allele_name[self.gene_symbol.len()..]
        } else {
            &self.allele_name
        }
    }

    /// Returns the core allele associated with this allele, if any.
    /// The PharmVar defs usually have the gene name on them, this just strips that after checking.
    pub fn core_allele(&self) -> Option<&str> {
        self.core_allele.as_deref().map(|s| {
            if s.starts_with(&self.gene_symbol) {
                &s[self.gene_symbol.len()..]
            } else {
                s
            }
        })
    }
}

/// This captures a full PharmVar variant definition, we only parse the elements we need though
#[derive(Clone, Debug, Deserialize)]
pub struct PharmvarVariantDefinition {
    /// The reference sequence, but in annoying form: "NC_000008.11"
    #[serde(alias = "referenceSequence")]
    pub ref_sequence: String,
    /// HGVS nomenclature variant, but the only info we have on REF/ALT sequence; we'll have to parse it
    pub hgvs: String,
    /// RS ID when available
    #[serde(alias = "rsId")]
    pub rs_id: Option<String>,
    /// Impact value, all core allele variants are expect to have these set
    pub impact: Option<String>,
    /// Looks like an internal variant ID number; stored as String but looks like an Integer
    #[serde(alias = "variantId")]
    pub variant_id: String,
    /// A more parseable version of HGVS; everything here will be 1-based
    pub position: String
}
