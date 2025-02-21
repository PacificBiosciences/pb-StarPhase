
use serde::Serialize;
use simple_error::bail;
use std::collections::BTreeMap;
use std::collections::btree_map::Entry::{Occupied, Vacant};

use crate::data_types::database::PgxMetadata;
use crate::data_types::normalized_variant::{NormalizedVariant, NormalizedGenotype};
use crate::hla::mapping::HlaMappingStats;

/// Intended to be serialized to JSON as the final result
#[derive(Debug, Serialize)]
pub struct PgxDiplotypes {
    /// Version of the tool that generated the calls
    pbstarphase_version: String,
    /// Metadata for the database
    database_metadata: PgxMetadata,
    /// Map from gene name to diplotype call
    gene_details: BTreeMap<String, PgxGeneDetails>
}

impl PgxDiplotypes {
    /// Basic constructor, will perform sanity checks if necessary
    pub fn new(database_metadata: PgxMetadata) -> PgxDiplotypes {
        PgxDiplotypes {
            pbstarphase_version: crate::cli::core::FULL_VERSION.to_string(),
            database_metadata,
            gene_details: Default::default()
        }
    }

    /// Simple wrapper for our diplotype insertion to make sure we do not double insert
    /// # Arguments
    /// * `gene` - the gene name we are saving the diplotype for
    /// * `diplotype` - the diplotype call getting saved
    pub fn insert(&mut self, gene: String, diplotype: PgxGeneDetails) -> Result<(), Box<dyn std::error::Error>> {
        match self.gene_details.entry(gene) {
            Vacant(entry) => entry.insert(diplotype),
            Occupied(entry) => bail!("Entry for {} is already occupied.", entry.key())
        };
        Ok(())
    }

    pub fn database_metadata(&self) -> &PgxMetadata {
        &self.database_metadata
    }

    pub fn gene_details(&self) -> &BTreeMap<String, PgxGeneDetails> {
        &self.gene_details
    }
}

/// Wrapper for all of the details for a single gene
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct PgxGeneDetails {
    /// Contains the list of exactly matching diplotypes
    diplotypes: Vec<Diplotype>,
    /// Contains an optional list of simplified diplotypes
    simple_diplotypes: Option<Vec<Diplotype>>,
    /// Contains the list of identified variants
    variant_details: Option<Vec<PgxVariantDetails>>,
    /// Contains the list of alignments analyzed
    mapping_details: Option<Vec<PgxMappingDetails>>,
    // Contains the list of multi-mapping alignments analyzed
    multi_mapping_details: Option<Vec<PgxMultiMappingDetails>>
}

impl PgxGeneDetails {
    /// Basic constructor for wrapping our details up into a nice bundle
    pub fn new(diplotypes: Vec<Diplotype>, simple_diplotypes: Option<Vec<Diplotype>>, variant_details: Vec<PgxVariantDetails>) -> Result<PgxGeneDetails, Box<dyn std::error::Error>> {
        if let Some(sd) = simple_diplotypes.as_ref() {
            if diplotypes.len() != sd.len() {
                bail!("diplotypes and simple_diplotypes must be the same length");
            }
        }
        Ok(PgxGeneDetails { 
            diplotypes,
            simple_diplotypes,
            variant_details: Some(variant_details),
            mapping_details: None,
            multi_mapping_details: None
        })
    }

    /// This is the one for HLA (currently, it will probably change to multi-mapping when we synchronize the approaches)
    pub fn new_from_mappings(diplotypes: Vec<Diplotype>, simple_diplotypes: Option<Vec<Diplotype>>, mapping_details: Vec<PgxMappingDetails>) -> Result<PgxGeneDetails, Box<dyn std::error::Error>> {
        if let Some(sd) = simple_diplotypes.as_ref() {
            if diplotypes.len() != sd.len() {
                bail!("diplotypes and simple_diplotypes must be the same length");
            }
        }
        Ok(PgxGeneDetails { 
            diplotypes,
            simple_diplotypes,
            variant_details: None,
            mapping_details: Some(mapping_details),
            multi_mapping_details: None
        })
    }

    /// This is the one for CYP2D6
    pub fn new_from_multi_mappings(diplotypes: Vec<Diplotype>, simple_diplotypes: Option<Vec<Diplotype>>, multi_mapping_details: Vec<PgxMultiMappingDetails>) -> Result<PgxGeneDetails, Box<dyn std::error::Error>> {
        if let Some(sd) = simple_diplotypes.as_ref() {
            if diplotypes.len() != sd.len() {
                bail!("diplotypes and simple_diplotypes must be the same length");
            }
        }
        Ok(PgxGeneDetails { 
            diplotypes,
            simple_diplotypes,
            variant_details: None,
            mapping_details: None,
            multi_mapping_details: Some(multi_mapping_details)
        })
    }

    /// Generic wrapper function to create a "NO_MATCH" result, usually from some algorithm failure.
    /// This does not include additional details.
    pub fn no_match() -> PgxGeneDetails {
        let diplotypes = vec![
            Diplotype::new("NO_MATCH", "NO_MATCH")
        ];
        PgxGeneDetails {
            diplotypes,
            simple_diplotypes: None,
            variant_details: None,
            mapping_details: None,
            multi_mapping_details: None
        }
    }

    pub fn diplotypes(&self) -> &[Diplotype] {
        &self.diplotypes
    }

    pub fn simple_diplotypes(&self) -> &[Diplotype] {
        if let Some(sd) = self.simple_diplotypes.as_ref() {
            sd
        } else {
            self.diplotypes()
        }
    }

    pub fn variant_details(&self) -> Option<&[PgxVariantDetails]> {
        self.variant_details.as_deref()
    }

    pub fn mapping_details(&self) -> Option<&[PgxMappingDetails]> {
        self.mapping_details.as_deref()
    }
}

/// Contains all the information related to a single gene's diplotype result
#[derive(Clone, Debug, Serialize)]
pub struct Diplotype {
    /// short string for haplotype 1
    hap1: String,
    /// short string for haplotype 2
    hap2: String,
    /// combination diplotype call
    diplotype: String
    // TODO: we will likely include deeper evidence information in the future, e.g. variant calls, counts, etc.
}

impl Diplotype {
    pub fn new(hap1: &str, hap2: &str) -> Diplotype {
        Diplotype {
            hap1: hap1.to_string(),
            hap2: hap2.to_string(),
            diplotype: format!("{}/{}", hap1, hap2)
        }
    }

    /// If homozygous, return the single haplotype
    pub fn homozygous_haplotype(&self) -> Option<&str> {
        if self.hap1 == self.hap2 {
            Some(&self.hap1)
        } else {
            None
        }
    }

    pub fn diplotype(&self) -> &str {
        &self.diplotype
    }

    /// Returns a PharmCAT formatted diplotype. 
    /// See https://pharmcat.org/using/Outside-Call-Format/#diplotypes for more details.
    pub fn pharmcat_diplotype(&self) -> String {
        let h1 = if self.hap1.contains('+') {
            format!("[{}]", self.hap1)
        } else {
            self.hap1.clone()
        };
        let h2 = if self.hap2.contains('+') {
            format!("[{}]", self.hap2)
        } else {
            self.hap2.clone()
        };
        format!("{h1}/{h2}")
    }
}

impl PartialEq for Diplotype {
    fn eq(&self, other: &Self) -> bool {
        // this allows for a swap in hap1/hap2 and we still report identity
        (self.hap1 == other.hap1 && self.hap2 == other.hap2) ||
            (self.hap1 == other.hap2 && self.hap2 == other.hap1)
    }
}

/// Contains all the details for a variant that was identified through our process
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct PgxVariantDetails {
    /// CPIC assigned variant ID
    cpic_variant_id: u64,
    /// CPIC assigned name
    cpic_name: String,
    /// DBSNP id when available
    dbsnp: Option<String>,
    /// The normalized variant we loaded
    normalized_variant: NormalizedVariant,
    /// The normalized genotype we loaded
    normalized_genotype: NormalizedGenotype
}

impl PgxVariantDetails {
    pub fn new(cpic_variant_id: u64, cpic_name: String, dbsnp: Option<String>, normalized_variant: NormalizedVariant, normalized_genotype: NormalizedGenotype) -> PgxVariantDetails {
        PgxVariantDetails { 
            cpic_variant_id,
            cpic_name,
            dbsnp, 
            normalized_variant, 
            normalized_genotype
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct PgxMappingDetails {
    /// Read mapping ID
    read_qname: String,
    /// The ID of the best matching star allele from the database
    best_hla_id: String,
    /// The star allele string for the best match
    best_star_allele: String,
    /// The mapping stats for the best match
    best_mapping_stats: HlaMappingStats,
    /// If true, this mapping was ignored due to high error rate
    is_ignored: bool
}

impl PgxMappingDetails {
    pub fn new(read_qname: String, best_hla_id: String, best_star_allele: String, best_mapping_stats: HlaMappingStats, is_ignored: bool) -> PgxMappingDetails {
        PgxMappingDetails {
            read_qname,
            best_hla_id,
            best_star_allele,
            best_mapping_stats,
            is_ignored
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct PgxMultiMappingDetails {
    /// Read mapping ID
    read_qname: String,
    /// Coordinates in the read corresponding to the extracted sequence
    read_position: std::ops::Range<usize>,
    /// The final consensus ID this was assigned to
    consensus_id: usize,
    /// The final star allele this was assigned to 
    consensus_star_allele: String
}

impl PgxMultiMappingDetails {
    pub fn new(
        read_qname: String, read_position: std::ops::Range<usize>, consensus_id: usize, consensus_star_allele: String
    ) -> PgxMultiMappingDetails {
        PgxMultiMappingDetails {
            read_qname, read_position, consensus_id, consensus_star_allele
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pgx_diplotypes() {
        let hap1 = "A";
        let hap2 = "B";
        let diplotype = vec![Diplotype::new(hap2, hap1)];

        let mut diplotypes = PgxDiplotypes::new(Default::default());
        let gene_details = PgxGeneDetails::new(diplotype, None, vec![]).unwrap();
        assert!(gene_details.mapping_details.is_none());
        diplotypes.insert("CACNA1S".to_string(), gene_details.clone()).unwrap();

        let map = &diplotypes.gene_details;
        assert_eq!(map.len(), 1);
        assert_eq!(map.get("CACNA1S").unwrap(), &gene_details);
    }

    #[test]
    #[should_panic]
    fn test_duplicate_diplotype() {
        let hap1 = "A";
        let hap2 = "B";
        let diplotype = vec![Diplotype::new(hap2, hap1)];

        let mut diplotypes = PgxDiplotypes::new(Default::default());
        let gene_details = PgxGeneDetails::new(diplotype, None, vec![]).unwrap();
        diplotypes.insert("CACN1S".to_string(), gene_details.clone()).unwrap();
        diplotypes.insert("CACN1S".to_string(), gene_details).unwrap();
    }

    #[test]
    fn test_new_from_mappings() {
        let hap1 = "A";
        let hap2 = "B";
        let diplotype = vec![Diplotype::new(hap2, hap1)];

        let mut diplotypes = PgxDiplotypes::new(Default::default());
        let gene_details = PgxGeneDetails::new_from_mappings(diplotype, None, vec![]).unwrap();
        assert!(gene_details.variant_details.is_none());
        diplotypes.insert("HLA-A".to_string(), gene_details.clone()).unwrap();

        let map = &diplotypes.gene_details;
        assert_eq!(map.len(), 1);
        assert_eq!(map.get("HLA-A").unwrap(), &gene_details);
    }

    #[test]
    fn test_diplotype() {
        // this is a basic test for now, will likely get more complicated over time
        let hap1 = "A";
        let hap2 = "B";
        let diplotype = Diplotype::new(hap2, hap1);
        assert_eq!(diplotype.diplotype(), "B/A");
    }

    #[test]
    fn test_pharmcat_diplotype() {
        let diplotype = Diplotype::new("*4", "*1");
        assert_eq!(diplotype.pharmcat_diplotype(), "*4/*1");

        let diplotype = Diplotype::new("*4x2", "*1");
        assert_eq!(diplotype.pharmcat_diplotype(), "*4x2/*1");

        let diplotype = Diplotype::new("*4 + *68", "*1");
        assert_eq!(diplotype.pharmcat_diplotype(), "[*4 + *68]/*1");
    }
}
