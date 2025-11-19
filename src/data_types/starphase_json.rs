
use serde::Serialize;
use simple_error::bail;
use std::collections::{BTreeMap, BTreeSet};
use std::collections::btree_map::Entry::{Occupied, Vacant};

use crate::data_types::normalized_variant::{NormalizedVariant, NormalizedGenotype};
use crate::data_types::pgx_diplotype::{Diplotype, InexactDiplotype};
use crate::database::pgx_database::PgxMetadata;
use crate::hla::mapping::HlaMappingStats;

/// Intended to be serialized to JSON as the final result
#[derive(Debug, Serialize)]
pub struct StarphaseJson {
    /// Version of the tool that generated the calls
    pbstarphase_version: String,
    /// Metadata for the database
    database_metadata: PgxMetadata,
    /// Map from gene name to diplotype call
    gene_details: BTreeMap<String, PgxGeneDetails>
}

impl StarphaseJson {
    /// Basic constructor, will perform sanity checks if necessary
    pub fn new(database_metadata: PgxMetadata) -> Self {
        Self {
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
    /// Contains the list of inexact matching diplotypes; only used if diplotypes is empty and for the non-mapping-based genes
    inexact_diplotypes: Option<Vec<InexactDiplotype>>,
    /// Contains the list of identified variants
    variant_details: Option<Vec<PgxVariantDetails>>,
    /// Contains the list of alignments analyzed
    mapping_details: Option<Vec<PgxMappingDetails>>,
    // Contains the list of multi-mapping alignments analyzed
    multi_mapping_details: Option<Vec<PgxMultiMappingDetails>>
}

impl PgxGeneDetails {
    /// Basic constructor for wrapping our details up into a nice bundle.
    /// This is for exact matches down to the sub-allele level.
    /// # Arguments
    /// * `diplotypes` - the list of exact matching diplotypes
    /// * `simple_diplotypes` - the list of simple diplotypes, or core alleles
    /// * `variant_details` - the list of variant details
    pub fn new_suballele_match(diplotypes: Vec<Diplotype>, simple_diplotypes: Option<Vec<Diplotype>>, variant_details: Vec<PgxVariantDetails>) -> Result<PgxGeneDetails, Box<dyn std::error::Error>> {
        if let Some(sd) = simple_diplotypes.as_ref() {
            if diplotypes.len() != sd.len() {
                bail!("diplotypes and simple_diplotypes must be the same length");
            }
        }
        Ok(PgxGeneDetails { 
            diplotypes,
            simple_diplotypes,
            inexact_diplotypes: None,
            variant_details: Some(variant_details),
            mapping_details: None,
            multi_mapping_details: None
        })
    }

    /// This is for exact matches down to the core allele level.
    /// We have to tracked extended diplotypes here because there is not a sub-allele match.
    /// # Arguments
    /// * `diplotypes` - the list of main diplotypes
    /// * `inexact_diplotypes` - the list of inexact matching diplotypes
    /// * `simple_diplotypes` - the list of simple diplotypes, or core alleles
    /// * `variant_details` - the list of variant details
    pub fn new_core_match(diplotypes: Vec<Diplotype>, inexact_diplotypes: Vec<InexactDiplotype>, simple_diplotypes: Option<Vec<Diplotype>>, variant_details: Vec<PgxVariantDetails>) -> Result<PgxGeneDetails, Box<dyn std::error::Error>> {
        if let Some(sd) = simple_diplotypes.as_ref() {
            if diplotypes.len() != sd.len() {
                bail!("diplotypes and simple_diplotypes must be the same length");
            }
        }
        if inexact_diplotypes.len() != diplotypes.len() {
            bail!("diplotypes and inexact_diplotypes must be the same length");
        }
        Ok(PgxGeneDetails { 
            diplotypes,
            simple_diplotypes,
            inexact_diplotypes: Some(inexact_diplotypes),
            variant_details: Some(variant_details),
            mapping_details: None,
            multi_mapping_details: None
        })
    }

    /// Constructor for when we have inexact matching diplotypes
    /// # Arguments
    /// * `inexact_diplotypes` - the list of inexact matching diplotypes
    /// * `variant_details` - the list of variant details
    pub fn new_inexact_diplotypes(inexact_diplotypes: Vec<InexactDiplotype>, variant_details: Vec<PgxVariantDetails>) -> Result<PgxGeneDetails, Box<dyn std::error::Error>> {
        // no exact matches, so we just use a NO_MATCH diplotype
        let diplotypes = vec![
            Diplotype::new("NO_MATCH", "NO_MATCH")
        ];

        Ok(PgxGeneDetails {
            diplotypes,
            simple_diplotypes: None,
            inexact_diplotypes: Some(inexact_diplotypes),
            variant_details: Some(variant_details),
            mapping_details: None,
            multi_mapping_details: None
        })
    }

    /// This is the one for HLA (currently, it will probably change to multi-mapping when we synchronize the approaches).
    /// # Arguments
    /// * `diplotypes` - the list of exact matching diplotypes
    /// * `simple_diplotypes` - the list of simple diplotypes, which is always None for HLA currently
    /// * `mapping_details` - the list of mapping details
    pub fn new_from_mappings(diplotypes: Vec<Diplotype>, simple_diplotypes: Option<Vec<Diplotype>>, mapping_details: Vec<PgxMappingDetails>) -> Result<PgxGeneDetails, Box<dyn std::error::Error>> {
        if let Some(sd) = simple_diplotypes.as_ref() {
            if diplotypes.len() != sd.len() {
                bail!("diplotypes and simple_diplotypes must be the same length");
            }
        }
        Ok(PgxGeneDetails { 
            diplotypes,
            simple_diplotypes,
            inexact_diplotypes: None,
            variant_details: None,
            mapping_details: Some(mapping_details),
            multi_mapping_details: None
        })
    }

    /// This is the one for CYP2D6. Here, we expect all of the diplotype formats to be present.
    /// # Arguments
    /// * `diplotypes` - the list of exact matching diplotypes
    /// * `simple_diplotypes` - the list of simple diplotypes
    /// * `inexact_diplotypes` - the list of inexact matching diplotypes
    /// * `multi_mapping_details` - the list of multi-mapping details
    pub fn new_from_multi_mappings(
        diplotypes: Vec<Diplotype>,
        simple_diplotypes: Option<Vec<Diplotype>>,
        inexact_diplotypes: Option<Vec<InexactDiplotype>>,
        multi_mapping_details: Vec<PgxMultiMappingDetails>
    ) -> Result<PgxGeneDetails, Box<dyn std::error::Error>> {
        if let Some(sd) = simple_diplotypes.as_ref() {
            if diplotypes.len() != sd.len() {
                bail!("diplotypes and simple_diplotypes must be the same length");
            }
        }
        Ok(PgxGeneDetails { 
            diplotypes,
            simple_diplotypes,
            inexact_diplotypes,
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
            inexact_diplotypes: None,
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

    /// This will deduplicate the simple diplotypes, which are always core alleles.
    /// This is useful because we may have multiple inexact diplotypes that resolve to the same core alleles.
    pub fn dedup_simple_diplotypes(&self) -> Vec<Diplotype> {
        let dedup_set: BTreeSet<Diplotype> = self.simple_diplotypes().iter()
            .cloned()
            .collect();
        dedup_set.into_iter().collect()
    }

    pub fn inexact_diplotypes(&self) -> Option<&[InexactDiplotype]> {
        self.inexact_diplotypes.as_deref()
    }

    pub fn variant_details(&self) -> Option<&[PgxVariantDetails]> {
        self.variant_details.as_deref()
    }

    pub fn mapping_details(&self) -> Option<&[PgxMappingDetails]> {
        self.mapping_details.as_deref()
    }
}

/// Contains all the details for a variant that was identified through our process
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct PgxVariantDetails {
    /// Assigned variant ID
    variant_id: u64,
    /// Assigned name - CPIC has names, but PharmVar we use RS IDs and default to variant ID if not available
    variant_name: String,
    /// DBSNP id when available
    dbsnp: Option<String>,
    /// The normalized variant we loaded
    normalized_variant: NormalizedVariant,
    /// The normalized genotype we loaded
    normalized_genotype: NormalizedGenotype,
    /// If true, this is a core variant
    is_core_variant: bool
}

impl PgxVariantDetails {
    pub fn new(
        variant_id: u64, variant_name: String, dbsnp: Option<String>,
        normalized_variant: NormalizedVariant, normalized_genotype: NormalizedGenotype, is_core_variant: bool
    ) -> PgxVariantDetails {
        PgxVariantDetails { 
            variant_id,
            variant_name,
            dbsnp, 
            normalized_variant, 
            normalized_genotype,
            is_core_variant
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

    use crate::data_types::pgx_diplotype::InexactHaplotype;
    use crate::data_types::region_variants::{RegionVariant, VariantAlleleRelationship};

    #[test]
    fn test_starphase_json() {
        let hap1 = "A";
        let hap2 = "B";
        let diplotype = vec![Diplotype::new(hap2, hap1)];

        let mut diplotypes = StarphaseJson::new(Default::default());
        let gene_details = PgxGeneDetails::new_suballele_match(diplotype, None, vec![]).unwrap();
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

        let mut diplotypes = StarphaseJson::new(Default::default());
        let gene_details = PgxGeneDetails::new_suballele_match(diplotype, None, vec![]).unwrap();
        diplotypes.insert("CACN1S".to_string(), gene_details.clone()).unwrap();
        diplotypes.insert("CACN1S".to_string(), gene_details).unwrap();
    }

    #[test]
    fn test_new_from_mappings() {
        let hap1 = "A";
        let hap2 = "B";
        let diplotype = vec![Diplotype::new(hap2, hap1)];

        let mut diplotypes = StarphaseJson::new(Default::default());
        let gene_details = PgxGeneDetails::new_from_mappings(diplotype, None, vec![]).unwrap();
        assert!(gene_details.variant_details.is_none());
        diplotypes.insert("HLA-A".to_string(), gene_details.clone()).unwrap();

        let map = &diplotypes.gene_details;
        assert_eq!(map.len(), 1);
        assert_eq!(map.get("HLA-A").unwrap(), &gene_details);
    }

    #[test]
    fn test_new_inexact_diplotypes() {
        // Create some test diplotypes for inexact matching
        // the exact values do not matter, we just need to make sure the constructor works
        let inexact_diplotypes = vec![
            InexactDiplotype::new(
                InexactHaplotype::new("*1".to_string(), BTreeSet::from([
                    RegionVariant::new("test_variant_1".to_string(), true, VariantAlleleRelationship::Match)
                ])),
                InexactHaplotype::new("*2".to_string(), BTreeSet::from([
                    RegionVariant::new("test_variant_2".to_string(), true, VariantAlleleRelationship::Match)
                ]))
            ),
            InexactDiplotype::new(
                InexactHaplotype::new("*3".to_string(), BTreeSet::from([
                    RegionVariant::new("test_variant_3".to_string(), true, VariantAlleleRelationship::Match)
                ])),
                InexactHaplotype::new("*4".to_string(), BTreeSet::from([
                    RegionVariant::new("test_variant_4".to_string(), true, VariantAlleleRelationship::Match)
                ]))
            )
        ];

        // Create some test variant details
        let variant_details = vec![
            PgxVariantDetails::new(
                12345,
                "test_variant_1".to_string(),
                Some("rs123456".to_string()),
                NormalizedVariant::new("chr1".to_string(), 1000, "A", "T", None).unwrap(),
                NormalizedGenotype::new(crate::data_types::normalized_variant::Genotype::HeterozygousUnphased, None),
                true
            ),
            PgxVariantDetails::new(
                67890,
                "test_variant_2".to_string(),
                None,
                NormalizedVariant::new("chr1".to_string(), 2000, "C", "G", None).unwrap(),
                NormalizedGenotype::new(crate::data_types::normalized_variant::Genotype::HomozygousAlternate, None),
                true
            )
        ];

        // Test the new_inexact_diplotypes constructor
        let gene_details = PgxGeneDetails::new_inexact_diplotypes(inexact_diplotypes.clone(), variant_details.clone()).unwrap();

        // Verify the structure
        assert_eq!(gene_details.diplotypes.len(), 1);
        assert_eq!(gene_details.diplotypes[0].diplotype(), "NO_MATCH/NO_MATCH");

        assert!(gene_details.simple_diplotypes.is_none());

        assert!(gene_details.inexact_diplotypes.is_some());
        assert_eq!(gene_details.inexact_diplotypes.as_ref().unwrap(), &inexact_diplotypes);

        assert!(gene_details.variant_details.is_some());
        assert_eq!(gene_details.variant_details.as_ref().unwrap(), &variant_details);

        assert!(gene_details.mapping_details.is_none());
        assert!(gene_details.multi_mapping_details.is_none());
    }

    #[test]
    fn test_new_core_match() {
        // Create test diplotypes for core matching
        let diplotypes = vec![
            Diplotype::new("*1", "*2")
        ];

        // Create test inexact diplotypes (should match length of diplotypes)
        let inexact_diplotypes = vec![
            InexactDiplotype::new(
                InexactHaplotype::new("*1".to_string(), BTreeSet::from([
                    RegionVariant::new("test_variant_1".to_string(), true, VariantAlleleRelationship::Match),
                    RegionVariant::new("test_variant_2".to_string(), true, VariantAlleleRelationship::Unexpected)
                ])),
                InexactHaplotype::new("*2".to_string(), BTreeSet::from([
                    RegionVariant::new("test_variant_2".to_string(), true, VariantAlleleRelationship::Match)
                ]))
            )
        ];

        // Create test variant details
        let variant_details = vec![
            PgxVariantDetails::new(
                12345,
                "test_variant_1".to_string(),
                Some("rs123456".to_string()),
                NormalizedVariant::new("chr1".to_string(), 1000, "A", "T", None).unwrap(),
                NormalizedGenotype::new(crate::data_types::normalized_variant::Genotype::HeterozygousUnphased, None),
                true
            ),
            PgxVariantDetails::new(
                67890,
                "test_variant_2".to_string(),
                None,
                NormalizedVariant::new("chr1".to_string(), 2000, "C", "G", None).unwrap(),
                NormalizedGenotype::new(crate::data_types::normalized_variant::Genotype::HomozygousAlternate, None),
                true
            )
        ];

        // Test the new_core_match constructor
        let gene_details = PgxGeneDetails::new_core_match(
            diplotypes.clone(),
            inexact_diplotypes.clone(),
            Some(diplotypes.clone()),
            variant_details.clone()
        ).unwrap();

        // Verify the structure
        assert_eq!(gene_details.diplotypes, diplotypes);
        assert_eq!(gene_details.simple_diplotypes, Some(diplotypes));
        assert_eq!(gene_details.inexact_diplotypes, Some(inexact_diplotypes));
        assert_eq!(gene_details.variant_details, Some(variant_details));
        assert!(gene_details.mapping_details.is_none());
        assert!(gene_details.multi_mapping_details.is_none());
    }
}
