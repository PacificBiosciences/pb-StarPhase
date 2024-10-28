
use serde::{Deserialize, Serialize};
use simple_error::{SimpleError, bail};
use std::collections::BTreeMap;

use crate::data_types::coordinates::Coordinates;

/// Contains the configuration for the HLA database.
/// These are generally coordinates that we do not expect to change except between reference builds.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct HlaConfig {
    /// High-level coordinates of the CYP2D regions
    hla_coordinates: BTreeMap<String, Coordinates>,
    /// Strand orientation for the gene relative to the reference genome
    #[serde(default="HlaConfig::default_strand")]
    hla_is_forward_strand: BTreeMap<String, bool>,
    /// Specific subregion, like exons
    hla_exons: BTreeMap<String, Vec<Coordinates>>,
}

impl HlaConfig {
    /// This function should be called after loading a config to verify that everything required to run the algorithms is present.
    pub fn validate_config(&self) -> Result<(), SimpleError> {
        // make sure all expected regions are defined
        let expected_hla_coordinates = [
            "HLA-A", "HLA-B"
        ];
        for &k in expected_hla_coordinates.iter() {
            if !self.hla_coordinates.contains_key(k) {
                bail!("Coordinates for \"{}\" were not found in provided hla_coordinates.", k);
            }
        }

        // make sure all exon regions are defined
        let expected_num_exons = 8;
        for &k in expected_hla_coordinates.iter() {
            if !self.hla_exons.contains_key(k) {
                bail!("Data for \"{}\" was not found in provided hla_exons.", k);
            }
            let exons = self.hla_exons.get(k).unwrap();
            if exons.len() != expected_num_exons {
                bail!("Found {} exons for \"{}\", expected {}.", exons.len(), k, expected_num_exons);
            }
        }

        Ok(())
    }

    // getters
    pub fn hla_coordinates(&self) -> &BTreeMap<String, Coordinates> {
        &self.hla_coordinates
    }

    pub fn hla_is_forward_strand(&self) -> &BTreeMap<String, bool> {
        &self.hla_is_forward_strand
    }

    pub fn hla_exons(&self) -> &BTreeMap<String, Vec<Coordinates>> {
        &self.hla_exons
    }

    // Defaults for strand orientations. Having this as a separate functions allows us to do an update without requiring a DB update.
    fn default_strand() -> BTreeMap<String, bool> {
        let mut hla_is_forward_strand: BTreeMap<String, bool> = Default::default();
        hla_is_forward_strand.insert("HLA-A".to_string(), true);
        hla_is_forward_strand.insert("HLA-B".to_string(), false);
        hla_is_forward_strand
    }
}

impl Default for HlaConfig {
    fn default() -> Self {
        let mut hla_coordinates: BTreeMap<String, Coordinates> = Default::default();
        let preshift = 1; // coordinate below are from UCSC and not 0-base shifted
        let postshift = 0;
        // HLA-A
        // chr6:29,942,532-29,945,870 from UCSC browser
        // 03:01:01:01 blats to chr6:29942254-29945755, partially updated below (old end was higher); TODO: systematically solve this?
        hla_coordinates.insert("HLA-A".to_string(), Coordinates::new("chr6".to_string(), 29942254 - preshift, 29945870 - postshift));
        // HLA-B
        // chr6:31,353,875-31,357,179 from UCSC browser
        // 08:01:01:01 blats to chr6:31353362-31357442, updated below; TODO: systematically solve this?
        hla_coordinates.insert("HLA-B".to_string(), Coordinates::new("chr6".to_string(), 31353362 - preshift, 31357442 - postshift));

        let hla_is_forward_strand: BTreeMap<String, bool> = HlaConfig::default_strand();

        /*
        // we were originally using this as a filter, but it's no longer necessary; preserving in case we need this in the future
        // the HLA lengths
        static ref HLA_CDNA_LENGTHS: HashMap<String, usize> = {
            let mut hla_lengths: HashMap<String, usize> = Default::default();
            hla_lengths.insert("HLA-A".to_string(), 1098);
            hla_lengths.insert("HLA-B".to_string(), 1089);
            hla_lengths
        };
        */

        // TODO: figure out how to automate this part in the future?
        let mut hla_exons: BTreeMap<String, Vec<Coordinates>> = Default::default();
        hla_exons.insert("HLA-A".to_string(), 
            vec![
                // 0-based, exclusive, sorted; copied from RefSeq file
                Coordinates::new("chr6".to_string(), 29942532 - preshift, 29942626 - postshift),
                Coordinates::new("chr6".to_string(), 29942757 - preshift, 29943026 - postshift),
                Coordinates::new("chr6".to_string(), 29943268 - preshift, 29943543 - postshift),
                Coordinates::new("chr6".to_string(), 29944122 - preshift, 29944397 - postshift),
                Coordinates::new("chr6".to_string(), 29944500 - preshift, 29944616 - postshift),
                Coordinates::new("chr6".to_string(), 29945059 - preshift, 29945091 - postshift),
                Coordinates::new("chr6".to_string(), 29945234 - preshift, 29945281 - postshift),
                Coordinates::new("chr6".to_string(), 29945451 - preshift, 29945870 - postshift)
            ]
        );
        hla_exons.insert("HLA-B".to_string(), 
            vec![
                // 0-based, exclusive, sorted; copied from RefSeq file
                Coordinates::new("chr6".to_string(), 31353875 - preshift, 31354296 - postshift),
                Coordinates::new("chr6".to_string(), 31354479 - preshift, 31354526 - postshift),
                Coordinates::new("chr6".to_string(), 31354633 - preshift, 31354665 - postshift),
                Coordinates::new("chr6".to_string(), 31355107 - preshift, 31355223 - postshift),
                Coordinates::new("chr6".to_string(), 31355317 - preshift, 31355592 - postshift),
                Coordinates::new("chr6".to_string(), 31356167 - preshift, 31356442 - postshift),
                Coordinates::new("chr6".to_string(), 31356688 - preshift, 31356957 - postshift),
                Coordinates::new("chr6".to_string(), 31357086 - preshift, 31357179 - postshift)
            ]
        );

        Self { 
            hla_coordinates, 
            hla_is_forward_strand,
            hla_exons
        }
    }
}

/// Wrapper for all the HLA allele informations
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct HlaAlleleDefinition {
    /// The identifier from IMGT
    hla_id: String,
    /// The gene name this goes to; e.g. HLA-A
    gene_name: String,
    /// The assigned star allele as a Vec; e.g. ["01", "01", "01"] is a three field allele
    star_allele: Vec<String>,
    /// The DNA sequence for the record
    dna_sequence: Option<String>,
    /// The cDNA sequence for the record
    cdna_sequence: String
}

impl HlaAlleleDefinition {
    /// Creates a new HlaAlleleDefinition and performs some checks along the way
    /// # Arguments
    /// * `hla_id` - the identifier, expected to be of form "HLA:HLA00001"
    /// * `description` - basically, the star allele, should be of form "A*01:01:01:01"
    /// * `dna_sequence` - the optional DNA sequence, should be ACGT symbols only
    /// * `cdna_sequence` - the cDNA sequence, should be ACGT symbols only
    pub fn new(hla_id: String, description: &str, dna_sequence: Option<String>, cdna_sequence: String) -> Result<HlaAlleleDefinition, Box<dyn std::error::Error>> {
        let star_split: Vec<&str> = description.split('*').collect();
        if star_split.len() != 2 {
            bail!("Star split length != 2 for allele description: {description}");
        }
        let gene_name: String = format!("HLA-{}", star_split[0]);

        let star_allele: Vec<String> = star_split[1].split(':').map(String::from).collect();
        if star_allele.len() > 4 {
            bail!("Unexpected number of fields for allele description: {description}");
        }

        let allowed_symbols = ['A', 'C', 'G', 'T'];
        if let Some(dna) = dna_sequence.as_ref() {
            if !dna.chars().all(|c| allowed_symbols.contains(&c)) {
                bail!("DNA sequence contains non-ACGT symbols.");
            }
        }
        if !cdna_sequence.chars().all(|c| allowed_symbols.contains(&c)) {
            bail!("cDNA sequence contains non-ACGT symbols.");
        }

        Ok(HlaAlleleDefinition {
            hla_id,
            gene_name,
            star_allele,
            dna_sequence,
            cdna_sequence
        })
    }

    pub fn hla_id(&self) -> &str {
        &self.hla_id
    }

    pub fn gene_name(&self) -> &str {
        &self.gene_name
    }

    pub fn star_allele(&self) -> &[String] {
        &self.star_allele
    }

    pub fn dna_sequence(&self) -> Option<&str> {
        self.dna_sequence.as_deref()
    }

    pub fn cdna_sequence(&self) -> &str {
        &self.cdna_sequence
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::path::PathBuf;

    use crate::util::file_io::load_json;

    #[test]
    fn test_config_full_length() {
        // full file
        let test_fn = PathBuf::from("test_data/HLA_configs/full_length.json");
        let config: HlaConfig = load_json(&test_fn).unwrap();
        assert!(config.validate_config().is_ok());
    }

    #[test]
    fn test_config_missing_regions() {
        // this one is missing a CYP2D6 coordinate
        let test_fn = PathBuf::from("test_data/HLA_configs/missing_regions.json");
        let config: HlaConfig = load_json(&test_fn).unwrap();
        assert!(config.validate_config().is_err());
    }

    #[test]
    fn test_config_missing_exons() {
        // this one is missing a CYP2D6 exon
        let test_fn = PathBuf::from("test_data/HLA_configs/missing_exons.json");
        let config: HlaConfig = load_json(&test_fn).unwrap();
        assert!(config.validate_config().is_err());
    }
    
    #[test]
    fn test_good_allele_def() {
        let test_name = "test_name".to_string();
        let test_gene = "A";
        let test_star = "01:01:01:01";
        let test_description = format!("{test_gene}*{test_star}");
        let test_dna = Some("ACGT".to_string());
        let test_cdna = "CG".to_string();
        let test_result = HlaAlleleDefinition::new(
            test_name.clone(),
            &test_description,
            test_dna.clone(),
            test_cdna.clone()
        ).unwrap();
        assert_eq!(test_result, HlaAlleleDefinition {
            hla_id: test_name,
            gene_name: "HLA-A".to_string(),
            star_allele: vec!["01".to_string(); 4],
            dna_sequence: test_dna,
            cdna_sequence: test_cdna
        });
    }

    #[test]
    fn test_bad_fields() {
        // not getting tested
        let test_name = "test_name".to_string();
        let test_dna = Some("ACGT".to_string());
        let test_cdna = "CG".to_string();

        // too many fields
        let test_gene = "A";
        let test_star = "01:01:01:01:01";
        let test_description = format!("{test_gene}*{test_star}");
        let test_result = HlaAlleleDefinition::new(
            test_name,
            &test_description,
            test_dna,
            test_cdna
        );
        assert!(test_result.is_err());
    }

    #[test]
    fn test_bad_alleles() {
        let test_name = "test_name".to_string();
        let test_gene = "A";
        let test_star = "01:01:01:01";
        let test_description = format!("{test_gene}*{test_star}");
        
        // bad dna
        let test_bad = "BOB".to_string();
        let test_good = "CG".to_string();
        let test_result = HlaAlleleDefinition::new(
            test_name.clone(),
            &test_description,
            Some(test_bad.clone()),
            test_good.clone()
        );
        assert!(test_result.is_err());

        // bad cnda
        let test_result = HlaAlleleDefinition::new(
            test_name.clone(),
            &test_description,
            Some(test_good),
            test_bad
        );
        assert!(test_result.is_err());
    }
}