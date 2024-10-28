
use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use simple_error::bail;

/// Wrapper for all the HLA allele informations
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct AlleleDefinition {
    /// The identifier from the upstream DB
    id: String,
    /// The gene name this is associated with
    gene_name: String,
    /// The assigned star allele
    star_allele: String,
    /// The variants defining this allele
    variants: Vec<VariantDefinition>
}

impl AlleleDefinition {
    /// Creates a new AlleleDefinition and performs some checks along the way
    /// # Arguments
    /// * `id` - the identifier, expected to be unique
    /// * `description` - basically, the star allele, should be of form "{gene}*{star_allele}", e.g.: "CYP2D6*1"
    /// * `dna_sequence` - the DNA sequence, should be ACGT symbols only
    pub fn new(opt_id: Option<String>, description: &str, variants: Vec<VariantDefinition>) -> Result<AlleleDefinition, Box<dyn std::error::Error>> {
        let star_split: Vec<&str> = description.split('*').collect();
        if star_split.len() != 2 {
            bail!("Star split length != 2 for allele description: {description}");
        }

        let gene_name: String = star_split[0].to_string();
        let star_allele: String = star_split[1].to_string();

        // if we do not have a special ID, just use the full star allele format
        let id = opt_id.unwrap_or(format!("{gene_name}*{star_allele}"));
        
        Ok(AlleleDefinition {
            id,
            gene_name,
            star_allele,
            variants,
        })
    }

    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn gene_name(&self) -> &str {
        &self.gene_name
    }

    pub fn star_allele(&self) -> &str {
        &self.star_allele
    }

    pub fn variants(&self) -> &[VariantDefinition] {
        &self.variants
    }
}

/// This corresponds to an individual variant or sequence replacement that is part of an allele definition
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct VariantDefinition {
    /// an external identifier for this variant
    id: Option<String>,
    /// the chromosome
    chrom: String,
    /// the 0-based position 
    position: usize,
    /// the reference sequence that is getting replaced
    reference: String,
    /// the alternate sequence replacing the reference
    alternate: String,
    /// additional data can go here, such as the VI field
    extras: BTreeMap<String, String>
}

impl VariantDefinition {
    pub fn new(id: Option<String>, chrom: String, position: usize, reference: String, alternate: String, extras: BTreeMap<String, String>) 
        -> Result<VariantDefinition, Box<dyn std::error::Error>> {
        // check reference and alternate for ACGT only
        let allowed_symbols = ['A', 'C', 'G', 'T'];
        if !reference.chars().all(|c| allowed_symbols.contains(&c)) {
            bail!("Reference sequence contains non-ACGT symbols: {reference}");
        }
        if !alternate.chars().all(|c| allowed_symbols.contains(&c)) {
            bail!("Reference sequence contains non-ACGT symbols: {alternate}");
        }

        Ok(VariantDefinition {
            id,
            chrom, 
            position,
            reference,
            alternate,
            extras
        })
    }

    pub fn position(&self) -> usize {
        self.position
    }

    pub fn reference(&self) -> &str {
        &self.reference
    }

    pub fn alternate(&self) -> &str {
        &self.alternate
    }

    pub fn extras(&self) -> &BTreeMap<String, String> {
        &self.extras
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_good_allele_def() {
        let test_name = "test_name".to_string();
        let test_gene = "CYP2D6";
        let test_star = "1";
        let test_description = format!("{test_gene}*{test_star}");
        let test_result = AlleleDefinition::new(
            Some(test_name.clone()),
            &test_description,
            vec![]
        ).unwrap();
        assert_eq!(test_result, AlleleDefinition {
            id: test_name,
            gene_name: test_gene.to_string(),
            star_allele: test_star.to_string(),
            variants: vec![]
        });
    }

    #[test]
    fn test_bad_allele_name() {
        let test_result = AlleleDefinition::new(
            None, "A bad name", vec![]
        );
        assert!(test_result.is_err());
    }

    #[test]
    fn test_variant_definition() {
        let id = Some("random_id".to_string());
        let chrom = "chr22".to_string();
        let position = 10;
        let reference = "A".to_string();
        let alternate = "C".to_string();
        let extras: BTreeMap<String, String> = Default::default();
        let vd = VariantDefinition::new(
            id.clone(), chrom.clone(), position, reference.clone(), alternate.clone(), extras.clone()
        ).unwrap();
        assert_eq!(vd, VariantDefinition {
            id,
            chrom,
            position,
            reference,
            alternate,
            extras
        });
    }

    #[test]
    fn test_bad_variant_definition() {
        let id = Some("random_id".to_string());
        let chrom = "chr22".to_string();
        let position = 10;
        let reference = "B".to_string();
        let alternate = "C".to_string();
        let extras: BTreeMap<String, String> = Default::default();
        let vd = VariantDefinition::new(
            id.clone(), chrom.clone(), position, reference.clone(), alternate.clone(), extras.clone()
        );
        assert!(vd.is_err());
    }
}