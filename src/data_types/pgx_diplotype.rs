
use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};

use crate::data_types::region_variants::{RegionVariant, VariantAlleleRelationship};

/// Contains all the information related to a single gene's diplotype result
#[derive(Clone, Debug, Deserialize, Eq, Serialize)]
pub struct Diplotype {
    /// short string for haplotype 1
    hap1: String,
    /// short string for haplotype 2
    hap2: String,
    /// combination diplotype call
    diplotype: String
}

impl Diplotype {
    pub fn new(hap1: &str, hap2: &str) -> Diplotype {
        Diplotype {
            hap1: hap1.to_string(),
            hap2: hap2.to_string(),
            diplotype: format!("{hap1}/{hap2}")
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

    // getters
    pub fn hap1(&self) -> &str {
        &self.hap1
    }

    pub fn hap2(&self) -> &str {
        &self.hap2
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

impl PartialOrd for Diplotype {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Diplotype {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // sort the haplotypes to make sure we are comparing the same ones
        let mut ordered_haps = vec![&self.hap1, &self.hap2];
        let mut ordered_other_haps = vec![&other.hap1, &other.hap2];
        ordered_haps.sort();
        ordered_other_haps.sort();
        ordered_haps.cmp(&ordered_other_haps)
    }
}

/// Intended to capture the haplotype details for an inexact diplotype match
#[derive(Clone, Debug, Deserialize, PartialEq,Serialize)]
pub struct InexactDiplotype {
    /// Wrapper for the string formats of the diplotype
    basic_diplotype: Diplotype,
    /// haplotype 1 - these can be None for high complexity genes like CYP2D6 where we cannot populate this simply
    haplotype_1: Option<InexactHaplotype>,
    /// haplotype 2
    haplotype_2: Option<InexactHaplotype>,
}

impl InexactDiplotype {
    /// Constructor takes the two inexact haplotypes and derives the diplotype from them
    /// # Arguments
    /// * `haplotype_1` - the first inexact haplotype
    /// * `haplotype_2` - the second inexact haplotype
    pub fn new(haplotype_1: InexactHaplotype, haplotype_2: InexactHaplotype) -> InexactDiplotype {
        let basic_diplotype = Diplotype::new(&haplotype_1.full_haplotype(), &haplotype_2.full_haplotype());
        InexactDiplotype {
            basic_diplotype,
            haplotype_1: Some(haplotype_1),
            haplotype_2: Some(haplotype_2)
        }
    }

    pub fn new_diplotype_only(diplotype: Diplotype) -> InexactDiplotype {
        InexactDiplotype {
            basic_diplotype: diplotype,
            haplotype_1: None,
            haplotype_2: None
        }
    }

    // getters
    pub fn basic_diplotype(&self) -> &Diplotype {
        &self.basic_diplotype
    }

    pub fn haplotype_1(&self) -> Option<&InexactHaplotype> {
        self.haplotype_1.as_ref()
    }

    pub fn haplotype_2(&self) -> Option<&InexactHaplotype> {
        self.haplotype_2.as_ref()
    }
}

/// Enum describing the type of match a haplotype is, which is basically a summary of the variant statuses
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
pub enum InexactMatchType {
    /// The match type is unknown, mostly a placeholder for future
    Unknown,
    /// Does not match a core allele (and therefor not a sub-allele)
    NoMatch,
    /// Matches a core allele, but NOT a sub-allele
    CoreMatch,
    /// Exact match to a sub-allele
    SubAlleleMatch
}

/// Wrapper for the details of an inexact haplotype match
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct InexactHaplotype {
    /// The base haplotype name
    base_haplotype: String,
    /// The type of match this haplotype is
    match_type: InexactMatchType,
    /// The list of variants and their statuses relative to the base haplotype
    variant_relationships: BTreeSet<RegionVariant>,
}

impl InexactHaplotype {
    /// Constructor takes the base haplotype name and the list of variants and their statuses relative to the base haplotype
    /// # Arguments
    /// * `base_haplotype` - the base haplotype name
    /// * `variant_relationships` - the list of variants and their statuses relative to the base haplotype
    pub fn new(base_haplotype: String, variant_relationships: BTreeSet<RegionVariant>) -> InexactHaplotype {
        // find any variants that are not a match
        let mut core_match = true;
        let mut suballele_match = true;
        for variant in variant_relationships.iter() {
            if variant.variant_state() != VariantAlleleRelationship::Match {
                suballele_match = false;
                if variant.is_vi() {
                    core_match = false;
                }
            }
        }

        // classify the match type after scanning the variants
        let match_type = if suballele_match {
            InexactMatchType::SubAlleleMatch
        } else if core_match {
            InexactMatchType::CoreMatch
        } else {
            InexactMatchType::NoMatch
        };

        InexactHaplotype {
            base_haplotype,
            match_type,
            variant_relationships
        }
    }

    /// Returns the full haplotype name, which includes the base haplotype name and the list of variants and their statuses relative to the base haplotype.
    pub fn full_haplotype(&self) -> String {
        let mut haplotype = self.base_haplotype.clone();
        let mut mod_made = false;

        // add any variants that are not a match
        for variant in self.variant_relationships.iter() {
            if variant.variant_state() != VariantAlleleRelationship::Match {
                haplotype.push_str(&format!(" {variant}"));
                mod_made = true;
            }
        }

        // if we made any modifications, then wrap the haplotype in parentheses
        if mod_made {
            format!("({haplotype})")
        } else {
            haplotype
        }
    }

    // getters
    pub fn base_haplotype(&self) -> &str {
        &self.base_haplotype
    }

    pub fn match_type(&self) -> InexactMatchType {
        self.match_type
    }

    pub fn variant_relationships(&self) -> &BTreeSet<RegionVariant> {
        &self.variant_relationships
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn test_inexact_diplotype() {
        // Test new_diplotype_only constructor
        let diplotype = Diplotype::new("*1", "*2");
        let inexact = InexactDiplotype::new_diplotype_only(diplotype.clone());
        
        assert_eq!(inexact.basic_diplotype(), &diplotype);
        assert_eq!(inexact.haplotype_1(), None);
        assert_eq!(inexact.haplotype_2(), None);

        // Test new constructor with haplotypes
        let hap1 = InexactHaplotype::new("*1".to_string(), Default::default());
        let hap2 = InexactHaplotype::new("*2".to_string(), Default::default());
        let inexact_with_haps = InexactDiplotype::new(hap1.clone(), hap2.clone());
        
        assert_eq!(inexact_with_haps.basic_diplotype().diplotype(), "*1/*2");
        assert_eq!(inexact_with_haps.haplotype_1(), Some(&hap1));
        assert_eq!(inexact_with_haps.haplotype_2(), Some(&hap2));
    }

    #[test]
    fn test_inexact_haplotype() {
        // Test SubAlleleMatch - all variants match
        let variants_match = BTreeSet::from([
            RegionVariant::new("rs123".to_string(), true, VariantAlleleRelationship::Match),
            RegionVariant::new("rs456".to_string(), false, VariantAlleleRelationship::Match),
        ]);
        let hap_suballele = InexactHaplotype::new("*1.001".to_string(), variants_match.clone());
        assert_eq!(hap_suballele.base_haplotype(), "*1.001");
        assert_eq!(hap_suballele.match_type(), InexactMatchType::SubAlleleMatch);
        assert_eq!(hap_suballele.variant_relationships(), &variants_match);
        assert_eq!(hap_suballele.full_haplotype(), "*1.001"); // No modifications, no parentheses

        // Test CoreMatch - non-VI variants don't match, but VI variants do
        let variants_core = BTreeSet::from([
            RegionVariant::new("rs123".to_string(), true, VariantAlleleRelationship::Match), // VI match
            RegionVariant::new("rs456".to_string(), false, VariantAlleleRelationship::Unexpected), // non-VI doesn't match
        ]);
        let hap_core = InexactHaplotype::new("*1.001".to_string(), variants_core.clone());
        assert_eq!(hap_core.base_haplotype(), "*1.001");
        assert_eq!(hap_core.match_type(), InexactMatchType::CoreMatch);
        assert_eq!(hap_core.variant_relationships(), &variants_core);
        assert_eq!(hap_core.full_haplotype(), "(*1.001 +rs456)"); // Modified, should be in parentheses

        // Test NoMatch - VI variant doesn't match
        let variants_no_match = BTreeSet::from([
            RegionVariant::new("rs123".to_string(), true, VariantAlleleRelationship::Missing), // VI doesn't match
            RegionVariant::new("rs456".to_string(), false, VariantAlleleRelationship::Unexpected),
        ]);
        let hap_no_match = InexactHaplotype::new("*1".to_string(), variants_no_match.clone());
        assert_eq!(hap_no_match.base_haplotype(), "*1");
        assert_eq!(hap_no_match.match_type(), InexactMatchType::NoMatch);
        assert_eq!(hap_no_match.variant_relationships(), &variants_no_match);
        assert_eq!(hap_no_match.full_haplotype(), "(*1 -rs123 +rs456)"); // Both variants modified

        // Test empty variants - should be SubAlleleMatch
        let hap_empty = InexactHaplotype::new("*3".to_string(), Default::default());
        assert_eq!(hap_empty.match_type(), InexactMatchType::SubAlleleMatch);
        assert_eq!(hap_empty.full_haplotype(), "*3"); // No modifications
    }
}
