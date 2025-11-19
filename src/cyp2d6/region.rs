
use std::fmt::Display;

use crate::cyp2d6::region_label::Cyp2d6RegionLabel;
use crate::data_types::region_variants::{RegionVariant, VariantAlleleRelationship};

/// Detail level for reporting an allelic region
#[derive(Clone, Copy,PartialEq)]
pub enum Cyp2d6DetailLevel {
    /// Only reports a core allele, this is what is relevant for PGx clinically (currently)
    CoreAlleles,
    /// Reports sub-alleles, which are non-coding variation
    SubAlleles,
    /// Same as sub-alleles, except it also include any detected differences encoded as +/- the variant; really only for power users and debugging
    DeepAlleles
}

pub struct Cyp2d6Region {
    /// Optional unique ID for a region, this is useful when two region labels are otherwise identical, but the sequences are different
    unique_id: Option<usize>,
    /// The high-level identifier for the sequence
    label: Cyp2d6RegionLabel,
    /// The specific variants that were identified, usually only for CYP2D6 derivatives
    variants: Option<Vec<RegionVariant>>
}

impl Cyp2d6Region {
    /// Constructor
    pub fn new(label: Cyp2d6RegionLabel, variants: Option<Vec<RegionVariant>>) -> Self {
        // TODO: do we want to enforce checks here?
        Self {
            unique_id: None,
            label,
            variants
        }
    }

    /// Enable a unique ID
    pub fn set_unique_id(&mut self, unique_id: usize) {
        self.unique_id = Some(unique_id);
    }

    /// Marks the label for this as a false allele, preserving other data.
    pub fn mark_false_allele(&mut self) {
        self.label.mark_false_allele();
    }

    /// Returns a label for an allele of the form {index}_{full_allele}.
    /// Importantly, the full_allele is not translated to human readable, but should be unique if the `unique_id` has been set.
    pub fn index_label(&self) -> String {
        if let Some(ud) = self.unique_id {
            format!("{ud}_{}", self.label.full_allele())
        } else {
            format!("X_{}", self.label.full_allele())
        }
    }

    /// This function will provide a deep haplotype label that annotates differences relative to the original label.
    /// It is based on a combination of the assigned `label` and also the state of the attached variants.
    pub fn deep_label(&self) -> String {
        let index_label = self.index_label();
        let mut join_set = vec![index_label];
        if let Some(v_list) = self.variants.as_ref() {
            // we have a variant list to scan
            for variant in v_list.iter() {
                match variant.variant_state() {
                    VariantAlleleRelationship::Match |
                    VariantAlleleRelationship::UnknownUnexpected => {
                        // Match is the ideal state, in which case no-op required OR
                        // it's unknown, but also was not expected
                    },
                    VariantAlleleRelationship::Unexpected => {
                        // this is an "extra" variant
                        join_set.push(format!("+{}", variant.label()));
                    },
                    VariantAlleleRelationship::Missing => {
                        // this is a variant that is missing
                        join_set.push(format!("-{}", variant.label()));
                    },
                    // for either ambiguous result, the implication is that something unusual is happening
                    VariantAlleleRelationship::AmbiguousUnexpected | // ambiguous result, but we did not expect the variant
                    VariantAlleleRelationship::AmbiguousMissing | // ambiguous result, but we expected the variant
                    VariantAlleleRelationship::Unknown | // this one should not happen
                    VariantAlleleRelationship::UnknownMissing => {
                        // It's not really clear what to do with any of these, maybe flag with ? for now
                        // we may just want to ignore them though
                        join_set.push(format!("?{}", variant.label()));
                    }
                };
            }
        }
        join_set.join(" ")
    }

    // getters
    pub fn label(&self) -> &Cyp2d6RegionLabel {
        &self.label
    }

    pub fn variants(&self) -> Option<&[RegionVariant]> {
        self.variants.as_deref()
    }
}

impl Display for Cyp2d6Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.label)?;
        Ok(())
    }
}
