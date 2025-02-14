
use serde::Serialize;

/// Enum describing the expected states that a variant can be in relative to some described haplotype
#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub enum VariantAlleleRelationship {
    /// Generic for anything undescribed
    Unknown,
    /// Indicates the variant was identified and expected; this is the 99% path
    Match,
    /// Indicates the variant was identified, but was not expected for the allele
    Unexpected,
    /// Indicated the variant was not identified, but was expected for the allele
    Missing,
    /// Indicates the variant was ambiguous in ref/alt, but was not expected for the allele
    AmbiguousUnexpected,
    /// Indicates the variant was ambiguous in ref/alt, but was expected for the allele
    AmbiguousMissing,
    /// Indicates the variant was unknown, but was not expected for the allele
    UnknownUnexpected,
    /// Indicates the variant was unknown, but was expected for the allele
    UnknownMissing
}

/// Wrapper for all the variants that were identified for the region
#[derive(Clone, Debug, Serialize)]
pub struct Cyp2d6RegionVariant {
    /// User-friendly label, this is expected to be cross-referenced to something else; e.g. rsID or {chr}:{pos} {ref}>{alt}
    label: String,
    /// True if flagged as VI
    is_vi: bool,
    /// The state of the variant relative to some allele
    variant_state: VariantAlleleRelationship
}

impl Cyp2d6RegionVariant {
    /// Constructor
    pub fn new(label: String, is_vi: bool, variant_state: VariantAlleleRelationship) -> Self {
        Self {
            label, is_vi, variant_state
        }
    }

    // getters
    pub fn label(&self) -> &str {
        &self.label
    }

    pub fn is_vi(&self) -> bool {
        self.is_vi
    }

    pub fn variant_state(&self) -> VariantAlleleRelationship {
        self.variant_state
    }
}