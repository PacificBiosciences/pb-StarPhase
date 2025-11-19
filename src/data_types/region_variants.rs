
use serde::{Deserialize, Serialize};

/// Enum describing the expected states that a variant can be in relative to some described haplotype
#[derive(Clone, Copy, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
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
#[derive(Clone, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
pub struct RegionVariant {
    /// User-friendly label, this is expected to be cross-referenced to something else; e.g. rsID or {chr}:{pos} {ref}>{alt}
    label: String,
    /// True if flagged as VI, this is also an alias for "core variants"
    is_vi: bool,
    /// The state of the variant relative to some allele
    variant_state: VariantAlleleRelationship
}

impl RegionVariant {
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

impl std::fmt::Display for RegionVariant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let precursor = match self.variant_state {
            VariantAlleleRelationship::Match => "=",
            VariantAlleleRelationship::Unexpected => "+",
            VariantAlleleRelationship::Missing => "-",
            // these are all ambiguous or unknown, so we report a ?
            VariantAlleleRelationship::AmbiguousUnexpected |
            VariantAlleleRelationship::AmbiguousMissing |
            VariantAlleleRelationship::UnknownUnexpected |
            VariantAlleleRelationship::UnknownMissing |
            VariantAlleleRelationship::Unknown => "?"
        };
        write!(f, "{}{}", precursor, self.label)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_region_variant() {
        // Test constructor and getters
        let variant = RegionVariant::new("rs123456".to_string(), true, VariantAlleleRelationship::Match);
        assert_eq!(variant.label(), "rs123456");
        assert_eq!(variant.is_vi(), true);
        assert_eq!(variant.variant_state(), VariantAlleleRelationship::Match);

        let variant_non_vi = RegionVariant::new("chr1:1000A>T".to_string(), false, VariantAlleleRelationship::Unexpected);
        assert_eq!(variant_non_vi.label(), "chr1:1000A>T");
        assert_eq!(variant_non_vi.is_vi(), false);
        assert_eq!(variant_non_vi.variant_state(), VariantAlleleRelationship::Unexpected);

        // Test Display implementation for Match
        let match_var = RegionVariant::new("rs123".to_string(), true, VariantAlleleRelationship::Match);
        assert_eq!(match_var.to_string(), "=rs123");

        // Test Display implementation for Unexpected
        let unexpected_var = RegionVariant::new("rs456".to_string(), false, VariantAlleleRelationship::Unexpected);
        assert_eq!(unexpected_var.to_string(), "+rs456");

        // Test Display implementation for Missing
        let missing_var = RegionVariant::new("rs789".to_string(), true, VariantAlleleRelationship::Missing);
        assert_eq!(missing_var.to_string(), "-rs789");

        // Test Display implementation for AmbiguousUnexpected
        let amb_unexpected = RegionVariant::new("rs012".to_string(), false, VariantAlleleRelationship::AmbiguousUnexpected);
        assert_eq!(amb_unexpected.to_string(), "?rs012");

        // Test Display implementation for AmbiguousMissing
        let amb_missing = RegionVariant::new("rs345".to_string(), true, VariantAlleleRelationship::AmbiguousMissing);
        assert_eq!(amb_missing.to_string(), "?rs345");

        // Test Display implementation for UnknownUnexpected
        let unk_unexpected = RegionVariant::new("rs678".to_string(), false, VariantAlleleRelationship::UnknownUnexpected);
        assert_eq!(unk_unexpected.to_string(), "?rs678");

        // Test Display implementation for UnknownMissing
        let unk_missing = RegionVariant::new("rs901".to_string(), true, VariantAlleleRelationship::UnknownMissing);
        assert_eq!(unk_missing.to_string(), "?rs901");

        // Test Display implementation for Unknown
        let unknown = RegionVariant::new("rs234".to_string(), false, VariantAlleleRelationship::Unknown);
        assert_eq!(unknown.to_string(), "?rs234");
    }
}