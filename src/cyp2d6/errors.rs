
/// Errors that can be produced by the CYP2D6 calling algorithm, these result in a failure state
#[derive(thiserror::Error, Debug, PartialEq)]
pub enum CallerError {
    #[error("no alleles were found that can start a chain")]
    NoChainingHead,
    #[error("no valid chains were identified")]
    NoChainsFound,
    #[error("no successful chain scoring pairs")]
    NoScorePairs
}

/// Warnings that can be produced by the CYP2D6 calling algorithm
#[derive(Debug, PartialEq)]
pub enum CallerWarning {
    /// Indicates an allele that was detected but not included as part of the final reported diplotype
    DanglingAllele { allele_name: String }
}
