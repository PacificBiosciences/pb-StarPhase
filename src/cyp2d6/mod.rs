
/// The entry function for diplotyping CYP2D6
pub mod caller;
/// Functionality that connects multiple regions together into a pair of full chains (i.e., the diplotype)
pub mod chaining;
/// Constants that control how the alleles are labeled and defined
pub mod definitions;
/// Errors and warnings that can come from the CYP2D6 algorithm that we want to cleanly handle
pub mod errors;
/// Functionality for extracting and haplotying CYP2D6 alleles in a sequence
pub mod haplotyper;
/// Wrapper for region labeling and constraining based on the labels
pub mod region_label;
/// Contains functionality for generating visualization of CYP2D6-related components
pub mod visualization;
