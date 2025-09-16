
use std::path::Path;

/// Gets a list of sample names from a given VCF file
/// # Arguments
/// * `filename` - the VCF file to load
/// # Errors
/// * if the filename fails to load as a VCF
/// * if the sample name fails to parse from utf8
pub fn get_vcf_samples(filename: &Path) -> Result<Vec<String>, Box<dyn std::error::Error>> {
    use rust_htslib::bcf;
    use rust_htslib::bcf::Read;
    let vcf_reader: bcf::IndexedReader = bcf::IndexedReader::from_path(filename)?;
    let vcf_header: bcf::header::HeaderView = vcf_reader.header().clone();
    let mut sample_names = vec![];
    for sv in vcf_header.samples().iter() {
        let vcf_sample_string: String = std::str::from_utf8(sv)?.to_string();
        sample_names.push(vcf_sample_string);
    }
    Ok(sample_names)
}
