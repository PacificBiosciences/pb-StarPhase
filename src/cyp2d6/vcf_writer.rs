
use log::debug;
use rust_htslib::bcf;
use rust_htslib::bcf::header::Header;
use rust_htslib::bcf::Writer;
use std::collections::BTreeMap;
use std::path::Path;

use crate::cyp2d6::haplotyper::LoadedVariants;
use crate::cyp2d6::region::Cyp2d6Region;
use crate::cyp2d6::region_label::Cyp2d6RegionType;
use crate::data_types::region_variants::VariantAlleleRelationship;

const GT_UNKNOWN: i32 = 0;
const GT_REF: i32 = 2; // encode "0"
const GT_ALT: i32 = 4; // encode "1"
const GT_END_SENTINEL: i32 = i32::MIN+1;

/// Writes a VCF file containing the CYP2D6 diplotypes
/// # Arguments
/// * `regions` - the regions to write, some of which we will skip because they are not CYP2D6 alleles
/// * `vcf_fn` - the path to the VCF file to write
/// * `loaded_variants` - the loaded variants to use for the VCF file
/// # Errors
/// * if the VCF file cannot be written
pub fn write_cyp2d6_vcf(regions: &[Cyp2d6Region], vcf_fn: &Path, loaded_variants: &LoadedVariants) -> Result<(), Box<dyn std::error::Error>> {
    debug!("Writing CYP2D6 alleles to {vcf_fn:?}...");

    // header information
    let ver: &str = crate::cli::core::FULL_VERSION.as_str(); // clippy gets weird about direct access
    let cli_version = format!("\"{ver}\"");
    let cli_string = format!("\"{}\"", std::env::args().collect::<Vec<String>>().join(" "));

    // Create VCF header
    let mut header = Header::new();

    // Add file format version
    header.push_record(b"##fileformat=VCFv4.2");
    
    // Add file date
    let date = chrono::Utc::now().format("%Y%m%d").to_string();
    header.push_record(format!("##fileDate={date}").as_bytes());
    
    // Add source program
    header.push_record(b"##source=pb-StarPhase");
    
    // Add reference genome
    header.push_record(b"##reference=GRCh38"); // TODO: do we want to pass in the file? I don't think it matters here
    
    // Add contig information for chromosome 22 (where CYP2D6 is located)
    header.push_record(b"##contig=<ID=chr22,length=50818468>");
    
    // Add INFO field definitions - nothing to do here yet, leaving a boiler plate for future use
    header.push_record(b"##INFO=<ID=VI,Number=1,Type=String,Description=\"Variant impact\">");
    
    // Add FORMAT field definitions
    header.push_record(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

    // Add program version and CLI string
    header.push_record(format!("##pbstarphase_version={cli_version}").as_bytes());
    header.push_record(format!("##pbstarphase_command={cli_string}").as_bytes());

    // Add sample column header
    let mut num_regions = 0;
    for region in regions.iter() {
        if region.label().region_type() == Cyp2d6RegionType::Cyp2d6 {
            let region_name = region.index_label();
            header.push_sample(region_name.as_bytes());
            num_regions += 1;
        }
    }
    
    // Create VCF writer
    let mut writer = Writer::from_path(vcf_fn, &header, false, bcf::Format::Vcf)?;

    // figure out which variants are getting written for each region
    // this is the index of the variant in the loaded variants, and then an optional assignment
    let mut variants_to_write: BTreeMap<usize, Vec<VariantAlleleRelationship>> = Default::default();
    let mut region_index = 0;
    for region in regions.iter() {
        if region.label().region_type() == Cyp2d6RegionType::Cyp2d6 {
            for variant in region.variants().unwrap().iter() {
                let label = variant.label();
                let variant_state = variant.variant_state();

                let idx = loaded_variants.index_label(label)?;
                let variant_entry = variants_to_write.entry(idx)
                    .or_insert(vec![VariantAlleleRelationship::Unknown; num_regions]);

                variant_entry[region_index] = variant_state;
            }

            // increment the region index
            region_index += 1;
        }
    }

    for (&variant_idx, variant_states) in variants_to_write.iter() {
        let variant = &loaded_variants.ordered_variants()[variant_idx];
        let variant_meta = &loaded_variants.variant_metadata()[variant_idx];

        // Create a variant record
        let mut record = writer.empty_record();

        // set the variant ID, chromosome, position, and alleles
        record.set_id(variant_meta.label.as_bytes())?;
        record.set_rid(Some(0)); // 0 is the index of the contig in the header, which is currently always 0 for us
        record.set_pos(variant.position()); // 0-based position of the variant
        record.set_alleles(&[
            variant.get_allele0(),
            variant.get_allele1()
        ])?;

        // Set quality - TODO: if we ever enable this, we will need to calculate it; I imagine we would have DP and other metrics at that point
        record.set_qual(bcf::record::Numeric::missing());

        // Set the "VI" tag which is just pass-through from the database
        if let Some(vi_label) = variant_meta.vi_label.as_ref() {
            record.push_info_string(b"VI", &[vi_label.as_bytes()])?;
        }

        // Set genotype for each sample
        let mut alleles = vec![];
        for variant_state in variant_states.iter() {
            let gt = match variant_state {
                // ambiguity, so report "."
                VariantAlleleRelationship::AmbiguousUnexpected |
                VariantAlleleRelationship::AmbiguousMissing => GT_UNKNOWN,
                // unknown or missing, so report "0" or REF allele
                VariantAlleleRelationship::UnknownUnexpected |
                VariantAlleleRelationship::UnknownMissing |
                VariantAlleleRelationship::Unknown |
                VariantAlleleRelationship::Missing => GT_REF,
                // unexpected or match, so report "1" or ALT allele
                VariantAlleleRelationship::Unexpected |
                VariantAlleleRelationship::Match => GT_ALT,
            };

            alleles.push(gt);
            alleles.push(GT_END_SENTINEL); // sentinel value for end
        }
        record.push_format_integer(b"GT", &alleles)?;

        // Write the record
        writer.write(&record)?;
    }

    Ok(())
}
