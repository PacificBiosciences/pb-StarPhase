

use serde::Serialize;
use std::collections::BTreeMap;

use crate::cyp2d6::caller::convert_chain_to_hap;
use crate::cyp2d6::region::{Cyp2d6DetailLevel, Cyp2d6Region};
use crate::cyp2d6::region_variants::Cyp2d6RegionVariant;

/// Primary object that gets converted into JSON
#[derive(Serialize)]
pub struct DeeplotypeDebug {
    /// hap1 various complete forms
    hap1: HeeplotypeDebug,
    /// hap2 various complete forms
    hap2: HeeplotypeDebug,
    /// allelic details
    alleles: BTreeMap<String, Vec<Cyp2d6RegionVariant>>
}

impl DeeplotypeDebug {
    /// Constructor from a provided set of haplotypes
    pub fn new(best_diplotype_indices: &[Vec<usize>], hap_regions: &[Cyp2d6Region], cyp_translate: &BTreeMap<String, String>) -> Self {
        assert!(best_diplotype_indices.len() == 2);

        let hap1 = HeeplotypeDebug::new(&best_diplotype_indices[0], hap_regions, cyp_translate);
        let hap2 = HeeplotypeDebug::new(&best_diplotype_indices[1], hap_regions, cyp_translate);

        let mut alleles: BTreeMap<String, Vec<Cyp2d6RegionVariant>> = Default::default();
        for region in hap_regions.iter() {
            if let Some(variants) = region.variants() {
                // this region has defined variants, so lets add it to the debug outputs
                let core_label = region.index_label();
                let vec_variants = variants.to_vec();
                assert!(alleles.insert(core_label, vec_variants).is_none());
            }
        }

        Self {
            hap1,
            hap2,
            alleles
        }
    }
}

/// Wrapper around a haplotype containing all the representations of the full sequence
#[derive(Serialize)]
struct HeeplotypeDebug {
    /// suballeles and also any differences
    deep_form: String,
    /// suballeles included
    suballele_form: String,
    /// just core alleles
    core_form: String
}

impl HeeplotypeDebug {
    /// Constructor
    fn new(hap_chain: &[usize], hap_regions: &[Cyp2d6Region], cyp_translate: &BTreeMap<String, String>) -> Self {
        let deep_form = convert_chain_to_hap(hap_chain, hap_regions, Cyp2d6DetailLevel::DeepAlleles, cyp_translate);
        let suballele_form = convert_chain_to_hap(hap_chain, hap_regions, Cyp2d6DetailLevel::SubAlleles, cyp_translate);
        let core_form = convert_chain_to_hap(hap_chain, hap_regions, Cyp2d6DetailLevel::CoreAlleles, cyp_translate);
        Self {
            deep_form,
            suballele_form,
            core_form
        }
    }
}
