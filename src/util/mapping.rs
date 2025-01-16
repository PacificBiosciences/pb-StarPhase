
use minimap2::{Aligner, PresetSet};

use crate::data_types::mapping::MappingStats;

/// This is our standard aligner type that we use in almost all places.
/// After calling this function, mapopt and idxopt can still be changed until a sequence/file is provided to index.
pub fn standard_hifi_aligner() -> Aligner<PresetSet> {
    let mut aligner = Aligner::builder()
        .map_hifi()
        .with_cigar();
    aligner.mapopt.best_n = 5; // temporary workaround caused by error injected in v0.1.21
    aligner
}

/// Given a list of mappings, pick the best one out
/// # Arguments
/// * `mappings` - the list to select from
/// * `unmapped_from_target` - if True, use unmapped relative to target instead of query
/// * `penalize_unmapped` - if True, then penalize unmapped bases
/// * `base_length_override` - if Some, used as baseline length instead of query/target
pub fn select_best_mapping<'a>(
    mappings: &'a [minimap2::Mapping],
    unmapped_from_target: bool,
    penalize_unmapped: bool,
    base_length_override: Option<usize>
) -> (Option<&'a minimap2::Mapping>, MappingStats) {
    // default has 100% mismatch rate
    let mut best_stats: MappingStats = MappingStats::new(
        base_length_override.unwrap_or(1),
        base_length_override.unwrap_or(1),
        0
    );
    let mut best_mapping: Option<&'a minimap2::Mapping> = None;
    for m in mappings.iter() {
        // scoring is based on the lowest edit distance, including unmapped
        let nm = m.alignment.as_ref().unwrap().nm as usize;

        // some sanity checks while we debug
        let (base_length, unmapped) = if unmapped_from_target {
            let bl = base_length_override.unwrap_or(m.target_len as usize);
            let um = bl - (m.target_end - m.target_start) as usize;
            (bl, um)
        } else {
            let bl = base_length_override.unwrap_or(m.query_len.unwrap().get() as usize);
            let um = bl - (m.query_end - m.query_start) as usize;
            (bl, um)
        };
        let stats = MappingStats::new(base_length, nm, unmapped);
        if stats.custom_score(penalize_unmapped).score() < best_stats.custom_score(penalize_unmapped).score() {
            best_stats = stats;
            best_mapping = Some(m);
        }
    }

    (best_mapping, best_stats)
}