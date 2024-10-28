
use log::debug;
use simple_error::{bail, SimpleError};

use crate::data_types::mapping::MappingStats;
use crate::hla::mapping::HlaMappingStats;

/// Container for a processing HLA mapping. Mostly helpful for wrapping the `is_better_match` function.
pub struct HlaProcessedMatch {
    /// Name for this haplotype
    haplotype: String,
    /// Collection of full mappings
    full_mappings: Vec<Option<minimap2::Mapping>>,
    /// Collection of end-to-end stats, mostly for tie-breaking
    full_mapping_stats: Vec<Option<MappingStats>>,
    /// Collection of processed cigars for comparing
    processed_cigars: Vec<Option<Vec<usize>>>,
    /// Ranges that were actually included in the processing
    processed_ranges: Vec<std::ops::Range<usize>>,
}

impl HlaProcessedMatch {
    // constructors
    /// Creates a new processed match with the given haplotype name
    pub fn new(haplotype: String) -> Result<Self, SimpleError> {
        if haplotype.is_empty() {
            bail!("Haplotype name cannot be empty")
        }
        Ok(Self {
            haplotype,
            full_mappings: vec![],
            full_mapping_stats: vec![],
            processed_cigars: vec![],
            processed_ranges: vec![]
        })
    }

    /// Creates a terrible default, everything should beat this
    pub fn worst_match(num_sequences: usize) -> Self {
        Self {
            haplotype: Default::default(),
            full_mappings: vec![None; num_sequences],
            full_mapping_stats: vec![None; num_sequences],
            processed_cigars: vec![None; num_sequences],
            processed_ranges: vec![0..0; num_sequences]
        }
    }
    
    /// Adds a mapping to this processed match; can be None to indicate missing cDNA/DNA comparators
    /// # Arguments 
    /// * `mapping` - the minimap2 full mapping
    /// * `seq_len` - baseline sequence length (`target`) from the alignment
    pub fn add_mapping(&mut self, mapping: Option<minimap2::Mapping>, seq_len: usize) -> Result<(), SimpleError> {
        let (opt_processed_cigar, opt_mapping_stats, processed_range) = if let Some(m) = mapping.as_ref() {
            // first, get the processed cigar
            if m.strand != minimap2::Strand::Forward {
                bail!("Reverse strand mappings are not supported by HlaProcessedMatch");
            }
            let cigar = m.alignment.as_ref().unwrap().cigar.as_ref().unwrap();
            let target_offset = m.target_start as usize;
            let target_len = m.target_len as usize;
            let clip_start = m.query_start as usize;
            let clip_end = seq_len - m.query_end as usize;
            let processed_cigar = process_mm_cigar(
                cigar, 
                target_offset, target_len,
                clip_start, clip_end
            )?;
            
            // regions that overlap from this processed segment
            let pc_start = target_offset - clip_start;
            let pc_end = (m.target_end as usize) + clip_end;

            // create the mapping also
            let nm = m.alignment.as_ref().unwrap().nm as usize;
            let unmapped = seq_len - (m.query_end - m.query_start) as usize;
            let mapping_stats = MappingStats::new(seq_len, nm, unmapped);
            
            // make sure our last count equals the NM + unmapped; otherwise we goofed while processing
            assert_eq!(processed_cigar.len(), target_len+1); // we should have a lookup before/after each base
            assert_eq!(*processed_cigar.last().unwrap(), nm + unmapped);
            (Some(processed_cigar), Some(mapping_stats), pc_start..pc_end)
        } else {
            (None, None, 0..0)
        };

        // add the relevant info
        self.full_mappings.push(mapping);
        self.full_mapping_stats.push(opt_mapping_stats);
        self.processed_cigars.push(opt_processed_cigar);
        self.processed_ranges.push(processed_range);

        Ok(())
    }

    /// Wrapper function for determining if this match is better than the other
    pub fn is_better_match(&self, rhs: &Self) -> Result<bool, SimpleError> {
        // sanity check
        if self.processed_cigars.len() != rhs.processed_cigars.len() {
            bail!("RHS has different processed cigar length");
        }

        // useful when we want to get at something specific in our debug comparisons, should be false in all commits
        let manual_debug = false; //self.haplotype == "HLA:HLA24487";

        // we need to compare these sequentially
        for (i, (opt_lhs_pc, opt_rhs_pc)) in self.processed_cigars.iter().zip(rhs.processed_cigars.iter()).enumerate() {
            if let (Some(lhs_pc), Some(rhs_pc)) = (opt_lhs_pc, opt_rhs_pc) {
                // get the two ranges and figure out the overlap
                let lhs_range = &self.processed_ranges[i];
                let rhs_range = &rhs.processed_ranges[i];
                let overlap_start = lhs_range.start.max(rhs_range.start);
                let overlap_end = lhs_range.end.min(rhs_range.end);
                assert!(overlap_start < overlap_end);

                let lhs_nm = lhs_pc[overlap_end] - lhs_pc[overlap_start];
                let rhs_nm = rhs_pc[overlap_end] - rhs_pc[overlap_start];

                if manual_debug {
                    debug!("DEBUG_MODE");
                    debug!("lhs: {}", self.haplotype);
                    debug!("rhs: {}", rhs.haplotype);
                    debug!("overlap: {overlap_start}..{overlap_end}");
                    debug!("lhs_nm: {lhs_nm}");
                    debug!("rhs_nm: {rhs_nm}");
                    debug!("lhs_m : {:?}", self.full_mappings[i]);
                    debug!("rhs_m : {:?}", rhs.full_mappings[i]);
                }

                match lhs_nm.cmp(&rhs_nm) {
                    std::cmp::Ordering::Less => {
                        // LHS is better
                        return Ok(true);
                    },
                    std::cmp::Ordering::Equal => {
                        // equal at this level, so iterate to next
                    },
                    std::cmp::Ordering::Greater => {
                        // RHS is better
                        return Ok(false);
                    },
                };
            } else if opt_lhs_pc.is_none() && opt_rhs_pc.is_none() {
                // both are absent, iterate
            } else if opt_lhs_pc.is_some() {
                // LHS has a value, but RHS does not; mark LHS as better
                return Ok(true);
            } else {
                // RHS has a value, but LHS does not; mark RHS as better
                assert!(opt_rhs_pc.is_some());
                return Ok(false);
            }
        }

        // the following is only programmed with cDNA & DNA in mind; verify that here
        assert_eq!(self.full_mapping_stats.len(), 2);
        assert_eq!(rhs.full_mapping_stats.len(), 2);

        // tie-break by comparing the scores from end-to-end mappings
        let lhs_stats = HlaMappingStats::from_mapping_stats(
            self.full_mapping_stats[0].clone(),
            self.full_mapping_stats[1].clone()
        );
        let rhs_stats = HlaMappingStats::from_mapping_stats(
            rhs.full_mapping_stats[0].clone(),
            rhs.full_mapping_stats[1].clone()
        );
        Ok(lhs_stats.mapping_score() < rhs_stats.mapping_score())
    }

    // getters
    pub fn haplotype(&self) -> &str {
        &self.haplotype
    }

    pub fn full_mappings(&self) -> &[Option<minimap2::Mapping>] {
        &self.full_mappings
    }

    pub fn full_mapping_stats(&self) -> &[Option<MappingStats>] {
        &self.full_mapping_stats
    }
}

/// This will load a cigar and output a Vec of length equal to the target sequence + 1 with a value equal to the number of edits before that position.
/// Nice little shortcut for making sure we are comparing equivalent regions.
/// The value at index "i" is the number of edits before position "i" in the string.
/// Thus, the first value in this vec is always 0, and the last should equal NM + unmapped.
fn process_mm_cigar(cigar: &[(u32, u8)], target_offset: usize, target_len: usize, clip_start: usize, clip_end: usize) -> Result<Vec<usize>, SimpleError> {
    // alignment starts target_offset into the vec, so everything before is a 0
    assert!(clip_start <= target_offset);
    let mut ret = vec![0; target_offset - clip_start + 1]; // there should always be at least 1 zero at the start
    let mut current_nm = 0;

    // if we have soft clipping at the start, add a region that is basically a bunch of mismatches
    for _i in 0..clip_start {
        current_nm += 1;
        ret.push(current_nm);
    }

    for &(length, cigar_type) in cigar.iter() {
        // cigar types: https://github.com/lh3/minimap2/blob/69e36299168d739dded1c5549f662793af10da83/minimap.h#L57
        match cigar_type {
            // I - insertion of bases, we should increase current NM but don't insert anything
            1 => { current_nm += length as usize; },
            // D | X - deletion of reference character or mismatch to reference; either way, we increase by one for each base
            2 | 8 => {
                for _i in 0..length {
                    current_nm += 1;
                    ret.push(current_nm);
                }
            },
            // = - matches, so just copy NM value at this point
            7 => ret.extend(std::iter::repeat(current_nm).take(length as usize)),
            // we should not have any of the others
            unexpected => bail!("Unexpected cigar type: {unexpected}")
        };
    }

    // if we have soft clipping at the end, add a region that is basically a bunch of mismatches
    for _i in 0..clip_end {
        current_nm += 1;
        ret.push(current_nm);
    }

    // extend to fill out whatever remains, with 1 extra so we have lookups after each position
    assert!(ret.len() <= target_len+1);
    ret.extend(std::iter::repeat(current_nm).take(target_len+1  - ret.len()));
    Ok(ret)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_process_mm_cigar() {

        // (length, type)
        let cigar = [
            (2, 7), // match
            (1, 8), // mismatch,
            (2, 7), // match
            (1, 1), // insertion
            (2, 7), // match
            (1, 2), // deletion
            (2, 7), // match
        ];

        // first try exact match overlap
        // ==X==I==D==
        let expected_pc = vec![0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 3];
        let target_offset = 0;
        let target_len = 10;
        let clip_start = 0;
        let clip_end = 0;
        let result = process_mm_cigar(&cigar, target_offset, target_len, clip_start, clip_end).unwrap();
        assert_eq!(expected_pc, result);

        // Now add some clipping and offsets to the test
        // -SS==X==I==D==SSS--
        let expected_pc = vec![0, 0, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 7, 8, 8, 8];
        let target_offset = 3;
        let target_len = 18;
        let clip_start = 2;
        let clip_end = 3;
        let result = process_mm_cigar(&cigar, target_offset, target_len, clip_start, clip_end).unwrap();
        assert_eq!(expected_pc, result);
    }

}