
use serde::Serialize;
use std::ops::AddAssign;

use crate::data_types::mapping::{MappingScore, MappingStats};

/// Wraps the mapping stats for read to an HLA locus
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct HlaMappingStats {
    /// the cDNA mapping stats
    cdna_stats: Option<MappingStats>,
    /// the DNA mapping stats
    dna_stats: Option<MappingStats>
}

impl HlaMappingStats {
    /// Basic constructor, but will verify some assumptions
    /// # Arguments
    /// * `cdna_len` - the total length of the cDNA fragment that was mapped to the read
    /// * `cdna_nm` - the edit distance (NM tag) of the mapping against the cDNA
    /// * `cdna_unmapped` - the number of unmapped bases in the cDNA
    /// * `dna_len` - the total length of the DNA fragment that was mapped to the read
    /// * `dna_nm` - the edit distance (NM tag) of the mapping against the DNA
    /// * `dna_unmapped` - the number of unmapped bases in the DNA
    /// # Panics
    /// * if both the cDNA and DNA are unset
    /// * if cDNA fields are only partially filled, i.e. all must be Some or None
    /// * if DNA fields are only partially filled, i.e. all must be Some or None
    pub fn new(
        cdna_len: Option<usize>, cdna_nm: Option<usize>, cdna_unmapped: Option<usize>,
        dna_len: Option<usize>, dna_nm: Option<usize>, dna_unmapped: Option<usize>
    ) -> HlaMappingStats {
        // make sure at least one is set
        assert!(cdna_len.is_some() || dna_len.is_some());
        // make sure all of cDNA is either Some or None
        assert!(cdna_len.is_some() == cdna_nm.is_some() && cdna_len.is_some() == cdna_unmapped.is_some());
        // make sure all of DNA is either Some or None
        assert!(dna_len.is_some() == dna_nm.is_some() && dna_len.is_some() == dna_unmapped.is_some());

        let cdna_stats = cdna_len.map(|l| {
            MappingStats::new(l, cdna_nm.unwrap(), cdna_unmapped.unwrap())
        });

        let dna_stats = dna_len.map(|l| {
            MappingStats::new(l, dna_nm.unwrap(), dna_unmapped.unwrap())
        });

        HlaMappingStats {
            cdna_stats, dna_stats
        }
    }

    /// Wrapper constructor for directly pulling in establishing mapping stats.
    /// Does not perform sanity checks on the input.
    pub fn from_mapping_stats(cdna_stats: Option<MappingStats>, dna_stats: Option<MappingStats>) -> Self {
        Self {
            cdna_stats,
            dna_stats
        }
    }

    /// Calculates the mapping scores for this mapping.
    /// Defaults to worst score when mappings are unavailable.
    /// # Panics
    /// * if a dev somehow creates one with incomplete cDNA or DNA entries, don't do that!
    pub fn mapping_score(&self) -> HlaMappingScore {
        let cdna_score = match self.cdna_stats.as_ref() {
            Some(cds) => {
                cds.mapping_score()
            },
            None => MappingScore::worst_score()
        };
        let dna_score = match self.dna_stats.as_ref() {
            Some(ds) => {
                ds.mapping_score()
            },
            None => MappingScore::worst_score()
        };
        HlaMappingScore::new(cdna_score, dna_score)
    }

    /// Wrapper for writing a scoring string to the screen
    /// # Panics
    /// * if a dev somehow creates one with incomplete cDNA or DNA entries, don't do that!
    pub fn score_string(&self) -> String {
        let cdna_str = match self.cdna_stats.as_ref() {
            Some(cds) => cds.score_string(),
            None => "None".to_string()
        };
        let dna_str = match self.dna_stats.as_ref() {
            Some(ds) => ds.score_string(),
            None => "None".to_string()
        };
        format!(
            "{}, {}", cdna_str, dna_str
        )
    }

    pub fn has_cdna(&self) -> bool {
        self.cdna_stats.is_some()
    }

    pub fn has_dna(&self) -> bool {
        self.dna_stats.is_some()
    }
}


/// Contains the score for aligning an HLA sequence against a read.
/// This is basically an error rate, defined as (edit_distance + unmapped) / seq_len
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct HlaMappingScore {
    /// the score for the cDNA alignment
    cdna_score: MappingScore,
    /// the score for the DNA alignment
    dna_score: MappingScore
}

impl HlaMappingScore {
    /// Simple constructor placeholder
    pub fn new(cdna_score: MappingScore, dna_score: MappingScore) -> HlaMappingScore {
        HlaMappingScore {
            cdna_score, dna_score
        }
    }

    /// Useful wrapper when we just want to directly pass values in
    pub fn from_values(cdna_value: f64, dna_value: f64) -> HlaMappingScore {
        let cdna_score = MappingScore::new(cdna_value);
        let dna_score = MappingScore::new(dna_value);
        HlaMappingScore {
            cdna_score, dna_score
        }
    }

    // Returns the worst possible score (cDNA + DNA) for a mapping
    pub fn worst_score() -> HlaMappingScore {
        HlaMappingScore {
            cdna_score: MappingScore::worst_score(),
            dna_score: MappingScore::worst_score()
        }
    }

    /// Returns an empty score so we can accumulate scores with AddAssign
    pub fn zero_score() -> HlaMappingScore {
        HlaMappingScore {
            cdna_score: MappingScore::zero_score(),
            dna_score: MappingScore::zero_score()
        }
    }

    /// Convenient comparator since Eq and Ord cannot be derived with f64
    pub fn min(self, other: HlaMappingScore) -> HlaMappingScore {
        if self <= other {
            self
        } else {
            other
        }
    }

    pub fn cdna_score(&self) -> MappingScore {
        self.cdna_score
    }

    pub fn dna_score(&self) -> MappingScore {
        self.dna_score
    }
}

// we are not really using this anymore, but it will not hurt to keep it
impl AddAssign for HlaMappingScore {
    fn add_assign(&mut self, rhs: Self) {
        // an individual score is at most 1.0, so we can just add them together
        self.cdna_score += rhs.cdna_score;
        self.dna_score += rhs.dna_score;
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mapping_stats() {
        let mapping_stats = HlaMappingStats::new(
            Some(10), Some(1), Some(0),
            Some(20), Some(0), Some(1)
        );
        assert_eq!(mapping_stats.mapping_score(), HlaMappingScore::from_values(0.1, 0.05));
    }

    #[test]
    #[should_panic]
    fn test_empty_stats() {
        let _mapping_stats = HlaMappingStats::new(
            None, None, None,
            None, None, None
        );
    }

    #[test]
    #[should_panic]
    fn test_partial_cdna_stats() {
        let _mapping_stats = HlaMappingStats::new(
            Some(10), None, None,
            None, None, None
        );
    }

    #[test]
    #[should_panic]
    fn test_partial_dna_stats() {
        let _mapping_stats = HlaMappingStats::new(
            None, None, None,
            Some(10), None, None
        );
    }

    #[test]
    fn test_score_min() {
        let s1 = HlaMappingScore::from_values(1.0, 0.5);
        let s2 = HlaMappingScore::from_values(0.9, 1.0);
        let s3 = HlaMappingScore::from_values(1.0, 0.2);

        assert_eq!(s1.min(s2), s2);
        assert_eq!(s1.min(s3), s3);
        assert_eq!(s2.min(s3), s2);
    }
}