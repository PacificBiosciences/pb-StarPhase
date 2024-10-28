
use serde::Serialize;
use std::ops::AddAssign;

/// Wraps the mapping stats for read to an HLA locus
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct MappingStats {
    /// the length of the sequence
    seq_len: usize,
    /// the NM tag for the mapping
    nm: usize,
    /// the number of unmapped based
    unmapped: usize,
    /// the number of clipped bases at the start
    clipped_start: Option<usize>,
    /// the number of clipped bases at the end
    clipped_end: Option<usize>
}

impl MappingStats {
    /// Basic constructor, but will verify some assumptions
    /// # Arguments
    /// * `seq_len` - the total length of the fragment that was mapped
    /// * `nm` - the edit distance (NM tag) of the mapping
    /// * `unmapped` - the number of unmapped bases
    pub fn new(
        seq_len: usize, nm: usize, unmapped: usize
    ) -> MappingStats {
        MappingStats {
            seq_len,
            nm, 
            unmapped, 
            clipped_start: None, 
            clipped_end: None
        }
    }

    /// Basic constructor, but will verify some assumptions
    /// # Arguments
    /// * `seq_len` - the total length of the fragment that was mapped
    /// * `nm` - the edit distance (NM tag) of the mapping
    /// * `unmapped` - the number of unmapped bases
    /// * `clipped_start` - the number of clipped bases at the start
    /// * `clipped_end` - the number of clipped bases at the end
    pub fn new_with_clippings(
        seq_len: usize, nm: usize, unmapped: usize,
        clipped_start: usize, clipped_end: usize
    ) -> MappingStats {
        MappingStats {
            seq_len,
            nm, 
            unmapped, 
            clipped_start: Some(clipped_start), 
            clipped_end: Some(clipped_end)
        }
    }

    /// Calculates the default mapping scores for this mapping.
    /// This approach penalizes based on both the number of mismatches and unmapped bases.
    pub fn mapping_score(&self) -> MappingScore {
        // let dna_score = MappingScore::score_value(self.seq_len, self.nm, self.unmapped);
        // MappingScore::new(dna_score)
        self.custom_score(true)
    }

    /// Calculates a parameterized mapping score.
    /// # Arguments
    /// * `penalize_unmapped` - if True, this will penalize the score based on unmapped reads
    pub fn custom_score(&self, penalize_unmapped: bool) -> MappingScore {
        let seq_len = if penalize_unmapped {
            self.seq_len
        } else {
            // remove unmapped from consideration for length
            self.seq_len - self.unmapped
        };
        let nm = self.nm;
        let unmapped = if penalize_unmapped {
            self.unmapped
        } else {
            // remove unmapped as a penalty
            0
        };
        let dna_score = MappingScore::score_value(seq_len, nm, unmapped);
        MappingScore::new(dna_score)
    }

    /// Wrapper for writing a scoring string to the screen
    pub fn score_string(&self) -> String {
        self.custom_score_string(true)
    }

    /// Wrapper for writing a custom scoring string to the screen.
    pub fn custom_score_string(&self, penalize_unmapped: bool) -> String {
        let score = self.custom_score(penalize_unmapped);
        if penalize_unmapped {
            // unmapped is a penalty
            format!("{:.5}=({}+{})/{}", score.score(), self.nm, self.unmapped, self.seq_len)
        } else {
            // unmapped is removed from counting
            format!("{:.5}={}/({}-{})", score.score(), self.nm, self.seq_len, self.unmapped)
        }
    }

    // getters
    pub fn clipped_start(&self) -> Option<usize> {
        self.clipped_start
    }

    pub fn clipped_end(&self) -> Option<usize> {
        self.clipped_end
    }

    pub fn nm(&self) -> usize {
        self.nm
    }

    pub fn unmapped(&self) -> usize {
        self.unmapped
    }
}


/// Contains the score for aligning an HLA sequence against a read.
/// This is basically an error rate, defined as (edit_distance + unmapped) / seq_len
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct MappingScore {
    /// the score for the alignment
    score: f64
}

impl MappingScore {
    /// Simple constructor placeholder
    pub fn new(score: f64) -> MappingScore {
        MappingScore {
            score
        }
    }

    /// Returns the worst possible value for a mapping
    pub fn worst_value() -> f64 {
        1.0
    }

    // Returns the worst possible score (cDNA + DNA) for a mapping
    pub fn worst_score() -> MappingScore {
        MappingScore {
            score: MappingScore::worst_value()
        }
    }

    /// Returns an empty score so we can accumulate scores with AddAssign
    pub fn zero_score() -> MappingScore {
        MappingScore {
            score: 0.0
        }
    }

    /// Convenient comparator since Eq and Ord cannot be derived with f64
    pub fn min(self, other: MappingScore) -> MappingScore {
        if self <= other {
            self
        } else {
            other
        }
    }

    /// Calculates the harmonic mean of a bunch of scores.
    /// # Arguments
    /// * `scores` - the values to derive the harmonic mean from
    /// # Panics
    /// * if any of the provided scores are 0.0
    pub fn harmonic_mean(scores: &[MappingScore]) -> MappingScore {
        let mut harmonic_sum: f64 = 0.0;
        
        for score in scores.iter() {
            assert!(score.score() > 0.0, "dna_score must be > 0.0");
            harmonic_sum += 1.0 / score.score();
        }

        MappingScore {
            score: if harmonic_sum > 0.0 { scores.len() as f64 / harmonic_sum } else { 0.0 }
        }
    }

    pub fn score(&self) -> f64 {
        self.score
    }

    /// Calculates the score of a set of values.
    /// If `map_nm + unmapped == 0.0`, this will add a partial error to keep the score above 0 (harmonic mean requirement).
    pub fn score_value(mapping_len: usize, map_nm: usize, unmapped: usize) -> f64 {
        let numerator = ((map_nm + unmapped) as f64).max(0.1);
        let denominator = (mapping_len) as f64;
        numerator / denominator
    }
}

// we are not really using this anymore, but it will not hurt to keep it
impl AddAssign for MappingScore {
    fn add_assign(&mut self, rhs: Self) {
        // an individual score is at most 1.0, so we can just add them together
        self.score += rhs.score();
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mapping_stats() {
        let mapping_stats = MappingStats::new(
            10, 1, 0,
        );
        assert_eq!(mapping_stats.mapping_score(), MappingScore::new(0.1));
    }

    #[test]
    fn test_score_min() {
        let s1 = MappingScore::new(1.0);
        let s2 = MappingScore::new(0.9);
        let s3 = MappingScore::new(0.2);

        assert_eq!(s1.min(s2), s2);
        assert_eq!(s1.min(s3), s3);
        assert_eq!(s2.min(s3), s3);
    }

    #[test]
    fn test_harmonic_mean() {
        let s1 = MappingScore::new(0.2);
        let s2 = MappingScore::new(0.4);
        let s3 = MappingScore::new(0.2);

        let arrayed = [s1, s2, s3];
        let expected = MappingScore::new(
            3.0 / (5.0+2.5+5.0), // 0.24
        );
        assert_eq!(expected, MappingScore::harmonic_mean(&arrayed));
    }
}