
use serde::{Deserialize, Serialize};

/// Wrapper for basic region coordinates
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Coordinates {
    /// Chromosome string
    chrom: String,
    /// 0-based start, inclusive
    start: u64,
    /// 0-based end, exclusive
    end: u64
}

impl Coordinates {
    /// Typical constructor with some verification
    pub fn new(chrom: String, start: u64, end: u64) -> Coordinates {
        assert!(start <= end);
        Coordinates {
            chrom, start, end
        }
    }

    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    pub fn start(&self) -> u64 {
        self.start
    }

    pub fn end(&self) -> u64 {
        self.end
    }

    /// Wrapper for sending to htslib fetch
    pub fn fetch_definition(&self) -> (&str, u64, u64) {
        (&self.chrom, self.start, self.end)
    }
}

impl std::fmt::Display for Coordinates {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // start and end are 0-based, so shift start+1
        write!(f, "{}:{}-{}", self.chrom, self.start+1, self.end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coordinates() {
        let chrom = "chr1".to_string();
        let start = 10;
        let end = 20;
        let coordinate = Coordinates::new(chrom.clone(), start, end);
        assert_eq!(coordinate.fetch_definition(), (chrom.as_str(), start, end));
    }

    #[test]
    #[should_panic]
    fn test_bad_coordinates() {
        let chrom = "chr1".to_string();
        let start = 10;
        let end = 5;
        let _coordinate = Coordinates::new(chrom.clone(), start, end);
    }
}