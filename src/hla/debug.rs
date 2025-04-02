
use serde::Serialize;
use simple_error::bail;
use std::collections::BTreeMap;
use std::collections::btree_map::Entry::{Occupied, Vacant};

/// Primary object that gets converted into JSON
#[derive(Default, Serialize)]
pub struct HlaDebug {
    /// Each gene & read has a detailed set of mappings
    read_mapping_stats: BTreeMap<String, BTreeMap<String, ReadMappingStats>>,
    /// Each gene has stats on the dual consensus stats
    dual_passing_stats: Option<BTreeMap<String, DualPassingStats>>
}

impl HlaDebug {
    /// Constructor, default works for now, but we may change in the future as things get added
    pub fn new() -> HlaDebug {
        Default::default()
    }

    /// Adds a read to the debug collection
    /// # Arguments
    /// * `gene` - the gene name we are inserting into
    /// * `qname` - the read name to insert
    /// * `stats` - the HLA sequence mapping stats for this read
    /// # Errors
    /// * if the read has already been added
    pub fn add_read(&mut self, gene: String, qname: String, stats: ReadMappingStats) -> Result<(), Box<dyn std::error::Error>> {
        let gene_entry = self.read_mapping_stats.entry(gene).or_default();
        match gene_entry.entry(qname.clone()) {
            Occupied(_entry) => {
                bail!("Entry {qname} is already occupied");
            },
            Vacant(entry) => {
                entry.insert(stats);
                Ok(())
            }
        }
    }

    /// Adds dual passing stats to the debug collection
    /// # Arguments
    /// * `gene` - the gene name we are inserting into
    /// * `stats` - the stats to add for the gene
    pub fn add_dual_passing_stats(&mut self, gene: String, stats: DualPassingStats) -> Result<(), Box<dyn std::error::Error>> {
        if self.dual_passing_stats.is_none() {
            self.dual_passing_stats = Some(Default::default());
        }

        let dps = self.dual_passing_stats.as_mut().unwrap();
        match dps.entry(gene) {
            Occupied(occupied_entry) => {
                bail!("Entry {} is already occupied", occupied_entry.key());
            },
            Vacant(vacant_entry) => {
                vacant_entry.insert(stats);
                Ok(())
            }
        }
    }
}

/// Wrapper for an individual read stats
#[derive(Clone, Debug, Default, Serialize)]
pub struct ReadMappingStats {
    /// The best matching ID
    best_match_id: Option<String>,
    /// The best matching star allele
    best_match_star: Option<String>,
    /// The details around mappings for the read
    mapping_stats: BTreeMap<String, PairedMappingStats>
}

impl ReadMappingStats {
    /// Generic constructor, may add some constraints later
    pub fn new() -> ReadMappingStats {
        Default::default()
    }

    /// sets the best match, this usually is not known when constructed
    pub fn set_best_match(&mut self, new_best_id: String, new_best_star: String) {
        self.best_match_id = Some(new_best_id);
        self.best_match_star = Some(new_best_star);
    }
    
    pub fn best_match_id(&self) -> Option<&str> {
        self.best_match_id.as_deref()
    }

    pub fn best_match_star(&self) -> Option<&str> {
        self.best_match_star.as_deref()
    }

    pub fn mapping_stats(&self) -> &BTreeMap<String, PairedMappingStats> {
        &self.mapping_stats
    }

    /// Adds a mapping to the stats for this read
    /// # Arguments
    /// * `hla_id` - the HLA ID for the mapping
    /// * `cdna_mm2` - the minimap2 mapping for the cDNA sequence, optional
    /// * `dna_mm2` - the minimap2 mapping for the DNA sequence, optional
    /// # Errors
    /// * if the HLA ID has already been inserted
    pub fn add_mapping(&mut self, hla_id: String, cdna_mm2: Option<&minimap2::Mapping>, dna_mm2: Option<&minimap2::Mapping>) 
        -> Result<(), Box<dyn std::error::Error>> {
        match self.mapping_stats.entry(hla_id.clone()) {
            Occupied(_entry) => {
                bail!("Entry {hla_id} is already occupied!");
            },
            Vacant(entry) => {
                let cdna_mapping: Option<DetailedMappingStats> = cdna_mm2.map(DetailedMappingStats::from_mapping);
                let dna_mapping: Option<DetailedMappingStats> = dna_mm2.map(DetailedMappingStats::from_mapping);
                let paired_mapping = PairedMappingStats {
                    cdna_mapping, dna_mapping
                };
                entry.insert(paired_mapping);
                Ok(())
            }
        }
    }
}

/// Wrapper for a read to allele comparison, both cDNA and DNA are optional
#[derive(Clone, Debug, Serialize)]
pub struct PairedMappingStats {
    /// the cDNA mapping stats
    cdna_mapping: Option<DetailedMappingStats>,
    /// the DNA mapping stats
    dna_mapping: Option<DetailedMappingStats>
}

/// Wrapper for an individual mapping statistic set
#[derive(Clone, Debug, Serialize)]
struct DetailedMappingStats {
    /// length of the sequence
    query_len: usize,
    /// length of the target
    target_len: usize,
    /// length of the match
    match_len: usize,
    /// number of mismatches
    nm: usize,
    /// number of unmapped bases relative to query
    query_unmapped: usize,
    /// number of unmapped bases relative to target
    target_unmapped: usize,
    /// Cigar string
    cigar: String,
    /// The incorrect matching bases
    md: String
}

impl DetailedMappingStats {
    /// Creates a detailed mapping stat we can serialize from the minimap2 mapping
    /// # Arguments
    /// * `mapping` - the minimap2 mapping, we will pull relevant info out
    /// # Panics
    /// * if the mapping is missing detailed alignment information
    fn from_mapping(mapping: &minimap2::Mapping) -> DetailedMappingStats {
        let query_len = mapping.query_len.unwrap().get() as usize;
        let target_len = mapping.target_len as usize;
        let match_len = mapping.match_len as usize;
        let alignment = mapping.alignment.as_ref().unwrap();
        let nm = alignment.nm as usize;
        let query_unmapped = query_len - (mapping.query_end - mapping.query_start) as usize;
        let target_unmapped = target_len - (mapping.target_end - mapping.target_start) as usize;
        let cigar = alignment.cigar_str.as_ref().unwrap().clone();
        let md = alignment.md.as_ref().unwrap().clone();

        DetailedMappingStats {
            query_len,
            target_len,
            match_len,
            nm,
            query_unmapped,
            target_unmapped,
            cigar,
            md
        }
    }
}

#[derive(Default, Serialize)]
pub struct DualPassingStats {
    /// If true, indicates that a dual consensus (heterozygous) will get reported
    is_passing: bool,
    /// If true, indicates that a dual consensus was found
    is_dual: bool,
    /// Counts assigned to consensus 1
    counts1: Option<usize>,
    /// Counts assigned to consensus 2
    counts2: Option<usize>,
    /// Minor count fraction
    maf: Option<f64>,
    /// CDF of the result given the input targets
    cdf: Option<f64>
}

impl DualPassingStats {
    /// Creates stats for a dual consensus result
    pub fn new_dual(is_passing: bool, counts1: usize, counts2: usize, maf: f64, cdf: f64) -> Self {
        Self {
            is_passing,
            is_dual: true,
            counts1: Some(counts1),
            counts2: Some(counts2),
            maf: Some(maf),
            cdf: Some(cdf)
        }
    }

    /// Creates stats for a non-dual consensus result
    pub fn new_non_dual() -> Self {
        Self {
            is_passing: false,
            is_dual: false,
            ..Default::default()
        }
    }

    pub fn is_passing(&self) -> bool {
        self.is_passing
    }
}