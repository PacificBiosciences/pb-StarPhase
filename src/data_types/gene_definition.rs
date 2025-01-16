
use bio::bio_types::strand::Strand;
use bio::io::gff;
use log::{debug, info, trace, warn};
use serde::{Deserialize, Serialize};
use simple_error::{bail, SimpleError};
use std::collections::{BTreeMap, BTreeSet};
use std::collections::btree_map::Entry;
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::path::Path;

use crate::data_types::coordinates::Coordinates;

/// RefSeq stores the lastest coordinate information here; a version is located inside with #!annotation-source
const REFSEQ_LATEST: &str = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz";

/// Useful wrapper for a single gene
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct GeneDefinition {
    /// The gene name
    gene_name: String,
    /// Coordinates for the gene
    coordinates: Coordinates,
    /// if True, this gene is on the forward strand; otherwise, reverse-complement
    is_forward_strand: bool,
    /// MANE transcript ID
    transcript_id: Option<String>,
    /// Ordered exon coordinates; these are *always* in forward order, regardless of the gene strand
    exons: Vec<Coordinates>,
    /// If true, then these gene is not guaranteed to be present in a sample; could be absent, hemizygous, or diploid
    #[serde(default)] // default to False if not present
    is_absent_capable: bool
}

impl GeneDefinition {
    /// Creates a new definition with basic information. Exons are added individually later.
    /// # Arguments
    /// * `gene_name` - the name associated with the definition; e.g. "HLA-A"
    /// * `coordinates` - the coordinates from RefSeq, 0-based
    /// * `is_forward_strand` - if True, this is on the reference strand
    pub fn new(gene_name: String, coordinates: Coordinates, is_forward_strand: bool) -> Self {
        Self {
            gene_name,
            coordinates, 
            is_forward_strand,
            transcript_id: None,
            exons: vec![],
            is_absent_capable: false
        }
    }

    /// Adds a transcript ID to the gene definition
    /// # Arguments
    /// * `transcript_id` - the MANE transcript name for this gene
    /// * `update_coordinates` - if Some, then the coordinates are replaced with these new ones
    /// # Errors
    /// * if the gene definition already has a transcript ID
    pub fn add_transcript_id(&mut self, transcript_id: String, update_coordinates: Option<Coordinates>) -> Result<(), SimpleError> {
        if self.transcript_id.is_some() {
            bail!("Transcript ID has already been added to GeneDefinition: {self:?}");
        } else {
            self.transcript_id = Some(transcript_id);
            if let Some(uc) = update_coordinates {
                self.coordinates = uc;
            }
            Ok(())
        }
    }

    /// Adds the next exon coordinate to the gene definition, must be added exon order.
    /// For '-' strand, this means in reverse.
    /// # Arguments
    /// * `new_coordinate` - the coordinates that are getting added, must be in order
    /// # Errors
    /// * if the chromosome provided does not match the gene definition
    /// * if the new coordinate is not sequential with the existing coordinates, no overlaps allowed
    pub fn add_exon(&mut self, new_coordinate: Coordinates) -> Result<(), SimpleError> {
        // verify same chromosome
        if self.coordinates.chrom() != new_coordinate.chrom() {
            bail!("Exon chromosome does not match gene chromosome.");
        }

        // in the file, they are in strand order; so if on the '-' strand, we need to insert in reverse
        if self.is_forward_strand {
            // make sure it's after the previous
            if let Some(last_exon) = self.exons.last() {
                if last_exon.end() > new_coordinate.start() {
                    bail!("New exon ({}) must start after the previous exon; previous: {}, new: {}", self.gene_name, last_exon, new_coordinate);
                }
            }

            // forward strand, add it to the end
            self.exons.push(new_coordinate);
        } else {
            // make sure it's before the first
            if let Some(first_exon) = self.exons().first() {
                if new_coordinate.end() > first_exon.start() {
                    bail!("New exon ({}) must start before the first exon; first: {}, new: {}", self.gene_name, first_exon, new_coordinate);
                }
            }

            // reverse strand, prepend it
            self.exons.insert(0, new_coordinate);
        }

        Ok(())
    }

    /// This will check two new coordinates to see if we need to extend our existing definition.
    /// Returns true if a change is made.
    pub fn extend_coordinates(&mut self, alt_start: u64, alt_end: u64) -> bool {
        let start_shifted = self.coordinates.extend_start(alt_start);
        let end_shifted = self.coordinates.extend_end(alt_end);
        start_shifted || end_shifted
    }

    /// Call if a gene can be absent
    pub fn set_absent_capable(&mut self) {
        self.is_absent_capable = true;
    }

    // getters
    pub fn gene_name(&self) -> &str {
        &self.gene_name
    }

    pub fn coordinates(&self) -> &Coordinates {
        &self.coordinates
    }

    pub fn is_forward_strand(&self) -> bool {
        self.is_forward_strand
    }

    pub fn transcript_id(&self) -> Option<&str> {
        self.transcript_id.as_deref()
    }

    pub fn exons(&self) -> &[Coordinates] {
        &self.exons
    }

    pub fn is_absent_capable(&self) -> bool {
        self.is_absent_capable
    }
}

/// Contains all the gene definitions from a source, typically RefSeq (either URL or file).
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct GeneCollection {
    /// Version string for the gene collection
    version: String,
    /// The gene lookup
    gene_dict: BTreeMap<String, GeneDefinition>
}

impl GeneCollection {
    /// General constructor
    pub fn new(version: String, gene_dict: BTreeMap<String, GeneDefinition>) -> Self {
        Self {
            version,
            gene_dict
        }
    }

    /// Loads the RefSeq collection via a URL
    /// # Arguments
    /// * `opt_url` - an optional URL to parse through, which is expected to be a .gff.gz file; if None, it will pull the latest
    /// * `gene_targets` - an optional set of gene targets; if None, it will load everything
    /// # Errors
    /// * if there is an error fetching the given URL
    /// * pretty much any parsing error or unexpected gene/mRNA/exon combination
    pub fn load_refseq_url(opt_url: Option<String>, gene_targets: Option<&BTreeSet<String>>) -> Result<Self, Box<dyn std::error::Error>> {
        // Use latest if nothing provided
        let refseq_url: String = opt_url.unwrap_or(REFSEQ_LATEST.to_string());
        info!("\tQuerying RefSeq via {refseq_url}");

        // hit the end point so we can parse it
        let result = reqwest::blocking::get(refseq_url)?;
        debug!("Response received.");

        // wrap it in a decoder and BufReader
        let gz_decoder = flate2::read::GzDecoder::new(result);
        let buf_reader = BufReader::new(gz_decoder);

        // now it's good to test
        Self::load_refseq(Box::new(buf_reader), gene_targets)
    }

    /// Loads RefSeq data from a given file path, mostly a wrapper for `load_refseq(...)`.
    /// # Arguments
    /// * `filename` - the file we want to parse through, which is expected to be a .gff.gz file
    /// * `gene_targets` - an optional set of gene targets; if None, it will load everything
    /// # Errors
    /// * pretty much any parsing error or unexpected gene/mRNA/exon combination
    pub fn load_refseq_file(filename: &Path, gene_targets: Option<&BTreeSet<String>>) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(filename).unwrap();
        let gz_decoder = flate2::read::GzDecoder::new(file);
        let buf_reader = BufReader::new(gz_decoder);
        Self::load_refseq(Box::new(buf_reader), gene_targets)
    }

    /// This does the bulk of the work, taking a generic boxed reader as input.
    /// # Arguments
    /// * `reader` - the file we want to parse through, any fetching or unzipping should happen before this
    /// * `gene_targets` - an optional set of gene targets; if None, it will load everything
    /// # Errors
    /// * pretty much any parsing error or unexpected gene/mRNA/exon combination
    fn load_refseq(mut reader: Box<dyn BufRead>, gene_targets: Option<&BTreeSet<String>>) -> Result<Self, Box<dyn std::error::Error>> {
        // this gets returned
        let mut gene_dict: BTreeMap<String, GeneDefinition> = Default::default();
        let mut opt_version: Option<String> = None;

        // this lets us track which transcript is the MANE one; key is transcript, value is gene
        let mut mane_transcripts: BTreeMap<String, String> = Default::default();

        // search for a version in the header; this does not come through the gff::Reader unfortunately
        let mut line: String = Default::default();
        while opt_version.is_none() {
            line.clear(); // need to clear buffer before reading
            let read_len = reader.read_line(&mut line)?;
            if read_len == 0 {
                bail!("Reached EOF before identifying RefSeq version");
            }

            if line.starts_with("##") {
                // no-op, header line
            } else if line.starts_with("#!") {
                let line_frags = line.trim().split(" ").collect::<Vec<&str>>();
                if line_frags[0] == "#!annotation-source" {
                    opt_version = Some(
                        line_frags[1..].join(" ")
                    );
                }
            } else {
                bail!("Reached end of header without finding a RefSeq version");
            }
        }

        let version = match opt_version {
            Some(v) => v,
            None => bail!("Version string was not found while parsing header.")
        };

        // parse through the file, which is a TSV
        let mut gff_reader = gff::Reader::new(reader, gff::GffType::GFF3);
        let mut chrom_dict: BTreeMap<String, String> = Default::default();
        
        for result in gff_reader.records() {
            let record = result?;
            
            // these are ordered in the way we expect to encounter them
            match record.feature_type() {
                "region" => {
                    // these can be chromosomes, so we need to parse them out
                    if let Some((chrom_id, chrom_str)) = parse_record_to_chromosome(&record) {
                        debug!("Chrom dict update: {chrom_id} -> {chrom_str}");
                        match chrom_dict.entry(chrom_id) {
                            Entry::Vacant(e) => e.insert(chrom_str),
                            Entry::Occupied(e) => {
                                let k = e.key();
                                bail!("Found duplicate chrom_id: \"{k}\"");
                            }
                        };
                    }
                },
                "gene" => {
                    // check if this is a passing gene entry
                    if let Some((gene_name, is_forward_strand)) = parse_record_gene_def(&record, gene_targets) {
                        // we will need full coordinates for later, parse them out here
                        let coordinates = parse_record_coordinates(&record, &chrom_dict)?;
                        let gene_def = GeneDefinition::new(gene_name, coordinates, is_forward_strand);
                        debug!("Gene definition update: {} -> {gene_def:?}", gene_def.gene_name());
                        match gene_dict.entry(gene_def.gene_name().to_string()) {
                            Entry::Vacant(e) => e.insert(gene_def),
                            Entry::Occupied(_e) => {
                                bail!("Found duplicate gene definition: \"{}\"", gene_def.gene_name());
                            }
                        };
                    }
                },
                "mRNA" => {
                    if let Some((transcript_id, gene_name, is_forward_strand)) = parse_record_mrna_def(&record, gene_targets) {
                        debug!("MANE trancript update: {transcript_id} -> {gene_name}");

                        // compare and update gene definition
                        if let Some(gene_def) = gene_dict.get_mut(&gene_name) {
                            if is_forward_strand != gene_def.is_forward_strand() {
                                bail!("Found MANE transcript with different strand orientation to gene definition: {} != {}", gene_name, transcript_id);
                            }
                            let mane_coordinates = parse_record_coordinates(&record, &chrom_dict)?;
                            debug!("MANE coordinates update: {mane_coordinates}");
                            gene_def.add_transcript_id(transcript_id.clone(), Some(mane_coordinates))?;
                        } else {
                            bail!("Found a MANE transcript but no core gene definition.");
                        }

                        // save the transcript -> gene entry for exon lookups
                        match mane_transcripts.entry(transcript_id) {
                            Entry::Vacant(e) => e.insert(gene_name),
                            Entry::Occupied(e) => {
                                bail!("Found multiple transcripts with same id: \"{}\"", e.key());
                            }
                        };
                    }
                },
                "exon" => {
                    if let Some((exon_id, transcript_id, is_forward_strand)) = parse_record_exon_def(&record, &mane_transcripts) {
                        // it should not be possible to panic from these unwraps due to upstream checks
                        let gene_name = mane_transcripts.get(&transcript_id).unwrap();
                        let gene_def = gene_dict.get_mut(gene_name).unwrap();

                        // make sure strands match
                        if is_forward_strand != gene_def.is_forward_strand() {
                            bail!("Found exon with different strand orientation to gene definition: {} != {}", gene_name, exon_id);
                        }
                        
                        // we need exon-level coordinates
                        let coordinates = parse_record_coordinates(&record, &chrom_dict)?;
                        debug!("Exon found for {gene_name}-{transcript_id}: {exon_id} => {coordinates}");
                        gene_def.add_exon(coordinates)?;
                    }
                },
                _ => {
                    // for debugging, but we can skip in general
                    // todo!("Unhandled feature type \"{ft}\": {record:?}")
                }
            }
        }

        Ok(Self {
            version,
            gene_dict
        })
    }

    /// Wrapper function that will go through our list of genes and create a copy of the relevant ones that are undefined.
    /// # Arguments
    /// * `copy_keys` - (copy_to_key, copy_from_key) indicating where to copy data to and from
    pub fn copy_missing_genes(&mut self, copy_keys: &BTreeMap<String, String>) -> Result<(), Box<dyn std::error::Error>> {
        for (copy_to_key, copy_from_key) in copy_keys.iter() {
            debug!("Copying gene definition for {copy_from_key} to {copy_to_key}...");
            let copy_value = self.gene_dict.get(copy_from_key)
                .ok_or(SimpleError::new(format!("Cannot copy definition from {copy_from_key} to {copy_to_key}; {copy_from_key} does not exist")))?
                .clone();
            match self.gene_dict.entry(copy_to_key.clone()) {
                Entry::Vacant(vacant_entry) => {
                    // we need to copy values into here
                    vacant_entry.insert(copy_value);
                },
                Entry::Occupied(_occupied_entry) => {
                    // this one has data already, so no need to copy
                },
            };
        }

        Ok(())
    }

    // getters
    pub fn version(&self) -> &str {
        &self.version
    }

    pub fn gene_dict(&self) -> &BTreeMap<String, GeneDefinition> {
        &self.gene_dict
    }

    pub fn gene_dict_mut(&mut self) -> &mut BTreeMap<String, GeneDefinition> {
        &mut self.gene_dict
    }
}

/// Parses a GFF record and returns an ID to chromosome pair if found.
/// # Arguments
/// * `record` - the record to parse
fn parse_record_to_chromosome(record: &gff::Record) -> Option<(String, String)> {
    if record.feature_type() == "region" && record.source() == "RefSeq" {
        let attributes = record.attributes();
        let genome = attributes.get("genome");
        if genome.is_some() && genome.unwrap() == "chromosome" {
            // we found a chrom mapping
            let chrom_key = record.seqname().to_string();
            let chrom_name = attributes.get("chromosome");
            match chrom_name {
                Some(cn) => {
                    // add the chr prefix if needed
                    let cn_mod = if cn.starts_with("chr") {
                        cn.clone()
                    } else {
                        format!("chr{cn}")
                    };
                    Some((chrom_key, cn_mod))
                },
                None => {
                    warn!("Found allowed gene record with no chromosome: {record:?}");
                    None
                }
            }
        } else {
            None
        }
    } else {
        None
    }
}

/// Parses a record and returns the gene name and strand if it matches one in our DB
/// # Arguments
/// * `record` - the record to parse
/// * `gene_targets` - an optional set of genes that are allowed; if None, then all genes are allowed
fn parse_record_gene_def(record: &gff::Record, gene_targets: Option<&BTreeSet<String>>) -> Option<(String, bool)> {
    // only gene records
    if record.feature_type() == "gene" &&
        // best ref seq is what we trust here; contains is because sometimes the entries have multiple like "BestRefSeq%2CGnomon"
        record.source().contains("BestRefSeq") &&
        // there are spurious records on the other contigs, so lets filter to the "NC_" prefixed for now
        record.seqname().starts_with("NC_") {
        if let Some(gene_name) = record.attributes().get("Name") {
            let allowed_gene = match gene_targets {
                Some(gt) => gt.contains(gene_name),
                None => true
            };

            if allowed_gene {
                match record.strand() {
                    Some(Strand::Forward) => Some((gene_name.clone(), true)),
                    Some(Strand::Reverse) => Some((gene_name.clone(), false)),
                    Some(Strand::Unknown) | None => {
                        warn!("Found allowed gene with no strand: {record:?}");
                        None
                    }
                }
            } else {
                None
            }
        } else {
            None
        }
    } else {
        None
    }
}

/// Parses a record and and returns the mRNA transcript ID, gene name, and orientation if it is the MANE transcript for one of our genes
/// # Arguments
/// * `record` - the record to parse
/// * `gene_targets` - an optional set of genes that are allowed; if None, then all genes are allowed
fn parse_record_mrna_def(record: &gff::Record, gene_targets: Option<&BTreeSet<String>>) -> Option<(String, String, bool)> {
    if record.feature_type() == "mRNA" &&
        record.source() == "BestRefSeq" &&
        record.seqname().starts_with("NC_") {
        if let Some(gene_name) = record.attributes().get("gene") {
            // check if it's a gene we care about
            let allowed_gene = gene_targets.is_none_or(|gt| gt.contains(gene_name));

            // also check if it's a main transcript
            let is_mane_transcript = record.attributes().get("tag").is_some_and(|tag| tag == "MANE Select");

            if allowed_gene && is_mane_transcript {
                // allowed and mane
                if let Some(transcript_id) = record.attributes().get("transcript_id") {
                    match record.strand() {
                        Some(Strand::Forward) => Some((transcript_id.clone(), gene_name.clone(), true)),
                        Some(Strand::Reverse) => Some((transcript_id.clone(), gene_name.clone(), false)),
                        Some(Strand::Unknown) | None => {
                            warn!("Found allowed MANE transcript record with no strand: {record:?}");
                            None
                        }
                    }
                } else {
                    warn!("Found allowed MANE transcript record with no transcript_id: {record:?}");
                    None
                }
            } else {
                None
            }
        } else {
            None
        }
    } else {
        None
    }
}

/// Parses a record and and returns the exon ID, mRNA transcript ID, and orientation if it's one we care about
/// # Arguments
/// * `record` - the record to parse
/// * `transcript_targets` - the set of transcripts that are targets so far; these *should* be in-line before any exons
fn parse_record_exon_def(record: &gff::Record, transcript_targets: &BTreeMap<String, String>) -> Option<(String, String, bool)> {
    if record.feature_type() == "exon" &&
        record.source() == "BestRefSeq" &&
        record.seqname().starts_with("NC_") {
        // check if the transcript_id is in the MANE target list
        if let Some(transcript_id) = record.attributes().get("transcript_id") {
            if transcript_targets.contains_key(transcript_id) {
                if let Some(exon_id) = record.attributes().get("ID") {
                    match record.strand() {
                        Some(Strand::Forward) => Some((exon_id.clone(), transcript_id.clone(), true)),
                        Some(Strand::Reverse) => Some((exon_id.clone(), transcript_id.clone(), false)),
                        Some(Strand::Unknown) | None => {
                            warn!("Found allowed exon record with no strand: {record:?}");
                            None
                        }
                    }
                } else {
                    warn!("Found allowed exon record with no ID: {record:?}");
                    None
                }
            } else {
                // not one we care about
                None
            }
        } else {
            // turns out this is pretty common, so put it in trace instead of warn for now
            trace!("Found exon record with no transcript ID: {record:?}");
            None
        }
    } else {
        None
    }
}

/// Parses a GFF record for the coordinates
/// # Arguments
/// * `record` - the record to parse
/// * `chrom_dict` - translation of RefSeq chromosome labels to normal labels
fn parse_record_coordinates(record: &gff::Record, chrom_dict: &BTreeMap<String, String>) -> Result<Coordinates, Box<dyn std::error::Error>> {
    let raw_chrom = record.seqname();
    let chrom = match chrom_dict.get(raw_chrom) {
        Some(c) => c.clone(),
        None => bail!("No chromosome definition found for \"{raw_chrom}\"")
    };
    // inputs are 1-based, so shift the start coordinate
    let start = *record.start() - 1;
    let end = *record.end();
    Ok(Coordinates::new(chrom, start, end))
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::path::PathBuf;

    /// Handy wrapper for loading the previously calculated gene definitions that we hard-coded
    fn load_expected_gene_defs() -> BTreeMap<String, GeneDefinition> {
        let mut gene_defs: BTreeMap<String, GeneDefinition> = Default::default();

        // set up HLA-A
        let shift = 1;
        let mut a_def = GeneDefinition::new(
            "HLA-A".to_string(),
            Coordinates::new("chr6".to_string(), 29942532 - shift, 29945870),
            true
        );
        a_def.add_transcript_id("NM_002116.8".to_string(), None).unwrap();
        let a_exons = [
            // 0-based, exclusive, sorted; copied from RefSeq file
            Coordinates::new("chr6".to_string(), 29942532 - shift, 29942626),
            Coordinates::new("chr6".to_string(), 29942757 - shift, 29943026),
            Coordinates::new("chr6".to_string(), 29943268 - shift, 29943543),
            Coordinates::new("chr6".to_string(), 29944122 - shift, 29944397),
            Coordinates::new("chr6".to_string(), 29944500 - shift, 29944616),
            Coordinates::new("chr6".to_string(), 29945059 - shift, 29945091),
            Coordinates::new("chr6".to_string(), 29945234 - shift, 29945281),
            Coordinates::new("chr6".to_string(), 29945451 - shift, 29945870)
        ];
        for exon in a_exons.into_iter() {
            a_def.add_exon(exon).unwrap();
        }
        gene_defs.insert("HLA-A".to_string(), a_def);

        // set up HLA-B
        let mut b_def = GeneDefinition::new(
            "HLA-B".to_string(),
            Coordinates::new("chr6".to_string(), 31353875 - shift, 31357179),
            false
        );
        b_def.add_transcript_id("NM_005514.8".to_string(), None).unwrap();
        let b_exons = [
            // 0-based, exclusive, sorted; copied from RefSeq file
            Coordinates::new("chr6".to_string(), 31353875 - shift, 31354296),
            Coordinates::new("chr6".to_string(), 31354479 - shift, 31354526),
            Coordinates::new("chr6".to_string(), 31354633 - shift, 31354665),
            Coordinates::new("chr6".to_string(), 31355107 - shift, 31355223),
            Coordinates::new("chr6".to_string(), 31355317 - shift, 31355592),
            Coordinates::new("chr6".to_string(), 31356167 - shift, 31356442),
            Coordinates::new("chr6".to_string(), 31356688 - shift, 31356957),
            Coordinates::new("chr6".to_string(), 31357086 - shift, 31357179)
        ];
        // these have to be inserted in reverse
        for exon in b_exons.into_iter().rev() {
            b_def.add_exon(exon).unwrap();
        }
        gene_defs.insert("HLA-B".to_string(), b_def);

        gene_defs
    }

    #[test]
    fn test_refseq_file_parsing() {
        // pull from a local smaller copy, prevents us from spamming RefSeq and having to parse the whole file in a test
        let filename = PathBuf::from("./test_data/refseq_faux/refseq_small.gff.gz");
        let mut targets: BTreeSet<String> = Default::default();
        targets.insert("HLA-A".to_string());
        targets.insert("HLA-B".to_string());

        let loaded_db = GeneCollection::load_refseq_file(&filename, Some(&targets)).unwrap();
        
        // check version parsing
        assert_eq!(loaded_db.version(), "NCBI RefSeq GCF_000001405.40-RS_2024_08");

        // check high level containment
        let gene_dict = loaded_db.gene_dict();
        assert_eq!(gene_dict.len(), targets.len());
        for target in targets.iter() {
            // check that only our targets are loaded
            assert!(gene_dict.contains_key(target));
        }

        let expected_results = load_expected_gene_defs();
        for (gene, exp_gene_def) in expected_results.iter() {
            println!("Checking {gene}...");

            let gene_def = gene_dict.get(gene).unwrap();
            assert_eq!(gene_def.gene_name(), gene);
            assert_eq!(gene_def.coordinates(), exp_gene_def.coordinates());
            assert_eq!(gene_def.is_forward_strand(), exp_gene_def.is_forward_strand());
            assert_eq!(gene_def.transcript_id(), exp_gene_def.transcript_id());
            assert_eq!(gene_def.exons(), exp_gene_def.exons());

            println!("{gene} matches.")
        }
    }

    #[test]
    #[ignore]
    fn test_refseq_url() {
        // pull from a local smaller copy, prevents us from spamming RefSeq and having to parse the whole file in a test
        let fixed_url = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/GCF_000001405.40-RS_2024_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz".to_string();
        let mut targets: BTreeSet<String> = Default::default();
        targets.insert("HLA-A".to_string());
        targets.insert("HLA-B".to_string());

        let loaded_db = GeneCollection::load_refseq_url(Some(fixed_url), Some(&targets)).unwrap();
        
        // check version parsing
        assert_eq!(loaded_db.version(), "NCBI RefSeq GCF_000001405.40-RS_2024_08");

        // check high level containment
        let gene_dict = loaded_db.gene_dict();
        assert_eq!(gene_dict.len(), targets.len());
        for target in targets.iter() {
            // check that only our targets are loaded
            assert!(gene_dict.contains_key(target));
        }

        let expected_results = load_expected_gene_defs();
        for (gene, exp_gene_def) in expected_results.iter() {
            println!("Checking {gene}...");

            let gene_def = gene_dict.get(gene).unwrap();
            assert_eq!(gene_def.gene_name(), gene);
            assert_eq!(gene_def.coordinates(), exp_gene_def.coordinates());
            assert_eq!(gene_def.is_forward_strand(), exp_gene_def.is_forward_strand());
            assert_eq!(gene_def.transcript_id(), exp_gene_def.transcript_id());
            assert_eq!(gene_def.exons(), exp_gene_def.exons());

            println!("{gene} matches.")
        }
    }
}