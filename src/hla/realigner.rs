
use log::{debug, trace, warn};
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use simple_error::{bail, SimpleError};
use std::collections::BTreeMap;
use std::ops::Deref;
use std::path::Path;
use tempfile::NamedTempFile;

use crate::data_types::coordinates::Coordinates;
use crate::data_types::database::PgxDatabase;
use crate::data_types::mapping::MappingStats;
use crate::data_types::pgx_diplotypes::PgxMappingDetails;
use crate::hla::debug::ReadMappingStats;
use crate::hla::mapping::HlaMappingStats;
use crate::util::file_io::save_fasta;
use crate::util::homopolymers::{hpc_bytes, hpc_pos};
use crate::util::mapping::{select_best_mapping, standard_hifi_aligner};
use crate::util::sequence::reverse_complement;

/// This contains the functionality responsible for comparing and then remapping any sequences in our dataset.
pub struct HlaRealigner<'a> {
    /// The source database
    database: &'a PgxDatabase,
    /// A temporary file containing all the database sequences
    tmp_db_fasta: NamedTempFile,
    /// The aligner to the database fasta file
    db_aligner: minimap2::Aligner<minimap2::Built>,
    /// Collection of aligners for the specific reference gene regions
    gene_aligners: HashMap<String, minimap2::Aligner<minimap2::Built>>,
    /// Collection of the sequences that are mapped inside `gene_aligners`
    gene_ref_sequence: HashMap<String, Vec<u8>>
}

impl<'a> HlaRealigner<'a> {
    /// Creates a new HLA realigner from a database
    /// # Arguments
    /// * `gene_list` - the genes we care about in the datbase
    /// * `database` - pre-loaded database containing all the HLA alleles
    /// * `reference_genome` - the input reference genome
    pub fn new(gene_list: &[String], database: &'a PgxDatabase, reference_genome: &ReferenceGenome) -> Result<Self, Box<dyn std::error::Error>> {
        // create a FASTA file that contains all the sequences from our database, which we will ultimately map against
        let tmp_db_fasta = tempfile::Builder::new()
            .prefix("starphase")
            .suffix(".fa")
            .tempfile()?;
        let tmp_db_fasta_fn = match tmp_db_fasta.path().as_os_str().to_str() {
            Some(filename) => filename,
            None => {
                bail!("Failed to generate NamedTempFile");
            }
        };

        debug!("Creating temporary FASTA file: {tmp_db_fasta_fn:?}");
        create_hla_fasta(gene_list, database, tmp_db_fasta.path())?;

        // prep our alignment to the full database set
        let db_aligner = standard_hifi_aligner()
            .with_index(tmp_db_fasta.path(), None)?;
        
        // create a gene aligner for each
        let mut gene_aligners: HashMap<String, minimap2::Aligner<minimap2::Built>> = Default::default();
        let mut gene_ref_sequence: HashMap<String, Vec<u8>> = Default::default();
        for gene in gene_list.iter() {
            // get the gene aligner, or make it if it does not exist
            let gene_def = match database.hla_config().gene_definition(gene) {
                Some(gd) => gd,
                None => bail!("Gene definition for {gene} not found.")
            };
            let region = gene_def.coordinates();

            // build a reference DNA aligner
            let buffer = 100;
            let reference_coordinates = Coordinates::new(region.chrom().to_string(), region.start() - buffer, region.end() + buffer);
            debug!("Buffered coordinates for {gene}: {reference_coordinates:?}");
            let reference_sequence = reference_genome.get_slice(reference_coordinates.chrom(), reference_coordinates.start() as usize, reference_coordinates.end() as usize);
            let g_aligner = standard_hifi_aligner()
                .with_seq(reference_sequence)?;
            gene_aligners.insert(gene.clone(), g_aligner);
            gene_ref_sequence.insert(gene.clone(), reference_sequence.to_vec());
        };

        Ok(Self {
            database,
            tmp_db_fasta,
            db_aligner,
            gene_aligners,
            gene_ref_sequence
        })
    }

    /// Given a record from an aligned BAM file, this function will re-align it against the HLA database.
    /// If it finds a mapping, it will also align to hg38 for the purpose of identifying offsets.
    /// It returns a re-aligned record containing the DNA, cDNA, and offsets for consensus.
    /// # Arguments 
    /// * `record` - the original aligned read
    pub fn realign_record(&self, record: &rust_htslib::bam::Record) -> Result<RealignmentResult, Box<dyn std::error::Error>> {
        // get the qname for debug output
        let qname = std::str::from_utf8(record.qname()).unwrap_or("UTF8_decode_error");

        // shared settings for mappings - we only need cigar and md for debugging
        let output_cigar: bool = true;
        let output_md: bool = true;
        let max_frag_len: Option<usize> = None;
        let extra_flags: Option<Vec<u64>> = None;

        // align the read to our database sequence
        let read_bytes = record.seq().as_bytes();
        let read_len = read_bytes.len();
        let d_mappings = self.db_aligner.map(
            &read_bytes,
            output_cigar, output_md, max_frag_len, extra_flags.as_deref(), None
        )?;

        // TODO: can we use the shared select_best_mapping function? trick here is we have some cut-offs encoded
        // find the single best mapping for the read
        let mut best_stats = MappingStats::new(read_len, read_len, 0);
        let mut best_mapping = None;
        let penalize_unmapped = false; // we want the lowest ED score as long as mappability is high
        for m in d_mappings.iter() {
            // scoring is based on the lowest edit distance, including unmapped
            let nm = m.alignment.as_ref().unwrap().nm as usize;

            // unmapped is going to be based on the target, not the query
            let unmapped = (m.target_len - (m.target_end - m.target_start)) as usize;
            let stats = MappingStats::new(m.target_len as usize, nm, unmapped);
            trace!("\t{stats:?}");

            // require a percentage to map, a low ED in what is mapped, and that the ED score is lower than current best
            let max_unmapped_frac = 0.5; // TODO: parameter? tuning?
            let max_ed_frac = 0.03; // TODO: parameter? tuning? cli_settings.max_error_rate <- currently 0.07, which is relative to GRCh38
            if stats.mapping_score().score() <= max_unmapped_frac &&
                stats.custom_score(penalize_unmapped).score() <= max_ed_frac &&
                stats.custom_score(penalize_unmapped).score() < best_stats.custom_score(penalize_unmapped).score() {
                // strand is variable here, but typically matches the gene strand; e.g. HLA-A will be Forward, HLA-B Reverse
                best_stats = stats;
                best_mapping = Some(m);
            }
        }

        // check if the mappings failed
        if best_mapping.is_none() {
            debug!("{qname} => None, no sufficiently good mappings to database, ignoring read");
            let mapping_stats = HlaMappingStats::from_mapping_stats(None, Some(best_stats));
            let ignored_details = PgxMappingDetails::new(
                qname.to_string(),
                "REFERENCE".to_string(),
                "REFERENCE".to_string(),
                mapping_stats,
                true
            );
            return Ok(RealignmentResult::failure(
                ignored_details
            ))
        }

        // we found a mapping, check if it's good or not
        let bm = best_mapping.unwrap();
        let target_name: String = bm.target_name.as_ref()
            .ok_or(SimpleError::new("Mapping failed to generate a target_name"))?
            .deref()
            .clone();

        // we have a target name also
        let best_def = self.database.hla_sequences().get(&target_name).unwrap();
        let best_gene = best_def.gene_name();
        let best_gene_def = self.database.hla_config().gene_definition(best_gene)
            .ok_or(SimpleError::new(format!("Could not find matching gene definition for {target_name}")))?;
        let best_star = best_def.star_allele().join(":");
        debug!("{qname} => {best_gene:?} {best_star:?} => {:?}, {}, {}", bm.strand, best_stats.score_string(), best_stats.custom_score_string(penalize_unmapped));
        
        // TODO: I think this happens with bad read mappings; we could theoretically handle it by rev-comping the read; it's rare though, so punt for now
        if bm.strand != minimap2::Strand::Forward {
            debug!("\t => Ignoring {:?} strand mapping", bm.strand);
            let mapping_stats = HlaMappingStats::from_mapping_stats(None, Some(best_stats));
            let ignored_details = PgxMappingDetails::new(
                qname.to_string(),
                target_name.clone(),
            format!("{best_gene}*{best_star}"),
                mapping_stats,
                true
            );
            return Ok(RealignmentResult::failure(
                ignored_details
            ))
        }

        // after this, we assume everything is forward
        assert!(bm.strand == minimap2::Strand::Forward);

        // add this best mapping to our read_debug stats
        let mut all_mapping_stats = ReadMappingStats::new();
        all_mapping_stats.add_mapping(target_name.clone(), None, Some(bm))?;
        all_mapping_stats.set_best_match(target_name.clone(), best_star.clone());
        
        // also save this record as matching something
        let mapping_stats = HlaMappingStats::from_mapping_stats(None, Some(best_stats));
        let details = PgxMappingDetails::new(
            qname.to_string(),
            target_name.clone(),
            format!("{best_gene}*{best_star}"),
            mapping_stats,
            false
        );

        // get the gene aligner, or make it if it does not exist
        let gene_aligner = self.gene_aligners.get(best_gene)
            .ok_or(SimpleError::new(format!("Failed to get gene_aligner for {best_gene}")))?;
        let reference_sequence = self.gene_ref_sequence.get(best_gene)
            .ok_or(SimpleError::new(format!("Failed to get reference_sequence for {best_gene}")))?;

        // pull out the segment that mapped best to the database allele
        let db_segment_start = bm.query_start as usize;
        let db_segment_end = bm.query_end as usize;
        debug!("\t{db_segment_start}..{db_segment_end} aligned to {best_star}:{}..{}", bm.target_start, bm.target_end);
        
        // add a buffer around the best mapping to see if we can match hg38 also
        // TODO: figure out ideal buffer size; seems like maybe we should make this quite a bit larger
        let buffer = 1000;
        let buffered_segment_start = db_segment_start.saturating_sub(buffer); // subtract a buffer, but not before index 0
        let buffered_segment_end = (db_segment_end+buffer).min(read_bytes.len()); // add the buffer, but not longer than the readlen

        // map just that segment to the hg38 gene
        let gene_mappings = gene_aligner.map(
            &read_bytes[buffered_segment_start..buffered_segment_end],
            output_cigar, output_md, max_frag_len, extra_flags.as_deref(), None
        )?;

        // we want the unmapped bases relative to the target (reference), and we want to penalize on unmapped also
        let unmapped_from_target = true;
        let penalize_unmapped = true;
        let (best_mapping, _best_stats) = select_best_mapping(
            &gene_mappings, unmapped_from_target, penalize_unmapped, None
        );

        let realigned_record = if let Some(ref_mapping) = best_mapping {
            // these are all read mappings (reference = forward) compared to reference; theoretically, all are Forward
            if ref_mapping.strand == minimap2::Strand::Forward {
                /*
                * At this point, we have our best read mapping against the database *DNA* sequences in `bm`; anything cDNA-only cannot be aligned to.
                * We also have that same segment cut out of the read and aligned to GRCh38 in `best_mapping`.
                * In the RealignedHlaRecord::new, only `best_mapping` is used, which means all offsets are implicitly connected to GRCh38.
                * With high diversity alleles, this means our offsets are often really crappy.
                *
                * Instead, we should calculate our offsets based on the original `bm` mapping, BUT
                *   we cannot just ignore GRCh38 because the database alleles are not equal length.
                * IF the mapping to `bm.target_start==0`, then we should re-align to GRCh38 and use min(bm.query_start, best_mapping.query_start) for the segment. Offset will be based on best_mapping.target_start
                * IF the mapping to `bm.target_start>0`, then we have an offset that needs to get calculated on both. If we map the allele to GRCh38, we get an initial offset. We then add that to `bm.target_start` to get total offset.
                */

                // these are the offsets into our read sequence that aligned to hg38
                // adjusted because we spliced a segment out of the read that started at `buffered_segment_start`
                let adjusted_hg38_segment_start = buffered_segment_start + ref_mapping.query_start as usize;
                let adjusted_hg38_segment_end = buffered_segment_start + ref_mapping.query_end as usize;

                debug!("\t{adjusted_hg38_segment_start}..{adjusted_hg38_segment_end} aligned to hg38:{}..{}", ref_mapping.target_start, ref_mapping.target_end);

                // best start/end extends as far as possible in both directions
                let optimal_segment_start = db_segment_start.min(adjusted_hg38_segment_start);
                let optimal_segment_end = db_segment_end.max(adjusted_hg38_segment_end);

                // now figure out the offset
                let (dna_offset, hpc_offset) = if adjusted_hg38_segment_start < db_segment_start {
                    // if our reference mapping worked, then it will have an earlier start; and that means we can just use the target_start from that mapping
                    let d = ref_mapping.target_start as usize;
                    let h = hpc_pos(reference_sequence, d);
                    (d, h)
                } else {
                    // ELSE our reference mappings started after the mapping to the database allele; we need to estimate our start/stop using two mappings
                    // our DB segment started partway through the allele, we need to align the DB allele and add the offsets together
                    let best_def_seq = best_def.dna_sequence().unwrap();
                    let revcomp_bytes;

                    // the database sequence might be rev-comp relative to hg38
                    let best_def_bytes = if best_gene_def.is_forward_strand() {
                        best_def_seq.as_bytes()
                    } else {
                        revcomp_bytes = reverse_complement(best_def_seq.as_bytes())?;
                        &revcomp_bytes
                    };

                    // now map it
                    let best_def_mappings = gene_aligner.map(
                        best_def_bytes,
                        output_cigar, output_md, max_frag_len, extra_flags.as_deref(), None
                    )?;

                    // make sure we skip anything on the reverse strand
                    let best_def_mappings = best_def_mappings.into_iter()
                        .filter(|m| m.strand == minimap2::Strand::Forward)
                        .collect::<Vec<minimap2::Mapping>>();

                    // we are mapping a database definition to hg38; unmapped comes from the query, and we do want to penalize
                    let unmapped_from_target = false;
                    let penalize_unmapped = true;
                    let (allele_mapping, _best_stats) = select_best_mapping(
                        &best_def_mappings, unmapped_from_target, penalize_unmapped, None
                    );

                    if let Some(am) = allele_mapping {
                        debug!("\t{best_star}:{}..{} aligned to hg38:{}..{}", am.query_start, am.query_end, am.target_start, am.target_end);
                        let added_offset = (am.target_start - am.query_start).max(0);
                        debug!("\t\tadded_offset={added_offset}");

                        // offset from aligning the allele to the reference PLUS
                        //      offset from aligning the read to the allele
                        // EQUALS total theoretical offset of the read into the consensus
                        let d = (added_offset + bm.target_start) as usize;
                        // for HPC version, we have to add the HPC offset of both also
                        let h = hpc_pos(reference_sequence, added_offset as usize) + hpc_pos(best_def_seq.as_bytes(), bm.target_start as usize);
                        (d, h)
                    } else {
                        warn!("Failed to map allele {best_gene:?} {best_star:?} to reference, ignoring offset adjustment");
                        let d = ref_mapping.target_start as usize;
                        let h = hpc_pos(reference_sequence, d);
                        (d, h)
                    }
                };

                // build the segment from what is selected
                debug!("\tSelecting {optimal_segment_start}..{optimal_segment_end} as segment with dna_offset={dna_offset}");
                Some(RealignedHlaRecord::new(
                    record, optimal_segment_start..optimal_segment_end, dna_offset, hpc_offset
                )?)
            } else {
                // TODO: we can probably handle this by reversing the read and re-mapping it, but tricky and rare
                // This seems to happen when we find a segment that is not actually mapped to the reference genome. For example,
                //   a soft-clipped part of the read contains the segment.
                warn!("Best remapping of {qname} was to Reverse strand, ignoring.");
                None
            }
        } else {
            warn!("Remapping of {qname} to reference failed, ignoring.");
            None
        };

        Ok(RealignmentResult::new(
            best_gene.to_string(),
            all_mapping_stats,
            details,
            realigned_record
        ))
    }

    // getters
    pub fn tmp_db_fasta(&self) -> &NamedTempFile {
        &self.tmp_db_fasta
    }
}

/// Wrapper for a realignment result, containing the realignment stats and the realigned record.
#[derive(Debug, Default)]
pub struct RealignmentResult {
    /// The gene this was remapped into
    gene_name: String,
    /// Contains the mapping stats from the read to the best match
    read_mapping_stats: ReadMappingStats,
    /// Contains the best match details
    mapping_details: PgxMappingDetails,
    /// Optional, only present if successful
    realigned_record: Option<RealignedHlaRecord>
}

impl RealignmentResult {
    /// Basic constructor
    pub fn new(
        gene_name: String, read_mapping_stats: ReadMappingStats, mapping_details: PgxMappingDetails, 
        realigned_record: Option<RealignedHlaRecord>
    ) -> Self {
        Self {
            gene_name,
            read_mapping_stats,
            mapping_details,
            realigned_record
        }
    }

    /// Generic constructor for failure case, but with some stats
    pub fn failure(mapping_details: PgxMappingDetails) -> Self {
        Self {
            mapping_details,
            ..Default::default()
        }
    }

    /// Indicates if a realigned record is present
    pub fn is_realigned(&self) -> bool {
        self.realigned_record.is_some()
    }

    // getters
    pub fn gene_name(&self) -> &str {
        &self.gene_name
    }

    pub fn read_mapping_stats(&self) -> &ReadMappingStats {
        &self.read_mapping_stats
    }

    pub fn mapping_details(&self) -> &PgxMappingDetails {
        &self.mapping_details
    }
    
    pub fn realigned_record(&self) -> Option<&RealignedHlaRecord> {
        self.realigned_record.as_ref()
    }
}


/// This is a wrapper return value from re-aligning a record to a specific genomic context.
/// Main purpose is to encapsulate the derived data we care about.
#[derive(Debug, Default, PartialEq)]
pub struct RealignedHlaRecord {
    /// Contains the record after re-aligning to a specific region; this will have the full position information
    realigned_record: rust_htslib::bam::Record,
    /// Contains the sequence that was successfully mapped to the target region
    dna_sequence: Vec<u8>,
    /// Contains the offset into the target region where this sequence starts mapping
    dna_offset: usize,
    /// Contains the HPC version of `dna_sequence`
    hpc_sequence: Vec<u8>,
    /// Contains the HPC offset, which is calculated by looking at the sequence predicted before the observation
    hpc_offset: usize
}

impl RealignedHlaRecord {
    /// Creates a new re-aligned HLA record
    /// # Arguments
    /// * `record` - the original alignment record
    /// * `dna_segment` - the segment of the record that needs to be spliced out
    /// * `dna_offset` - the offset assigned to that segment
    /// * `hpc_offset` - the offset for the HPC version of the sequence
    pub fn new(
        // gene_name: &str,
        record: &rust_htslib::bam::Record,
        dna_segment: std::ops::Range<usize>,
        dna_offset: usize,
        hpc_offset: usize
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let realigned_record = record.clone();
        let read_bytes = record.seq().as_bytes();
        let dna_sequence = read_bytes[dna_segment].to_vec();
        let hpc_sequence = hpc_bytes(&dna_sequence);
        let qname = record.qname();
        
        trace!("Created realigned record:");
        trace!("\tqname:      {}", std::str::from_utf8(qname)?);
        trace!("\tdna_len:    {}", dna_sequence.len());
        trace!("\tdna_offset: {dna_offset}");
        trace!("\thpc_len:    {}", hpc_sequence.len());
        trace!("\thpc_offset: {hpc_offset}");

        Ok(Self {
            realigned_record,
            dna_sequence,
            dna_offset,
            hpc_sequence,
            hpc_offset
        })
    }

    // getters
    pub fn realigned_record(&self) -> &rust_htslib::bam::Record {
        &self.realigned_record
    }

    pub fn dna_sequence(&self) -> &[u8] {
        &self.dna_sequence
    }

    pub fn dna_offset(&self) -> usize {
        self.dna_offset
    }

    pub fn hpc_sequence(&self) -> &[u8] {
        &self.hpc_sequence
    }

    pub fn hpc_offset(&self) -> usize {
        self.hpc_offset
    }
}

/// Wrapper function that creates the mapping from our database to a basic BTreeMap and then saves it using save_fasta.
/// In this instance, we want everything on hg38 orientation, so we flip any reverse genes to handle that.
/// # Arguments
/// * `gene_list` - the list of genes we care about, will ignore all others
/// * `database` - the full database, which we pull HLA sequences from for mapping
/// * `filename` - the Path to the file we want to write
fn create_hla_fasta(
    gene_list: &[String], database: &PgxDatabase, filename: &Path
) -> Result<(), Box<dyn std::error::Error>> {
    // collect into set for easier lookups, also need to deref it for "contains" query
    let gene_set: HashSet<&str> = gene_list.iter()
        .map(|s| s.as_str())
        .collect();

    // add the definition if it's in our gene list AND has DNA sequence
    let mut seq_dict: BTreeMap<String , String> = Default::default();
    for (hla_allele, hla_def) in database.hla_sequences().iter() {
        let gene_name = hla_def.gene_name();
        if gene_set.contains(gene_name) {
            let gene_def = database.hla_config().gene_definition(gene_name)
                .ok_or(SimpleError::new(format!("Expected {gene_name} in gene_definitions")))?;
            if let Some(dna_seq) = hla_def.dna_sequence() {
                let fw_dna_seq = if gene_def.is_forward_strand() {
                    dna_seq.to_string()
                } else {
                    String::from_utf8(reverse_complement(dna_seq.as_bytes())?)?
                };
                assert!(!seq_dict.contains_key(hla_allele));
                seq_dict.insert(hla_allele.to_string(), fw_dna_seq);
            }
        }
    }

    // now just save it to the FASTA
    save_fasta(&seq_dict, filename)
}

// TODO: tests - mostly on RealignedHlaRecord, maybe we can do some on HlaRealigner?
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_realigned_record() {
        // make a basic record to test on
        let qname = "test_qname";
        let full_sequence = "AACCGGTTAACCGGTTAACCGGTT";
        let qual = vec![255; full_sequence.len()];
        let mut record = rust_htslib::bam::Record::new();
        record.set(qname.as_bytes(), None, full_sequence.as_bytes(), &qual);

        // now fill in the details
        let dna_segment = 4..10;
        let dna_offset = 4;
        let hpc_offset = 2;

        let realigned_record = RealignedHlaRecord::new(&record, dna_segment, dna_offset, hpc_offset).unwrap();
        let expected_realigned_record = RealignedHlaRecord {
            realigned_record: record,
            dna_sequence: b"GGTTAA".to_vec(),
            dna_offset,
            hpc_sequence: b"GTA".to_vec(),
            hpc_offset,
        };
        assert_eq!(realigned_record, expected_realigned_record);
    }

    // possible TODO: do we want to unit test the HLA re-aligner? it will be complicated to set up and end-to-end tests should suffice
}