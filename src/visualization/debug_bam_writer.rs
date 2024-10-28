
use log::{debug, info, warn};
use minimap2::Aligner;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use simple_error::bail;
use std::path::PathBuf;
use std::collections::{BTreeMap, HashMap};

use crate::data_types::coordinates::Coordinates;

/// Creates an unmapped record with the minimal content we might need
pub fn unmapped_record(qname: &str, sequence: &str, tags: &BTreeMap<String, String>) -> Result<rust_htslib::bam::Record, Box<dyn std::error::Error>> {
    // we need to create a new record that is unmapped
    let mut record = rust_htslib::bam::Record::new();
    let qual = vec![255_u8; sequence.len()]; // this method has no provided quals
    let cigar = None; // cigar comes after we align
    
    // set the basics for the record
    record.set(
        qname.as_bytes(), 
        cigar,
        sequence.as_bytes(), 
        &qual
    );

    // set the chromosome and position
    record.set_tid(-1);
    record.set_pos(-1);
    record.set_mapq(20); // do we care about adjusting this?
    for (key, value) in tags.iter() {
        record.push_aux(key.as_bytes(), rust_htslib::bam::record::Aux::String(value))?;
    }

    Ok(record)
}

/// This is a debug BAM writer, all records are kept in memory until `write_all_records(...)` is called.
/// This is not meant to write a full BAM, but just a few small sections.
pub struct DebugBamWriter<'a> {
    /// Output path
    out_fn: PathBuf,
    /// Reference genome
    reference_genome: &'a ReferenceGenome,
    /// The actual record writer
    writer: Option<rust_htslib::bam::Writer>,
    /// Contains all records that will eventually get written to the BAM
    records: HashMap<String, Vec<rust_htslib::bam::Record>>
}

impl<'a> DebugBamWriter<'a> {
    /// Creates a new writer by building the header and prepping the writer for later.
    /// # Arguments
    /// * `out_fn` - the output filename to save everything to
    /// * `reference_genome` - need to build out the header
    pub fn new(out_fn: PathBuf, reference_genome: &'a ReferenceGenome) -> Result<Self, rust_htslib::errors::Error> {
        // create a default header
        let mut header = rust_htslib::bam::Header::new();
        let mut header_record = rust_htslib::bam::header::HeaderRecord::new(b"HD");
        header_record.push_tag(b"VN", "1.5");
        header_record.push_tag(b"SO", "coordinate");
        header.push_record(&header_record);

        // for each chromosome in our output, we need to add it to the header
        for chromosome in reference_genome.contig_keys().iter() {
            // @SQ	SN:chr22	LN:50818468
            let mut header_record = rust_htslib::bam::header::HeaderRecord::new(b"SQ");
            let target_length = reference_genome.get_full_chromosome(chromosome).len();
            header_record.push_tag(b"SN", chromosome);
            header_record.push_tag(b"LN", target_length);
            header.push_record(&header_record);
        }

        // finally, init the writer for later
        let writer = Some(rust_htslib::bam::Writer::from_path(&out_fn, &header, rust_htslib::bam::Format::Bam)?);
        Ok(Self {
            out_fn, 
            reference_genome,
            writer,
            records: Default::default()
        })
    }

    /// For a given collection of unmapped records, this will attempt to map them to a specified region and add them to our collection.
    /// # Arguments
    /// * `unmapped_records` - the record information to map and save as a BAM record
    /// * `target_region` - the region targeted for alignment
    pub fn map_records_to_region(&mut self, unmapped_records: &[rust_htslib::bam::Record], target_region: &Coordinates) -> Result<(), Box<dyn std::error::Error>> {
        if self.writer.is_none() {
            bail!("This writer has already written everything!");
        }
        debug!("Generating records for {target_region}...");
        
        // build a mapper for the region
        let region_sequence = self.reference_genome.get_slice(target_region.chrom(), target_region.start() as usize, target_region.end() as usize);
        let dna_aligner: Aligner = Aligner::builder()
            .map_hifi()
            .with_cigar()
            .with_seq(region_sequence)?;

        // we only need cigar and md for debugging
        // other settings for mapping
        let output_cigar: bool = true;
        let output_md: bool = true;
        let max_frag_len: Option<usize> = None;
        let extra_flags = None;

        // get the tid which will be shared for all records in this region
        let tid = match self.writer.as_ref().unwrap().header().tid(target_region.chrom().as_bytes()) {
            Some(t) => t as i32,
            None => bail!("Could not find chromosome \"{}\" in reference genome.", target_region.chrom())
        };

        // get the vec we populate with records on this chromosome
        let record_vec = self.records.entry(target_region.chrom().to_string()).or_default();

        for umr in unmapped_records.iter() {
            // we will need to pull all of these out for setting later
            let qname = umr.qname();
            let sequence = umr.seq().as_bytes();
            let qual = umr.qual();

            // first, map the sequence
            let mappings = dna_aligner.map(
                &sequence,
                output_cigar, output_md, max_frag_len, extra_flags.clone()
            )?;
            
            // pick the mapping with the largest match_len of those created
            if mappings.is_empty() {
                // technically can happen, but it shouldn't; regardless we don't want to panic here
                debug!("Failed to map unmapped record: {:?}", std::str::from_utf8(umr.qname()));
                continue;
            }

            // find the max based on (match_length - edit_distance)
            // in the event of a tie, pick the earlier position
            let core_mapping = mappings.iter()
                .max_by_key(|&m| (m.match_len - m.alignment.as_ref().unwrap().nm, std::cmp::Reverse(m.target_start)))
                .unwrap();
            
            let mut cigar = CigarString(
                core_mapping.alignment.as_ref().unwrap()
                    .cigar.as_ref().unwrap().iter()
                    .map(|&(l, t)| {
                        match t {
                            0 => Cigar::Match(l),
                            1 => Cigar::Ins(l),
                            2 => Cigar::Del(l),
                            _ => panic!("unhandled cigar type: {t}")
                        }
                    })
                    .collect()
            );
    
            let start_delta = core_mapping.query_start;
            if start_delta > 0 {
                cigar.insert(0, Cigar::SoftClip(start_delta as u32));
            }
            
            let end_delta = sequence.len() - core_mapping.query_end as usize;
            if end_delta > 0 {
                cigar.push(Cigar::SoftClip(end_delta as u32));
            }
            
            // these should come from the mapping
            let cigar = Some(&cigar);
            let pos = target_region.start() as i64 + core_mapping.target_start as i64;
    
            // initialize the mapped record by just cloning the other record and then overwrite the particulars
            let mut record = umr.clone();
            record.set(
                qname, 
                cigar, 
                &sequence, 
                qual
            );
    
            // set the chromosome and position
            record.set_tid(tid);
            record.set_pos(pos);
            record.set_mapq(20); // do we care about adjusting this?

            // push it into our "to-save" list
            record_vec.push(record);
        }

        // done with all records, should be good!
        Ok(())
    }

    /// Writes all records out that we've given to this struct
    pub fn write_all_records(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        if let Some(writer) = self.writer.as_mut() {
            info!("Writing all records to {:?}...", self.out_fn);
            for chromosome in self.reference_genome.contig_keys().iter() {
                if let Some(record_vec) = self.records.get_mut(chromosome) {
                    // we have records to write on this chromosome
                    // first, we need to order by position
                    record_vec.sort_by_key(|a| a.pos());

                    // now send it to file
                    for record in record_vec.iter() {
                        writer.write(record)?;
                    }
                }
            }
        } else {
            bail!("This writer has already written everything!");
        }

        // delete our BAM file handle, flushing it and closing the file prior to indexing
        self.writer = None;

        // build the index
        info!("Building index...");
        let idx_type = rust_htslib::bam::index::Type::Bai;
        match rust_htslib::bam::index::build(
            &self.out_fn,
            None,
            idx_type,
            1
        ) {
            Ok(()) => {
                info!("Finished building index for {:?}", self.out_fn);
            },
            Err(e) => {
                warn!("Error while building index for {:?}: {}", self.out_fn, e);
                warn!("Continuing with other processing...");
            }
        };

        // finish up
        Ok(())
    }
}