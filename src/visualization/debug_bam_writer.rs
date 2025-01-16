
use log::{debug, error, info, warn};
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use simple_error::bail;
use std::path::PathBuf;
use std::collections::{BTreeMap, HashMap};

use crate::data_types::coordinates::Coordinates;
use crate::util::mapping::{select_best_mapping, standard_hifi_aligner};
use crate::util::sequence::reverse_complement;

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

/// Clears an existing record of mapping information and adds new tags if provided
pub fn clear_record(original_record: &rust_htslib::bam::Record, tags: &BTreeMap<String, String>) -> Result<rust_htslib::bam::Record, Box<dyn std::error::Error>> {
    let mut record = original_record.clone(); // this *should* copy any relevant tags that already exist
    let qname = original_record.qname();
    let cigar = None;
    let sequence = original_record.seq().as_bytes();
    let qual = original_record.qual();

    // this should just clear the cigar string
    record.set(qname, cigar, &sequence, qual);

    // set the chromosome and position
    record.set_tid(-1);
    record.set_pos(-1);
    record.set_mapq(20); // do we care about adjusting this?
    for (key, value) in tags.iter() {
        match record.remove_aux(key.as_bytes()) {
            // either it was remove or was not there to begin with
            Ok(()) |
            Err(rust_htslib::errors::Error::BamAuxTagNotFound) => {},
            // some other crazy error
            Err(e) => { return Err(Box::new(e)); }
        }
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
        let dna_aligner = standard_hifi_aligner()
            .with_seq(region_sequence)?;
        
        /*
        // TODO: tried to figure out how to prevent alignments from spanning our intentional gaps between haplotypes, but failed
        //       lets revisit this in the future
        let dna_aligner = Aligner {
                mapopt: minimap2::MapOpt {
                    max_gap: 500,
                    ..Aligner::builder().map_hifi().mapopt
                },
                ..Aligner::builder().map_hifi()
            }.with_cigar()
            .with_seq(region_sequence)?;
        debug!("mapopt: {:?}", dna_aligner.mapopt);
        */
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
            let qname_str = std::str::from_utf8(qname).unwrap_or("UTF8_FAILURE");
            let mut sequence = umr.seq().as_bytes();
            let mut qual = umr.qual();
            let reverse_qual: Vec<u8>; // only if needed

            if sequence.is_empty() {
                let qname_str = std::str::from_utf8(qname).unwrap_or("UTF8_ERROR");
                error!("Failed to map {qname_str}, empty sequence; ignoring.");
                continue;
            }

            // first, map the sequence
            let mappings = dna_aligner.map(
                &sequence,
                output_cigar, output_md, max_frag_len, extra_flags, None
            )?;
            
            // pick the mapping with the largest match_len of those created
            if mappings.is_empty() {
                // technically can happen, but it shouldn't; regardless we don't want to panic here
                debug!("Failed to map {qname_str} to {target_region}");
                continue;
            }

            // find the best mapping using our shared sub-routine
            let (best_mapping, _best_stats) = select_best_mapping(&mappings, false, false, None);
            let core_mapping = best_mapping.unwrap_or_else(|| {
                warn!("No good re-mapping found for {qname_str} to {target_region}, selecting first.");
                &mappings[0]
            });

            /*
            TODO: explore other scoring systems 
                we used to do (match_len - nm) to select the best mapping, but this can generate errors with similar sequences
                another option is to go with the first returned result, which is what minimap2 would select; this also
                has issues where you get a lot of mis-mapped reads due to minimap2 scoring (i.e. length tends to beat mismatches)
                our current approach is select the mapping with lowest ratio of edits / match, but this can technically also generate sub-par results
                we could also just put all the mappings out and label the secondaries / supplements accordingly, but that's some non-trivial effort
                overall, select_best_mapping is a good balance between minimap2 scoring and our sensitivity for now
             */
            
            // these are almost always on the forward strand
            let mut set_reverse = false;
            let mut set_forward = false;
            if core_mapping.strand == minimap2::Strand::Reverse {
                if umr.is_reverse() {
                    // this is a double reverse, so set to forward; but we still need to rev-comp and rev-qual
                    set_forward = true;
                } else {
                    // it is not currently set to reverse, so set it
                    set_reverse = true;
                }

                // looks like we need to revcomp the sequence
                sequence = reverse_complement(&sequence)?;

                // should also swap qual I believe
                reverse_qual = qual.iter().cloned().rev().collect();
                qual = &reverse_qual;
            }
            // assert!(core_mapping.strand == minimap2::Strand::Forward);

            // convert the cigar into BAM ready format
            let cigar = convert_mapping_to_cigar(core_mapping, None, None);
            
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

            // check for flag adjustments
            if set_reverse {
                record.set_reverse();
            } else if set_forward {
                record.unset_reverse();
            }
    
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

/// Simple wrapper that converts a minimap2 mapping to a CigarString for the rust_htslib API
/// # Arguments
/// * `mapping` - the minimap2 Mapping we want to convert
/// * `clipped_start_extra` - additional bases that were pre-clipped at the start
/// * `clipped_end_extra` - additional bases that were pre-clipped at the end
pub fn convert_mapping_to_cigar(mapping: &minimap2::Mapping, clipped_start_extra: Option<usize>, clipped_end_extra: Option<usize>) -> CigarString {
    let mut cigar = CigarString(
        mapping.alignment.as_ref().unwrap()
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

    let start_extra = clipped_start_extra.unwrap_or_default() as i32;
    let end_extra = clipped_end_extra.unwrap_or_default() as i32;

    // if it's reverse, we need to swap the clippings here
    let is_reverse = mapping.strand == minimap2::Strand::Reverse;
    let start_delta = mapping.query_start + start_extra;
    let end_delta = mapping.query_len.unwrap().get() - mapping.query_end + end_extra;
    let (start_delta, end_delta) = if is_reverse { (end_delta, start_delta) } else { (start_delta, end_delta) };

    if start_delta > 0 {
        cigar.insert(0, Cigar::SoftClip(start_delta as u32));
    }
    
    if end_delta > 0 {
        cigar.push(Cigar::SoftClip(end_delta as u32));
    }
    cigar
}