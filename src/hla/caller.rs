
use bio::bio_types::genome::AbstractInterval;
use log::{debug, error, info, warn};
use rust_htslib::bam::Read;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use simple_error::bail;
use statrs::distribution::{Binomial, Continuous, Discrete, DiscreteCDF, Normal};
use std::collections::BTreeMap;
use std::path::PathBuf;
use waffle_con::cdwfa_config::{CdwfaConfig, CdwfaConfigBuilder, CdwfaConfigBuilderError};
use waffle_con::consensus::{Consensus, ConsensusDWFA};
use waffle_con::dual_consensus::{DualConsensus, DualConsensusDWFA};

use crate::cli::diplotype::DiplotypeSettings;
use crate::data_types::coordinates::Coordinates;
use crate::data_types::database::PgxDatabase;
use crate::data_types::pgx_diplotypes::{Diplotype, PgxGeneDetails, PgxMappingDetails};
use crate::hla::alleles::{HlaAlleleDefinition, NORMALIZING_HLA_GENES};
use crate::hla::debug::{HlaDebug, ReadMappingStats};
use crate::hla::mapping::HlaMappingStats;
use crate::hla::processed_match::HlaProcessedMatch;
use crate::hla::realigner::{HlaRealigner, RealignmentResult};
use crate::util::file_io::save_fasta;
use crate::util::mapping::{select_best_mapping, standard_hifi_aligner};
use crate::util::sequence::reverse_complement;
use crate::visualization::debug_bam_writer::{clear_record, convert_mapping_to_cigar, unmapped_record, DebugBamWriter};
use crate::visualization::igv_session_writer::{BUFFER_LEN, CONTIG_POSTFIX, IgvSessionWriter};

/// when a hap cannot be described for some reason
const UNKNOWN_HAP: &str = "UNKNOWN";
/// Intentionally absent haplotypes; e.g., hemizygous or absent haps
const NO_CALL_HAP: &str = ".";
/// Label when we expect a call, but had no data; i.e. it should *not* be absent, but it is
const NO_READS_HAP: &str = "NO_READS";

/// This is the original function to call for HLA diplotyping from a BAM file, which was built for HLA-A and HLA-B.
/// HLA-C also seems to work well, although it has not been deep tested.
/// # Arguments
/// * `gene_list` - the list of HLA genes to diplotype
/// * `database` - the pre-loaded database
/// * `bam_filenames` - list of BAM files containing reads to scan
/// * `reference_filename` - reference genome file path
/// * `cli_settings` - settings for diplotyping
/// # Errors
/// * if a gene is provided that we do not support
/// * if we cannot open or parse a BAM file correctly
pub fn diplotype_hla(
    gene_list: &[String], database: &PgxDatabase,
    bam_filenames: &[PathBuf], reference_genome: &ReferenceGenome,
    mut debug_bam_writer: Option<&mut DebugBamWriter>,
    cli_settings: &DiplotypeSettings
) -> Result<HashMap<String, PgxGeneDetails>, Box<dyn std::error::Error>> {
    // if we have disabled cDNA scoring AND the DNA requirement is NOT enabled; then we can cause errors later due to lack of comparator sequence
    if cli_settings.disable_cdna_scoring && !cli_settings.hla_require_dna {
        bail!("If cDNA scoring is disabled, require HLA DNA must be enabled");
    }

    // prep all the bam readers
    let mut bam_readers: Vec<rust_htslib::bam::IndexedReader> = vec![];
    for bam_fn in  bam_filenames.iter() {
        let mut b = rust_htslib::bam::IndexedReader::from_path(bam_fn)?;
        b.set_reference(reference_genome.filename())?;
        bam_readers.push(b);
    }

    //set up job configuration
    let mut ret: HashMap<String, PgxGeneDetails> = Default::default();
    let mut debug_stats = HlaDebug::new();
    
    for gene_name in gene_list.iter() {
        info!("Solving {gene_name}...");

        // get the coordinates for this gene
        let gene_definition = match database.hla_config().gene_definition(gene_name) {
            Some(c) => c,
            None => bail!("No gene definition for {gene_name}")
        };
        let gene_coordinates = gene_definition.coordinates();
        debug!("Gene coordinates: {gene_coordinates:?}");
        
        // build a reference DNA aligner
        let buffer = 100;
        let reference_coordinates = Coordinates::new(gene_coordinates.chrom().to_string(), gene_coordinates.start() - buffer, gene_coordinates.end() + buffer);
        debug!("Buffered coordinates: {reference_coordinates:?}");
        let reference_sequence = reference_genome.get_slice(reference_coordinates.chrom(), reference_coordinates.start() as usize, reference_coordinates.end() as usize);
        let dna_aligner = standard_hifi_aligner()
            .with_seq(reference_sequence)?;

        // we only need cigar and md for debugging
        // other settings for mapping
        let output_cigar: bool = true;
        let output_md: bool = true;
        let max_frag_len: Option<usize> = None;
        let extra_flags = None;

        // this is where we collect all the DNA segments; using BTreeMap for defined iteration order
        let mut read_segments: BTreeMap<String, String> = Default::default();
        let mut spliced_segments: BTreeMap<String, String> = Default::default();

        // store all the relevant mapping details
        let mut mapping_details: Vec<PgxMappingDetails> = Default::default();
        
        // iterate over each bam, and fetch the reads
        for (bam_index, bam) in bam_readers.iter_mut().enumerate() {
            match bam.fetch(reference_coordinates.fetch_definition()) {
                Ok(()) => {},
                Err(e) => {
                    let filename = &bam_filenames[bam_index];
                    warn!("Received error \"{e}\" while fetching {reference_coordinates} in {filename:?}, assuming no reads for region.");
                    continue;
                }
            };
            
            for read_entry in bam.records() {
                let mut read = read_entry.unwrap();
                
                //build out the cigar info
                read.cache_cigar();
                
                let qname: String = std::str::from_utf8(read.qname())?.to_string();
                let full_range = read.range();
                
                if full_range.start > reference_coordinates.start() || full_range.end < reference_coordinates.end() {
                    continue;
                }

                let read_bytes = read.seq().as_bytes();
                let d_mappings = dna_aligner.map(
                    &read_bytes,
                    output_cigar, output_md, max_frag_len, extra_flags, None
                )?;

                // we want to score relative to the target (reference) and penalize unmapped also
                let unmapped_from_target = true;
                let penalize_unmapped = true;
                let (best_mapping, best_stats) = select_best_mapping(
                    &d_mappings, unmapped_from_target, penalize_unmapped, None
                );

                if best_stats.mapping_score().score() > cli_settings.max_error_rate {
                    debug!("Best score for {qname} was {}, ignoring read.", best_stats.mapping_score().score());
                    
                    // this one is ignored, so add the stats here to make it clear that we intentionally ignored it later
                    let mapping_stats = HlaMappingStats::from_mapping_stats(None, Some(best_stats));
                    let ignored_details = PgxMappingDetails::new(
                        qname,
                        "REFERENCE".to_string(),
                        "REFERENCE".to_string(),
                        mapping_stats,
                        true
                    );
                    mapping_details.push(ignored_details);
                } else if let Some(bm) = best_mapping {
                    let is_reverse = bm.strand == minimap2::Strand::Reverse;
                    if is_reverse {
                        debug!("Ignoring reverse best mapping for {qname}: {}", best_stats.score_string());
                    } else {
                        debug!("Best score for {qname}: {}", best_stats.score_string());
                        let start = bm.query_start as usize;
                        let end = bm.query_end as usize;
                        let read_segment = std::str::from_utf8(&read_bytes[start..end])?.to_string();
                        read_segments.insert(qname.clone(), read_segment);
                        let (spliced_sequence, _offset) = splice_read(&mut read, database, gene_name)?;
                        let spliced_read = String::from_utf8(spliced_sequence)?;
                        spliced_segments.insert(qname, spliced_read);
                    }
                } else {
                    debug!("No mappings found for {qname}, ignoring read.");
                }
            }
        }

        let best_result = if read_segments.is_empty() {
            // we did not find any reads, definitely no way to get a consensus from that
            (NO_READS_HAP.to_string(), NO_READS_HAP.to_string())
        } else {
            let mut consensus_map: BTreeMap<String, String> = Default::default();
            let consensus = run_dual_consensus(&spliced_segments, cli_settings)?;
            let dual_passes = is_passing_dual(&consensus, cli_settings);
            
            // check if we found two using the cDNA
            let consensus = if dual_passes {
                // we did, continue on
                debug!("cDNA dual consensus successful.");
                consensus
            } else {
                // we did not, run dual on DNA
                debug!("cDNA dual consensus was homozygous, attempting dual consensus on DNA.");
                run_dual_consensus(&read_segments, cli_settings)?
            };

            // output debug records
            let mut debug_records = vec![];

            // re-run consensus on the groupings; this is required because we might have solved it via cDNA and need the DNA now
            let mut consensus_dwfa1 = ConsensusDWFA::with_config(dwfa_config_from_cli(cli_settings, false)?)?;
            let mut consensus_dwfa2 = ConsensusDWFA::with_config(dwfa_config_from_cli(cli_settings, false)?)?;
            for ((_qname, read_segment), &is_consensus1) in read_segments.iter().zip(consensus.is_consensus1().iter()) {
                if is_consensus1 {
                    consensus_dwfa1.add_sequence(read_segment.as_bytes())?;
                } else {
                    consensus_dwfa2.add_sequence(read_segment.as_bytes())?;
                }
            }

            let is_forward_strand = gene_definition.is_forward_strand();
            let c1_list = consensus_dwfa1.consensus()?;
            let con1 = std::str::from_utf8(c1_list[0].sequence())?.to_string();
            let con1_label = format!("consensus1_{gene_name}");
            if is_forward_strand {
                consensus_map.insert(con1_label, con1.clone());
            } else {
                consensus_map.insert(con1_label, String::from_utf8(reverse_complement(con1.as_bytes())?)?);
            }

            let (_map_results1, best_map1) = score_consensus(
                &dna_aligner, &reference_coordinates, &con1,
                database, gene_name, cli_settings
            )?;
            debug!("best_map1: {:?} {:?}", best_map1.best_match_id(), best_map1.best_match_star());

            // add the consensus to our debug records
            let tags = [("HP".to_string(), format!("1_consensus1_{gene_name}"))].into_iter().collect();
            match unmapped_record(
                &format!("consensus1_{gene_name}"),
                &con1,
                &tags
            ) {
                Ok(umr) => {
                    debug_records.push(umr);
                },
                Err(e) => {
                    error!("Error while creating unmapped record: {e}");
                }
            };
            
            // also add the matching haplotype for comparison
            if let Some(best_id) = best_map1.best_match_id() {
                let star_allele = best_map1.best_match_star().unwrap_or(UNKNOWN_HAP).to_string();
                let dna_sequence = database.hla_sequences().get(best_id).unwrap().dna_sequence();
                if let Some(sequence) = dna_sequence {
                    let sequence = if is_forward_strand {
                        sequence.to_string()
                    } else {
                        String::from_utf8(reverse_complement(sequence.as_bytes())?)?
                    };
                    let tags = [("HP".to_string(), format!("2_{gene_name}*{star_allele}"))].into_iter().collect();

                    match unmapped_record(
                        &star_allele,
                        &sequence,
                        &tags
                    ) {
                        Ok(umr) => {
                            debug_records.push(umr);
                        },
                        Err(e) => {
                            error!("Error while creating unmapped record: {e}");
                        }
                    };
                } else {
                    // I don't think we want a warning here since this is fairly common and outside user control
                }
            }

            // save the ID
            let best_id1 = best_map1.best_match_id().unwrap_or(UNKNOWN_HAP).to_string();
            debug_stats.add_read(gene_name.clone(), "consensus1".to_string(), best_map1)?;

            let best_dual_result = if consensus.is_dual() {
                // we have a dual consensus, type the second one also
                let c2_list = consensus_dwfa2.consensus()?;
                let con2 = std::str::from_utf8(c2_list[0].sequence())?.to_string();
                let con2_label = format!("consensus2_{gene_name}");
                if is_forward_strand {
                    consensus_map.insert(con2_label, con2.clone());
                } else {
                    consensus_map.insert(con2_label, String::from_utf8(reverse_complement(con2.as_bytes())?)?);
                }

                let (_map_results2, best_map2) = score_consensus(
                    &dna_aligner, &reference_coordinates, &con2,
                    database, gene_name, cli_settings
                )?;
                debug!("best_map2: {:?} {:?}", best_map2.best_match_id(), best_map2.best_match_star());

                // add the consensus to our debug records
                let tags = [("HP".to_string(), format!("4_consensus2_{gene_name}"))].into_iter().collect();
                match unmapped_record(
                    &format!("consensus2_{gene_name}"),
                    &con2,
                    &tags
                ) {
                    Ok(umr) => {
                        debug_records.push(umr);
                    },
                    Err(e) => {
                        error!("Error while creating unmapped record: {e}");
                    }
                };

                // also add the matching haplotype for comparison
                if let Some(best_id) = best_map2.best_match_id() {
                    // get the star allele and sequence
                    let star_allele = best_map2.best_match_star().unwrap_or(UNKNOWN_HAP).to_string();
                    let dna_sequence = database.hla_sequences().get(best_id).unwrap().dna_sequence();
                    if let Some(sequence) = dna_sequence {
                        let sequence = if is_forward_strand {
                            sequence.to_string()
                        } else {
                            String::from_utf8(reverse_complement(sequence.as_bytes())?)?
                        };

                        let tags = [("HP".to_string(), format!("5_{gene_name}*{star_allele}"))].into_iter().collect();
                        match unmapped_record(
                            &star_allele,
                            &sequence,
                            &tags
                        ) {
                            Ok(umr) => {
                                debug_records.push(umr);
                            },
                            Err(e) => {
                                error!("Error while creating unmapped record: {e}");
                            }
                        };
                    } else {
                        // I don't think we want a warning here since this is fairly common and outside user control
                    }
                }

                // save the id
                let best_id2 = best_map2.best_match_id().unwrap_or(UNKNOWN_HAP).to_string();
                debug_stats.add_read(gene_name.clone(), "consensus2".to_string(), best_map2)?;

                // we found a dual consensus, first check if the CDF and MAF are passing
                let total_count = read_segments.len();
                let counts1 = consensus.is_consensus1().iter().filter(|&&b| b).count();
                let counts2 = total_count - counts1;
                
                let dual_passed = is_passing_dual(&consensus, cli_settings);
                if dual_passed {
                    // MAF and CDF likelihood is above cutoff, so assume heterozygous is correct
                    (best_id1.clone(), best_id2)
                } else {
                    // MAF or CDF likelihood is below cutoff, so report homozygous for the dominant allele
                    debug!("MAF or CDF failed, returning homozygous result");
                    if counts1 > counts2 {
                        (best_id1.clone(), best_id1.clone())
                    } else {
                        (best_id2.clone(), best_id2)
                    }
                }
            } else {
                debug!("best_map2: No second consensus, homozygous result");
                (best_id1.clone(), best_id1.clone())
            };

            if let Some(debug_folder) = cli_settings.debug_folder.as_ref() {
                // save the consensus sequences
                let extension = format!("consensus_{gene_name}.fa");
                let consensus_fn = debug_folder.join(extension);
                debug!("Saving consensus for {gene_name} to {consensus_fn:?}");
                save_fasta(&consensus_map, &consensus_fn)?;
            }

            // pattern here is a little weird because we need the mutable reference multiple times in the loop
            if let Some(dbw) = debug_bam_writer.as_mut() {
                // add each extra target
                for eid in cli_settings.debug_hla_targets.iter() {
                    if let Some(hap_def) = database.hla_sequences().get(eid) {
                        // make sure it's a gene match, otherwise ignore for now
                        if hap_def.gene_name() != gene_name {
                            continue;
                        }
                        let dna_sequence = hap_def.dna_sequence();
                        if let Some(sequence) = dna_sequence {
                            // parse the star allele and get the correct sequence orientation
                            let star_allele = hap_def.star_allele().join(":");
                            let sequence = if is_forward_strand {
                                sequence.to_string()
                            } else {
                                String::from_utf8(reverse_complement(sequence.as_bytes())?)?
                            };

                            // save the record
                            let tags = [("HP".to_string(), format!("0_debug-target_{eid}_{gene_name}*{star_allele}"))].into_iter().collect();
                            match unmapped_record(
                                &star_allele,
                                &sequence,
                                &tags
                            ) {
                                Ok(umr) => {
                                    debug_records.push(umr);
                                },
                                Err(e) => {
                                    error!("Error while creating unmapped record: {e}");
                                }
                            };
                        } else {
                            warn!("Debug target \"{eid}\" does not have a DNA sequence in the database, ignoring it in output BAM");
                        }
                    } else {
                        warn!("Debug target \"{eid}\" was not found in the database, ignoring it in output BAM");
                    }
                }

                // we have a debug writer, so add all of the reads to this list also
                for ((qname, read_segment), &is_consensus1) in read_segments.iter().zip(consensus.is_consensus1().iter()) {
                    let con_type = if is_consensus1 { 1 } else { 2 };
                    let order_index = if is_consensus1 { 3 } else { 6 };
                    let phase_label = format!("{order_index}_consensus{con_type}_sequence");
                    let tags = [("HP".to_string(), phase_label)].into_iter().collect();
                    match unmapped_record(
                        qname,
                        read_segment,
                        &tags
                    ) {
                        Ok(umr) => {
                            debug_records.push(umr);
                        },
                        Err(e) => {
                            error!("Error while creating unmapped record: {e}");
                        }
                    };
                }

                // now map each record
                match dbw.map_records_to_region(&debug_records, &reference_coordinates) {
                    Ok(()) => {},
                    Err(e) => {
                        error!("Error while mappings records to debug BAM: {e}");
                    }
                };
            }

            best_dual_result
        };
        
        // For now, we only return one result, but lets leave the mechanisms for multiple in the future
        let best_combination = [best_result];

        // collect the diplotypes
        let diplotypes: Vec<Diplotype> = best_combination.iter().map(|(k1, k2)| {
            let star1 = match database.hla_sequences().get(k1) {
                Some(allele_def) => {
                    let s1 = allele_def.star_allele();
                    format!("*{}", s1.join(":"))
                },
                None => k1.clone()
            };
            let star2 = match database.hla_sequences().get(k2) {
                Some(allele_def) => {
                    let s2 = allele_def.star_allele();
                    format!("*{}", s2.join(":"))
                },
                None => k2.clone()
            };
            Diplotype::new(
                &star1,
                &star2
            )
        }).collect();

        debug!("Diplotype for {gene_name} => {:?}", diplotypes.iter().map(|d| d.diplotype()).collect::<Vec<&str>>());
        ret.insert(gene_name.clone(), PgxGeneDetails::new_from_mappings(
            diplotypes,
            None,
            mapping_details
        )?);
    }

    // if we have a debug output file, we can write it now
    // if let Some(debug_fn) = cli_settings.debug_hla_filename.as_ref() {
    if let Some(debug_folder) = cli_settings.debug_folder.as_ref() {
        let debug_fn = debug_folder.join("hla_debug.json");
        debug!("Saving HLA debug to {:?}", debug_fn);
        crate::util::file_io::save_json(&debug_stats, &debug_fn)?;
    }

    Ok(ret)
}

/// This is the updated functionality for batch HLA diplotyping from a BAM file.
/// This searches through all coordinates for mis-mapped sequences, which happens semi-frequently with the non-Class I alleles.
/// After assigning the reads to a gene bucket, it then runs DualConsensus like with `diplotype_hla`.
/// The extended HLA list has high diversity, so much so that we enabled a custom assembly writer with this functionality, like with D6.
/// # Arguments
/// * `gene_list` - the list of HLA genes to diplotype
/// * `database` - the pre-loaded database
/// * `bam_filenames` - list of BAM files containing reads to scan
/// * `reference_genome` - reference genome loaded into memory
/// * `debug_bam_writer` - optional BAM writer for debugging in the reference genome
/// * `debug_custom_writer` - optional BAM writer for custom visuals
/// * `cli_settings` - settings for diplotyping
/// # Errors
/// * if a gene is provided that we do not support
/// * if we cannot open or parse a BAM file correctly
pub fn diplotype_hla_batch(
    gene_list: &[String], database: &PgxDatabase,
    bam_filenames: &[PathBuf], reference_genome: &ReferenceGenome,
    mut debug_bam_writer: Option<&mut DebugBamWriter>,
    mut debug_custom_writer: Option<&mut IgvSessionWriter>,
    cli_settings: &DiplotypeSettings
) -> Result<HashMap<String, PgxGeneDetails>, Box<dyn std::error::Error>> {
    // if we have disabled cDNA scoring AND the DNA requirement is NOT enabled; then we can cause errors later due to lack of comparator sequence
    if cli_settings.disable_cdna_scoring && !cli_settings.hla_require_dna {
        bail!("If cDNA scoring is disabled, require HLA DNA must be enabled");
    }

    // prep all the bam readers
    let mut bam_readers: Vec<rust_htslib::bam::IndexedReader> = vec![];
    for bam_fn in  bam_filenames.iter() {
        let mut b = rust_htslib::bam::IndexedReader::from_path(bam_fn)?;
        b.set_reference(reference_genome.filename())?;
        bam_readers.push(b);
    }

    // stores the gene buckets for processing later
    let mut gene_buckets: BTreeMap<String, BTreeMap<String, RealignmentResult>> = Default::default();
    let mut qnames_checked: HashSet<String> = Default::default();

    // store all the relevant mapping details; this is really the best matching allele for a record in a vacuum
    let mut mapping_details: BTreeMap<String, Vec<PgxMappingDetails>> = Default::default();
    let mut read_debug = HlaDebug::new();

    // build our realigner from the database
    let realigner = HlaRealigner::new(
        gene_list, database, reference_genome
    )?;

    info!("Beginning batch HLA read search...");
    for gene_name in gene_list.iter() {
        // get the coordinates for this gene
        let gene_definition = match database.hla_config().gene_definition(gene_name) {
            Some(c) => c,
            None => bail!("No gene definition for {gene_name}")
        };
        let gene_coordinates = gene_definition.coordinates();
        debug!("Searching {gene_name} region ({} bp): {gene_coordinates}", gene_coordinates.len());

        // iterate over each bam, and fetch the reads
        for (bam_index, bam) in bam_readers.iter_mut().enumerate() {
            match bam.fetch(gene_coordinates.fetch_definition()) {
                Ok(()) => {},
                Err(e) => {
                    let filename = &bam_filenames[bam_index];
                    warn!("Received error \"{e}\" while fetching {gene_coordinates} in {filename:?}, assuming no reads for region.");
                    continue;
                }
            };

            for read_entry in bam.records() {
                // get the record and qname
                let record = read_entry.unwrap();
                let qname: String = std::str::from_utf8(record.qname())?.to_string();
                if !qnames_checked.insert(qname.clone()) {
                    // we have already checked this read, ignore it
                    continue;
                }

                let realigned_record = realigner.realign_record(&record)?;
                let best_gene = realigned_record.gene_name();
                if realigned_record.is_realigned() {
                    // add stat information
                    read_debug.add_read(best_gene.to_string(), qname.clone(), realigned_record.read_mapping_stats().clone())?;
                    
                    // more stats
                    let md = mapping_details.entry(best_gene.to_string()).or_default();
                    md.push(realigned_record.mapping_details().clone());

                    // get the gene bucket and add this record to the bucket, we will remap and splice within the bucket
                    let gene_bucket = gene_buckets.entry(best_gene.to_string()).or_default();
                    if let Some(old_record) = gene_bucket.insert(qname.clone(), realigned_record) {
                        // this should not happen with our previous check, but lets encode it to be sure
                        bail!("Found duplicate record in gene bucket: {old_record:?}");
                    }
                } else {
                    // if it failed, this is the only data we can expect to see; put it in the search bucket
                    let md = mapping_details.entry(gene_name.clone()).or_default();
                    md.push(realigned_record.mapping_details().clone());
                }
            }
        }
    }

    // normalize the coverage values for hemi-zygous identification later
    let mut read_total = 0;
    let mut hap_total = 0;
    for gene_name in NORMALIZING_HLA_GENES.iter() {
        if let Some(gb) = gene_buckets.get(gene_name) {
            read_total += gb.len();
            hap_total += 2;
        } else {
            warn!("No coverage detected for normalizing gene: {gene_name}");
        }
    }

    let normalized_coverage = if hap_total > 0 {
        let norm_coverage = read_total as f64 / hap_total as f64;
        debug!("Normalized haplotype coverage = {norm_coverage:.3}");
        Some(norm_coverage)
    } else {
        warn!("No coverage found for any normalizing genes.");
        None
    };

    debug!("Total reads found:");
    for gene_name in gene_list.iter() {
        let count = if let Some(gb) = gene_buckets.get(gene_name) {
            gb.len()
        } else {
            0
        };
        let norm_count = normalized_coverage.map(|nc| count as f64 / nc);
        debug!("\t{gene_name} => {count}, norm_count = {norm_count:?}");
    }

    // if we have a debug output file, we can write our read debug file out
    if let Some(debug_folder) = cli_settings.debug_folder.as_ref() {
        let debug_fn = debug_folder.join("read_debug.json");
        debug!("Saving read debug to {:?}", debug_fn);
        crate::util::file_io::save_json(&read_debug, &debug_fn)?;
    }

    //set up job configuration
    let mut ret: HashMap<String, PgxGeneDetails> = Default::default();
    let mut debug_stats = HlaDebug::new();

    // iterate over each gene bucket and solve it
    for gene_name in gene_list.iter() {
        info!("Solving {gene_name}...");
        
        // output debug records - these get re-mapped in BOTH the hg38 and the custom output BAM
        let mut debug_records = vec![];

        // get the bucket, we have one
        let gene_def = database.hla_config().gene_definition(gene_name).unwrap();
        let is_forward_strand = gene_def.is_forward_strand();
        let region = gene_def.coordinates();
        let region_sequence = reference_genome.get_slice(region.chrom(), region.start() as usize, region.end() as usize);
        let region_aligner = standard_hifi_aligner()
            .with_seq(region_sequence)?;

        // get whatever records we found
        let opt_bucket = gene_buckets.get(gene_name);
        let empty_bucket: BTreeMap<String, RealignmentResult> = Default::default();
        let realigned_records = opt_bucket.unwrap_or(&empty_bucket);
        debug!("Found {} realigned records for {gene_name}.", realigned_records.len());

        let best_result = if realigned_records.is_empty() {
            // we did not find any reads, definitely no way to get a consensus from that
            if gene_def.is_absent_capable() {
                (NO_CALL_HAP.to_string(), NO_CALL_HAP.to_string())
            } else {
                (NO_READS_HAP.to_string(), NO_READS_HAP.to_string())
            }
        } else {
            // collects the sequences for the FASTA output
            let mut consensus_map: BTreeMap<String, String> = Default::default();
            
            // run the consensus as if it's a dual
            let preliminary_consensus = run_dual_consensus_with_offsets(realigned_records, cli_settings)?;

            // check for hemizygosity based on the result
            let is_hemizygous = if gene_def.is_absent_capable() {
                // this one is capable of having absent haplotypes and the coverage is sufficiently low to be haploid
                // calculate a score for each and see how they compare
                debug!("Checking for hemizygosity...");
                is_hemizygous_better(&preliminary_consensus, cli_settings, normalized_coverage)?
            } else {
                // not capable of haploid OR looks diploid from coverage
                false
            };

            // now set the consensus based on whether we are hemizgyous result or not
            let consensus = if is_hemizygous {
                // return a boiler plate consensus that is not dual, none of the other details matter here since we re-consensus later
                DualConsensus::new(
                    Consensus::new(
                        vec![],
                        waffle_con::cdwfa_config::ConsensusCost::L1Distance,
                        vec![0; realigned_records.len()]
                    ),
                    None,
                    vec![true; realigned_records.len()],
                    vec![Some(0); realigned_records.len()],
                    vec![None; realigned_records.len()]
                )?
            } else {
                // this one should have two, so just return the initial DC result
                preliminary_consensus
            };

            // re-run consensus on the groupings; this is required because we might have solved it via cDNA and need the DNA now
            // the half-window is added to each offset to "center" the window
            let config = dwfa_config_from_cli(cli_settings, true)?;
            let half_window = config.offset_window / 2;
            
            let mut min_offset1 = usize::MAX;
            let mut min_offset2 = usize::MAX;

            // first, figure out the minimum offsets in each
            for ((_qname, realignment_result), &is_consensus1) in realigned_records.iter().zip(consensus.is_consensus1().iter()) {
                let record = realignment_result.realigned_record().unwrap();
                if is_consensus1 {
                    min_offset1 = min_offset1.min(record.dna_offset());
                } else {
                    min_offset2 = min_offset2.min(record.dna_offset());
                }
            }
            debug!("Minimum offset1: {min_offset1}");
            debug!("Minimum offset2: {min_offset2}");

            // now assign the sequences with appropriately shifted offsets
            let mut consensus_dwfa1 = ConsensusDWFA::with_config(config.clone())?;
            let mut consensus_dwfa2 = ConsensusDWFA::with_config(config)?;
            for ((_qname, realignment_result), &is_consensus1) in realigned_records.iter().zip(consensus.is_consensus1().iter()) {
                let record = realignment_result.realigned_record().unwrap();
                if is_consensus1 { // if the dual failed, everything should go into con1
                    let offset = if record.dna_offset() == min_offset1 { None } else { Some(record.dna_offset() - min_offset1 + half_window) };
                    consensus_dwfa1.add_sequence_offset(record.dna_sequence(), offset)?;
                } else {
                    let offset = if record.dna_offset() == min_offset2 { None } else { Some(record.dna_offset() - min_offset2 + half_window) };
                    consensus_dwfa2.add_sequence_offset(record.dna_sequence(), offset)?;
                }
            }

            // get consensus 1 and type it
            let c1_list = match consensus_dwfa1.consensus() {
                Ok(l) => l,
                Err(e) => {
                    error!("Failed to generate a consensus: {e}");
                    vec![]
                }
            };
            let con1 = if c1_list.is_empty() {
                // failed to generate a consensus, return empty string which gets tagged as unknown
                String::new()
            } else {
                let c1 = &c1_list[0];
                debug!("c1 scores: {:?}", c1.scores());    
                std::str::from_utf8(c1.sequence())?.to_string()
            };

            debug!("Scoring consensus1 of length {} bp...", con1.len());
            let (_map_results1, best_map1) = score_consensus(
                &region_aligner, region, &con1,
                database, gene_name, cli_settings
            )?;
            debug!("best_map1: {:?} {:?}", best_map1.best_match_id(), best_map1.best_match_star());

            // save the sequence into our FASTA lookup
            let con1_label = format!("consensus1_{gene_name}");
            if is_forward_strand {
                consensus_map.insert(con1_label, con1.clone());
            } else {
                consensus_map.insert(con1_label, String::from_utf8(reverse_complement(con1.as_bytes())?)?);
            }

            // consensus1 as an unmapped read
            let tags = [("HP".to_string(), format!("1_consensus1_{gene_name}"))].into_iter().collect();
            match unmapped_record(&format!("consensus1_{gene_name}"), &con1, &tags) {
                Ok(umr) => { debug_records.push(umr); },
                Err(e) => { error!("Error while creating unmapped record: {e}"); }
            };

            // also add the matching haplotype for comparison
            if let Some(best_id) = best_map1.best_match_id() {
                let star_allele = best_map1.best_match_star().unwrap_or(UNKNOWN_HAP).to_string();
                let dna_sequence = match database.hla_sequences().get(best_id) {
                    Some(ad) => ad.dna_sequence(),
                    None => None // like from UNKNOWN type
                };
                if let Some(sequence) = dna_sequence {
                    let sequence = if is_forward_strand {
                        sequence.to_string()
                    } else {
                        String::from_utf8(reverse_complement(sequence.as_bytes())?)?
                    };
                    let tags = [("HP".to_string(), format!("2_{gene_name}*{star_allele}"))].into_iter().collect();

                    match unmapped_record(&star_allele, &sequence, &tags) {
                        Ok(umr) => { debug_records.push(umr); },
                        Err(e) => { error!("Error while creating unmapped record: {e}"); }
                    };
                } else {
                    // I don't think we want a warning here since this is fairly common and outside user control
                }
            }
            
            // save the ID
            let best_id1 = best_map1.best_match_id().unwrap_or(UNKNOWN_HAP).to_string();
            debug_stats.add_read(gene_name.clone(), "consensus1".to_string(), best_map1.clone())?;

            // these are only populated if consensus 2 exists
            let opt_con2: Option<String>;
            let opt_best_map2: Option<ReadMappingStats>;

            let best_dual_result = if consensus.is_dual() {
                // we have a dual consensus, type the second one also
                let c2_list = match consensus_dwfa2.consensus() {
                    Ok(l) => l,
                    Err(e) => {
                        error!("Failed to generate a consensus: {e}");
                        vec![]
                    }
                };
                let con2 = if c2_list.is_empty() {
                    // failed to generate a consensus, return empty string which gets tagged as unknown
                    String::new()
                } else {
                    let c2 = &c2_list[0];
                    debug!("c2 scores: {:?}", c2.scores());
                    std::str::from_utf8(c2_list[0].sequence())?.to_string()
                };
                debug!("Scoring consensus2 of length {} bp...", con2.len());
                let (_map_results2, best_map2) = score_consensus(
                    &region_aligner, region, &con2,
                    database, gene_name, cli_settings
                )?;
                debug!("best_map2: {:?} {:?}", best_map2.best_match_id(), best_map2.best_match_star());

                // save the sequence into our FASTA lookup
                let con2_label = format!("consensus2_{gene_name}");
                if is_forward_strand {
                    consensus_map.insert(con2_label, con2.clone());
                } else {
                    consensus_map.insert(con2_label, String::from_utf8(reverse_complement(con2.as_bytes())?)?);
                }

                // save the second consensus to our debug unmapped reads
                let tags = [("HP".to_string(), format!("4_consensus2_{gene_name}"))].into_iter().collect();
                match unmapped_record(&format!("consensus2_{gene_name}"), &con2, &tags) {
                    Ok(umr) => { debug_records.push(umr); },
                    Err(e) => { error!("Error while creating unmapped record: {e}"); }
                };

                // also add the matching haplotype for comparison
                if let Some(best_id) = best_map2.best_match_id() {
                    // get the star allele and sequence
                    let star_allele = best_map2.best_match_star().unwrap_or(UNKNOWN_HAP).to_string();
                    let dna_sequence = match database.hla_sequences().get(best_id) {
                        Some(ad) => ad.dna_sequence(),
                        None => None // like from UNKNOWN type
                    };
                    if let Some(sequence) = dna_sequence {
                        let sequence = if is_forward_strand {
                            sequence.to_string()
                        } else {
                            String::from_utf8(reverse_complement(sequence.as_bytes())?)?
                        };

                        let tags = [("HP".to_string(), format!("5_{gene_name}*{star_allele}"))].into_iter().collect();
                        match unmapped_record(&star_allele, &sequence, &tags) {
                            Ok(umr) => { debug_records.push(umr); },
                            Err(e) => { error!("Error while creating unmapped record: {e}"); }
                        };
                    } else {
                        // I don't think we want a warning here since this is fairly common and outside user control
                    }
                }

                // save the ID
                let best_id2 = best_map2.best_match_id().unwrap_or(UNKNOWN_HAP).to_string();
                debug_stats.add_read(gene_name.clone(), "consensus2".to_string(), best_map2.clone())?;

                // lastly, check if we actually want to report this, or if it's a dual that doesn't pass our config
                let total_count = realigned_records.len();
                let counts1 = consensus.is_consensus1().iter().filter(|&&b| b).count();
                let counts2 = total_count - counts1;

                // populate the extras we need to pull out
                opt_con2 = Some(con2);
                opt_best_map2 = Some(best_map2);

                // check if the final consensus passes our cutoffs
                let dual_passed = is_passing_dual(&consensus, cli_settings);
                if dual_passed {
                    // MAF and CDF likelihood is above cutoff, so assume heterozygous is correct
                    (best_id1.clone(), best_id2)
                } else {
                    // MAF or CDF likelihood is below cutoff, so report homozygous for the dominant allele
                    debug!("MAF or CDF failed, returning homozygous result");
                    if counts1 > counts2 {
                        (best_id1.clone(), best_id1)
                    } else {
                        (best_id2.clone(), best_id2)
                    }
                }
            } else {
                debug!("best_map2: No second consensus, homozygous result");
                opt_con2 = None;
                opt_best_map2 = None;
                (best_id1.clone(), best_id1.clone())
            };

            // check if we have a hemizygous result 
            let best_dual_result = if is_hemizygous {
                assert_eq!(best_dual_result.0, best_dual_result.1); // sanity check; we should have a "homozygous" result from above
                debug!("Converting to hemizygous result...");
                (NO_CALL_HAP.to_string(), best_dual_result.0)
            } else {
                // heterozygous call, just return both for now
                best_dual_result
            };

            if let Some(debug_folder) = cli_settings.debug_folder.as_ref() {
                // save the consensus sequences
                let extension = format!("consensus_{gene_name}.fa");
                let consensus_fn = debug_folder.join(extension);
                debug!("Saving consensus for {gene_name} to {consensus_fn:?}");
                save_fasta(&consensus_map, &consensus_fn)?;
            }

            // pattern here is a little weird because we need the mutable reference multiple times in the loop
            if let Some(dbw) = debug_bam_writer.as_mut() {
                // add each extra target
                for eid in cli_settings.debug_hla_targets.iter() {
                    if let Some(hap_def) = database.hla_sequences().get(eid) {
                        // make sure it's a gene match, otherwise ignore for now
                        if hap_def.gene_name() != gene_name {
                            continue;
                        }
                        let dna_sequence = hap_def.dna_sequence();
                        if let Some(sequence) = dna_sequence {
                            // parse the star allele and get the correct sequence orientation
                            let star_allele = hap_def.star_allele().join(":");
                            let sequence = if is_forward_strand {
                                sequence.to_string()
                            } else {
                                String::from_utf8(reverse_complement(sequence.as_bytes())?)?
                            };

                            // save the record
                            let tags = [("HP".to_string(), format!("0_debug-target_{eid}_{gene_name}*{star_allele}"))].into_iter().collect();
                            match unmapped_record(&star_allele, &sequence, &tags) {
                                Ok(umr) => { debug_records.push(umr); },
                                Err(e) => { error!("Error while creating unmapped record: {e}"); }
                            };
                        } else {
                            warn!("Debug target \"{eid}\" does not have a DNA sequence in the database, ignoring it in output BAM");
                        }
                    } else {
                        warn!("Debug target \"{eid}\" was not found in the database, ignoring it in output BAM");
                    }
                }

                // we have a debug writer, so add all of the reads to this list also
                for ((_qname, realignment_result), &is_consensus1) in realigned_records.iter().zip(consensus.is_consensus1().iter()) {
                    let realigned_record = realignment_result.realigned_record().unwrap();
                    let con_type = if is_consensus1 { 1 } else { 2 };
                    let order_index = if is_consensus1 { 3 } else { 6 };
                    let phase_label = format!("{order_index}_consensus{con_type}_sequence");
                    let tags = [("HP".to_string(), phase_label)].into_iter().collect();

                    match clear_record(realigned_record.realigned_record(), &tags) {
                        Ok(umr) => { debug_records.push(umr); },
                        Err(e) => { error!("Error while creating unmapped record: {e}"); }
                    };
                }

                // now map each record
                match dbw.map_records_to_region(&debug_records, region) {
                    Ok(()) => {},
                    Err(e) => {
                        error!("Error while mappings records to debug BAM: {e}");
                    }
                };
            }

            // build our custom visualization for this region also
            if let Some(dcw) = debug_custom_writer.as_mut() {
                let region_name = format!("{gene_name}_{CONTIG_POSTFIX}");
                
                // generic buffer between regions
                let buffer_sequence: String = "N".repeat(BUFFER_LEN);
                
                // start with a buffer
                let mut region_sequence = buffer_sequence.clone();
                let mut regions: Vec<(Coordinates, String)> = vec![];

                // add consensus 1
                let start_index = region_sequence.len() as u64;
                region_sequence.push_str(&con1);
                let end_index = region_sequence.len() as u64;
                region_sequence.push_str(&buffer_sequence);
                regions.push((
                    Coordinates::new(region_name.clone(), start_index, end_index),
                    format!("0_{gene_name}*{}", best_map1.best_match_star().unwrap_or(UNKNOWN_HAP))
                ));

                /*
                  * TODO: right now, this will always add con2, even if we report homozygous con1 (or hom con2); is this the desired behavior or should we limit to only what is reported?
                  * I think it is given that this is for debugging
                  * if we change it, we have to handle all the scenarios for homozygous (and hom swap), which will be difficult
                  */

                // add consensus 2 if we have one
                if consensus.is_dual() {
                    let con2 = opt_con2.unwrap();
                    let best_map2 = opt_best_map2.unwrap();

                    let start_index = region_sequence.len() as u64;
                    region_sequence.push_str(&con2);
                    let end_index = region_sequence.len() as u64;
                    region_sequence.push_str(&buffer_sequence);
                    regions.push((
                        Coordinates::new(region_name.clone(), start_index, end_index),
                        format!("1_{gene_name}*{}", best_map2.best_match_star().unwrap_or(UNKNOWN_HAP))
                    ));
                }
                
                if let Err(e) = dcw.add_custom_region(
                    region_name, &region_sequence, &regions, debug_records
                ) {
                    error!("Error while adding custom region: {e}");
                }
            }

            best_dual_result
        };

        // For now, we only return one result, but lets leave the mechanisms for multiple in the future
        let best_combination = [best_result];

        // collect the diplotypes
        let diplotypes: Vec<Diplotype> = best_combination.iter().map(|(k1, k2)| {
            let star1 = match database.hla_sequences().get(k1) {
                Some(allele_def) => {
                    let s1 = allele_def.star_allele();
                    format!("*{}", s1.join(":"))
                },
                None => k1.clone()
            };
            let star2 = match database.hla_sequences().get(k2) {
                Some(allele_def) => {
                    let s2 = allele_def.star_allele();
                    format!("*{}", s2.join(":"))
                },
                None => k2.clone()
            };
            Diplotype::new(
                &star1,
                &star2
            )
        }).collect();

        debug!("Diplotype for {gene_name} => {:?}", diplotypes.iter().map(|d| d.diplotype()).collect::<Vec<&str>>());
        ret.insert(gene_name.clone(), PgxGeneDetails::new_from_mappings(
            diplotypes,
            None,
            mapping_details.remove(gene_name).unwrap_or_default()
        )?);
    }

    // if we have a debug output file, we can write it now
    if let Some(debug_folder) = cli_settings.debug_folder.as_ref() {
        let debug_fn = debug_folder.join("hla_debug.json");
        debug!("Saving HLA debug to {:?}", debug_fn);
        crate::util::file_io::save_json(&debug_stats, &debug_fn)?;
    }

    Ok(ret)
}

/// Wrapper function for if we allow an allele definition
/// # Arguments
/// * `hla_allele_def` - the definition in question
/// * `gene_name` - the name of the gene we are matching against
/// * `cli_settings` - any settings from the CLI; some control the behavior of this function
fn is_allowed_allele_def(hla_allele_def: &HlaAlleleDefinition, gene_name: &str, cli_settings: &DiplotypeSettings) -> bool {
    // require the gene name to match AND
    hla_allele_def.gene_name() == gene_name && 
        // require that either we have DNA sequence OR that the DNA sequence is not required
        (hla_allele_def.dna_sequence().is_some() || !cli_settings.hla_require_dna)
}

/// Wrapper script for setting up the consensus configuration from our CLI.
/// # Arguments
/// * `cli_settings` - the provided CLI options
/// * `allow_early_termination` - if True, this indicates the sequences may not fully span the full region
/// # Errors
/// * if the CLI options are not allowed by our config
fn dwfa_config_from_cli(cli_settings: &DiplotypeSettings, allow_early_termination: bool) -> Result<CdwfaConfig, CdwfaConfigBuilderError> {
    CdwfaConfigBuilder::default()
        .min_count(cli_settings.min_consensus_count)
        .min_af(cli_settings.min_consensus_fraction)
        .dual_max_ed_delta(cli_settings.dual_max_ed_delta)
        .allow_early_termination(allow_early_termination)
        .weighted_by_ed(false) // currently, I'm not convinced on either approach
        .consensus_cost(waffle_con::cdwfa_config::ConsensusCost::L1Distance)
        .max_queue_size(20)
        .max_capacity_per_size(10)
        .offset_window(400) // ideally, we would have accurate estimates, but this is tricky with the HLA genes
        .offset_compare_length(50) // 50 is default
        .build()
}

/// Runs dual consensus on an ordered collection of sequences. Useful because we may solve via cDNA or DNA depending on dataset.
/// If multiple equal options are found, this only returns the first one.
/// # Arguments
/// * `segments` - a map from segment ID to a segment sequence; traversed in order
/// * `cli_settings` - contains controls for the consensus step
/// # Errors
/// * if the cli_settings are not valid for consensus
/// * if the consensus itself has errors while running
fn run_dual_consensus(segments: &BTreeMap<String, String>, cli_settings: &DiplotypeSettings) -> Result<DualConsensus, Box<dyn std::error::Error>> {
    // now prep the priority consensus runner - we have HPC as priority, then full length
    let mut consensus_dwfa = DualConsensusDWFA::with_config(dwfa_config_from_cli(cli_settings, false)?)?;
    for (_qname, read_segment) in segments.iter() {
        consensus_dwfa.add_sequence(read_segment.as_bytes())?;
    }

    let mut consensus_list = consensus_dwfa.consensus()?;
    if consensus_list.len() > 1 {
        warn!("Found multiple solutions, selecting first.");
    }
    Ok(consensus_list.remove(0))
}

/// Runs dual consensus on an ordered collection of realigned HLA sequences.
/// This handles the partial alignments from our longer HLA loci (i.e., sequences that do not fully span the locus).
/// This will run it first on the cDNA and then the DNA if the cDNA does not split it appropriately
/// If multiple equal options are found, this only returns the first one.
/// # Arguments
/// * `segments` - a map from segment ID to a segment sequence; traversed in order
/// * `cli_settings` - contains controls for the consensus step
/// # Errors
/// * if the cli_settings are not valid for consensus
/// * if the consensus itself has errors while running
fn run_dual_consensus_with_offsets(segments: &BTreeMap<String, RealignmentResult>, cli_settings: &DiplotypeSettings) -> Result<DualConsensus, Box<dyn std::error::Error>> {
    // get the config, and also the half-window since we have to deal with offsets
    let config = dwfa_config_from_cli(cli_settings, true)?;
    // the half-window is added to each offset to "center" the window
    let half_window = config.offset_window / 2;

    // we need to adjust the offsets ourselves for accuracy
    let min_hpc_offset = segments.values()
        .map(|r| r.realigned_record().unwrap().hpc_offset())
        .min().unwrap();
    debug!("Min HPC offset: {min_hpc_offset}");

    // first, run on the HPC sequence
    let mut hpc_consensus_dwfa = DualConsensusDWFA::with_config(config.clone())?;
    for (_qname, realignment_result) in segments.iter() {
        let record = realignment_result.realigned_record().unwrap();
        let offset = if record.hpc_offset() == min_hpc_offset { None } else { Some(record.hpc_offset() - min_hpc_offset + half_window) };
        hpc_consensus_dwfa.add_sequence_offset(record.hpc_sequence(), offset)?;
    }
    let mut consensus_list = hpc_consensus_dwfa.consensus()?;
    if consensus_list.len() > 1 {
        warn!("Found multiple solutions, selecting first.");
    }
    let hpc_consensus = consensus_list.remove(0);
    debug!("DNA scores1: {:?}", hpc_consensus.scores1());
    debug!("DNA scores2: {:?}", hpc_consensus.scores2());
    if is_passing_dual(&hpc_consensus, cli_settings) {
        debug!("HPC consensus passed.");
        /*
        // just for debugging
        for ((qname, realignment_result), (s1, s2)) in segments.iter().zip(hpc_consensus.scores1().iter().zip(hpc_consensus.scores2().iter())) {
            println!("{qname}\t{}\t{s1:?}\t{s2:?}", realignment_result.realigned_record().unwrap().hpc_offset());
        }
        */
        return Ok(hpc_consensus);
    } else {
        debug!("HPC consensus did not pass, trying full-length DNA sequences...");
    }

    // HPC did not find a difference, so try full length DNA
    let min_dna_offset = segments.values()
        .map(|r| r.realigned_record().unwrap().dna_offset())
        .min().unwrap();
    debug!("Min DNA offset: {min_dna_offset}");

    // now add sequences with adjust offset
    let mut consensus_dwfa = DualConsensusDWFA::with_config(config)?;
    for (_qname, realignment_result) in segments.iter() {
        let record = realignment_result.realigned_record().unwrap();
        let offset = if record.dna_offset() == min_dna_offset { None } else { Some(record.dna_offset() - min_dna_offset + half_window) };
        consensus_dwfa.add_sequence_offset(record.dna_sequence(), offset)?;
    }
    let mut consensus_list = consensus_dwfa.consensus()?;
    if consensus_list.len() > 1 {
        warn!("Found multiple solutions, selecting first.");
    }
    let dna_consensus = consensus_list.remove(0);
    debug!("DNA scores1: {:?}", dna_consensus.scores1());
    debug!("DNA scores2: {:?}", dna_consensus.scores2());

    /*
    // debug zipper that prints qname, dna_offset, score1, score2 
    for ((qname, realignment_result), (s1, s2)) in segments.iter().zip(dna_consensus.scores1().iter().zip(dna_consensus.scores2().iter())) {
        println!("{qname}\t{}\t{s1:?}\t{s2:?}", realignment_result.realigned_record().unwrap().dna_offset());
    }
    */

    Ok(dna_consensus)
}

/// Checks a DualConsensus solution to see if it passes our MAF and CDF cutoffs from the CLI.
/// # Arguments
/// * `dual_consensus` - the consensus to check; if not dual, this will always return false
/// * `cli_settings` - contains the cutoffs we compare against
fn is_passing_dual(dual_consensus: &DualConsensus, cli_settings: &DiplotypeSettings) -> bool {
    if dual_consensus.is_dual() {
        // we found a dual consensus, first check if the CDF and MAF are passing
        let maf_cutoff = cli_settings.min_consensus_fraction;
        let total_count = dual_consensus.is_consensus1().len();
        let counts1 = dual_consensus.is_consensus1().iter().filter(|&&b| b).count();
        let counts2 = total_count - counts1;
        let minor_count = counts1.min(counts2) as f64;
        let maf = minor_count / (total_count as f64);
        
        // binomial distribution based decision making for het/hom
        // calculate based on the expected MAF; if maf > expected MAF, this will always pass our default threshold
        let distro = Binomial::new(cli_settings.expected_maf, total_count as u64).unwrap();
        let cdf = distro.cdf(minor_count as u64);

        let cdf_cutoff = cli_settings.min_cdf;
        let is_passing = maf >= maf_cutoff && cdf >= cdf_cutoff;
        debug!("DualConsensus detected: counts1={counts1}, counts2={counts2}, MAF={maf:.5}, CDF={cdf:.5}; is_passing={is_passing}");
        is_passing
    } else {
        false
    }
}

/// Given a consensus sequence from a particular region, this will score it against all the entries in the database by comparing both the cDNA and DNA sequences.
/// cDNA takes priority, and then DNA is mostly a tie-breaker (4th-field).
/// # Arguments
/// * `ref_aligner` - a preconstructed aligner to the reference sequence, which allows us to align and then splice out the cDNA
/// * `ref_coordinates` - coordinates of the reference sequence, which are need to create an artificial read
/// * `consensus` - the consensus sequence we are comparing to the DB; this must always be on the forward strand relative to hg38
/// * `database` - contains all the HLA entries
/// * `gene_name` - gene we are trying to label
/// * `cli_settings` - parameters that control this step
fn score_consensus(
    ref_aligner: &minimap2::Aligner<minimap2::Built>, ref_coordinates: &Coordinates, consensus: &str,
    database: &PgxDatabase, gene_name: &str, cli_settings: &DiplotypeSettings
) -> Result<(HashMap<String, HlaMappingStats>, ReadMappingStats), Box<dyn std::error::Error>> {
    // this can happen if the consensus fails
    if consensus.is_empty() {
        let empty_map = Default::default();
        let unknown_best = ReadMappingStats::new();
        return Ok((empty_map, unknown_best));
    }

    // we need cigar here
    // other settings for mapping
    let output_cigar: bool = true;
    let output_md: bool = true;
    let max_frag_len: Option<usize> = None;
    let extra_flags = None;

    // target is reference; query is the consensus
    let d_mappings = ref_aligner.map(
        consensus.as_bytes(),
        output_cigar, output_md, max_frag_len, extra_flags, None
    )?;
    
    if d_mappings.is_empty() {
        error!("Failed to align consensus to reference genome, returning empty result.");
        let empty_map = Default::default();
        let unknown_best = ReadMappingStats::new();
        return Ok((empty_map, unknown_best));
    }

    // we want to score relative to the reference, and penalize anything that is unmapped
    let unmapped_from_target = true;
    let penalize_unmapped = true;
    let (best_mapping, _best_score) = select_best_mapping(
        &d_mappings, unmapped_from_target, penalize_unmapped, None
    );

    // we know there should be Some
    let d_map = best_mapping.unwrap();
    assert!(d_map.strand == minimap2::Strand::Forward);
    
    // we need to include the whole consensus since for DRB1 (and others), it does not map well to reference
    let cigar_string = convert_mapping_to_cigar(d_map, None, None);

    // create an artificial bam record that we can feed to our score_read function
    let mut consensus_record = rust_htslib::bam::Record::new();
    let record_name = b"consensus";
    let qual = vec![255; consensus.len()];
    let cigar = Some(&cigar_string);
    consensus_record.set(
        record_name,
        cigar,
        consensus.as_bytes(),
        &qual
    );
    let start_align = ref_coordinates.start() as i64 + d_map.target_start as i64;
    consensus_record.set_pos(start_align);

    // score it as if we pulled it directly from a BAM file
    score_read(consensus_record, database, gene_name, cli_settings)
}

/// Work-horse function for scoring every allele in the HLA database against a single read.
/// Originally, this was intended to be called in parallel, once for each read sequence.
/// Now, it is called up to twice per gene, once for each consensus sequence after re-mapping it back to the reference.
/// Returns a hash map of all the scores as well as the single best ID.
/// # Arguments
/// * `read` - the read we are aligning against, which is always on the forward strand
/// * `database` - the full PGx database, including HLA sequences
/// * `gene_name` - the gene we want to compare against
/// # Errors
/// * if there are errors parsing the read record
/// * if there are errors when aligning the sequences
fn score_read(mut read: rust_htslib::bam::Record, database: &PgxDatabase, gene_name: &str, cli_settings: &DiplotypeSettings) 
    -> Result<(HashMap<String, HlaMappingStats>, ReadMappingStats), Box<dyn std::error::Error>> {
    // result is a map of HLA_ID -> score
    let mut ret: HashMap<String, HlaMappingStats> = Default::default();

    let is_forward_strand = database.hla_config().gene_definition(gene_name).unwrap().is_forward_strand();

    // get the read sequence for alignment
    let qname: String = std::str::from_utf8(read.qname())?.to_string();
    let read_sequence: Vec<u8> = read.seq().as_bytes();

    // put it on the correct strand
    let read_sequence: Vec<u8> = if is_forward_strand {
        // read is always on forward strand, and the gene is also; no-op
        read_sequence
    } else {
        // read is always on forward strand, but the gene is on reverse; make the read match the gene
        reverse_complement(&read_sequence)?
    };

    // figure out the cDNA sequence
    let prespliced_sequence = if !cli_settings.disable_cdna_scoring {
        let (fw_spliced, _fw_offset) = splice_read(&mut read, database, gene_name)?;
        if fw_spliced.is_empty() {
            warn!("Pre-splicing of sequence found no exons, setting cDNA to \"N\"");
            // "N".to_string()
            vec![b'N']
        } else if is_forward_strand {
            fw_spliced
        } else {
            reverse_complement(&fw_spliced)?
        }
    } else {
        // cDNA scoring is disabled, this doesn't matter
        vec![b'N']
    };

    // cDNA aligner
    let mut cdna_aligner = standard_hifi_aligner();
    cdna_aligner.mapopt.a = 5; // see notes below
    let cdna_aligner = cdna_aligner
        .with_seq(&prespliced_sequence)?;

    // DNA aligner is much more straight-forward, just use full sequence
    let mut dna_aligner = standard_hifi_aligner();
    dna_aligner.mapopt.a = 5; // see notes below
    let dna_aligner = dna_aligner
        .with_seq(&read_sequence)?;

    // default options: a: 1, b: 4, q: 6, e: 2, q2: 26, e2: 1
    // descriptions from minimap2 here: https://github.com/lh3/minimap2/blob/69e36299168d739dded1c5549f662793af10da83/minimap.h#L157
    // John suggested increasing "a" from 1 to some higher value to remove clipping
    // equality does not seem to err on the side of inclusion; we had two identical + 1 mismatch and a=2 did not work (2 + 2 - 4); a=3 did work
    // so in theory, if we have one matching base separated by a mismatch, then a > 4 to catch it
    // cdna_aligner.mapopt.a = 5;
    // dna_aligner.mapopt.a = 5;

    // we only need cigar and md for debugging
    // other settings for mapping
    let output_cigar: bool = true;
    let output_md: bool = true;
    let max_frag_len: Option<usize> = None;
    // see https://github.com/lh3/minimap2/blob/69e36299168d739dded1c5549f662793af10da83/minimap.h#L36 for more flag options
    let extra_flags = Some(vec![0x4000000]); // enables the X/= cigar tags

    // used to be controled by targets when we had one entry per read, now we can just assume this is fine
    let all_hla_targets: bool = true; // cli_settings.debug_hla_targets.contains(&"all".to_string());
    let mut read_mapping_stats = ReadMappingStats::new();

    /*
     * Note for future Matt: we tried using HPC DNA aligner as well. It sort of worked, but the problem is when true variation
     *      gets masked by the HPC. Also, if the HPC alleles are different lengths, it always picks the longer one.
     *      We need something smarter, and I don't think HPC is it due to variant masking...
     */
    let aligners = [
        cdna_aligner, // cDNA
        dna_aligner // DNA
    ];

    // this is all of the best results thus far
    let mut best_match: HlaProcessedMatch = HlaProcessedMatch::worst_match(aligners.len());
    for (hla_id, hla_allele_def) in database.hla_sequences().iter() {
        // only check the sequences for this particular gene
        if !is_allowed_allele_def(hla_allele_def, gene_name, cli_settings) {
            continue;
        }

        // these are the ordered rankings of our sequence mappings; cDNA -> HPC DNA -> full DNA
        let sequences = [
            // cDNA sequence - 2nd/3rd field
            if cli_settings.disable_cdna_scoring { None } else { Some(hla_allele_def.cdna_sequence()) },
            // lastly DNA sequence - 4th field, HPC tie-break
            hla_allele_def.dna_sequence()
        ];

        // sanity check for dev
        assert_eq!(aligners.len(), sequences.len());

        // these should be the same length as sequences at the end
        let mut current_match = HlaProcessedMatch::new(hla_id.clone())?;

        for (aligner, opt_sequence) in aligners.iter().zip(sequences.iter()) {
            let best_mapping = if let Some(sequence) = opt_sequence {
                // cDNA mapper first
                let mappings = aligner.map(
                    sequence.as_bytes(),
                    output_cigar, output_md, max_frag_len, extra_flags.as_deref(), None
                )?;

                // here, we only allow Forward mappings; reads have been put onto the correct strand, so anything on Reverse is non-sense
                // usually, this means the forward mapping is also pretty terrible (if it exists)
                let mappings: Vec<minimap2::Mapping> = mappings.into_iter()
                    .filter(|m| m.strand == minimap2::Strand::Forward)
                    .collect();

                // we are scoring relative to the query (read) with an unmapped penalty
                let unmapped_from_target = false;
                let penalize_unmapped = true;
                let (best_mapping, _bs) = select_best_mapping(
                    &mappings, unmapped_from_target, penalize_unmapped, None
                );

                // return the best mapping
                best_mapping.cloned()
            } else {
                None
            };

            // add what we found for this haplotype aligner
            current_match.add_mapping(best_mapping)?;
        }

        // first, do all the debug and tracking of the comparisons
        let hla_mapping_stats = HlaMappingStats::from_mapping_stats(
            current_match.full_mapping_stats()[0].clone(),
            current_match.full_mapping_stats()[1].clone()
        );
        let match_mapping = current_match.full_mappings()[0].as_ref();
        let d_mapping = current_match.full_mappings()[1].as_ref();

        // check if we want to store the full info for this one
        let hla_star_allele = database.hla_sequences().get(hla_id).unwrap().star_allele().join(":");
        if all_hla_targets || cli_settings.debug_hla_targets.contains(&hla_star_allele) {
            // save this info using the star allele notation for simplicity
            read_mapping_stats.add_mapping(hla_star_allele, match_mapping, d_mapping)?;
        } else if cli_settings.debug_hla_targets.contains(hla_id) {
            // we need to save this info using the HLA ID
            read_mapping_stats.add_mapping(hla_id.clone(), match_mapping, d_mapping)?;
        }

        // now do the comparison to the current best, and overwrite if better
        if current_match.is_better_match(&best_match)? {
            // pull out cigars for debug
            let match_cigar = match_mapping.map(|bm| bm.alignment.as_ref().unwrap().cigar_str.as_ref().unwrap());
            let d_cigar = d_mapping.map(|bm| bm.alignment.as_ref().unwrap().cigar_str.as_ref().unwrap());

            debug!("{qname} {hla_id} -> new best {:?}, {}", hla_allele_def.star_allele(), hla_mapping_stats.score_string());
            debug!("\tcDNA cigar   -> {:?}", match_cigar);
            debug!("\tcDNA mapping -> {:?}", match_mapping);
            debug!("\tDNA cigar    -> {:?}", d_cigar);
            debug!("\tDNA mapping  -> {:?}", d_mapping);

            // the current is better, time to replace things
            best_match = current_match;
        }

        // now save the result for this ID
        ret.insert(hla_id.clone(), hla_mapping_stats);
    }

    // set the best match before returning
    let best_hla_id = best_match.haplotype().to_string();
    if !best_hla_id.is_empty() {
        let best_hla_star = database.hla_sequences().get(&best_hla_id).unwrap().star_allele().join(":");
        read_mapping_stats.set_best_match(best_hla_id, best_hla_star);
    } else {
        // default sets it to "" and ""; so no operations needed here
    }
    Ok((ret, read_mapping_stats))
}

/// Accepts a bam record and extracts just the parts that align to the cDNA for a gene. Additionally returns an offset for any prefix bases that are not found.
/// # Arguments
/// * `read` - the record to parse; sequence must be on forward strand relative to reference genome
/// * `database` - contains coordinates for splicing
/// * `gene_name` - the gene we are splicing
pub fn splice_read(read: &mut rust_htslib::bam::Record, database: &PgxDatabase, gene_name: &str) -> Result<(Vec<u8>, usize), Box<dyn std::error::Error>> {
    // get the DNA string out
    let read_sequence: Vec<u8> = read.seq().as_bytes();
    
    // we should have already called this, but do it again to be safe
    read.cache_cigar();
            
    //build a lookup from reference coordinate -> sequence coordinate
    let mut coordinate_lookup: HashMap<usize, usize> = Default::default();
    for bp in read.aligned_pairs() {
        let segment_index = bp[0] as usize;
        let ref_index = bp[1] as usize;
        coordinate_lookup.insert(ref_index, segment_index);
    }

    // splice in the coordinates from the lookup
    let mut offset = 0;
    let mut splice_segments: Vec<std::ops::Range<usize>> = vec![];
    for exon_coordinates in database.hla_config().gene_definition(gene_name).unwrap().exons().iter() {
        let mut first_exon_base = exon_coordinates.start() as usize;
        let mut last_exon_base = (exon_coordinates.end() - 1) as usize;
        
        // find the first base in this exon that the read maps to
        while !coordinate_lookup.contains_key(&first_exon_base) && first_exon_base <= last_exon_base {
            first_exon_base += 1;
        }
        
        // find the last base (inclusive) in this exon that the read maps to
        while !coordinate_lookup.contains_key(&last_exon_base) && first_exon_base <= last_exon_base {
            last_exon_base -= 1;
        }

        // add in the offset if we have not found any bases yet
        if splice_segments.is_empty() {
            offset += first_exon_base - exon_coordinates.start() as usize;
        }

        // if there is a non-zero number of bases, we will have a splice segment to add
        if first_exon_base <= last_exon_base {
            // these are both inclusive
            let first_read_base = *coordinate_lookup.get(&first_exon_base).unwrap();
            let last_read_base = *coordinate_lookup.get(&last_exon_base).unwrap();

            // add one to make this a normal range
            splice_segments.push(first_read_base..(last_read_base+1));
        }
    }

    // we don't want to oversplice, so set the first one to the read start and the last one to the read end
    // splice_segments[0].0 = 0;
    // splice_segments.last_mut().unwrap().1 = read_string.len();

    // now create the pre-spliced read using the segments we found
    let mut prespliced_read = vec![];
    for segment in splice_segments.into_iter() {
        prespliced_read.extend_from_slice(&read_sequence[segment]);
    }
    Ok((prespliced_read, offset))
}

/// This will look at a dual consensus result and return True if a hemizygous result is better.
/// # Arguments
/// * `preliminary_consensus` - the DualConsensus to check
/// * `cli_settings` - contains settings for scoring
/// * `normalized_coverage` - optional normalized coverage; if absent, the scoring here will not really be super helpful
fn is_hemizygous_better(preliminary_consensus: &DualConsensus, cli_settings: &DiplotypeSettings, normalized_coverage: Option<f64>) -> Result<bool, Box<dyn std::error::Error>> {
    // we do not want to run Consensus in its entirety, so lets use our DC result and see if we can make a good heuristic
    // - add the delta assuming EITHER con1 or con2 is the correct consensus
    // - look at the balance of the consensuses also

    // figure out the cost if either c1 or c2 was the single consensus
    let read_count = preliminary_consensus.is_consensus1().len();
    let min_ed = if preliminary_consensus.is_dual() {
        let mut c1_cost = 0;
        let mut c2_cost = 0;
        for (o1, o2) in preliminary_consensus.scores1().iter().zip(preliminary_consensus.scores2().iter()) {
            // at least one of these should be set
            assert!(o1.is_some() || o2.is_some());

            // get the scores; if no score is present, we can assume that it reached the dual_max_ed_delta + the other score
            let s1 = o1.unwrap_or(o2.unwrap_or_default() + cli_settings.dual_max_ed_delta);
            let s2 = o2.unwrap_or(o1.unwrap_or_default() + cli_settings.dual_max_ed_delta);
            let min_cost = s1.min(s2);

            // add the delta to each
            c1_cost += s1 - min_cost;
            c2_cost += s2 - min_cost;
        }

        // the haploid cost is basically the ED from each read to one of the two consensuses
        c1_cost.min(c2_cost)
    } else {
        // we did not find a dual consensus result, this is going to come down to coverage comparison
        0
    };
    let ln_ed_penalty = 2.0;
    let haploid_ed_cost = ln_ed_penalty * min_ed as f64;
    
    let nc_hap = normalized_coverage.unwrap_or(read_count as f64);
    let nc_dev = nc_hap * 0.1;
    let haploid_norm = Normal::new(nc_hap, nc_dev)?;
    let haploid_norm_cost = haploid_norm.ln_pdf(read_count as f64).abs();
    let haploid_cost = haploid_ed_cost + haploid_norm_cost;

    // diploid cost is based on our ratios; expecting 50:50; count the number on consensus 1
    let obs1 = preliminary_consensus.is_consensus1().iter()
        .filter(|&&x| x)
        .count();
    let obs2 = read_count - obs1;
    let diploid_binomial = Binomial::new(0.5, read_count as u64)?;
    let diploid_balance_penalty = 2.0;
    let diploid_balance_cost = if preliminary_consensus.is_dual() {
        diploid_balance_penalty * diploid_binomial.ln_pmf(obs1 as u64).abs()
    } else {
        // if the result is not dual, then this should be 0.0
        0.0
    };

    // also factor in the coverage expectation
    let nc_dip = 2.0 * normalized_coverage.unwrap_or(read_count as f64);
    let diploid_norm = Normal::new(nc_dip, nc_dev)?;
    let diploid_norm_cost = diploid_norm.ln_pdf(read_count as f64).abs();
    let diploid_cost = diploid_balance_cost + diploid_norm_cost;

    debug!("\tPreliminary is_dual: {}", preliminary_consensus.is_dual());
    debug!("\tHaploid cost: {haploid_cost} = {haploid_ed_cost} ({min_ed} ED) + {haploid_norm_cost} ({read_count} v. {nc_hap})");
    debug!("\tDiploid cost: {diploid_cost} = {diploid_balance_cost} ({obs1}/{obs2} reads) + {diploid_norm_cost} ({read_count} v. {nc_dip})");

    if haploid_cost < diploid_cost {
        // cheaper cost if haploid
        Ok(true)
    } else {
        // diploid cost is lower, so use it
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam::record::{Cigar, CigarString};
    use std::str::FromStr;

    use super::*;
    use crate::data_types::mapping::MappingScore;
    use crate::util::file_io::load_json;

    /*
    // writing a test for this will be complicated and may be better served with our end-to-end testing; punting for now
    #[test]
    fn test_diplotype_hla() {
        panic!("no impl");
    }
    */

    #[test]
    fn test_is_allowed_allele_def() {
        // base case, should be allowed
        let mut cli_settings = Default::default();
        let gene_name = "HLA-A";
        let hla_allele_def = HlaAlleleDefinition::new(
            "HLA1".to_string(), "A*01", Some("ACGT".to_string()), "AG".to_string()
        ).unwrap();
        assert!(is_allowed_allele_def(&hla_allele_def, gene_name, &cli_settings));

        // wrong gene, should not be allowed
        let hla_allele_def = HlaAlleleDefinition::new(
            "HLA1".to_string(), "B*01", Some("ACGT".to_string()), "AG".to_string()
        ).unwrap();
        assert!(!is_allowed_allele_def(&hla_allele_def, gene_name, &cli_settings));

        // missing DNA, should not be allowed if require DNA is enabled
        cli_settings.hla_require_dna = true;
        let hla_allele_def = HlaAlleleDefinition::new(
            "HLA1".to_string(), "A*01", None, "AG".to_string()
        ).unwrap();
        assert!(!is_allowed_allele_def(&hla_allele_def, gene_name, &cli_settings));

        // remove the requirement and it's allowed again
        cli_settings.hla_require_dna = false;
        let hla_allele_def = HlaAlleleDefinition::new(
            "HLA1".to_string(), "A*01", None, "AG".to_string()
        ).unwrap();
        assert!(is_allowed_allele_def(&hla_allele_def, gene_name, &cli_settings));
    }

    fn load_default_database() -> PgxDatabase {
        let database_fn = PathBuf::from_str("./test_data/HLA-faux/database.json").unwrap();
        let database = load_json(&database_fn).unwrap();
        database
    }

    #[test]
    fn test_reference_alleles() {
        // these are the reference alleles for each, so they should exactly match
        let genes = vec!["HLA-A".to_string(), "HLA-B".to_string()];
        let hla_ids = vec!["HLA:HLA00037".to_string(), "HLA:HLA00132".to_string()];
        let hla_stars = vec!["03:01:01:01".to_string(), "07:02:01:01".to_string()];
        let mapping_positions = vec![29942254, 31353362];
        let is_revcomp = vec![false, true];

        // load our proxy database that just has the reference alleles
        let database = load_default_database();
        let mut cli_settings: DiplotypeSettings = Default::default();
        cli_settings.hla_require_dna = false;
        cli_settings.disable_cdna_scoring = false;
        cli_settings.debug_hla_targets = hla_ids.clone();

        for (i, gene_name) in genes.iter().enumerate() {
            let hla_key = &hla_ids[i];
            let hla_star = &hla_stars[i];
            let mut test_sequence = database.hla_sequences().get(hla_key).expect("the HLA key is present")
                .dna_sequence().expect("the DNA sequence is present")
                .to_string();

            // HLA-B is on the opposite strand, so we need to rev-comp it to match the reference correctly; otherwise splicing gets jacked
            if is_revcomp[i] {
                test_sequence = test_sequence.chars()
                    .rev()
                    .map(|c| match c {
                        'A' => 'T',
                        'C' => 'G',
                        'G' => 'C',
                        'T' => 'A',
                        _ => panic!("unexpected {c}")
                    })
                    .collect::<String>();
            }
            
            // create a read that exactly matches the test sequence
            let mut read = rust_htslib::bam::Record::new();
            let cigar_string = CigarString(vec![Cigar::Match(test_sequence.len() as u32); 1]);
            read.set(
                "read_name".as_bytes(),
                Some(&cigar_string),
                test_sequence.as_bytes(),
                &vec![20; test_sequence.len()]
            );
            read.set_pos(mapping_positions[i] - 1);

            // score it
            let (all_scores, mapping_stats) = score_read(read, &database, gene_name, &cli_settings).unwrap();
            let best_score = mapping_stats.best_match_id().unwrap();
            assert_eq!(hla_key, best_score);
            let best_score_hla = mapping_stats.best_match_star().unwrap();
            assert_eq!(hla_star, best_score_hla);

            let key_score = all_scores.get(hla_key).unwrap();
            let cdna_len = database.hla_sequences().get(hla_key).unwrap().cdna_sequence().len();
            let dna_len = test_sequence.len();
            assert_eq!(key_score, &HlaMappingStats::new(
                // figure out the oddity here, it's that HLA-A*03:01:01:01 is the reference allele, this should work now
                Some(cdna_len), Some(0), Some(0),
                Some(dna_len), Some(0), Some(0)
            ));
        }
    }

    // make the settings match the default from clap
    fn get_default_cli_settings() -> DiplotypeSettings {
        let mut default_settings: DiplotypeSettings = Default::default();
        default_settings.max_error_rate = 0.05;
        default_settings.min_cdf = 0.001;
        default_settings
    }

    #[test]
    fn test_score_bad_read() {
        // nothing should match this read
        let mut record = rust_htslib::bam::record::Record::new();
        record.set("test".as_bytes(), None, "ACGT".as_bytes(), &[255; 4]);

        // verify that we can successfully run this read, even though it's going to get thrown away
        let database = load_default_database();
        let gene_name = "HLA-A";
        let mut cli_settings = get_default_cli_settings();
        cli_settings.disable_cdna_scoring = true; // this read will never work with cDNA
        let (hash_scores, best_result) = score_read(record, &database, gene_name, &cli_settings).unwrap();
        
        // make sure there is no best score
        assert!(best_result.best_match_id().is_none());
        assert!(best_result.best_match_star().is_none());
        
        // we get one mapping stat per entry in the DB
        assert_eq!(best_result.mapping_stats().len(), 1);

        // finally, make sure all voting is marked as terrible, this will lead to getting discarded later
        for (_key, score) in hash_scores.iter() {
            // make sure both scores are the worst
            assert_eq!(score.mapping_score().cdna_score(), MappingScore::worst_score());
            assert_eq!(score.mapping_score().dna_score(), MappingScore::worst_score());
        }
    }

    /// Wrapper for a running test on whether a DualConsensus passes just based on counts.
    fn run_passing_test(c1: usize, c2: usize) -> bool {
        use waffle_con::consensus::Consensus;
        let mut cli_settings: DiplotypeSettings = Default::default();
        cli_settings.min_cdf = 0.001;
        cli_settings.min_consensus_fraction = 0.10;
        cli_settings.expected_maf = 0.5;
        
        let total = c1+c2;
        let mut is_consensus1 = vec![true; c1];
        is_consensus1.extend(vec![false; c2]);

        let dual_consensus = DualConsensus::new(
            // these do not matter
            Consensus::new(vec![], waffle_con::cdwfa_config::ConsensusCost::L1Distance, vec![]),
            Some(Consensus::new(vec![], waffle_con::cdwfa_config::ConsensusCost::L1Distance, vec![])),
            // this is all that matters for the test
            is_consensus1,
            // also don't matter
            vec![None; total],
            vec![None; total],
        ).unwrap();
        is_passing_dual(&dual_consensus, &cli_settings)
    }

    #[test]
    fn test_is_passing_dual() {
        // imbalanced, should fail
        assert!(!run_passing_test(3, 20));
        assert!(!run_passing_test(20, 3));

        // close enough to pass
        assert!(run_passing_test(10, 20));
        assert!(run_passing_test(20, 10));
    }

    fn run_is_hemi_better_test(c1: usize, c2: usize, normalized_coverage: f64, score_delta: usize) -> Result<bool, Box<dyn std::error::Error>> {
        // for the settings, we just need to set the dual max ED
        let mut cli_settings: DiplotypeSettings = Default::default();
        cli_settings.dual_max_ed_delta = 20;

        // dual consensus if we have counts for both
        assert!(c1 != 0);
        let is_dual = c2 != 0;

        // make this vec based on consensus counts provided
        let mut is_consensus1 = vec![true; c1];
        is_consensus1.extend(vec![false; c2]);

        // scores1 has exact match for consensus 1, and then delta for consensus 2
        let mut scores1 = vec![Some(0); c1];
        scores1.extend(vec![Some(score_delta); c2]);

        // scores 2 is the reverse of scores1
        let mut scores2 = vec![Some(score_delta); c1];
        scores2.extend(vec![Some(0); c2]);

        // this is the trickier part, building a DC based on provided params
        // the consensuses do not matter here
        let consensus1: Consensus = Default::default();
        let consensus2: Option<Consensus> = if is_dual { Some(Default::default()) } else { None };
        let consensus = DualConsensus::new(
            consensus1,
            consensus2,
            is_consensus1,
            scores1,
            scores2
        )?;

        is_hemizygous_better(&consensus, &cli_settings, Some(normalized_coverage))
    }

    #[test]
    fn test_is_hemizygous_better() {
        // homozygous result that should be hemizygous
        assert!(run_is_hemi_better_test(20, 0, 20.0, 1).unwrap());

        // homozygous result that should be homozygous
        assert!(!run_is_hemi_better_test(40, 0, 20.0, 1).unwrap());

        // heterozygous result that should be hemizygous; it is skewed with low delta
        assert!(run_is_hemi_better_test(18, 2, 20.0, 1).unwrap());

        // heterozygous result that should be heterozygous
        assert!(!run_is_hemi_better_test(18, 17, 20.0, 1).unwrap());

        // lower quality heterozygous with heavy skew, but signal is still strong enough to keep it
        assert!(!run_is_hemi_better_test(15, 6, 20.0, 20).unwrap());
    }
}