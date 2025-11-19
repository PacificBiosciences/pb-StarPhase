
use log::{debug, error, info, warn};
use itertools::Itertools;
use rust_htslib::bam::Read;
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use simple_error::bail;
use std::collections::BTreeMap;
use std::path::PathBuf;
use waffle_con::cdwfa_config::{CdwfaConfig, CdwfaConfigBuilder};
use waffle_con::consensus::{Consensus, ConsensusDWFA};
use waffle_con::multi_consensus::MultiConsensus;
use waffle_con::priority_consensus::{PriorityConsensus, PriorityConsensusDWFA};

use crate::cli::diplotype::DiplotypeSettings;
use crate::cyp2d6::chaining::{find_best_chain_pair, weight_sequence, ChainPenalties, SequenceWeights};
use crate::cyp2d6::debug::DeeplotypeDebug;
use crate::cyp2d6::haplotyper::{AlleleMapping, Cyp2d6Extractor};
use crate::cyp2d6::region::{Cyp2d6DetailLevel, Cyp2d6Region};
use crate::cyp2d6::region_label::{Cyp2d6RegionLabel, Cyp2d6RegionType};
use crate::cyp2d6::vcf_writer::write_cyp2d6_vcf;
use crate::cyp2d6::visualization::create_custom_cyp2d6_reference;
use crate::data_types::pgx_diplotype::Diplotype;
use crate::data_types::starphase_json::{PgxGeneDetails, PgxMultiMappingDetails};
use crate::database::pgx_database::PgxDatabase;
use crate::util::file_io::save_fasta;
use crate::util::homopolymers::hpc_with_guide;
use crate::visualization::debug_bam_writer::{unmapped_record, DebugBamWriter};
use crate::visualization::igv_session_writer::IgvSessionWriter;

/// This is the main function to call for CYP2D6 diplotyping from a BAM file.
/// # Arguments
/// * `database` - the pre-loaded database
/// * `bam_filenames` - list of BAM files containing reads to scan
/// * `reference_genome` - already loaded reference genome
/// * `cli_settings` - settings for diplotyping
/// # Errors
/// * if we cannot open or parse a BAM file correctly
pub fn diplotype_cyp2d6(
    database: &PgxDatabase,
    bam_filenames: &[PathBuf], reference_genome: &ReferenceGenome,
    debug_bam_writer: Option<&mut DebugBamWriter>,
    cli_settings: &DiplotypeSettings
) -> Result<PgxGeneDetails, Box<dyn std::error::Error>> {
    info!("Solving CYP2D6...");

    // load the D6 typing engine, it's complicated so we moved it into cyp2d6_typer.rs
    let d6_typer = Cyp2d6Extractor::new(
        database,
        reference_genome
    )?;

    // 4 - extract all the reads from our BAM file
    // this will store the records of interest
    let mut read_collection: HashMap<String, rust_htslib::bam::Record> = Default::default();
    // these are the corresponding sequences, which need to be cached until we are done
    let mut read_sequences: HashMap<String, String> = Default::default();

    // prep all the bam readers
    let mut bam_readers: Vec<rust_htslib::bam::IndexedReader> = vec![];
    for bam_fn in  bam_filenames.iter() {
        let mut b = rust_htslib::bam::IndexedReader::from_path(bam_fn)?;
        b.set_reference(reference_genome.filename())?;
        bam_readers.push(b);
    }

    // get the extraction region from CYP2D6
    let bam_region = database.cyp2d6_config().extraction_region();
    debug!("Parsing reads in region: {bam_region:?}");

    // iterate over each bam, and fetch the reads
    for (bam_index, bam) in bam_readers.iter_mut().enumerate() {
        // bam.fetch(coordinate.fetch_definition())?;
        match bam.fetch(bam_region.fetch_definition()) {
            Ok(()) => {},
            Err(e) => {
                let filename = &bam_filenames[bam_index];
                warn!("Received error \"{e}\" while fetching {bam_region} in {filename:?}, assuming no reads for region.");
                continue;
            }
        };
        
        for read_entry in bam.records() {
            let mut read = read_entry.unwrap();
            
            // make sure we do not do the same read twice
            let qname: String = std::str::from_utf8(read.qname())?.to_string();
            if read_collection.contains_key(&qname) {
                continue;
            }
            
            /*
            // TODO: thus far, we haven't applied this filter; leaving it in place in case we want to think about this in the future
            // make sure we care about the alignment
            if filter_out_alignment_record(&read, min_mapq) {
                continue;
            }
            */
            
            //build out the cigar info
            read.cache_cigar();

            let sequence = String::from_utf8(read.seq().as_bytes())?;
            
            // insert, make sure it was not already inserted
            assert!(read_collection.insert(qname.clone(), read).is_none());
            // assert!(read_sequences.insert(qname, filtered_sequence).is_none());
            assert!(read_sequences.insert(qname, sequence).is_none());
        }
    }

    // These are candidate for CLI parameters, but I think keep them hidden until we have an obvious reason not to do so
    // constants for removing matches that are too short from chaining and/or consensus steps
    let min_chain_frac = 0.5; // requires this fraction to go into the chaining step
    let min_consensus_frac = 0.5; // requires this fraction to go into the consensus step; we usually want this to be relatively high
    let min_typing_frac = 0.9; // requires high fraction of the allele to get a type; otherwise assigned to "UNKNOWN"

    // derive these, which are what is actually used below
    let max_missing_chain_frac = 1.0 - min_chain_frac;
    let max_missing_consensus_frac = 1.0 - min_consensus_frac;
    let max_missing_typing_frac = 1.0 - min_typing_frac;
    assert!(max_missing_chain_frac >= max_missing_consensus_frac); // make sure we never break this assumption

    // 4? - identify all putative D6, D7, hybrid, and deletion regions
    // we want these region results in a BTree because we need to traverse them consistently (i.e., sorted order) so we can map the downstream results.
    let mut regions_of_interest: BTreeMap<String, Vec<AlleleMapping>> = Default::default();
    for (read_id, record) in read_collection.iter() {
        debug!("Searching {read_id} at {}", record.pos());
        // let read_sequence = String::from_utf8(record.seq().as_bytes())?;
        let read_sequence = read_sequences.get(read_id).unwrap(); // we want the filtered sequence, not the original
        let penalize_unmapped = false; // penalize is fine here, we just want the locations at this point
        let initial_regions = d6_typer.find_base_type_in_sequence(
            read_sequence,
            penalize_unmapped,
            max_missing_chain_frac
        )?;
        debug!("Found {} regions of interest.", initial_regions.len());
        regions_of_interest.insert(read_id.clone(), initial_regions);
    }

    // TODO: many of these should likely be CLI options - depth and AF will be relevant when we get to targeted
    //       relevant question: what does `min_count` mean when our reads don't fully span a region anymore, should we weight by fraction covered?
    //           I think this will be more important for chaining that consensus 
    // prep the consensus algo
    let offset_compare_length = 100;
    let offset_window = 50; // we want it to look +-50 bp, but the config only lets us look before; so we have to shift things
    let config_offset_window = 2 * offset_window;
    let consensus_config = CdwfaConfigBuilder::default()
        .wildcard(Some(b'*'))
        .min_count(cli_settings.min_consensus_count)
        .min_af(cli_settings.min_consensus_fraction)
        .dual_max_ed_delta(cli_settings.dual_max_ed_delta)
        .allow_early_termination(true)
        .weighted_by_ed(false) // currently, I'm not convinced on either approach
        .consensus_cost(waffle_con::cdwfa_config::ConsensusCost::L1Distance)
        .max_queue_size(20)
        .max_capacity_per_size(10)
        .offset_window(config_offset_window)
        .offset_compare_length(offset_compare_length)
        .build()?;

    // now prep the priority consensus runner - we have HPC as priority, then full length
    let mut consensus_dwfa = PriorityConsensusDWFA::with_config(consensus_config.clone())?;
    
    // load all the corresponding sequences for the regions of interest into the consensus generator; these are in a fixed order
    let mut raw_sequences: Vec<&str> = vec![]; // stores just the part of the read matching the region
    let mut shifted_sequences: Vec<String> = vec![]; // stores the part of the read matching the region AND any prefix wildcards to shift it
    let mut hpc_shifted_sequences: Vec<String> = vec![]; // stores the HPC string AND any prefix wildcards to HPC shift it; note, this prefix will be shorter than above due to re-anchoring
    let mut base_offsets: Vec<usize> = vec![];
    let mut hpc_offsets: Vec<usize> = vec![];
    
    let mut sequence_ids: Vec<String> = vec![]; // an ID to reference the read, region, and allele name
    let mut flattened_regions_of_interest: Vec<(String, AlleleMapping)> = Default::default(); // flattened version of the regions_of_interest for quick lookup
    
    for (read_id, regions) in regions_of_interest.iter() {
        let read_sequence = read_sequences.get(read_id).unwrap().as_bytes();
        for region in regions.iter() {
            if region.mapping_stats().custom_score(true).score() > max_missing_consensus_frac {
                debug!("Ignoring {read_id}-{:?} for consensus generation: {}", region.region(), region.mapping_stats().custom_score_string(true));
                continue;
            }

            // get the prefix and postfix components
            let prefix = String::from_utf8(vec![b'*'; region.mapping_stats().clipped_start().unwrap()])?;
            // let postfix = String::from_utf8(vec![b'*'; region.mapping_stats.clipped_end().unwrap()])?;
            let seq = std::str::from_utf8(&read_sequence[region.region().clone()])?;
            
            // this is just the raw matching sequence
            raw_sequences.push(seq);

            // generate the full target sequence
            let full_sequence = prefix.clone() + seq;// + &postfix;
            shifted_sequences.push(full_sequence);
            base_offsets.push( if prefix.is_empty() {
                0    
            } else {
                prefix.len() + offset_window
            }); // add in the offset window buffer since we want it to look +- the approx offset

            let guide_id = region.allele_label();
            let guide_seq = d6_typer.get_allele(guide_id).unwrap();
            let (hpc_sequence, prefix_offset) = hpc_with_guide(seq, guide_seq, prefix.len())?;
            hpc_shifted_sequences.push(hpc_sequence);
            hpc_offsets.push(if prefix_offset == 0 {
                0
            } else {
                prefix_offset + offset_window
            }); // add in the offset window buffer since we want it to look +- the approx offset

            // add a sequence ID for easy tracking
            sequence_ids.push(format!("{read_id}_{}_{}_{}", region.region().start, region.region().end, region.allele_label()));

            // flatten the regions also
            flattened_regions_of_interest.push((read_id.clone(), region.clone()));
        }
    }

    debug!("sequence_ids: {sequence_ids:?}");

    // now add all of them to the DWFA
    // for ((hpc_ss, ss), (_fr_read_id, fr_region)) in (hpc_shifted_sequences.iter().zip(shifted_sequences.iter())).zip(flattened_regions_of_interest.iter()) {
    for (seq_index, hpc_ss) in hpc_shifted_sequences.iter().enumerate() {
        // let ss = &shifted_sequences[seq_index];
        let raw_seq = raw_sequences[seq_index];
        let fr_region = &flattened_regions_of_interest[seq_index].1;

        // this basically forces all *5 into a separate grouping from the get-go
        // we needed this because have three hyper-diverse alleles in one led to death spiral in one sample
        let seed = match fr_region.allele_label().region_type() {
            Cyp2d6RegionType::Cyp2d6Deletion => Some(0),
            Cyp2d6RegionType::Rep6 => Some(1),
            Cyp2d6RegionType::Rep7 => Some(2),
            Cyp2d6RegionType::Spacer => Some(3),
            Cyp2d6RegionType::LinkRegion => Some(4),
            _ => None
        };

        // we chain first by HPC and second by the full sequence
        let sequence_chain = vec![
            hpc_ss.as_bytes(),
            raw_seq.as_bytes()
        ];

        // build the offset chains, set anything with offset 0 to just be auto-start
        let hpc_off = if hpc_offsets[seq_index] == 0 {
            None
        } else {
            Some(hpc_offsets[seq_index])
        };
        let base_off = if base_offsets[seq_index] == 0 {
            None
        } else {
            Some(base_offsets[seq_index])
        };

        let offset_chain = vec![
            hpc_off, base_off
        ];
        consensus_dwfa.add_seeded_sequence_chain(sequence_chain, offset_chain, seed)?;
    }

    // make sure we found some sequence
    if consensus_dwfa.sequences().is_empty() {
        warn!("No reads found for CYP2D6 consensus generation.");

        // finally lets build our results
        let diplotypes = vec![Diplotype::new("NO_READS", "NO_READS")];
        debug!("Full diplotype for CYP2D6 => \"{}\"", diplotypes[0].diplotype());

        let pgx_gene_details = PgxGeneDetails::new_from_multi_mappings(
            diplotypes, None, None, vec![]
        )?;
        return Ok(pgx_gene_details);
    }

    // set to true if you need to debug and print out some sequences in a fasta-like system
    let debug_sequences: bool = false;

    // now solve the core consensus
    let raw_consensus_result = consensus_dwfa.consensus()?;
    debug!("Found {} raw consensus sequences", raw_consensus_result.consensuses().len());
    if debug_sequences {
        for (i, c) in raw_consensus_result.consensuses().iter().enumerate() {
            println!(">raw_con_{i}");
            println!("{}", std::str::from_utf8(c[0].sequence()).unwrap());
        }
    }

    debug!("Found {} expanded consensus sequence", raw_consensus_result.consensuses().len());
    if debug_sequences {
        for (i, c) in raw_consensus_result.consensuses().iter().enumerate() {
            println!(">expanded_con_{i}");
            println!("{}", std::str::from_utf8(c[1].sequence()).unwrap());
        }
    }

    // nowe we need to collapse those that are identical at HPC and map to same allele
    let consensus_result = merge_consensus_results(
        // &shifted_sequences, 
        &raw_sequences,
        &base_offsets,
        &consensus_config,
        &raw_consensus_result,
        &d6_typer,
        max_missing_typing_frac
    )?;

    debug!("Found {} final consensus sequences", consensus_result.consensuses().len());
    if debug_sequences {
        for (i, c) in consensus_result.consensuses().iter().enumerate() {
            println!(">final_con_{i}");
            println!("{}", std::str::from_utf8(c.sequence()).unwrap());
        }
    }

    /*
    // this is a test block that will compare HPC sequences
    for (i, c) in consensus_result.consensuses().iter().enumerate() {
        let hpc1 = hpc(std::str::from_utf8(c.sequence()).unwrap())?;
        for (j, c2) in consensus_result.consensuses().iter().enumerate().skip(i+1) {
            let ed = waffle_con::sequence_alignment::wfa_ed(c.sequence(), c2.sequence());
            let hpc2 = hpc(std::str::from_utf8(c2.sequence()).unwrap())?;
            let hpc_ed = waffle_con::sequence_alignment::wfa_ed(hpc1.as_bytes(), hpc2.as_bytes());
            println!("{i} {j} => {ed} {hpc_ed}");
        }
    }
    todo!("inspect above");
    */

    debug!("Sequence to consensus: {:?}", consensus_result.sequence_indices());
    assert_eq!(raw_sequences.len(), consensus_result.sequence_indices().len());

    // figure out what each consensus haplotype is
    let mut hap_regions = vec![];
    let mut sequences_labeled: HashSet<String> = Default::default();
    for (i, c) in consensus_result.consensuses().iter().enumerate() {
        let matches = consensus_result.sequence_indices().iter()
            .filter(|&&con_index| con_index == i)
            .count();

        let sequence_to_type = std::str::from_utf8(c.sequence()).unwrap()
            .trim_matches('*');

        debug!("Typing consensus #{i} with {matches} matches, {} wildcards trimmed", c.sequence().len() - sequence_to_type.len());

        // we are *expecting* full length sequences in our consensus
        let force_assignment = true; // at this point, we need a label even if there is some ambiguity
        let mut hap_region = match d6_typer.find_full_type_in_sequence(sequence_to_type, max_missing_typing_frac, force_assignment) {
            Ok(hl) => hl,
            Err(e) => {
                let unknown = Cyp2d6RegionLabel::new_unknown();
                error!("Error while typing consensus #{i}, setting to {unknown}.");
                error!("Typing error: {e}");
                Cyp2d6Region::new(unknown, None)
            }
        };

        // do a check to make sure we do not already have this sequence
        if sequences_labeled.contains(sequence_to_type) {
            // we somehow created two identical alleles and split them into two groups
            // label this second one as a FalseAllele so we ignore it later
            debug!("Detected duplicate allele in consensus {i}, marking as FalseAllele");
            hap_region.mark_false_allele();
        } else {
            // new sequence, no changes needed
            sequences_labeled.insert(sequence_to_type.to_string());
        };

        // set the region unique ID right before we push it into the list
        hap_region.set_unique_id(hap_regions.len());

        debug!("hap_label = \"{}\"", hap_region);
        debug!("deep_label = {:?}", hap_region.deep_label());
        hap_regions.push(hap_region);
    }

    // make the output BAM if requested
    if let Some(dbw) = debug_bam_writer {
        // create an unmapped record for each sequence that went into consensus
        let mut unmapped_records = vec![];
        for (seq_id, (raw_seq, &phase_id)) in raw_sequences.iter().zip(consensus_result.sequence_indices().iter()).enumerate() {
            // TODO: long-term, we probably want to trace this to a read if possible
            let qname = format!("seq_{seq_id}");
            let sequence = raw_seq.to_string(); // rev-comp not necessary because these are sourced from BAM records
            let tags = [
                ("HP".to_string(), hap_regions[phase_id].index_label())
            ].into_iter().collect();

            // add the record
            match unmapped_record(&qname, &sequence, &tags) {
                Ok(umr) => {
                    unmapped_records.push(umr);
                },
                Err(e) => {
                    error!("Error while creating unmapped record: {e}");
                }
            };
        }

        // this does all the alignment work for us
        match dbw.map_records_to_region(&unmapped_records, &bam_region) {
            Ok(()) => {},
            Err(e) => {
                error!("Error while mappings records to debug BAM: {e}");
            }
        };
    }

    // write the VCF file if requested
    if let Some(debug_folder) = cli_settings.debug_folder.as_ref() {
        let vcf_fn = debug_folder.join("cyp2d6_alleles.vcf.gz");
        debug!("Writing CYP2D6 diplotypes to {vcf_fn:?}");
        match write_cyp2d6_vcf(&hap_regions, &vcf_fn, d6_typer.loaded_variants()) {
            Ok(()) => {},
            Err(e) => {
                bail!("Error while writing CYP2D6 VCF: {e}");
            }
        };

        // save the index file also (.tbi)
        match hiphase::writers::vcf_util::build_bcf_index(vcf_fn, None, 1, true) {
            Ok(()) => {},
            Err(e) => {
                bail!("Error while building CYP2D6 VCF index: {e}");
            }
        };
    }

    // build up all the chains
    let mut qname_chains: BTreeMap<String, Vec<Vec<usize>>> = Default::default();
    let mut unique_chains: HashSet<Vec<usize>> = Default::default();
    let mut best_allele_mapping_counts: Vec<u64> = vec![0; hap_regions.len()];
    let mut qname_chain_scores: BTreeMap<String, Vec<SequenceWeights>> = Default::default();
    let mut multi_mapping_details = vec![];
    for (read_id, regions) in regions_of_interest.iter() {
        if regions.is_empty() {
            continue;
        }
        debug!("Labeling chains in {read_id}...");
        let read_sequence = read_sequences.get(read_id).unwrap().as_bytes();
        let mut putative_chains: Vec<Vec<usize>> = vec![vec![]];
        let mut weighted_chains: Vec<SequenceWeights> = vec![];
        for (region_index, region) in regions.iter().enumerate() {
            debug!("\tScanning {:?}", region.region());

            // pull the sequence
            let seq = std::str::from_utf8(&read_sequence[region.region().clone()])?;
            
            // score this sequence against each consensus
            let weighted_scores = weight_sequence(seq, &consensus_result, &hap_regions)?;
            if weighted_scores.is_empty() {
                // all mappings were bad, we should skip this one
                // this should only happen at the edges, print a warning if it doesn't
                if region_index != 0 && region_index != regions.len()-1 {
                    warn!("\tRemoved putative allele mid-read due to no good matches, chaining may be impacted.");
                } else {
                    debug!("\tRemoved putative allele at start/end due to no good matches.");
                }
                continue;
            }

            let min_ed = weighted_scores.iter()
                .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
                .unwrap().0;
            let num_minimum = weighted_scores.iter()
                .filter(|a| a.0 == min_ed)
                .count();

            // generate new chains by extending the current ones
            let mut new_pc = vec![];
            for pc in putative_chains.into_iter() {
                for (ci, score) in weighted_scores.iter().enumerate() {
                    if score.0 == min_ed {
                        let mut new_chain = pc.clone();
                        new_chain.push(ci);
                        new_pc.push(new_chain);

                        if num_minimum == 1 {
                            // this is a unique best match
                            best_allele_mapping_counts[ci] += 1;
                        }
                    }
                }
            }

            // now overwrite the original
            putative_chains = new_pc;
            
            // save the scores
            weighted_chains.push(weighted_scores);
        }

        if putative_chains.is_empty() || (putative_chains.len() == 1 && putative_chains[0].is_empty()) {
            debug!("\tNo chains found.");
        } else {
            if putative_chains.len() == 1 {
                if putative_chains[0].len() > 1 {
                    // we only add those with at least one pair to this list
                    // TODO: do we need to parse out if there are 3+ alleles into pairs?
                    unique_chains.insert(putative_chains[0].clone());

                    if putative_chains[0].len() > 2 {
                        for i in 0..(putative_chains[0].len()-1) {
                            let pair = vec![putative_chains[0][i], putative_chains[0][i+1]];
                            unique_chains.insert(pair);
                        }
                    }
                }
                debug!("\tMost likely chain: {:?}", putative_chains[0]);
            } else {
                debug!("\tAmbig chains found: {:?}", putative_chains);
            }

            qname_chains.insert(read_id.clone(), putative_chains.clone());
            qname_chain_scores.insert(read_id.clone(), weighted_chains);
        }
    }

    // remove any chains with non-unique sub-alleles
    for (_qname, chain_set) in qname_chains.iter_mut() {
        let new_chain_set: Vec<Vec<usize>> = chain_set.iter()
            .filter(|chain| {
                chain.iter()
                    .all(|&c_index| best_allele_mapping_counts[c_index] > 0)
            })
            .cloned()
            .collect();
        
        if new_chain_set.is_empty() {
            panic!("chain collapse: {chain_set:?} => {new_chain_set:?}");
        }

        *chain_set = new_chain_set;

        // make sure we didn't somehow remove all options
        assert!(!chain_set.is_empty());
    }

    // count the chain frequencies
    let mut single_frequency: BTreeMap<usize, f64> = Default::default();
    let mut chain_frequency: BTreeMap<Vec<usize>, f64> = Default::default();

    // also count the ambiguous ones, but these are down-weighted by the number of candidates
    for (qname, chain_set) in qname_chains.iter() {
        let weight = 1.0 / chain_set.len() as f64;
        for chain in chain_set.iter() {
            // increment the chain
            let entry = chain_frequency.entry(chain.clone()).or_default();
            *entry += weight;

            // increment each singleton
            for &c in chain.iter() {
                let entry = single_frequency.entry(c).or_default();
                *entry += weight;
            }
        }

        if chain_set.len() == 1 {
            for (&consensus_index, region) in chain_set[0].iter().zip(regions_of_interest.get(qname).unwrap().iter()) {
                multi_mapping_details.push(PgxMultiMappingDetails::new(
                    qname.clone(),
                    region.region().clone(),
                    consensus_index,
                    hap_regions[consensus_index].index_label()
                ));
            }
        }
    }

    // debug for unique counts
    debug!("Uniquely assigned table:");
    for (con_index, &unique_count) in best_allele_mapping_counts.iter().enumerate() {
        let label = hap_regions[con_index].label();

        // if there are no unique reads; THEN we mark this as a FalseAllele
        if unique_count == 0 &&
            label.region_type() != Cyp2d6RegionType::Unknown &&
            label.region_type() != Cyp2d6RegionType::FalseAllele {
            // this is a false allele, nothing is uniquely mapping to it; retain the original label as the subtype
            hap_regions[con_index].mark_false_allele();
        }
        debug!("\t{} => {unique_count}", hap_regions[con_index].index_label());
    }

    // debug output for the table
    debug!("Allele count table:");
    for (&con_index, count) in single_frequency.iter() {
        debug!("\t{con_index}_{} => {count}", hap_regions[con_index]);
    }

    // print the uniquely assigned chains as well
    debug!("Unique chains:");
    for chain in unique_chains.iter().sorted() {
        let string_form: Vec<String> = chain.iter()
            .map(|&c_index| format!("{c_index}_{}", hap_regions[c_index]))
            .collect();
        debug!("\t{string_form:?}");
    }

    // debug output for the table
    debug!("Chain count table:");
    for (chain, count) in chain_frequency.iter() {
        let string_form: Vec<String> = chain.iter()
            .map(|&c_index| hap_regions[c_index].index_label())
            .collect();
        debug!("\t{string_form:?} => {count}")
    }

    if let Some(debug_folder) = cli_settings.debug_folder.as_ref() {
        // make the output link graph if we have a debug folder
        let out_graph_fn = debug_folder.join("cyp2d6_link_graph.svg");
        debug!("Generating CYP2D6 graph at {out_graph_fn:?}");
        if let Err(e) = crate::cyp2d6::visualization::generate_debug_graph(&hap_regions, &chain_frequency, &out_graph_fn) {
            error!("Error while generating CYP2D6 debug graph: {e}");
        }

        // make the fasta output as well
        let consensus_fn = debug_folder.join("consensus_CYP2D6.fa");
        debug!("Saving consensus for CYP2D6 to {consensus_fn:?}");
        let mut consensus_map: BTreeMap<String, String> = Default::default();
        for (region, consensus) in hap_regions.iter().zip(consensus_result.consensuses().iter()) {
            let k = region.index_label();
            let v = std::str::from_utf8(consensus.sequence())?.to_string();
            consensus_map.insert(k, v);
        }
        save_fasta(&consensus_map, &consensus_fn)?;
    }

    // parameters that control chaining
    let infer_connections = cli_settings.infer_connections;
    let normalize_all_alleles = !cli_settings.normalize_d6_only;
    let penalties: ChainPenalties = Default::default();
    let ignore_chain_label_limits = false; // this should always be false in prod; true is just for testing
    let (best_result, chain_warnings) = find_best_chain_pair(
        database.cyp2d6_config(),
        &qname_chains, &qname_chain_scores, &hap_regions,
        infer_connections, normalize_all_alleles, penalties, ignore_chain_label_limits
    )?;
    if !chain_warnings.is_empty() {
        warn!("Chain warnings: {chain_warnings:?}");
    }
    
    if best_result.len() != 2 {
        bail!("best_result has non-2 length: {best_result:?}");
    }

    debug!("Best_results:");
    for chain in best_result.iter() {
        let string_form = chain.iter()
            .map(|&c_index| hap_regions[c_index].index_label())
            .join(" -> ");
        debug!("\t{string_form}");
    }

    // used in the return generation, but possible in the debug folder construction also
    let cyp_translate = d6_typer.cyp2d6_config().cyp_translate();

    if let Some(debug_folder) = cli_settings.debug_folder.as_ref() {
        // we have all the data to build a custom session file now; first, we need to build our custom reference genome
        match create_custom_cyp2d6_reference(
            reference_genome, database,
            &consensus_result, &hap_regions, &best_result
        ) {
            Ok(cust_ref) => {
                // collect all our reads for re-mapping into the special D6 regions
                let all_records = read_collection.into_iter()
                    .filter_map(|(k, v)| {
                        if !regions_of_interest.get(&k).unwrap().is_empty() {
                            // only keep a read if we found something meaningful inside it
                            Some(v)
                        } else {
                            None
                        }
                    })
                    .collect();

                // finally, put it all together
                let session_folder = debug_folder.join("cyp2d6_igv_custom");
                let mut session_writer = IgvSessionWriter::new(session_folder, false);
                match session_writer.add_custom_region(
                    cust_ref.contig_name.clone(),
                    &cust_ref.sequence,
                    &cust_ref.regions,
                    all_records
                ) {
                    Ok(()) => {
                        // TODO: move this outside of here so it is shared with HLA
                        // region added fine, write the session
                        if let Err(e) = session_writer.write_session() {
                            error!("Error while writing custom session file: {e}");
                        }
                    },
                    Err(e) => {
                        error!("Error while adding a CYP2D6 custom region: {e}");
                    }
                }
            },
            Err(e) => {
                error!("Error while creating custom CYP2D6 reference file: {e}");
            }
        }

        // also create a CYP2D6 deep haplotype output JSON
        let allele_fn = debug_folder.join("cyp2d6_alleles.json");
        debug!("Saving CYP2D6 alleles to {:?}", allele_fn);
        let allele_stats = DeeplotypeDebug::new(&best_result, &hap_regions, cyp_translate);
        crate::util::file_io::save_json(&allele_stats, &allele_fn)?;
    }

    // finally lets build our results
    let hap1 = convert_chain_to_hap(&best_result[0], &hap_regions, Cyp2d6DetailLevel::DeepAlleles, cyp_translate);
    let hap2 = convert_chain_to_hap(&best_result[1], &hap_regions, Cyp2d6DetailLevel::DeepAlleles, cyp_translate);
    let deeplotype = crate::data_types::pgx_diplotype::InexactDiplotype::new_diplotype_only(
        Diplotype::new(&hap1, &hap2)
    );
    debug!("Full inexact diplotype for CYP2D6 => \"{}\"", deeplotype.basic_diplotype().diplotype());

    // finally lets build our results
    let hap1 = convert_chain_to_hap(&best_result[0], &hap_regions, Cyp2d6DetailLevel::SubAlleles, cyp_translate);
    let hap2 = convert_chain_to_hap(&best_result[1], &hap_regions, Cyp2d6DetailLevel::SubAlleles, cyp_translate);
    let diplotypes = vec![Diplotype::new(&hap1, &hap2)];
    debug!("Full diplotype for CYP2D6 => \"{}\"", diplotypes[0].diplotype());

    let hap1_collapsed = convert_chain_to_hap(&best_result[0], &hap_regions, Cyp2d6DetailLevel::CoreAlleles, cyp_translate);
    let hap2_collapsed = convert_chain_to_hap(&best_result[1], &hap_regions, Cyp2d6DetailLevel::CoreAlleles, cyp_translate);
    let diplotypes_collapsed = vec![Diplotype::new(&hap1_collapsed, &hap2_collapsed)];
    debug!("Simple diplotype for CYP2D6 => \"{}\"", diplotypes_collapsed[0].diplotype());

    // build the PGx details for D6
    // TODO: this is currently only capturing mapping info at the consensus stage; ultimately we would want to know some details from the chaining
    //       and also probably the full length D6 sequences somehow
    //       additionally, the ambig_chains is getting split out above; we can avoid the split by checking for .len() == 1 and then keep reads in order
    //       ideally, this would let us follow a read through the *whole* process and create a final BAM as well instead of the current intermediate
    let pgx_gene_details = PgxGeneDetails::new_from_multi_mappings(
        diplotypes,
        Some(diplotypes_collapsed),
        Some(vec![deeplotype]),
        multi_mapping_details
    )?;
    Ok(pgx_gene_details)
}

/// This will take a priority based consensus and collapse HPC & assignment identical alleles into single entries.
/// # Arguments
/// * `sequences` - the full length sequences that built the consensus, should correspond to the second entries in the PriorityConsensus
/// * `cdwfa_config` - configuration used to determine consensus
/// * `raw_consensus_result` - the original consensus
/// * `d6_typer` - utility for assigning labels to alleles
/// * `max_missing_consensus_frac` - the maximum allowed missing from a consensus for identification
fn merge_consensus_results(
    sequences: &[&str],
    offsets: &[usize],
    cdwfa_config: &CdwfaConfig, 
    raw_consensus_result: &PriorityConsensus,
    d6_typer: &Cyp2d6Extractor,
    max_missing_consensus_frac: f64
) -> Result<MultiConsensus, Box<dyn std::error::Error>> {
    let mut consensus_set: BTreeMap<(String, String), Vec<usize>> = Default::default();
    let mut unknown_set: BTreeMap<String, Vec<usize>> = Default::default();
    for (i, consensus) in raw_consensus_result.consensuses().iter().enumerate() {
        // pull out the two consensus sequences
        let hpc_consensus: String = std::str::from_utf8(consensus[0].sequence())?.to_string();
        let full_consensus: &str = std::str::from_utf8(consensus[1].sequence())?;

        // figure out the label for this consensus
        let sequence_to_type = full_consensus.trim_matches('*');
        debug!("Typing consensus #{i}, {} wildcards trimmed", full_consensus.len() - sequence_to_type.len());
        
        // we are *expecting* full length sequences in our consensus
        let force_assignment = false; // if we have label ambiguity, then the allele is incomplete and we want it to get merged if possible
        let allele = match d6_typer.find_full_type_in_sequence(sequence_to_type, max_missing_consensus_frac, force_assignment) {
            Ok(al) => al,
            Err(e) => {
                let unknown = Cyp2d6RegionLabel::new_unknown();
                error!("Error while typing consensus #{i}, setting to {unknown}.");
                error!("Typing error: {e}");
                Cyp2d6Region::new(unknown, None)
            }
        };
        let allele_label = allele.label();

        // reduce the label for merging
        let detailed = true; // true means we keep sub-alleles such as "*4.001"
        let reduced_label = allele_label.simplify_allele(detailed, d6_typer.cyp2d6_config().cyp_translate());
        debug!("Reduced {allele_label} to {reduced_label} for merging.");

        if !allele_label.is_allowed_label() {
            // this one is not allowed, lets see if we can merge it into the parent later
            let entry = unknown_set.entry(hpc_consensus).or_default();
            entry.push(i);
        }
        else {
            // now save it to key (hpc sequence, allele label)
            let entry = consensus_set.entry((hpc_consensus, reduced_label)).or_default();
            entry.push(i);
        }   
    }

    // see if we can collapse the unknowns into a parent
    let mut intentional_ignore: HashSet<(String, String)> = Default::default();
    for (hpc_consensus, entries) in unknown_set.into_iter() {
        // figure out how many matches we have
        let mut other_keys = vec![];
        for (merge_key, _merge_set) in consensus_set.iter() {
            if hpc_consensus == merge_key.0 {
                other_keys.push(merge_key.clone());
            }
        }

        let hpc_count = other_keys.len();

        // TODO: handling of multiple or no matches may need to be explored further in the future
        match hpc_count {
            0 => {
                // we have no HPC neighbors, so just collapse all of them into one big UNKNOWN pile
                let unknown = Cyp2d6RegionLabel::new_unknown();
                assert!(consensus_set.insert((hpc_consensus, unknown.full_allele()), entries).is_none());
            },
            1 => {
                // we found exactly one parent option we can merge with, so lets do it!
                let other_key = other_keys.pop().unwrap();
                let entry = consensus_set.get_mut(&other_key).unwrap();
                debug!("Collapsing entries {entries:?} into HPC relative {} ({entry:?})", other_key.1);
                entry.extend(entries);
            }
            _n => {
                // there are multiple, and we don't have logic to resolve that currently...
                // for now, we will completely ignore the read set
                // TODO: there may come a point where we want to somehow compare these groupings to the HPC relatives and force them into one
                debug!("Multiple collapse options detected for entries {entries:?}, ignoring.");
                let unknown = Cyp2d6RegionLabel::new_unknown();
                intentional_ignore.insert((hpc_consensus.clone(), unknown.full_allele()));
                assert!(consensus_set.insert((hpc_consensus, unknown.full_allele()), entries).is_none());
            }
        };
    }

    let mut consensuses = vec![];
    let mut sequence_indices = vec![usize::MAX; raw_consensus_result.sequence_indices().len()];
    for (masked_consensus, con_indices) in consensus_set.iter() {
        let con_index = consensuses.len();
        debug!("Collapsing {con_indices:?} into {con_index}");
        let consensus: Consensus = if intentional_ignore.contains(masked_consensus) {
            // fill in the si
            let mut num_scored = 0;
            for (i, si) in raw_consensus_result.sequence_indices().iter().enumerate() {
                if con_indices.contains(si) {
                    sequence_indices[i] = con_index;
                    num_scored += 1;
                }
            }

            // this was a group we made but that we should really ignore
            Consensus::new(vec![], waffle_con::cdwfa_config::ConsensusCost::L1Distance, vec![0; num_scored])
        } else if con_indices.len() == 1 {
            // we can just copy the result, but we do still need to propagate the sequence index information
            for (i, &si) in raw_consensus_result.sequence_indices().iter().enumerate() {
                if si == con_indices[0] {
                    sequence_indices[i] = con_index;
                }
            }

            // return the copy of the full consensus (index 1)
            raw_consensus_result.consensuses()[con_indices[0]][1].clone()
        } else {
            // we have a merge situation
            let mut combined_consensus = ConsensusDWFA::with_config(cdwfa_config.clone())?;
            for (si, (seq, offset)) in (sequences.iter().zip(offsets.iter())).enumerate() {
                let assigned_consensus = raw_consensus_result.sequence_indices()[si];
                if con_indices.contains(&assigned_consensus) {
                    let opt_offset = if *offset == 0 { None } else { Some(*offset) };

                    // add the sequence
                    combined_consensus.add_sequence_offset(seq.as_bytes(), opt_offset)?;
                    // at the same time, mark this sequence index
                    sequence_indices[si] = con_index;
                }
            }
            
            // now run the consensus
            let mut new_consensus_set = combined_consensus.consensus()?;
            assert!(!new_consensus_set.is_empty());
            if new_consensus_set.len() > 1 {
                warn!("Multiple consensuses found during collapse, picking first.");
            }
            new_consensus_set.remove(0)
        };
        consensuses.push(consensus);
    }

    // make sure all sequences have been re-assigned to a valid index
    assert!(sequence_indices.iter().all(|&v| v < consensuses.len()));

    Ok(MultiConsensus::new(
        consensuses,
        sequence_indices
    ))
}

/// Given a particular chain, this will construct the user friendly visual of that chain.
/// E.g. [0, 0, 1] will become "*4x2 + *10"
/// # Arguments
/// * `chain` - the chains of alleles to tie together
/// * `hap_regions` - the haplotype labels for the internals, these can get converted to alleles where appropriate
/// * `detail_level` - sets the corresponding level of detail on the reported alleles
/// * `cyp_translate` - a map from internal CYP name to user friendly name
pub fn convert_chain_to_hap(chain: &[usize], hap_regions: &[Cyp2d6Region], detail_level: Cyp2d6DetailLevel, cyp_translate: &BTreeMap<String, String>) -> String {
    // track the number of non-deletion alleles we identify
    // this is robust to a potential *5x2 situation (which may also have reporting issue, but that's a future problem)
    let mut num_non_deletion = 0;
    
    // first identify the reportable indices in the chain
    let reportable_indices: Vec<usize> = chain.iter()
        .rev()
        .filter(|&&c_index| {
            // only keep CYP2D alleles that are not CYP2D7
            let hap_label = hap_regions[c_index].label();
            let keep_allele = hap_label.is_cyp2d() && hap_label.region_type() != Cyp2d6RegionType::Cyp2d7;
            if keep_allele && hap_label.region_type() != Cyp2d6RegionType::Cyp2d6Deletion {
                num_non_deletion += 1;
            }
            keep_allele
        })
        .copied()
        .collect();
    
    // secondary filtering and translation into a string
    reportable_indices.iter()
        .filter_map(|&c_index| {
            let hap_label = hap_regions[c_index].label();
            if hap_label.region_type() == Cyp2d6RegionType::Cyp2d6Deletion && num_non_deletion > 0 {
                // this will filter out *5 (deletion) alleles if something else is on the same chain; i.e. *5+*10 is not an allowed report (although it may happen biologically)
                None
            } else {
                // passed all above secondary filtering
                // convert the allele index into a human readable string
                let string_label = match detail_level {
                    Cyp2d6DetailLevel::CoreAlleles => hap_label.simplify_allele(false, cyp_translate),
                    Cyp2d6DetailLevel::SubAlleles => hap_label.simplify_allele(true, cyp_translate),
                    Cyp2d6DetailLevel::DeepAlleles => format!("({})", hap_regions[c_index].deep_label()),
                };
                Some(string_label)
            } 
        })
        .group_by(|v| v.clone()) // group them by the ID so we get adjacent counts
        .into_iter()
        .map(|(string_label, group)| {
            let group_len = group.count();
            // add in any xN values
            if group_len > 1 {
                format!("{string_label}x{group_len}")
            } else {
                string_label
            }
        })
        .join(" + ")
}

#[cfg(test)]
mod tests {
    use crate::cyp2d6::definitions::Cyp2d6Config;

    use super::*;

    /*
    TODO: tests that ideally exist but are difficult to implement:
    - diplotype_cyp2d6 - this is basically end-to-end test; for now, we will rely on the pipeline since this is difficult to encode
    - merge_consensus_results - also difficult to do, we need to have a D6 typer built for this, which means loading one from somewhere; pipeline will suffice for now
     */

    #[test]
    fn test_convert_chain_to_hap() {
        let hap_labels = vec![
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d7, None),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("1.001".to_string())),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("10".to_string())),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("1.002".to_string())),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("1.002".to_string()))
        ];
        let hap_regions: Vec<Cyp2d6Region> = hap_labels.into_iter().map(|l| Cyp2d6Region::new(l, None)).collect();

        // get the translator for testing
        let cyp2d6_config = Cyp2d6Config::default();
        let cyp_translate = cyp2d6_config.cyp_translate();
        
        // basic example
        let chain = vec![2, 2, 1, 0];
        let hap = convert_chain_to_hap(&chain, &hap_regions, Cyp2d6DetailLevel::SubAlleles, cyp_translate);
        assert_eq!(&hap, "*1.001 + *10x2"); // remember these get reversed

        // eventually this will get collapsed instead of separate
        let chain = vec![3, 1, 0];
        let hap = convert_chain_to_hap(&chain, &hap_regions, Cyp2d6DetailLevel::SubAlleles, cyp_translate);
        assert_eq!(&hap, "*1.001 + *1.002");

        // test the collapsed version
        let hap = convert_chain_to_hap(&chain, &hap_regions, Cyp2d6DetailLevel::CoreAlleles, cyp_translate);
        assert_eq!(&hap, "*1x2");

        // tests if there are two near-identical alleles with the same name get collapsed
        let chain = vec![3, 4];
        let hap = convert_chain_to_hap(&chain, &hap_regions, Cyp2d6DetailLevel::SubAlleles, cyp_translate);
        assert_eq!(&hap, "*1.002x2");

        // TODO: do we need to test Cyp2d6DetailLevel::DeepAlleles? it's in the debug, not core outputs
    }
}
