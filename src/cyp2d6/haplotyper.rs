
use hiphase::data_types::variants::Variant;
use hiphase::wfa_graph::{NodeAlleleMap, WFAGraph, WFAResult};
use itertools::Itertools;
use log::{debug, trace};
use minimap2::Aligner;
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use simple_error::bail;
use std::collections::BTreeMap;

use crate::cyp2d6::definitions::{Cyp2d6Config, generate_cyp_hybrids};
use crate::cyp2d6::region_label::{Cyp2d6RegionLabel, Cyp2d6RegionType};
use crate::data_types::database::PgxDatabase;
use crate::data_types::mapping::MappingStats;

/// The primary interface for identifying regions of interest within a sequence.
/// Can be used to find all base type regions (e.g. CYP2D6, link_region, etc.) within a longer sequence, or to find specific full-length star-alleles (e.g., CYP2D6*4.001).
pub struct Cyp2d6Extractor<'a> {
    /// Set of loaded variants required for WFA typing
    loaded_variants: LoadedVariants,
    /// Contains a map from a haplotype the set of included variants
    haplotype_lookup: BTreeMap<Cyp2d6RegionLabel, Vec<u8>>,
    /// These are the baseline sequences to go from a label to a sequence
    hybrid_sequences: HashMap<Cyp2d6RegionLabel, String>,
    /// The set of hybrids with deep typing (e.g. CYP2D6 and some hybrids)
    mapped_hybrids: HashSet<Cyp2d6RegionLabel>,
    /// Contains a reference to the reference genome we are using
    reference_genome: &'a ReferenceGenome,
    /// Contains a reference to the CYP2D6 config we are using
    cyp2d6_config: &'a Cyp2d6Config
}

impl<'a> Cyp2d6Extractor<'a> {
    /// Creates a new CYP2D6 region extractor based on a load PGx database and the provided reference genome sequences.
    /// # Arguments
    /// * `database` - a pre-loaded PGx database; should have information on the CYP2D6 allele definitions
    /// * `reference_genome` - a pre-loaded reference genome, which is used to stitch the exon/intron sequences together
    /// # Errors
    /// * if the variants do not load correctly from the database
    /// * if a variant in the database cannot be found (usually a database error/corruption)
    /// * if hybrid generation fails
    pub fn new(database: &'a PgxDatabase, reference_genome: &'a ReferenceGenome) -> Result<Cyp2d6Extractor<'a>, Box<dyn std::error::Error>> {
        // first, just load the variants in
        let loaded_variants = load_variant_database(database)?;
        
        // this is a map from an allele ID to the vector of 0s and 1s for the variant set
        let mut haplotype_lookup: BTreeMap<Cyp2d6RegionLabel, Vec<u8>> = Default::default();
        let variant_list = loaded_variants.ordered_variants();
        let num_variants = variant_list.len();
        for (_allele_id, allele_def) in database.cyp2d6_gene_def().iter() {
            assert_eq!(allele_def.gene_name(), "CYP2D6");
            let mut allele_assignments: Vec<u8> = vec![0; num_variants];
            for variant_def in allele_def.variants().iter() {
                // build the variant key
                let var_pos = variant_def.position();
                let var_ref = variant_def.reference();
                let var_alt = variant_def.alternate();

                // now set that bit to a 1
                let var_index = loaded_variants.index_variant(var_pos, var_ref, var_alt)?;
                allele_assignments[var_index] = 1;
            }
            haplotype_lookup.insert(
                Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some(allele_def.star_allele().to_string())),
                allele_assignments
            );
        }

        // cache this reference
        let cyp2d6_config = database.cyp2d6_config();

        // here is the core region of interest
        let cyp_coordinates = cyp2d6_config.cyp_coordinates();
        let full_d6_region = cyp_coordinates.get("CYP2D6").unwrap().clone();
        
        /*
        // 2 - construct a full D6 region that includes a little buffer around the first and last variants in our list
        //     we will use this to anchor our graph later
        //     this step was originally when we were doing some sanity checks, it's just stored in CYP_COORDINATES now
        let first_variant_pos = loaded_variants.first_variant_pos() as u64;
        let last_variant_pos = loaded_variants.last_variant_pos() as u64;
        let buffer = 50;
        
        // the upstream region (remember, rev-comp) will be from the end of exon 1 to just past the last variant position
        let first_exon = CYP_REGIONS.get("CYP2D6").unwrap().get("exon1").unwrap();
        let upstream_d6 = Coordinates::new("chr22".to_string(), first_exon.end(), last_variant_pos+buffer);

        // the downstream region (rev-comp) will be from just before the first variant position to the start of the final exon
        let last_exon = CYP_REGIONS.get("CYP2D6").unwrap().get("exon9").unwrap();
        let downstream_d6 = Coordinates::new("chr22".to_string(), first_variant_pos-buffer, last_exon.start());

        // this is the full target region for CYP2D6 variant calling; we will search for this in our samples
        let derived_d6_region = Coordinates::new("chr22".to_string(), downstream_d6.start(), upstream_d6.end());
        // let full_d7_region = CYP_COORDINATES.get("CYP2D7").unwrap().clone();
        let full_d6_sequence = String::from_utf8(reference_genome.get_slice(full_d6_region.chrom(), full_d6_region.start() as usize, full_d6_region.end() as usize).to_vec())?;
        println!(">full_d6_sequence");
        println!("{full_d6_sequence}");
        
        // this is a sanity check until we auto-derive the D6/D7 coordinates
        // NOTE: if this fails, then we need to check `extraction_region()` when we update.
        if derived_d6_region != full_d6_region {
            warn!("Hard-coded CYP2D6 coordinates do not match derived search space.");
        }
        */

        // make sure our designated coordinates is larger than the variant list set
        assert!(full_d6_region.start() <= variant_list[0].position() as u64);
        assert!(full_d6_region.end() >= variant_list.last().unwrap().position() as u64);

        // 3 - contruct search sequences that match the D6, D7, hybrid, and deletion regions
        let hybrid_sequences: HashMap<Cyp2d6RegionLabel, String> = generate_cyp_hybrids(reference_genome, cyp2d6_config)?;

        // also enumerate the ones that go through deep genotyping
        let mapped_hybrids: HashSet<Cyp2d6RegionLabel> = vec![
                Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, None),
                // hybrids currently in the database
                Cyp2d6RegionLabel::new(Cyp2d6RegionType::Hybrid, Some("CYP2D6::CYP2D7::exon9".to_string())) // *36 typically
            ].into_iter()
            .collect();
        
        Ok(Cyp2d6Extractor {
            loaded_variants,
            haplotype_lookup,
            hybrid_sequences,
            mapped_hybrids,
            reference_genome,
            cyp2d6_config
        })
    }

    /// This will search a given sequence for each of the D6 targets and select the best one when overlaps are detected.
    /// This is fairly common since D6 and D7 are similar, and obviously the fusions are as well.
    /// # Arguments
    /// * `search_sequence` - sequence we are searching through for alleles
    /// * `penalized_unmapped` - pass-through option on whether unmapped bases add to numerator (penalized) or subtract from denominator (non-penalized)
    /// * `max_missing_frac` - the maximum missing fraction that is allowed; primarily enforced when `penalized_unmapped` is true
    /// # Errors
    /// * if the aligner fails to build
    pub fn find_base_type_in_sequence(
        &self,
        search_sequence: &str, 
        penalize_unmapped: bool,
        max_missing_frac: f64
    ) -> Result<Vec<AlleleMapping>, Box<dyn std::error::Error>> {
        // first, handle the stupid case
        if search_sequence.is_empty() {
            return Ok(vec![]);
        }

        // these get converted into our AlleleMapping when finished
        let mut region_mappings: Vec<(Cyp2d6RegionLabel, std::ops::Range<usize>, MappingStats)> = vec![];
        let dna_aligner: Aligner = Aligner::builder()
            .map_hifi()
            .with_cigar()
            .with_seq(search_sequence.as_bytes())?;

        // this is the maximum edit penalty we allow
        // if `penalized_unmapped` is True, then this is (NM + unmapped) / total; otherwise, NM / (total - unmapped)
        let max_ed_frac = 0.05;
        
        // if we are penalizing unmapped bases, then this threshold will get tighter
        let max_penalized_frac = max_missing_frac;
        let penalize_during_search = false; // we don't want to penalize during the initial search, but we will at the end if enabled
        
        // we only need cigar and md for debugging
        // other settings for mapping
        let output_cigar: bool = true;
        let output_md: bool = true;
        let max_frag_len: Option<usize> = None;
        let extra_flags = None;

        // lets go through the keys in sorted order for output purposes
        let mut uncollapsed_regions: Vec<(std::ops::Range<usize>, MappingStats, Cyp2d6RegionLabel)> = vec![];
        let key_order: Vec<&Cyp2d6RegionLabel> = self.hybrid_sequences.keys()
            // .sorted()
            .sorted_by(|a, b| {
                // TODO: this is preserving the previous sort order (which apparently can matter given our current setup); ideally, we just used .sorted()
                //       current causes an issue in one sample when we change that; we should figure out a strategy that resolves
                //       the underlying issue at some point; for refactoring, lets just preserve the order by comparing the full_allele labels
                a.full_allele().cmp(&b.full_allele())
            })
            .collect();

        let penalized_types = [
            Cyp2d6RegionType::Cyp2d6Deletion, // this one needs the full region because it has to anchor two ends together
            Cyp2d6RegionType::Rep6, Cyp2d6RegionType::Rep7, // these regions are too similar and get clipped in ways that lead to incorrect answers
            // TODO: if we add an explicit clustering step, we may be able to relax this
            //       Xiao says that REP6 and REP7 are near identical except a few bases at the end; I bet that is the source of the issue
            //       if we cut those out, we could probably just have a "REP" region and push the deltas to the flanks
        ];
        
        for &target_id in key_order.iter() {
            // pull the actual target sequence
            let target_seq = self.hybrid_sequences.get(target_id).unwrap();

            // now find all of the mappings of this sequence
            let mappings = dna_aligner.map(
                target_seq.as_bytes(),
                output_cigar, output_md, max_frag_len, extra_flags.clone()
            )?;
            
            if mappings.is_empty() {
                trace!("\t{target_id}: None");
            }

            // it's possible to get multiple mappings
            for m in mappings.iter() {
                // all results are relative to the `target_seq`; aka, the D6/D7 allele
                let nm = m.alignment.as_ref().unwrap().nm as usize;
                let seq_len = target_seq.len();
                let unmapped = seq_len - (m.query_end - m.query_start) as usize;

                // the amount clipped at the start is the amount into query that we start
                let clipped_start = m.query_start as usize;
                // the amount clipped at the end is the length minus the query end point
                let clipped_end = target_seq.len() - m.query_end as usize;

                let mapping_stats = MappingStats::new_with_clippings(
                    seq_len, nm, unmapped,
                    clipped_start, clipped_end
                );

                // for *5 specifically, we need really tight controls because it is relatively precise signature so we don't want to miss much
                let penalize_star5 = penalize_during_search ||
                    penalized_types.contains(&target_id.region_type());
                let custom_score = mapping_stats.custom_score(penalize_star5);
                
                if custom_score.score() > max_ed_frac { // && target_id == "CYP2D6*5" {
                    // ignore this one
                    debug!("\tIgnoring {target_id}: {}-{} => {}", m.target_start, m.target_end, mapping_stats.custom_score_string(penalize_star5));
                } else {
                    debug!("\t{target_id}: {}-{} => {}", m.target_start, m.target_end, mapping_stats.custom_score_string(penalize_during_search));
                    uncollapsed_regions.push((
                        m.target_start as usize..m.target_end as usize,
                        mapping_stats,
                        target_id.clone(),
                    ));
                }
            }
        }

        // order the regions by the range they are a part of
        uncollapsed_regions.sort_by(|v1, v2| {
            (v1.0.start, v1.0.end).cmp(&(v2.0.start, v2.0.end))
        });

        // now collapse the overlapping regions
        // NOTE: this has the potential to fail for some weird edge cases, mainly if alignment sizes are vastly different, but this should be rare if it happens at all
        // NOTE: for the purpose of initial pullout, we really just want the larger region (observation suggest it was always the "best" as well); we can try that if this ever becomes an issue
        // the last bool here is if anything overlapping had the minimum mappings requirement
        let mut current_region: Option<(std::ops::Range<usize>, MappingStats, Cyp2d6RegionLabel)> = None;
        for (ucr_range, ucr_score, ucr_id) in uncollapsed_regions.into_iter() {
            match current_region.as_ref() {
                None => {
                    current_region = Some((ucr_range, ucr_score, ucr_id));
                }
                Some(cr) => {
                    if overlap_score(&ucr_range, &current_region.as_ref().unwrap().0) > 0.9 {
                        // these overlap enough that we compare them
                        // if it's a *5 involved, we need to compare the full lengths
                        let star5_pairing = penalized_types.contains(&ucr_id.region_type()) || penalized_types.contains(&cr.2.region_type());
                        let penalized_scoring = if star5_pairing { true } else { penalize_during_search };
                        let ucr_priority = get_allele_priority(&ucr_id);
                        let cr_priority = get_allele_priority(&cr.2);

                        // if our score is lower AND we have the same or better priority 
                        // OR if our priority is better
                        if (ucr_score.custom_score(penalized_scoring) < cr.1.custom_score(penalized_scoring) &&
                            ucr_priority >= cr_priority) || ucr_priority > cr_priority {
                            // the score of this new region is better than what we were looking at, so replace it
                            current_region = Some((ucr_range, ucr_score, ucr_id));
                        }
                    } else {
                        // not an overlap, save the current region and push this one
                        let unwrapped = current_region.unwrap();
                        region_mappings.push((unwrapped.2, unwrapped.0, unwrapped.1));
                        current_region = Some((ucr_range, ucr_score, ucr_id));
                    }
                }
            };
        }

        if let Some(unwrapped) = current_region {
            // handle the last one
            region_mappings.push((unwrapped.2, unwrapped.0, unwrapped.1));
        }

        // final results go here
        let mut ret: Vec<AlleleMapping> = vec![];
        
        debug!("Collapsed calls:");
        for (region_label, mapping_region, mapping_stats) in region_mappings.into_iter() {
            // 3b - the D7 alleles and most hybrids are done at this point
            let penalized_score = mapping_stats.custom_score(true);
            if penalized_score.score() > max_penalized_frac {
                debug!("\tIgnoring {region_label} at {mapping_region:?}, too short: {}", mapping_stats.custom_score_string(true));
            } else {
                debug!("\t{region_label} at {mapping_region:?}: {}", mapping_stats.custom_score_string(penalize_unmapped));
                ret.push(AlleleMapping::new(
                    region_label, mapping_region, mapping_stats
                ));
            }
        }

        Ok(ret)
    }

    /// This will search a given sequence for D6 targets and select the best one when overlaps are detected.
    /// Additionally, it will go deeper than just "CYP2D6" and will attempt a full typing of the allele.
    /// # Arguments
    /// * `search_sequence` - the sequence to search; unmapped bases are penalized by default
    /// * `max_missing_frac` - the maximum missing fraction that is allowed; primarily enforced when `penalized_unmapped` is true
    /// * `force_assignment` - if True, then ambiguous CYP2D6 results will be arbitrarily assigned one of the equal values
    /// # Errors
    /// * if finding base types fails
    /// * if assigning haplotype labels to a base type fails
    pub fn find_full_type_in_sequence(
        &self, search_sequence: &str, max_missing_frac: f64, force_assignment: bool
    ) -> Result<Cyp2d6RegionLabel, Box<dyn std::error::Error>> {
        // for each sequence, figure out what it best matches in our targets
        // this should be full length matches UNLESS the consensus is incomplete
        // I think this means we want to use a penalized version
        let penalize_unmapped = true;
        let best_matches = self.find_base_type_in_sequence(
            search_sequence,
            penalize_unmapped,
            max_missing_frac
        )?;

        if best_matches.is_empty() {
            bail!("no matches found");
        }
    
        // sometimes we get multiples, usually when wildcards are involved; get the one with the lowest score (score = mismatches)
        let best_match = best_matches.iter()
            .min_by(|a, b| 
                a.mapping_stats.custom_score(penalize_unmapped).partial_cmp(&b.mapping_stats.custom_score(penalize_unmapped)).unwrap()
            )
            .unwrap();
    
        let final_type = if self.mapped_hybrids.contains(best_match.allele_label()) {
            debug!("\tConverting {} to full allele definition...", best_match.allele_label());
            self.assign_haplotype(
                search_sequence.as_bytes(),
                force_assignment
            )?
        } else {
            best_match.allele_label().clone()
        };
        Ok(final_type)
    }

    /// This will take a chosen region from a sequence and deep genotype it in CYP2D6 using our WFA-based variant system.
    /// # Arguments
    /// * `sequence` - presumably matches a D6 allele that we want to find the best match for.
    /// * `force_assignment` - if True, then ambiguous results will be arbitrarily assigned one of the equal values
    /// # Errors
    /// * if graph backbone construction fails
    /// * if aligner construction fails
    /// * if WFA graph construction or traversal fails
    fn assign_haplotype(
        &self,
        sequence: &[u8],
        force_assignment: bool
    ) -> Result<Cyp2d6RegionLabel, Box<dyn std::error::Error>> {
        let cyp_coordinates = self.cyp2d6_config.cyp_coordinates();
        
        // get the relevant sequences from the reference
        let backbone_coordinates = cyp_coordinates.get("CYP2D6_wfa_backbone").unwrap();
        let chrom_seq = self.reference_genome.get_full_chromosome(backbone_coordinates.chrom());
        let graph_backbone = String::from_utf8(self.reference_genome.get_slice(backbone_coordinates.chrom(), backbone_coordinates.start() as usize, backbone_coordinates.end() as usize).to_vec())?;

        // we only need cigar and md for debugging
        // other settings for mapping
        let output_cigar: bool = true;
        let output_md: bool = true;
        let max_frag_len: Option<usize> = None;
        let extra_flags = None;

        // now build out a mapper so we make sure we have the right coordinates
        let dna_aligner: Aligner = Aligner::builder()
            .map_hifi()
            .with_cigar()
            .with_seq(graph_backbone.as_bytes())?;

        // now find all of the mappings of our consensus onto this backbone
        let mappings = dna_aligner.map(
            sequence,
            output_cigar, output_md, max_frag_len, extra_flags.clone()
        )?;

        // we *should* only find one, but that does not always happen
        assert!(!mappings.is_empty());
        let core_mapping = if mappings.len() == 1 {
            &mappings[0]
        } else {
            // we somehow have multiple mappings, pick the longest one
            let longest_index = mappings.iter()
                .enumerate()
                .map(|(i, m)| (m.block_len, i))
                .max().unwrap()
                .1;
            &mappings[longest_index]
        };

        // figure out which part of the backbone was actually use by pulling out the target coordinates and adding to the backbone start
        let aligned_start = backbone_coordinates.start() as usize + core_mapping.target_start as usize;
        let aligned_end = backbone_coordinates.start() as usize + core_mapping.target_end as usize;
        
        // figure out which part of the sequence is relevant
        let sub_sequence_start = core_mapping.query_start as usize;
        let sub_sequence_end = core_mapping.query_end as usize;
        let sub_sequence = &sequence[sub_sequence_start..sub_sequence_end];

        /*
        println!(">sub_sequence");
        println!("{}", std::str::from_utf8(sub_sequence).unwrap());
        println!(">graph_backbone");
        println!("{}", graph_backbone);
        */

        // 4 - extract just that end-to-end region for GraphWFA
        let start_time = std::time::Instant::now();
        let (wfa_graph, node_to_alleles): (WFAGraph, NodeAlleleMap) = 
            WFAGraph::from_reference_variants(
                chrom_seq, 
                self.loaded_variants.ordered_variants(), // these are both range style indices
                // use `target_region` to the do the full thing, even if the mapping is only partial
                // target_region.start() as usize, 
                // target_region.end() as usize
                // otherwise, let's fix to just the relevant aligned backbone region
                aligned_start,
                aligned_end
        )?;

        // we can probably make this smaller eventually, but this is a safe distance for now
        let wfa_prune_distance = 1000;
        let wfa_result: WFAResult = wfa_graph.edit_distance_with_pruning(sub_sequence, wfa_prune_distance)?;
        debug!(
            "\t{} WFAGraph result ({}) => num_nodes: {}, read_len: {}, edit_distance: {}", 
            "", start_time.elapsed().as_secs_f32(), wfa_graph.get_num_nodes(), sub_sequence.len(), wfa_result.score()
        );

        // 5 - assign the alleles and label it
        // we will populate these with the variant level info
        let num_variants = self.loaded_variants.ordered_variants().len();
        let mut alleles: Vec<u8> = vec![3; num_variants];
        let first_overlap = 0; // we are force matching everything

        // this loop will set the alleles for what was aligned, allowing us to compare to the pre-defined alleles for the haplotypes
        for traversed_index in wfa_result.traversed_nodes().iter() {
            for &(var_index, allele_assignment) in node_to_alleles.get(traversed_index).unwrap_or(&vec![]).iter() {
                let correct_index: usize = first_overlap+var_index;
                if alleles[correct_index] == 3 {
                    alleles[correct_index] = allele_assignment;
                } else if alleles[correct_index] != allele_assignment {
                    alleles[correct_index] = 2;
                }
            }
        }

        // 3 - use that information to score based on matching VI alleles, and then matching total alleles
        let mut best_id_set: HashSet<Cyp2d6RegionLabel> = Default::default();
        best_id_set.insert(Cyp2d6RegionLabel::new(Cyp2d6RegionType::Unknown, None));
        let mut best_score: (usize, usize) = (0, 0);

        for (allele_id, haplotype_vec) in self.haplotype_lookup.iter() {
            // these start as fully equal, and we remove things that do NOT match
            let mut vi_match = 0;
            let mut all_match = 0;

            // go through the pairs of alleles and add the ones that match together
            // let haplotype_vec = haplotype_lookup.get(allele_id).unwrap();
            assert_eq!(alleles.len(), haplotype_vec.len());

            for (i, (&seq_value, &hap_value)) in alleles.iter().zip(haplotype_vec.iter()).enumerate() {
                // hap_value will *always* be 0 or 1
                assert!(hap_value == 0 || hap_value == 1);

                // seq_value can be 0, 1, 2 (ambiguous), or 3 (unset, usually meaning 0)
                let is_match = match seq_value {
                    // both are clearly set, so compare for equality
                    0 | 1 => hap_value == seq_value,
                    // it could be either 0 or 1, so count it
                    2 => true,
                    // it was not set, this could be either unassigned OR a multi-allelic site; in either case, not counting is best
                    3 => false,
                    v => panic!("Unexpected seq_value={v}")
                };

                if is_match {
                    // these were a match, count the appropriate fields
                    all_match += 1;
                    if self.loaded_variants.is_vi(i) {
                        vi_match += 1;
                    }
                }
            }

            let combined_score = (vi_match, all_match);
            match combined_score.cmp(&best_score) {
                std::cmp::Ordering::Greater => {
                    // new best, clear out the old and add the new
                    trace!("new best: {allele_id} = {combined_score:?}");
                    best_id_set.clear();
                    best_id_set.insert(allele_id.clone());
                    best_score = combined_score;
                },
                std::cmp::Ordering::Equal => {
                    // equal result, add to the existing set
                    trace!("new equi: {allele_id} = {combined_score:?}");
                    best_id_set.insert(allele_id.clone());
                },
                std::cmp::Ordering::Less => {}, // score is less, do nothing
            };
        }

        let best_id = if best_id_set.len() == 1 {
            // only one result, drain it off
            best_id_set.drain().next().unwrap()
        } else {
            // sort the candidates by allele label (which is ~numerically for tie-breaking)
            let ordered_candidates: Vec<Cyp2d6RegionLabel> = best_id_set.into_iter()
                .sorted_by(|a, b| a.full_allele().cmp(&b.full_allele()))
                .collect();
            let ordered_labels: Vec<String> = ordered_candidates.iter().map(|l| l.to_string()).collect();

            if force_assignment {
                debug!("\tAmbiguous result detected, selecting first candidate; candidates: {ordered_labels:?}");
                ordered_candidates[0].clone()
            } else {
                debug!("\tAmbiguous result detected, setting to unknown; candidates: {ordered_labels:?}");
                Cyp2d6RegionLabel::new(Cyp2d6RegionType::Unknown, None)
            }
        };

        debug!("\t{} {best_id} -> {best_score:?}, ({:.4}, {:.4})", "", best_score.0 as f64 / self.loaded_variants.num_vi() as f64, best_score.1 as f64 / num_variants as f64);

        Ok(best_id)
    }

    /// Returns the sequence for a given allele
    /// # Arguments
    /// * `allele_name` - the allele to look up
    pub fn get_allele(&self, allele_name: &Cyp2d6RegionLabel) -> Option<&String> {
        self.hybrid_sequences.get(allele_name)
    }

    // getters
    pub fn loaded_variants(&self) -> &LoadedVariants {
        &self.loaded_variants
    }

    pub fn cyp2d6_config(&self) -> &Cyp2d6Config {
        self.cyp2d6_config
    }
}


/// Wrapper for variants that have been loaded in preparation for genotyping.
pub struct LoadedVariants {
    /// variants ordered by position and sequence
    ordered_variants: Vec<Variant>,
    /// a lookup from (position, REF, ALT) to index in `ordered_variants`
    variant_lookup: HashMap<(usize, String, String), usize>,
    /// one-to-one with `ordered_variants`; if True, this variant had a VI flag
    is_vi: Vec<bool>
}

impl LoadedVariants {
    /// Constructor for the loaded variants, see parameters.
    /// # Arguments
    /// * `ordered_variants` - variants ordered by position and sequence
    /// * `variant_lookup` - a lookup from (position, REF, ALT) to index in `ordered_variants`
    /// * `is_vi` - one-to-one with `ordered_variants`; if True, this variant had a VI flag indicating that it distinguishes the core allele
    /// # Errors
    /// * if `ordered_variants`, `is_vi`, and `variant_lookup` do not all have the same number of entriesa
    pub fn new(ordered_variants: Vec<Variant>, variant_lookup: HashMap<(usize, String, String), usize>, is_vi: Vec<bool>) -> Result<LoadedVariants, Box<dyn std::error::Error>> {
        if ordered_variants.len() != is_vi.len() {
            bail!("ordered_variants and is_vi must be same length");
        }
        if ordered_variants.len() != variant_lookup.len() {
            bail!("ordered_variants and variant_lookup must be same length");
        }
        Ok(LoadedVariants {
            ordered_variants, variant_lookup, is_vi
        })
    }

    pub fn ordered_variants(&self) -> &[Variant] {
        &self.ordered_variants
    }

    /// Searches for a given variant in our loaded variants and return the index
    /// # Arguments
    /// * `position` - the coordinate to find, 0-based
    /// * `reference` - the reference allele
    /// * `alternate` - the alternate allele
    /// # Errors
    /// * if the allele is not found
    pub fn index_variant(&self, position: usize, reference: &str, alternate: &str) -> Result<usize, Box<dyn std::error::Error>> {
        match self.variant_lookup.get(&(position, reference.to_string(), alternate.to_string())) {
            Some(&idx) => Ok(idx),
            None => bail!("({position}, {reference}, {alternate}) not found")
        }
    }

    /// Retrieves the position of the first variant in the dataset
    pub fn first_variant_pos(&self) -> i64 {
        self.ordered_variants[0].position()
    }

    /// Retrieves the position of the lastt variant in the dataset
    pub fn last_variant_pos(&self) -> i64 {
        self.ordered_variants.last().unwrap().position()
    }

    /// Returns true if the variant at the given index is VI
    pub fn is_vi(&self, index: usize) -> bool {
        self.is_vi[index]
    }

    /// Returns the total number of variants marked as VI
    pub fn num_vi(&self) -> usize {
        self.is_vi.iter().filter(|&x| *x).count()
    }
}

/// This will parse the relevant CYP2D6 variants in preparation for using in GraphWFA
/// # Arguments
/// * `database` - the PGx database pre-loaded, we need to translate this for GraphWFA
/// # Errors
/// * if allele translation to UTF8 fails
fn load_variant_database(database: &PgxDatabase) -> Result<LoadedVariants, Box<dyn std::error::Error>> {
    let mut inserted_variants: HashSet<(usize, String, String)> = Default::default();
    let mut variant_list: Vec<Variant> = vec![];
    let mut vi_set: HashSet<(usize, String, String)> = Default::default();
    let mut all_set: HashSet<(usize, String, String)> = Default::default();
    for (_allele_id, allele_def) in database.cyp2d6_gene_def().iter() {
        for variant_def in allele_def.variants().iter() {
            let var_pos = variant_def.position();
            let var_ref = variant_def.reference().to_string();
            let var_alt = variant_def.alternate().to_string();
            let is_vi = variant_def.extras().get("VI").is_some();
            let var_key = (var_pos, var_ref.clone(), var_alt.clone());
            
            // mark if this is a VI variant
            if is_vi {
                vi_set.insert(var_key.clone());
            }
            all_set.insert(var_key.clone());
            
            if inserted_variants.contains(&var_key) {
                // do nothing, we already inserted this one
            } else {
                // create the variant and insert it
                let variant = if var_ref.len() == 1 {
                    if var_alt.len() == 1 {
                        Variant::new_snv(
                            0, var_pos as i64,
                            var_ref.into_bytes(), var_alt.into_bytes(),
                            0, 1)
                    } else {
                        Variant::new_insertion(
                            0, var_pos as i64,
                            var_ref.into_bytes(), var_alt.into_bytes(),
                            0, 1
                        )
                    }
                } else if var_alt.len() == 1 {
                    Variant::new_deletion(
                        0, var_pos as i64,
                        var_ref.len(), var_ref.into_bytes(), var_alt.into_bytes(),
                        0, 1)
                } else {
                    Variant::new_indel(
                        0, var_pos as i64,
                        var_ref.len(), var_ref.into_bytes(), var_alt.into_bytes(),
                        0, 1
                    )
                };
                variant_list.push(variant);

                // mark this key as inserted
                inserted_variants.insert(var_key);
            }
        }
    }

    // sort the variants so we can slice it up when we need to later
    variant_list.sort_by_key(|v| v.position());
    /*
    for v in variant_list.iter() {
        println!("{v:?}");
    }
    */
    let num_variants = variant_list.len();
    let first_variant_pos = variant_list[0].position() as u64;
    let last_variant_pos = variant_list.last().unwrap().position() as u64;
    debug!("Found {} unique variants for GraphWFA from chr22:{}-{}", num_variants, first_variant_pos+1, last_variant_pos+1);

    // build a lookup table for each of the variants and also the haplotypes
    let mut variant_lookup: HashMap<(usize, String, String), usize> = Default::default();
    // quick lookup from variant index to know if it is VI or not
    let mut is_vi_lookup: Vec<bool> = vec![false; num_variants];
    for (i, variant) in variant_list.iter().enumerate() {
        let var_key = (
            variant.position() as usize,
            String::from_utf8(variant.get_allele0().to_vec())?,
            String::from_utf8(variant.get_allele1().to_vec())?
        );

        // if this is VI, label it as such
        if vi_set.contains(&var_key) {
            is_vi_lookup[i] = true;
        }

        // now save this in our hashmap index
        variant_lookup.insert(var_key, i);
    }

    LoadedVariants::new(variant_list, variant_lookup, is_vi_lookup)
}

/// Contains an allele name, the region it was found, and the mapping stats.
#[derive(Clone, Debug)]
pub struct AlleleMapping {
    /// The allele that is mapped
    allele_label: Cyp2d6RegionLabel,
    /// The coordinate range inside it mapped to
    region: std::ops::Range<usize>,
    /// The score for the mapping
    mapping_stats: MappingStats
}

impl AlleleMapping {
    /// Constructor
    pub fn new(allele_label: Cyp2d6RegionLabel, region: std::ops::Range<usize>, mapping_stats: MappingStats) -> AlleleMapping {
        AlleleMapping {
            allele_label,
            region,
            mapping_stats
        }
    }

    // getters
    pub fn allele_label(&self) -> &Cyp2d6RegionLabel {
        &self.allele_label
    }

    pub fn region(&self) -> &std::ops::Range<usize> {
        &self.region
    }

    pub fn mapping_stats(&self) -> &MappingStats {
        &self.mapping_stats
    }
}

/// Given two ranges, this will report the overlap score of the ranges computed as shared / min(l1, l2).
/// By this definition, any non-overlapping get a score a 0.0, and then overlaps are scored based on the smaller region.
/// This means that if either region is _fully_ contained in the other, then the score will be 1.0.
/// # Arguments
/// * `r1` - the first range
/// * `r2` - the second range
fn overlap_score(r1: &std::ops::Range<usize>, r2: &std::ops::Range<usize>) -> f64 {
    let min_end = r1.end.min(r2.end);
    let max_start = r1.start.max(r2.start);
    if max_start >= min_end {
        0.0
    } else {
        let l1 = r1.len() as f64;
        let l2 = r2.len() as f64;
        let shared = (min_end - max_start) as f64;
        /*
        // old reciprocal overlap score, which is not the best here
        2.0 * shared / (l1+l2)
        */
        shared / l1.min(l2) 
    }
}

/// Returns an allele priority where higher values take precedence over others when the overlaps are high.
/// # Arguments
/// * `allele_id` - the allele ID we want to get priority for
fn get_allele_priority(allele_id: &Cyp2d6RegionLabel) -> usize {
    match allele_id.region_type() {
        Cyp2d6RegionType::Cyp2d6Deletion => 1, // *5 often has high overlap with other components; if detected, we should definitely prioritize it
        _ => 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /*
    TODO: tests that ideally exist but are difficult to implement:
    - Cyp2d6Extractor - this relies on a DB and a reference, can we encode this? the reference is the sticking point currently; big on disk, and long to load
        - find_base_type_in_sequence - same issue
        - find_full_type_in_sequence - same issue
        - assign_haplotype - same issue
    - generate_cyp_hybrids - requires reference, this is the root issue from above
    */

    #[test]
    fn test_load_variant_database() {
        // can be test on our real DB file
        let test_db_fn = std::path::PathBuf::from("./data/v0.9.0/cpic_20240404.json.gz");
        let database: PgxDatabase = crate::util::file_io::load_json(&test_db_fn).unwrap();
        let vcb = load_variant_database(&database).unwrap();

        // check all the high level, easy-to-verify stats
        assert_eq!(vcb.first_variant_pos(), 42126309);
        assert_eq!(vcb.last_variant_pos(), 42132374);
        assert_eq!(vcb.ordered_variants().len(), 387);
        assert_eq!(vcb.num_vi(), 144);
    }

    #[test]
    fn test_overlap_score() {
        assert_eq!(overlap_score(&(0..1), &(1..2)), 0.0); // no overlap
        assert_eq!(overlap_score(&(0..10), &(1..5)), 1.0); // fully contained
        assert_eq!(overlap_score(&(0..10), &(5..100)), 0.5); // half shared of first
        assert_eq!(overlap_score(&(15..100), &(0..20)), 0.25); // quarter shared of second
    }
}