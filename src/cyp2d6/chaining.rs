
use log::{debug, log_enabled, trace};
use minimap2::Aligner;
use simple_error::bail;
use std::collections::BTreeMap;
use waffle_con::multi_consensus::MultiConsensus;

use crate::cyp2d6::caller::convert_chain_to_hap;
use crate::cyp2d6::definitions::Cyp2d6Config;
use crate::cyp2d6::errors::{CallerError, CallerWarning};
use crate::cyp2d6::region_label::{Cyp2d6RegionLabel, Cyp2d6RegionType};
use crate::data_types::mapping::MappingStats;
use crate::util::stats::multinomial_ln_pmf;

/// This is a wrapper that will be the same length as the number of identified consensus sequences.
/// Additionally, the scores are edit distance and the overlap score
pub type SequenceWeights = Vec<(usize, f64)>;

/// Given a sequence and a multi-consensus, this will compare the sequence to each consensus and score it based on edits and overlaps.
/// This assumes that the full `sequence` will be used, so penalties for unmapped bases from `sequence` are included when scoring.
/// If the best found mapping has a score > 0.05, then all mappings are marked as equally bad.
/// This step is necessary to generate SequenceWeights prior to finding the best chain solution.
/// # Arguments
/// * `sequence` - the sequence to compare to the consensuses, we expect this to be fully represented
/// * `consensus` - the multi-consensus containing all allowed consensuses
/// * `con_labels` - human readable labels for the consensuses (e.g. star-alleles)
pub fn weight_sequence(sequence: &str, consensus: &MultiConsensus, con_labels: &[Cyp2d6RegionLabel]) -> Result<SequenceWeights, Box<dyn std::error::Error>> {
    let dna_aligner: Aligner = Aligner::builder()
        .map_hifi()
        .with_cigar()
        .with_seq(sequence.as_bytes())?;
    let seq_len = sequence.len();

    // we only need cigar and md for debugging
    // other settings for mapping
    let output_cigar: bool = true;
    let output_md: bool = true;
    let max_frag_len: Option<usize> = None;
    let extra_flags = None;
    
    // default is the length of the sequence getting "deleted" with 0 overlap
    let mut ret = vec![(seq_len, 0.0); con_labels.len()];

    // track the minimum observed NM-only ED; if it's larger than the maximum then this read does not map well to anything
    let penalize_unmapped = true; // we are mapping consensuses against a sequence from a read, penalize if it cannot fill out that sequence space
    let maximum_allowed_ed: f64 = 0.05; // TODO: I have a feeling we _should_ lower this, maybe down to 2-3%?
    let mut min_ed_frac: f64 = 1.0;
    
    for (con_index, (con, label)) in consensus.consensuses().iter().zip(con_labels.iter()).enumerate() {
        let con_seq = con.sequence();

        if !label.is_allowed_label() {
            // we ignore all Unknown and FalseAlleles
            continue;
        }
        
        // first, map the sequence
        let mappings = dna_aligner.map(
            con_seq,
            output_cigar, output_md, max_frag_len, extra_flags.clone()
        )?;

        for m in mappings.iter() {
            let con_len = con_seq.len();
            let nm = m.alignment.as_ref().unwrap().nm as usize;
            let unmapped = seq_len - (m.target_end - m.target_start) as usize;

            // the amount clipped at the start is the amount into query that we start
            let clipped_start = m.query_start as usize;
            // the amount clipped at the end is the length minus the query end point
            let clipped_end = con_len - m.query_end as usize;

            let mapping_stats = MappingStats::new_with_clippings(
                seq_len, nm, unmapped,
                clipped_start, clipped_end
            );

            // we are mapping consensuses against a sequence from a read, penalize if it cannot fill out that sequence space
            // this means we count both the nm() and unmapped() *against* each mapping
            let match_score = mapping_stats.nm() + mapping_stats.unmapped();
            let overlap_score = 1.0 - (clipped_start + clipped_end) as f64 / con_len as f64;
            let mapping_score = (match_score, overlap_score);
            
            debug!("\t\t{con_index}_{label} => ({:.4}, {:.4}) => {}", mapping_score.0, mapping_score.1, mapping_stats.custom_score_string(penalize_unmapped));
            
            // if the edit distance is less OR it's equal but the overlap score is higher
            if mapping_score.0 < ret[con_index].0 || (
                mapping_score.0 == ret[con_index].0 && mapping_score.1 > ret[con_index].1
            ) {
                ret[con_index] = mapping_score;
                min_ed_frac = min_ed_frac.min(mapping_stats.custom_score(penalize_unmapped).score());
            }
        }
    }

    if min_ed_frac <= maximum_allowed_ed {
        // we found something that is seemingly matching at least one sequence
        Ok(ret)
    } else {
        // all of the mappings are bad, send back empty vec for ignoring
        Ok(vec![])
    }
}

/// Contains the costs associated with each type of penalty
#[derive(Clone, Debug)]
pub struct ChainPenalties {
    /// this is the penalty for duplicating an allele, increasing this may lead to under-estimating of CN
    pub lasso_penalty: f64, 
    /// the log penalty for each edit
    pub ln_ed_penalty: f64,
    /// a penalty applied for an unexpected chain combination
    pub unexpected_chain_penalty: f64,
    /// penalty applied for each inferred edge
    pub inferred_edge_penalty: f64
}

impl Default for ChainPenalties {
    fn default() -> Self {
        Self { 
            lasso_penalty: 4.0,
            ln_ed_penalty: 2.0, // -(0.01_f64.ln())
            unexpected_chain_penalty: 10.0,
            inferred_edge_penalty: 2.0
        }
    }
}

impl ChainPenalties {
    /// Creates a new collection of chaining penalties with the given penalty values
    pub fn new(lasso_penalty: f64, ln_ed_penalty: f64, unexpected_chain_penalty: f64, inferred_edge_penalty: f64) -> ChainPenalties {
        ChainPenalties { 
            lasso_penalty,
            ln_ed_penalty,
            unexpected_chain_penalty,
            inferred_edge_penalty
        }
    }
}

/// Wrapper for chain scoring, switched to struct for ease of use and to prevent bugs.
/// API is lite because this is private.
struct ChainScore {
    /// The lasso cost of having too many or too few alleles; i.e. raw dups and dels have a cost
    pub allele_expected_penalty: f64,
    /// The number of best observations that are not met by the chain pair
    pub unmet_observations: u64,
    /// The number of edits that are required to make all observations match
    pub edit_distance: u64,
    /// The log penalty for the edit distance value
    pub ln_ed_penalty: f64,
    /// Multi-nomial likelihood penalty
    pub mn_llh_penalty: f64,
    /// Penalty for any chain events that are unexpected
    pub unexpected_chain_penalty: f64,
    /// Penalty for including inferred connections
    pub inferred_chain_penalty: f64,
    /// For debug, the allele IDs for the corresponding reduced probs and reduced coverage
    pub reduced_alleles: Vec<usize>,
    /// For debug, the probability of observing each allele
    pub reduced_probs: Vec<f64>,
    /// For debug, the actual observed counts (rounded)
    pub reduced_coverage: Vec<u64>,
    /// Associated chain index 1
    pub chain_index1: usize,
    /// Associated chain index 2
    pub chain_index2: usize,
}

impl ChainScore {
    /// The primary score metric for this ChainScore, represented as log likelihood.
    fn primary_score(&self) -> f64 {
        self.ln_ed_penalty + self.mn_llh_penalty + self.allele_expected_penalty + self.unexpected_chain_penalty + self.inferred_chain_penalty
    }

    /// Handy score string for debugging
    fn primary_score_string(&self) -> String {
        format!("{:.2} ({} ed / {} reads) + {:.2} (MNLLH) + {:.2} (CN) + {:.2} (unexp) + {:.2} (infer)",
            self.ln_ed_penalty, self.edit_distance, self.unmet_observations,
            self.mn_llh_penalty,
            self.allele_expected_penalty,
            self.unexpected_chain_penalty,
            self.inferred_chain_penalty
        )
    }

    /// Simplified compare function now.
    fn compare(&self, other: &ChainScore) -> std::cmp::Ordering {
        // all of them are within the boundaries, so just return whichever has the best cumulative score
        self.primary_score().partial_cmp(&other.primary_score()).unwrap()
    }
}

/// This will search through all the observed chaining counts and return the best combination of chains, as well as warnings if we encounter any dangling chains.
/// # Arguments
/// * `cyp2d6_config` - a generally static collection for the CYP2D6 configuration
/// * `obs_chains` - a map from sequence ID to a list of equally possible chains
/// * `chain_scores` - a map from sequence ID to scores for each consensus allele in the chain; this is a tuple of form (match_score, overlap_score) where 0.0 is worst, 1.0 is best
/// * `hap_labels` - a set of labels, primarily for debug output
/// * `infer_connections` - if True, then this will infer connections between alleles that do not have direct observations
/// * `normalize_all_alleles` - if True, all alleles are used to normalize coverage; otherwise, just those passing `is_normalizing_allele` will get used
/// * `penalities` - the set of penalties that describe the calculations we do for tie-breaking
/// * `ignore_chain_label_limits` - in prod, this should be false; but we use true for testing simplified chains
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
pub fn find_best_chain_pair(
    cyp2d6_config: &Cyp2d6Config,
    obs_chains: &BTreeMap<String, Vec<Vec<usize>>>, chain_scores: &BTreeMap<String, Vec<SequenceWeights>>,
    hap_labels: &[Cyp2d6RegionLabel],
    infer_connections: bool, normalize_all_alleles: bool,
    penalties: ChainPenalties,
    ignore_chain_label_limits: bool
) -> Result<(Vec<Vec<usize>>, Vec<CallerWarning>), Box<dyn std::error::Error>> {
    let mut caller_warnings: Vec<CallerWarning> = vec![];

    if penalties.lasso_penalty < 0.0 {
        bail!("Lasso penalty must be >= 0.0");
    }

    // first, identify all connections
    let num_haps = hap_labels.len();
    let mut downstream_possible: Vec<Vec<bool>> =  vec![vec![false; num_haps]; num_haps]; // 2D [upstream][downstream] -> edge present
    for (_qname, putative_chains) in obs_chains.iter() {
        for chain in putative_chains.iter() {
            if chain.len() > 1 {
                for i in 1..chain.len() {
                    let upstream = chain[i-1];
                    let downstream = chain[i];

                    // add in the downstream possibility for the upstream
                    if hap_labels[upstream].is_allowed_label() && hap_labels[downstream].is_allowed_label() {
                        if ignore_chain_label_limits || hap_labels[upstream].is_allowed_label_pair(&hap_labels[downstream]) {
                            downstream_possible[upstream][downstream] = true;
                        } else {
                            debug!("Ignoring observed chain: {upstream}_{} -> {downstream}_{}", hap_labels[upstream], hap_labels[downstream])
                        }
                    }
                }
            }
        }
    }

    // add any inferred connections
    let mut inferred_possible: Vec<Vec<bool>> =  vec![vec![false; num_haps]; num_haps]; // 2D [upstream][downstream] -> edge present
    let cyp_translate = cyp2d6_config.cyp_translate();
    if infer_connections {
        let detailed_inference = false; // we do not care if it is *4.001 or *4.002, just *4 works
        
        debug!("Inferred population connections:");
        let mut found_inference = false;
        for (i, h1) in hap_labels.iter().enumerate() {
            // old method goes from D6 -> D6 -> D7; we want to use the chain inferrences now though
            let h1_mod = h1.simplify_allele(detailed_inference, cyp_translate);
            
            // Note: This was a version that we tested for fixing the erroneous connection of *3 + *68 / *4.
            //       The problem here is that it creates situations where copy number can spiral due to a lack of checking the core alleles.
            //       The final solution is likely a combination of this with some more intelligent parsing.
            //       In the short term, this is an extreme edge case that is partly correct, we will push this to a later patch.
            let downstream_no_link = !downstream_possible[i].iter().any(|&b| b);

            for (j, h2) in hap_labels.iter().enumerate() {
                let upstream_no_link = !downstream_possible.iter().any(|v| v[j]);

                // simpler check, if the pairing is allowed, then we will infer it in the absence of any outbound edges
                if (downstream_no_link || upstream_no_link) && // make sure at least one end of this has no observed connections
                    !downstream_possible[i][j] && // probably redundant check
                    hap_labels[i].is_allowed_label() && 
                    hap_labels[j].is_allowed_label() && 
                    hap_labels[i].is_allowed_label_pair(&hap_labels[j]) {
                    // the labels are allowed to be a pair, so we will infer it as possible
                    let h2_mod = h2.simplify_allele(detailed_inference, cyp_translate);
                    debug!("\t{i}_{h1} ({h1_mod}) => {j}_{h2} ({h2_mod}) = inferred");
                    inferred_possible[i][j] = true;
                    found_inference = true;
                }
            }
        }

        if !found_inference {
            debug!("\tNone");
        }
    }

    // figure out which regions can start a chain
    let head_indices: Vec<usize> = hap_labels.iter().enumerate()
        .filter_map(|(i, label)| {
            if ignore_chain_label_limits || // if we are ignoring labels, then allow each thing to be a head index
                label.is_candidate_chain_head(normalize_all_alleles) { // this is a valid candidate head
                // add it to the list
                Some(i)
            } else {
                // not valid, so filter it out
                None
            }
        }).collect();
    
    // make sure we have chain starts, otherwise everything else will fail
    if head_indices.is_empty() {
        return Err(Box::new(CallerError::NoChainingHead));
    }

    // build all possible chains that can come from the head indices
    let mut remaining_chains: Vec<Vec<usize>> = head_indices.iter().map(|&start_index| vec![start_index]).collect();
    let mut possible_chains = vec![];
    let max_copy_number = 3; // TODO: this *could* be a CLI parameter in the future, we just need to be careful given that downstream application will balk
    while let Some(current_chain) = remaining_chains.pop() {
        // push the chain IF it ends with a D6/D7 allele (we expect to always end in D7 unless we have dropout)
        // let last_label = hap_labels[*current_chain.last().unwrap()].as_str();
        // if ignore_chain_label_limits || last_label.starts_with("CYP2D") {
        // TODO: can we add restrictions back in at some point?

        let (is_allowed_inferrence, is_allowed_candidate) = check_chain_inferrences(cyp2d6_config, &current_chain, hap_labels, &inferred_possible);
        if !is_allowed_inferrence {
            // there is an inferrence happening that is not allowed, discard this candidate entirely
            continue;
        }
        
        let simplified_chain = convert_chain_to_hap(&current_chain, hap_labels, true, cyp_translate);
        if ignore_chain_label_limits || (!simplified_chain.is_empty() && is_allowed_candidate) {
            // only add a chain if we are ignore labels OR 
            //     (if the chain produces a non-empty haplotype AND
            //      it is not over-inferring something)
            possible_chains.push(current_chain.clone());
        }
        
        let current_index = *current_chain.last().unwrap();
        for (extension_index, &extension_possible) in downstream_possible[current_index].iter().enumerate() {
            if !extension_possible {
                // extension with this isn't even allowed, so skip it
                continue;
            }
        
            let extension_count = current_chain.iter().filter(|&&v| v == extension_index).count();
            if extension_count >= max_copy_number {
                // we already have the maximum allowed copies of this (somewhere), don't allow any more in the chain
                // this allows for non-adjacency
                continue;
            }

            // these checks prevent infinite loops by restricting it to 2
            let mut new_chain = current_chain.clone();
            new_chain.push(extension_index);
            remaining_chains.push(new_chain);
        }

        if infer_connections {
            // we didn't found a direct downstream possibility, lets see if we can find an inferred connection
            // add an extension for each inferrence, they are filtered later
            for (extension_index, &extension_possible) in inferred_possible[current_index].iter().enumerate() {
                if !extension_possible {
                    // extension with this isn't even allowed, so skip it
                    continue;
                }
                
                let extension_count = current_chain.iter().filter(|&&v| v == extension_index).count();
                if extension_count >= 3 {
                    // we already have three copies of this (somewhere), don't allow a third copy at this time
                    // this allows for non-adjacency
                    continue;
                }
    
                // these checks prevent infinite loops by restricting it to 2
                let mut new_chain = current_chain.clone();
                new_chain.push(extension_index);
                remaining_chains.push(new_chain);
            }
        }
    }

    if possible_chains.is_empty() {
        // we did not find any valid chains, likely due to low coverage
        return Err(Box::new(CallerError::NoChainsFound));
    }

    debug!("Possible chains (N={}):", possible_chains.len());
    for chain in possible_chains.iter() {
        debug!("\t{chain:?}");
    }

    trace!("Multiple chain pairs detected, scoring them:");
    let mut score_sets: Vec<ChainScore> = vec![];
    // figure out the best by what is explained
    for i in 0..possible_chains.len() {
        for j in i..possible_chains.len() {
            let mut read_combined_ed: usize = 0;
            let mut hap_weights = vec![0.0; hap_labels.len()];

            for (_qname, chain_weights) in chain_scores.iter() {
                let (score, chain_match) = containment_score(&possible_chains[i], &possible_chains[j], chain_weights);
                // squared for RMS calculation - saturating is really only necessary for unit tests
                read_combined_ed = read_combined_ed.saturating_add(score);
                
                // divide support between the matches
                let split_frac = 1.0 / chain_match.len() as f64;
                for (chain_offset, &con_index) in chain_match.iter()
                    .flat_map(|chain| chain.iter().enumerate()) { // convert each chain into offset + consensus index iterator
                    // the split fraction scales the chain weight for this part of the read (chain_offset) when compared to the given consensus index
                    // then we select the overlap fraction which is the second value
                    hap_weights[con_index] += split_frac * chain_weights[chain_offset][con_index].1;
                }
            }

            // at this point we have the hap weights assigned, now divide by the allele counts
            let mut hap_counts = vec![0; hap_labels.len()];
            for &con_index in possible_chains[i].iter().chain(possible_chains[j].iter()) {
                hap_counts[con_index] += 1;
            }
            
            // this applies a penalty factor for each deviation from 1.0 copies
            let allele_expected_penalty: f64 = penalties.lasso_penalty * hap_labels.iter().zip(hap_counts.iter())
                .filter_map(|(label, hc)| {
                    if label.is_allowed_label() && // make sure the label is allowed
                        (ignore_chain_label_limits || // if we are ignoring label (usually debug only) OR
                            label.is_normalizing_allele(normalize_all_alleles) || // if the allele is for normalizing OR
                            label.is_reported_allele() // the allele is going to appear in the output
                        ) {
                        if *hc > 0 {
                            Some((*hc - 1) as f64) // we are expecting 1 copy, so add the delta
                        } else {
                            // this is an allowed allele, but we do not have any copies in our hap_count, so ignore it
                            // we used to count this as Some(1.0), indicating that an allele was expected but is missing
                            // now, we just let the edit penalty handle it
                            None
                        }
                    } else {
                        // this is not an allowed label, and we are not using it in normalization
                        None
                    }
                })
                .sum::<f64>();

            // count the number of chains that we cannot exactly find in this diplotype
            let mut unmet_observations = 0;
            for (_qname, putative_chains) in obs_chains.iter() {
                let supported = putative_chains.iter()
                    .any(|chain| {
                        // returns true if the putative chain is a sub-chain of either possible chain we're currently looking at
                        is_sub(&possible_chains[i], chain) ||
                        is_sub(&possible_chains[j], chain)
                    });
                
                // if nothing support the putative chains, then we have a read that is an unmet observation
                if !supported {
                    unmet_observations += 1;
                }
            }

            let edit_distance = read_combined_ed as u64;
            let ln_ed_penalty = (read_combined_ed as f64) * penalties.ln_ed_penalty;
            
            // now the multinomial calculation
            let mut reduced_alleles: Vec<usize> = vec![];
            let mut reduced_counts: Vec<i32> = vec![];
            let mut reduced_coverage: Vec<u64> = vec![];

            // one loop to populate all these vecs
            for (hap_index, hl) in hap_labels.iter().enumerate() {
                let hap_count = hap_counts[hap_index];
                if hap_count > 0 && (ignore_chain_label_limits || hl.is_normalizing_allele(normalize_all_alleles)) {
                    reduced_alleles.push(hap_index);
                    reduced_counts.push(hap_count);
                    let hap_weight = hap_weights[hap_index].round() as u64;
                    reduced_coverage.push(hap_weight);
                }
            }

            // now normalize the probabilities to 1.0
            let total_hap_count: i32 = reduced_counts.iter().sum();
            let reduced_probs: Vec<f64> = reduced_counts.into_iter()
                .map(|c| (c as f64) / (total_hap_count as f64))
                .collect();

            if reduced_probs.is_empty() || reduced_coverage.iter().sum::<u64>() == 0 {
                // this happens when we have a haplotype with no actual D6 alleles, it's not a valid result
                continue;
            }
            assert_eq!(reduced_probs.len(), reduced_coverage.len());
            
            /*
            // turns out Multinomial is not using all ln format under the hood so it overflows, we need a custom one
            let total_coverage: u64 = reduced_coverage.iter().sum();
            let multinomial = Multinomial::new(&reduced_probs, total_coverage)?;
            let mn_llh_penalty = -multinomial.ln_pmf(&reduced_coverage);
            */
            let mn_llh_penalty = multinomial_ln_pmf(&reduced_probs, &reduced_coverage).abs();
            
            // check if we have any unexpected chain pairs at the CYP2D level
            let expectation_mismatch = if ignore_chain_label_limits { 0 } 
                else { unexpected_count(&possible_chains[i], hap_labels, cyp2d6_config) + unexpected_count(&possible_chains[j], hap_labels, cyp2d6_config) };
            let unexpected_chain_penalty = (expectation_mismatch as f64) * penalties.unexpected_chain_penalty;

            // count up the number of edges that are from inferrence
            let num_inferred_edges = if infer_connections {
                let mut count = 0;
                for chain in [&possible_chains[i], &possible_chains[j]] {
                    for window in chain.windows(2) {
                        let h1 = window[0];
                        let h2 = window[1];
                        if inferred_possible[h1][h2] {
                            count += 1;
                        }
                    }
                }
                count
            } else {
                0
            };
            let inferred_chain_penalty = (num_inferred_edges as f64) * penalties.inferred_edge_penalty;
            
            // trace!("\t{:?} {:?} => {:?} + {} + {} + {} => {:?}", possible_chains[i], possible_chains[j], combined_explained, allele_delta_penalty, allele_expected_penalty, full_chain_penalty, allele_counts);
            score_sets.push(ChainScore {
                allele_expected_penalty,
                unmet_observations,
                edit_distance,
                ln_ed_penalty,
                mn_llh_penalty,
                unexpected_chain_penalty,
                inferred_chain_penalty,
                reduced_alleles,
                reduced_probs,
                reduced_coverage,
                chain_index1: i,
                chain_index2: j,
            });
        }
    }

    score_sets.sort_by(|a, b| {
        // check if the combined explanations are "close enough" for equality
        // TODO: now that we have a unified score, we could probably convert this into proper compare with PartialCmp and Cmp; low priority
        a.compare(b)
    });
    score_sets.reverse();

    if log_enabled!(log::Level::Debug) {
        let num_shown = 50;
        debug!("Scored chain pairs (best {num_shown}):");
        let skip_count = if score_sets.len() > num_shown { score_sets.len() - num_shown } else { 0 };
        for chain_score in score_sets.iter().skip(skip_count) {
            debug!("\t{:?} {:?} => {:.4} ({})",
                possible_chains[chain_score.chain_index1], 
                possible_chains[chain_score.chain_index2], 
                chain_score.primary_score(),
                chain_score.primary_score_string()
            );
            debug!("\t\talleles: {:?}, probs: {:?}; obs: {:?}",
                chain_score.reduced_alleles,
                chain_score.reduced_probs,
                chain_score.reduced_coverage
            );
        }
    }

    if score_sets.is_empty() {
        // this shouldn't be possible anymore with the prior checks, but lets put an error catch here just in case
        return Err(Box::new(CallerError::NoScorePairs));
    }

    // we sorted them above from worst to best
    let highest_score = score_sets.last().unwrap();
    let best_pair = (highest_score.chain_index1, highest_score.chain_index2);

    let mut best_chain_pair = vec![
        possible_chains[best_pair.0].clone(),
        possible_chains[best_pair.1].clone()
    ];
    best_chain_pair.sort();
    let best_chains = best_chain_pair;

    // build all the chains for each head index
    let mut index_used: Vec<bool> = vec![false; num_haps];
    for chain in best_chains.iter() {
        for &i in chain.iter() {
            index_used[i] = true;
        }
    }

    // now go through the used indices and report any that were not used as dangling alleles (i.e., unchained to the result)
    for (i, &b) in index_used.iter().enumerate() {
        if !b {
            // this one was not used
            caller_warnings.push(CallerWarning::DanglingAllele { allele_name: format!("{}_{}", i, hap_labels[i]) });
        }
    }

    Ok((best_chains, caller_warnings))
}

/// Checks a putative haplotype chain for correct inferrences.
/// Returns a boolean of the form `(is_allowed_inferrence, is_allowed_candidate)` where
/// `is_allowed_inferrence` is only true if there are no inferrence issues in the most recent D6 pairing and
/// `is_allowed_candidate` is only true if it ends with an allowed CYP2D allele OR there are no inferrences since the most recent CYP2D allele.
/// # Arguments
/// * `cyp2d6_config` - a generally static collection for the CYP2D6 configuration
/// * `chain` - the putative chain to check
/// * `hap_labels` - the haplotype labels we need for interpreting the chain
/// * `inferred_connections` - contains the list of inferred chain links
fn check_chain_inferrences(cyp2d6_config: &Cyp2d6Config, chain: &[usize], hap_labels: &[Cyp2d6RegionLabel], inferred_connections: &[Vec<bool>]) -> (bool, bool) {
    // sanity check
    assert!(!chain.is_empty());

    // config items
    let cyp_translate = cyp2d6_config.cyp_translate();

    // check if the last link is a CYP2D allele or not
    let last_hap_index = *chain.last().unwrap();
    let last_is_cyp2d = hap_labels[last_hap_index].is_cyp2d();

    // get the most recent CYP2D allele before last
    let mut opt_index = None;
    for (chain_index, &hap_index) in chain.iter().enumerate().rev().skip(1) {
        if hap_labels[hap_index].is_cyp2d() {
            opt_index = Some(chain_index);
            break;
        }
    }

    // finally, check if we have an inferrence in between
    let mut inferrence_detected = false;
    for window in chain[opt_index.unwrap_or(0)..(chain.len())].windows(2) {
        let h1 = window[0];
        let h2 = window[1];
        if inferred_connections[h1][h2] {
            // this one was inferred
            inferrence_detected = true;
        }
    }

    if inferrence_detected {
        // we found an inferrence, we need to verify that it's an okay one
        if last_is_cyp2d {
            // the last allele in the chain is a D6 allele, verify that the previous D6 allele is valid
            if let Some(previous_d6_index) = opt_index {
                let detailed_inference = false; // we do not care if it is *4.001 or *4.002, just *4 works
                let infer_d7_tail = true; // if true, then we infer connections from D6 to D7

                // build out the allele labels
                let previous_hap_index = chain[previous_d6_index];
                let h1 = &hap_labels[previous_hap_index];
                let h2 = &hap_labels[last_hap_index];
                let h1_mod = h1.simplify_allele(detailed_inference, cyp_translate);
                let h2_mod = h2.simplify_allele(detailed_inference, cyp_translate);

                let connected = previous_hap_index != last_hap_index && // we do not allow inferred exact duplications
                    cyp2d6_config.inferred_connections().contains(&(h1_mod.clone(), h2_mod.clone())); // check if it is otherwise known
                let d7_tail_connection = infer_d7_tail && 
                    // must be a D7 link
                    h2.region_type() == Cyp2d6RegionType::Cyp2d7 && 
                    // must be a non-D7 CYP2D allele
                    h1.region_type() != Cyp2d6RegionType::Cyp2d7 &&
                    h1.region_type().is_cyp2d();

                // allowed IF either it's a known inferred connection OR it's a D7 tail
                let allowed = connected || d7_tail_connection;
                (allowed, allowed)
            } else {
                // there is no previous D6 allele, so we should be okay on both counts
                (true, true)
            }
        } else {
            // we do not end with CYP2D allele and we have an inferrence
            // in this case, it is not breaking any rules, but it also should not be allowed to be a full chain for scoring
            (true, false)
        }
    } else {
        // we did not detect an inferrence, so we can return true on both counts
        (true, true)
    }
}

/// This return the best possible score for a collection of putative `chain_weights` when compared to two candidate chain sets.
/// It also returns a list of the chains that can generate that best score.
/// This score is shifted relative to the absolute best score (if we ignored chaining), such that exactly matching the best will return a "0" score.
/// # Arguments
/// * `chain_set1` - the first chain set to search
/// * `chain_set2` - the second chain set to search
/// * `chain_weights` - an observed chain and the collected scores; each allele in the chain has a Vec of score tuples (edit distance, overlap); edit distance is raw integer and overlap is 1.0 max
fn containment_score(chain_set1: &[usize], chain_set2: &[usize], chain_weights: &[SequenceWeights]) -> (usize, Vec<Vec<usize>>) {
    // calculate the *best* possible score
    let optimum_weight: usize = chain_weights.iter()
        .map(|scores| {
            scores.iter()
                .map(|(w, _o)| *w)
                .min().unwrap()
        }).sum();

    // calculate the *worst* as well
    let worst_weight: usize = chain_weights.iter()
        .map(|scores| {
            scores.iter()
                .map(|(w, _o)| *w)
                .max().unwrap()
        }).sum();

    // initialize to double the worst score for the sake of distinguishing
    let mut best_score = 2*worst_weight;
    let mut best_chains = vec![];
    let weight_len = chain_weights.len();

    for other in [chain_set1, chain_set2] {
        if other.len() < weight_len {
            // this one is too short to ever represent the full chain
            continue;
        }

        for start_index in 0..(other.len() - weight_len+1) {
            let total_weight: usize = other[start_index..(start_index+weight_len)].iter().zip(chain_weights.iter())
                .map(|(&i, weight_dict)| {
                    weight_dict[i].0
                })
                .sum();
            
            if total_weight < best_score {
                best_score = total_weight;
                best_chains.clear();
            }
            if total_weight == best_score {
                best_chains.push(other[start_index..(start_index+weight_len)].to_vec());
            }
        }
    }

    assert!(best_score >= optimum_weight);

    (best_score - optimum_weight, best_chains)
}

/// This will look at a chain with haplotype labels and report out the total number of chain pairs that are not expected.
/// Checks for chains that are empty or that do not start with CYP2D6.
/// Also checks for any chain links that are not a part of our inferred set (meaning, we probably are not expecting it).
/// # Arguments
/// * `chain` - the chain to scan
/// * `hap_labels` - labels from index to String for each haplotype
fn unexpected_count(chain: &[usize], hap_labels: &[Cyp2d6RegionLabel], cyp2d6_config: &Cyp2d6Config) -> u32 {
    let reduced_chain: Vec<String> = chain.iter()
        .filter_map(|&c_index| {
            // remove all D7 alleles
            if hap_labels[c_index].is_cyp2d() && hap_labels[c_index].region_type() != Cyp2d6RegionType::Cyp2d7 {
                // convert the allele index into a human readable string
                let string_label = hap_labels[c_index].simplify_allele(false, cyp2d6_config.cyp_translate());
                Some(string_label)
            } else {
                None
            }
        }).collect();
    
    let mut errors_detected = 0;
    if reduced_chain.is_empty() || // we do not expect something to be empty, that's for sure
        !reduced_chain[0].starts_with('*') {  // we expect a reduced chain to start with D6
        errors_detected += 1;
    } else {
        // chain is not empty AND the chain starts with a star-allele, so should be okay on this front
    };

    // check if there is an unexpected singleton
    if reduced_chain.len() == 1 && cyp2d6_config.unexpected_singletons().contains(&reduced_chain[0]) {
        // this is a single allele that we do not expect to find as a single allele under most circumstances
        errors_detected += 1;
    }

    // check for any aberrant chain pairs
    for chain_pair in reduced_chain.windows(2) {
        if !cyp2d6_config.inferred_connections().contains(&(chain_pair[0].clone(), chain_pair[1].clone())) {
            // this is not a connection we would normally infer, so it's unexpected
            errors_detected += 1;
        }
    }

    errors_detected
}

/// From here: https://stackoverflow.com/questions/47043167/does-rust-contain-a-way-to-directly-check-whether-or-not-one-vector-is-a-substr
/// This is a simple function that will check if one slice (the haystack) contains another (the needle)
/// # Arguments
/// * `haystack` - the slice to search
/// * `needle` - the sub-slice to search for in the haystack
fn is_sub<T: PartialEq>(haystack: &[T], needle: &[T]) -> bool {
    haystack.windows(needle.len()).any(|c| c == needle)
}

#[cfg(test)]
mod tests {
    use super::*;

    use waffle_con::cdwfa_config::ConsensusCost;
    use waffle_con::consensus::Consensus;

    /// Wrapper function that will generate CYP2D6 chains for us as "reads" that always span just two regions of interest.
    /// User provides a list of haplotypes and the observed chains, then scores are auto-generated to match.
    /// This only generates constant values where the "match" is a perfect match and all others have a constant ED.
    /// # Arguments
    /// * `hap_labels` - the labels for the haplotypes
    /// * `chains` - all observed chains, these can just be expected diplotypes for ease of use
    fn create_pairwise_chains(hap_labels: &[Cyp2d6RegionLabel], chains: &[Vec<usize>]) -> (BTreeMap<String, Vec<Vec<usize>>>, BTreeMap<String, Vec<SequenceWeights>>) {
        // contains the best observed chains for each
        let mut obs_chains: BTreeMap<String, Vec<Vec<usize>>> = Default::default();
        // contains scores against all alleles
        let mut chain_scores: BTreeMap<String, Vec<SequenceWeights>> = Default::default();

        let mut read_index: usize = 0;
        for chain in chains.iter() {
            assert!(chain.len() >= 2);
            // create a read for each pair
            for window in chain.windows(2) {
                let seq_name = format!("read_{read_index}");
                obs_chains.insert(seq_name.clone(), vec![window.to_vec()]);

                // for each link in the chain, create penalties for all other comparators
                let mut weights = vec![];
                for &hap_index in chain.iter() {
                    // set all scores to "bad"; 100 ED and fully overlapping
                    let mut all_scores = vec![(100, 1.0); hap_labels.len()];
                    // set the one match to "good"; 0 ED
                    all_scores[hap_index] = (0, 1.0);
                    weights.push(all_scores);
                }
                chain_scores.insert(seq_name, weights);
                read_index += 1;
            }
        }

        (obs_chains, chain_scores)
    }

    #[test]
    fn test_find_best_chain_pair() {
        let hap_labels = vec![
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("A".to_string())),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("B".to_string())),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("C".to_string())),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("D".to_string()))
        ];
        let obs_chains = [
            // 0 -> 1 -> 2
            ("seq_1".to_string(), vec![vec![0, 2]]),
            ("seq_2".to_string(), vec![vec![1, 1]])
        ].into_iter().collect();
        let chain_scores = [
            ("seq_1".to_string(), vec![
                vec![(0, 1.0), (1, 1.0), (1, 1.0), (1, 1.0)],
                vec![(1, 1.0), (1, 1.0), (0, 1.0), (1, 1.0)]
            ]),
            ("seq_2".to_string(), vec![
                vec![(1, 1.0), (0, 1.0), (1, 1.0), (1, 1.0)],
                vec![(1, 1.0), (0, 1.0), (1, 1.0), (1, 1.0)]
            ])
        ].into_iter().collect();
        let penalties = Default::default();
        let cyp2d6_config = Cyp2d6Config::default();
        let (chain_result, danglers) = find_best_chain_pair(&cyp2d6_config, &obs_chains, &chain_scores, &hap_labels, false, true, penalties, true).unwrap();
        assert_eq!(chain_result, vec![
            vec![0, 2],
            vec![1, 1]
        ]);
        assert_eq!(danglers, vec![CallerWarning::DanglingAllele { allele_name: "3_CYP2D6*D".to_string() }]);
    }

    #[test]
    fn test_ambiguous_find_best_chain_pair() {
        let hap_labels = vec![
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("A".to_string())),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("B".to_string()))
        ];
        let obs_chains = [
            // B -> A -> A
            ("seq_0".to_string(), vec![vec![1]]),
            ("seq_1".to_string(), vec![vec![1, 0]]),
            ("seq_2".to_string(), vec![vec![0, 0]]),
            ("seq_3".to_string(), vec![vec![0]]),
            // B -> A
            ("seq_4".to_string(), vec![vec![1]]),
            ("seq_5".to_string(), vec![vec![1, 0]]),
            ("seq_6".to_string(), vec![vec![0]])
        ].into_iter().collect();
        
        // scores are (edit distance, overlap score)
        let chain_scores = [
            ("seq_0".to_string(), vec![
                vec![(10, 1.0), (0, 1.0)]
            ]),
            ("seq_1".to_string(), vec![
                vec![(10, 1.0), (0, 1.0)],
                vec![(0, 1.0), (10, 1.0)]
            ]),
            ("seq_2".to_string(), vec![
                vec![(0, 1.0), (10, 1.0)],
                vec![(0, 1.0), (10, 1.0)]
            ]),
            ("seq_3".to_string(), vec![
                vec![(0, 1.0), (10, 1.0)]
            ]),
            ("seq_4".to_string(), vec![
                vec![(10, 1.0), (0, 1.0)]
            ]),
            ("seq_5".to_string(), vec![
                vec![(10, 1.0), (0, 1.0)],
                vec![(0, 1.0), (10, 1.0)]
            ]),
            ("seq_6".to_string(), vec![
                vec![(0, 1.0), (10, 1.0)]
            ])
        ].into_iter().collect();

        // no lasso penalty
        let penalties = ChainPenalties::new(0.0, -(0.01_f64.ln()), 0.0, 2.0);
        let cyp2d6_config = Cyp2d6Config::default();
        let (chain_result, danglers) = find_best_chain_pair(&cyp2d6_config, &obs_chains, &chain_scores, &hap_labels, false, true, penalties, true).unwrap();
        assert_eq!(chain_result, vec![
            vec![1],
            vec![1, 0, 0, 0] // without lasso, it will just greedily add these
        ]);
        assert_eq!(danglers, vec![]);

        // now test with lasso penalty, it should drop the first 0
        let penalties = ChainPenalties::new(3.0, -(0.01_f64.ln()), 0.0, 2.0);
        let (chain_result, danglers) = find_best_chain_pair(&cyp2d6_config, &obs_chains, &chain_scores, &hap_labels, false, true, penalties, true).unwrap();
        assert_eq!(chain_result, vec![
            vec![1],
            vec![1, 0, 0] // with lasso, it restrict to what is observed
        ]);
        assert_eq!(danglers, vec![]);
    }

    #[test]
    fn test_weight_sequence() {
        // I think this is doable, but we will see
        let consensus = MultiConsensus::new(
            vec![
                Consensus::new(b"AGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGAACATGGGATTATGGGCAAGGGTAACAGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGGACATGGGATTATGGGCAAGGGTAAC".to_vec(), ConsensusCost::L1Distance, vec![0]),
                Consensus::new(b"AGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGCACATGGGATTATGGGCAAGGGTAACAGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGGACATGGGATTATGGGCAAGGGTAAC".to_vec(), ConsensusCost::L1Distance, vec![0]),
                Consensus::new(b"AGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGGACATGGGATTATGGGCAAGGGTAACAGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGGACATGGGATTATGGGCAAGGGTAAC".to_vec(), ConsensusCost::L1Distance, vec![0])
                //                                                                                                  ^change here
            ],
            vec![0] // these don't matter for the test
        );
        let con_labels = vec![
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("A".to_string())),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("C".to_string())),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("G".to_string())),
        ];

        // test mostly match, but best is the first entry
        let sequence = "AGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGAACATGGGATTATGGGCAAGGGTAACAGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGGACATGGGATTATGGGCAAGGGTAAC";
        let score = weight_sequence(&sequence, &consensus, &con_labels).unwrap();
        assert_eq!(score.iter().min_by(|a, b| {
            a.partial_cmp(b).unwrap()
        }).unwrap(), &score[0]);

        // test equal match (G -> T) to each
        let sequence = "AGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGNACATGGGATTATGGGCAAGGGTAACAGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGGACATGGGATTATGGGCAAGGGTAAC";
        let score = weight_sequence(&sequence, &consensus, &con_labels).unwrap();
        assert_eq!(score[0], score[1]);
        assert_eq!(score[0], score[2]);
    }

    #[test]
    fn test_inferred_alleles() {
        let hap_labels = [
            // *3 -> D7
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("3".to_string())),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::LinkRegion, None),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Rep7, None),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Spacer, None),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d7, None),
            // *4 -> *68, but through the same link types
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, Some("4".to_string())),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Hybrid, Some("CYP2D6::CYP2D7::exon2".to_string()))
        ];
        let chains = [
            // two components for *3 -> link_region; <BREAK> Rep7 -> D7
            vec![0, 1],
            vec![2, 3, 4],
            // two components for *4 -> link_region; <BREAK> Rep7 -> *68
            // inferrence will be required to connect them correctly
            vec![5, 1],
            vec![2, 3, 6]
        ];

        // create pairwise chains from the above
        let (obs_chains, chain_scores) = create_pairwise_chains(&hap_labels, &chains);

        // first, test with no inferrence, which should just find *3 / *4
        let infer = false;
        let penalties: ChainPenalties = Default::default();
        let cyp2d6_config = Cyp2d6Config::default();
        let (chain_result, danglers) = find_best_chain_pair(&cyp2d6_config, &obs_chains, &chain_scores, &hap_labels, infer, true, penalties.clone(), false).unwrap();
        assert_eq!(chain_result, vec![
            vec![0, 1], // *3
            vec![5, 1] // *4
        ]);
        // no inferrence leads to a bunch of danglers
        assert_eq!(danglers, vec![
            CallerWarning::DanglingAllele { allele_name: "2_REP7".to_string() },
            CallerWarning::DanglingAllele { allele_name: "3_spacer".to_string() },
            CallerWarning::DanglingAllele { allele_name: "4_CYP2D7".to_string() },
            CallerWarning::DanglingAllele { allele_name: "6_CYP2D6::CYP2D7::exon2".to_string() }
        ]);

        // now add in inferrence and re-test
        // first, test with no inferrence, which should just find *3 / *4
        let cyp2d6_config = Cyp2d6Config::default();
        let infer = true;
        let (chain_result, danglers) = find_best_chain_pair(&cyp2d6_config, &obs_chains, &chain_scores, &hap_labels, infer, true, penalties, false).unwrap();
        assert_eq!(chain_result, vec![
            vec![0, 1, 2, 3, 4], // *3 + D7
            vec![5, 1, 2, 3, 6] // *4 + *68, this should get inferred as the correct solution here
        ]);
        assert!(danglers.is_empty()); // should be no danglers now
    }

    #[test]
    fn test_chaining_errors() {
        // defaults that do not matter for this test
        let cyp2d6_config = Cyp2d6Config::default();
        let obs_chains = Default::default();
        let chain_scores = Default::default();
        let infer = false;
        let penalties = Default::default();
        
        // simple case where we just give the algorithm a bunch of non-starting hap labels
        let hap_labels = vec![
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d7, None),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::LinkRegion, None),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Spacer, None),
            Cyp2d6RegionLabel::new(Cyp2d6RegionType::Unknown, None)
        ];
        let result = find_best_chain_pair(&cyp2d6_config, &obs_chains, &chain_scores, &hap_labels, infer, true, penalties, false);
        assert!(result.is_err());
        assert_eq!(result.err().unwrap().downcast::<CallerError>().unwrap().as_ref(), &CallerError::NoChainingHead);

        // possible future TODO: we have two other error types, I don't think we can actually reach them though
        //                       we will lazily add if a user manages to do it somehow
    }
}