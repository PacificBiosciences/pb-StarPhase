
use itertools::Itertools;
use log::warn;
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use simple_error::bail;
use std::collections::BTreeMap;
use std::path::Path;
use waffle_con::multi_consensus::MultiConsensus;

use crate::cyp2d6::region::Cyp2d6Region;
use crate::cyp2d6::region_label::Cyp2d6RegionType;
use crate::data_types::coordinates::Coordinates;
use crate::data_types::database::PgxDatabase;
use crate::util::mapping::standard_hifi_aligner;
use crate::visualization::igv_session_writer::{BUFFER_LEN, CONTIG_POSTFIX};

/// Accepts information about available alleles as well as the reads spanning those alleles and converts it into a visual graph.
/// # Arguments
/// * `hap_regions` - a list of haplotype region labels we can translate into visual Strings
/// * `chain_frequency` - this is the frequency of each observed chain from the data
/// * `filename` - the output file path to save the results
pub fn generate_debug_graph(hap_regions: &[Cyp2d6Region], chain_frequency: &BTreeMap<Vec<usize>, f64>, filename: &Path) -> Result<(), Box<dyn std::error::Error>> {
    use layout::backends::svg::SVGWriter;
    use layout::core::base::Orientation;
    use layout::core::color::Color;
    use layout::core::geometry::Point;
    use layout::core::style::{LineStyleKind, StyleAttr};
    use layout::core::utils::save_to_file;
    use layout::std_shapes::shapes::{Arrow, Element, LineEndKind, ShapeKind, RecordDef};
    use layout::topo::layout::VisualGraph;

    // first lets part the chain frequencies into totals for edges and node
    let mut single_counts = vec![0.0; hap_regions.len()];
    let mut pair_counts: BTreeMap<(usize, usize), f64> = Default::default();

    for (chain, frequency) in chain_frequency.iter() {
        for index in chain.iter() {
            single_counts[*index] += *frequency;
        }

        for window in chain.windows(2) {
            let entry = pair_counts.entry((window[0], window[1])).or_default();
            *entry += *frequency;
        }
    }

    // the graph we're going to fill in
    let mut vg = VisualGraph::new(Orientation::LeftToRight);

    // create all the nodes with a count
    let mut node_handles = vec![];
    for (i, hr) in hap_regions.iter().enumerate() {
        if hr.label().is_allowed_label() {
            let shape = ShapeKind::Record(
                RecordDef::Array(vec![
                    RecordDef::new_text(&i.to_string()),
                    RecordDef::new_text(&hr.label().full_allele()),
                    RecordDef::new_text(&format!("{:.2}", single_counts[i]))
                ])
            );
            let look = StyleAttr::simple();
            let sz = Point::new(175.0, 100.0);

            let node = Element::create(shape, look, Orientation::TopToBottom, sz);
            let handle = vg.add_node(node);
            node_handles.push(Some(handle));
        } else {
            node_handles.push(None);
        }
    }

    let min_arrow_width = 2;
    let max_arrow_width = 5;
    let min_edge_size = *pair_counts.values().min_by(|a, b| a.total_cmp(b)).unwrap_or(&0.0);
    let max_edge_size = (*pair_counts.values().max_by(|a, b| a.total_cmp(b)).unwrap_or(&0.0)) // get the maximum
        .max(min_edge_size+1.0); // this makes sure the maximum > minimum by at least 1.0

    // create the edges between them
    for (chain_pair, frequency) in pair_counts.iter() {
        // figure out a scaled width
        let arrow_fraction = (*frequency - min_edge_size) / (max_edge_size - min_edge_size);
        let inferred_width = min_arrow_width + ((max_arrow_width - min_arrow_width) as f64 * arrow_fraction).round() as usize;

        // figure out the total color on the heatmap space
        let red_component = ((arrow_fraction * 255.0).floor() as u32) << 16;
        let blue_component = ((1.0 - arrow_fraction) * 255.0).floor() as u32;
        let inferred_color = red_component + blue_component;
        
        // for some reason, we have to shift and add 0xff
        let color = Color::new((inferred_color << 8) + 0xff);
        
        // Add an edge between the nodes.
        let arrow = Arrow::new(
            LineEndKind::None,
            LineEndKind::Arrow,
            LineStyleKind::Normal,
            &format!("{frequency:.2}"),
            &StyleAttr::new(
                color,
                inferred_width,
                Option::Some(Color::fast("white")),
                0,
                15
            ),
            &None,
            &None,
        );
        vg.add_edge(arrow, node_handles[chain_pair.0].unwrap(), node_handles[chain_pair.1].unwrap());
    }

    // there is a rogue panic in the layout-rs crate that seems rare, unclear how to reproduce it yet
    // similar issue submitted to layout-rs in 2022, but no response from devs
    // TODO: catch_unwind is clearly a stop-gap, we don't want this in here forever
    //       I like this crate for ease-of-use, but we may need something more robust with active maintainers
    #[allow(clippy::blocks_in_conditions)]
    match std::panic::catch_unwind(|| {
        let mut vg = vg; // need to move it inside this scope for rust to be happy (panic unwinding)
        
        // Render the nodes to some rendering backend.
        let mut svg = SVGWriter::new();
        vg.do_it(false, false, false, &mut svg);
        
        // Save the output.
        save_to_file(filename.as_os_str().to_str().unwrap(), &svg.finalize())
    }) {
        Ok(_v) => Ok(()),
        Err(e) => bail!("Received panic while writing SVG file: {:?}", e.downcast_ref::<&str>())
    }
}

/// This is just a wrapper for the output from the next function, mainly to make clippy happy.
pub struct CustomReference {
    /// name for the contig
    pub contig_name: String,
    /// the full custom sequence
    pub sequence: String,
    /// list of coordinates / label pairs
    pub regions: Vec<(Coordinates, String)>
}

/// Creates a customized reference genome sequence from our consensus, which we can output to file and use for IGV visuals.
/// This function is specific to CYP2D6 and the database config for it.
/// # Arguments
/// * `reference_genome` - the actual reference genome data (GRCh38)
/// * `database` - the config we loaded
/// * `consensus` - the multi-consensus result, which contains a bunch of sub-units from the D6 region
/// * `hap_regions` - the assigned region labels for each consensus
/// * `best_result` - the best chain pair; should be two Vecs of unknown length
/// # Errors
/// * if there are UTF-8 parsing errors
pub fn create_custom_cyp2d6_reference(
    reference_genome: &ReferenceGenome, database: &PgxDatabase,
    consensus: &MultiConsensus, hap_regions: &[Cyp2d6Region], best_result: &[Vec<usize>]
) -> Result<CustomReference, Box<dyn std::error::Error>> {
    // gene name followed by a post-fix for easier searching
    let contig_name = format!("CYP2D6_{CONTIG_POSTFIX}");
    
    // generic buffer between regions
    let buffer_sequence: String = "N".repeat(BUFFER_LEN);
    
    // start with a buffer
    let mut ret = buffer_sequence.clone();
    let mut regions: Vec<(Coordinates, String)> = vec![];

    for haplotype_chain in best_result.iter() {
        // process the first element
        let hap_index = haplotype_chain[0];
        let hap_sequence = std::str::from_utf8(consensus.consensuses()[hap_index].sequence())?;

        // add the coordinates of what was added
        let coordinates = Coordinates::new(contig_name.clone(),
            ret.len() as u64,
            (ret.len() + hap_sequence.len()) as u64
        );
        let region_name = format!("{hap_index}_{}", hap_regions[hap_index]);
        regions.push((coordinates, region_name));

        // now extend our region
        ret.push_str(hap_sequence);

        // handle everything else in pairs so we can check for gaps to fill
        for (&prev_index, &hap_index) in haplotype_chain.iter().tuple_windows() {
            // get the types
            let prev_type = hap_regions[prev_index].label().region_type();
            let hap_type = hap_regions[hap_index].label().region_type();
            let mut overlap_len = 0; // for almost every case, there is no overlap
            if prev_type.is_rep() && hap_type.is_cyp2d() {
                // we have a gap here to fill from the end of REP6 to the start of D6
                let chrom = "chr22";
                let rep6_end = database.cyp2d6_config().cyp_coordinates().get("REP6").unwrap().end() as usize;
                let d6_start = database.cyp2d6_config().cyp_coordinates().get("CYP2D6").unwrap().start() as usize;
                let gap_sequence = std::str::from_utf8(reference_genome.get_slice(chrom, rep6_end, d6_start))?;
                ret.push_str(gap_sequence);
            } else if prev_type == Cyp2d6RegionType::Spacer && hap_type.is_cyp2d() {
                // we have a gap here to fill from the end of REP7 to the start of D7
                let chrom = "chr22";
                let rep7_end = database.cyp2d6_config().cyp_coordinates().get("spacer").unwrap().end() as usize;
                let d7_start = database.cyp2d6_config().cyp_coordinates().get("CYP2D7").unwrap().start() as usize;
                let gap_sequence = std::str::from_utf8(reference_genome.get_slice(chrom, rep7_end, d7_start))?;
                ret.push_str(gap_sequence);
            } else if prev_type == Cyp2d6RegionType::Cyp2d6Deletion && hap_type == Cyp2d6RegionType::Spacer {
                // figure out the overlap
                let align_window = 500;
                let align_start = ret.len().saturating_sub(align_window); // shouldn't go below 0, but this protects us
                let align_sequence = &ret[align_start..];
                
                // build the aligner to our smaller tail region
                let dna_aligner = standard_hifi_aligner()
                    .with_seq(align_sequence.as_bytes())?;

                // we only need cigar and md for debugging
                // other settings for mapping
                let output_cigar: bool = true;
                let output_md: bool = true;
                let max_frag_len: Option<usize> = None;
                let extra_flags = None;

                // build the target
                let full_query = consensus.consensuses()[hap_index].sequence();
                let query_end = full_query.len().min(align_window);
                let query_sequence = &full_query[..query_end];

                // first, map the sequence
                let mappings = dna_aligner.map(
                    query_sequence,
                    output_cigar, output_md, max_frag_len, extra_flags, None
                )?;

                // finally search for a matching overlap and save it if we find one
                let bp_thresh = 5; // lets give ourselves a small wiggle room
                for m in mappings.iter() {
                    // we want a mapping that starts at 0 in our query and end at the end of the target
                    // expected length is ~230 bp or so
                    if m.query_start <= bp_thresh && (m.target_len - m.target_end) <= bp_thresh {
                        // this one is a pretty good match, save the overlap and skip out
                        overlap_len = m.query_end as usize;
                        break;
                    }
                }

                if overlap_len == 0 {
                    // we didn't find an overlap
                    warn!("No overlap found between adjacent *5 and spacer region, output reference may have unexplained gaps.");
                }
            } else {
                // TODO: other types to handle?
            }

            // first pass, naively extend by each sequence in the chain
            let hap_sequence = std::str::from_utf8(&consensus.consensuses()[hap_index].sequence()[overlap_len..])?;

            // add the coordinates of what was added
            let coordinates = Coordinates::new(contig_name.to_string(),
                (ret.len() - overlap_len) as u64,
                (ret.len() + hap_sequence.len()) as u64
            );
            
            // this mirrors what shows up in debug mode; I think that's fine to not have translations like *68 here for now
            let region_name = format!("{hap_index}_{}", hap_regions[hap_index]);
            regions.push((coordinates, region_name));

            // now extend our region
            ret.push_str(hap_sequence);
        }

        // add a buffer between each one
        ret.push_str(&buffer_sequence);
    }

    // add a buffer at the end also
    ret.push_str(&buffer_sequence);
    Ok(CustomReference {
        contig_name,
        sequence: ret,
        regions
    })
}