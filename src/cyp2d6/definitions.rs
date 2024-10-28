
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use rustc_hash::FxHashMap as HashMap;
use std::collections::{BTreeMap, BTreeSet};
use serde::{Deserialize, Serialize};
use simple_error::{SimpleError, bail};

use crate::cyp2d6::region_label::{Cyp2d6RegionLabel, Cyp2d6RegionType};
use crate::data_types::coordinates::Coordinates;

// these are the fixed buffers around the *5 region that we search for
// in the future we may increase the pre-buffer, but this is basically the *minimum* to accurately find it in WGS
static STAR5_PRE_BUFFER: usize = 500;
static STAR5_POST_BUFFER: usize = 3000;

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Cyp2d6Config {
    /// High-level coordinates of the CYP2D regions
    cyp_coordinates: BTreeMap<String, Coordinates>,
    /// Specific subregion, like exons
    cyp_regions: BTreeMap<String, BTreeMap<String, Coordinates>>,
    /// Coordinates for *5 deleted region
    cyp2d6_star5_del: Coordinates,
    /// Translation from certain alleles to a known star-allele
    cyp_translate: BTreeMap<String, String>,
    /// Inferred connections we allow
    inferred_connections: BTreeSet<(String, String)>,
    /// List of alleles we do not expect to find alone
    unexpected_singletons: BTreeSet<String>
}

impl Cyp2d6Config {
    /// This function should be called after loading a config to verify that everything required to run the algorithms is present.
    pub fn validate_config(&self) -> Result<(), SimpleError> {
        // make sure all expected regions are defined
        let expected_cyp_coordinates = [
            "CYP2D6", "CYP2D7", "REP6", "REP7", "spacer", "link_region", "CYP2D6_wfa_backbone"
        ];
        for &ecc in expected_cyp_coordinates.iter() {
            if !self.cyp_coordinates.contains_key(ecc) {
                bail!("Coordinates for \"{}\" were not found in provided cyp_coordinates.", ecc);
            }
        }

        // make sure all exon regions are defined for our hybrids
        let expected_cyp_region_keys = [
            "CYP2D6", "CYP2D7"
        ];
        for &ecrk in expected_cyp_region_keys.iter() {
            if !self.cyp_regions.contains_key(ecrk) {
                bail!("Data for \"{}\" was not found in provided cyp_regions.", ecrk);
            }
            let cr = self.cyp_regions.get(ecrk).unwrap();
            let expected_cyp_region_exons = 1..10;
            for ecre in expected_cyp_region_exons {
                let exon_label = format!("exon{ecre}");
                if !cr.contains_key(&exon_label) {
                    bail!("Data for \"{}\" is missing coordinates for \"{}\" in cyp_regions.", ecrk, exon_label);
                }
            }
        }

        // cyp2d6_star5_del - just needs a definition
        // cyp_translate - nothing required here to prevent crash
        // inferred_connections - nothing required here to prevent crash
        // unexpected_singletons - nothing required here to prevent crash

        Ok(())
    }

    /// This will return the target extraction region for CYP2D6 and CYP2D7 based on the database coordinates.
    pub fn extraction_region(&self) -> Coordinates {
        let cyp_coordinates = self.cyp_coordinates();
        let cyp2d6_star5_del = self.cyp2d6_star5_del();

        // figure out the BAM mapping coordinates
        let full_d6_region = cyp_coordinates.get("CYP2D6").unwrap();
        let full_d7_region = cyp_coordinates.get("CYP2D7").unwrap();
        let rep6_region = cyp_coordinates.get("REP6").unwrap();
        let rep7_region = cyp_coordinates.get("REP7").unwrap();
        // link_region and spacer are not necessary here
        let bam_region = Coordinates::new(
            full_d6_region.chrom().to_string(),
            [
                full_d6_region.start(), 
                cyp2d6_star5_del.start() - STAR5_PRE_BUFFER as u64,
                full_d7_region.start(),
                rep6_region.start(),
                rep7_region.start()
            ].into_iter().min().unwrap(),
            [
                full_d6_region.end(), 
                cyp2d6_star5_del.end() + STAR5_POST_BUFFER as u64,
                full_d7_region.end(),
                rep6_region.end(),
                rep7_region.end()
            ].into_iter().max().unwrap()
        );
        bam_region
    }

    // getters
    pub fn cyp_coordinates(&self) -> &BTreeMap<String, Coordinates> {
        &self.cyp_coordinates
    }

    pub fn cyp_regions(&self) -> &BTreeMap<String, BTreeMap<String, Coordinates>> {
        &self.cyp_regions
    }

    pub fn cyp2d6_star5_del(&self) -> &Coordinates {
        &self.cyp2d6_star5_del
    }

    pub fn cyp_translate(&self) -> &BTreeMap<String, String> {
        &self.cyp_translate
    }

    pub fn inferred_connections(&self) -> &BTreeSet<(String, String)> {
        &self.inferred_connections
    }

    pub fn unexpected_singletons(&self) -> &BTreeSet<String> {
        &self.unexpected_singletons
    }
}

impl Default for Cyp2d6Config {
    fn default() -> Self {
        let mut cyp_coordinates: BTreeMap<String, Coordinates> = Default::default();
        let preshift = 1;
        let postshift = 0;
        
        // These were our original coordinates that DID NOT include the REP regions
        // CYP2D6 
        // these coordinates were generated by using a +-50bp buffer around the variants in our DB, we assert! this below
        let d6_start = 42126260 - preshift;
        let d6_end = 42132424 - postshift;
        cyp_coordinates.insert("CYP2D6".to_string(), Coordinates::new("chr22".to_string(), d6_start, d6_end));
        // CYP2D7
        // the coordinates were generated by mapping the above and looking at where it lands in D7
        let d7_start = 42139966 - preshift;
        let d7_end = 42145903 - postshift;
        cyp_coordinates.insert("CYP2D7".to_string(), Coordinates::new("chr22".to_string(), d7_start, d7_end));
        
        /*
        // These are new coordinates, D6 is recommended by Xiao (chr22:42123192-42132193) and then we used BLAT to find the equivalent in D7
        // end coordinate adjusted due to a variant in D6
        cyp_coordinates.insert("CYP2D6".to_string(), Coordinates::new("chr22".to_string(), 42123192 - preshift, 42132424 - postshift));
        cyp_coordinates.insert("CYP2D7".to_string(), Coordinates::new("chr22".to_string(), 42135344 - preshift, 42145903 - postshift));
        */

        // regions upstream of D6
        let rep6_start = 42123192 - preshift;
        let rep6_end = 42125963 - postshift;
        cyp_coordinates.insert("REP6".to_string(), Coordinates::new("chr22".to_string(), rep6_start, rep6_end));
        // cyp_coordinates.insert("spacer_CYP2D6".to_string(), Coordinates::new("chr22".to_string(), 42125963 - preshift, 42125965 - postshift)); // spacer does not exist
        // cyp_coordinates.insert("fiveprime_CYP2D6".to_string(), Coordinates::new("chr22".to_string(), 42125965 - preshift, 42126260 - postshift)); // region between spacer and start is very small

        // regions upstream of D7 and after link region; spacer starts where REP7 ends
        let rep7_start = 42135344 - preshift;
        let rep7_end = 42138115 - postshift;
        cyp_coordinates.insert("REP7".to_string(), Coordinates::new("chr22".to_string(), rep7_start, rep7_end));
        let spacer_end = 42139679 - postshift;
        cyp_coordinates.insert("spacer".to_string(), Coordinates::new("chr22".to_string(), rep7_end, spacer_end));
        // cyp_coordinates.insert("fiveprime_CYP2D7_spacer".to_string(), Coordinates::new("chr22".to_string(), 42139679 - preshift, 42139966 - postshift)); // region between spacer and start is very small
        
        // region between D6 and REP7 (typically); starts at the end of D6 and goes to start of REP7
        cyp_coordinates.insert("link_region".to_string(), Coordinates::new("chr22".to_string(), d6_end, rep7_start));

        // these are the coordinates used for WFA realignment
        cyp_coordinates.insert("CYP2D6_wfa_backbone".to_string(), Coordinates::new("chr22".to_string(), d6_start, d6_end));
        
        // now save the exon-level information
        let mut cyp_regions: BTreeMap<String, BTreeMap<String, Coordinates>> = Default::default();
        cyp_regions.insert("CYP2D6".to_string(), 
            {
                // taken from RefSeq, remember these are on reverse strand
                let mut regions: BTreeMap<String, Coordinates> = Default::default();
                regions.insert("exon1".to_string(), Coordinates::new("chr22".to_string(), 42130612 - preshift, 42130810 - postshift));
                regions.insert("exon2".to_string(), Coordinates::new("chr22".to_string(), 42129738 - preshift, 42129909 - postshift));
                regions.insert("exon3".to_string(), Coordinates::new("chr22".to_string(), 42129033 - preshift, 42129185 - postshift));                
                regions.insert("exon4".to_string(), Coordinates::new("chr22".to_string(), 42128784 - preshift, 42128944 - postshift));
                regions.insert("exon5".to_string(), Coordinates::new("chr22".to_string(), 42128174 - preshift, 42128350 - postshift));
                regions.insert("exon6".to_string(), Coordinates::new("chr22".to_string(), 42127842 - preshift, 42127983 - postshift));
                regions.insert("exon7".to_string(), Coordinates::new("chr22".to_string(), 42127447 - preshift, 42127634 - postshift));
                regions.insert("exon8".to_string(), Coordinates::new("chr22".to_string(), 42126851 - preshift, 42126992 - postshift));
                regions.insert("exon9".to_string(), Coordinates::new("chr22".to_string(), 42126499 - preshift, 42126752 - postshift));

                // introns could be derived if we need them
                // spacer region - should be "GGT" which is the dup ACC in the spacer ends; unclear if we need this annotated currently
                // regions.insert("spacer".to_string(), Coordinates::new("chr22".to_string(), 42125963 - preshift, 42125965 - postshift));

                regions
            }
        );
        cyp_regions.insert("CYP2D7".to_string(),
            {
                let mut regions: BTreeMap<String, Coordinates> = Default::default();
                regions.insert("exon1".to_string(), Coordinates::new("chr22".to_string(), 42144284 - preshift, 42144483 - postshift));
                regions.insert("exon2".to_string(), Coordinates::new("chr22".to_string(), 42143410 - preshift, 42143581 - postshift));
                regions.insert("exon3".to_string(), Coordinates::new("chr22".to_string(), 42142728 - preshift, 42142880 - postshift));
                regions.insert("exon4".to_string(), Coordinates::new("chr22".to_string(), 42142479 - preshift, 42142639 - postshift));
                regions.insert("exon5".to_string(), Coordinates::new("chr22".to_string(), 42141868 - preshift, 42142044 - postshift));
                regions.insert("exon6".to_string(), Coordinates::new("chr22".to_string(), 42141534 - preshift, 42141675 - postshift));
                regions.insert("exon7".to_string(), Coordinates::new("chr22".to_string(), 42141152 - preshift, 42141339 - postshift));
                regions.insert("exon8".to_string(), Coordinates::new("chr22".to_string(), 42140555 - preshift, 42140696 - postshift));
                // for some reason, RefSeq has a really big exon 9 for CYP2D7, we changed it to be same size as D6 using UCSC
                // regions.insert("exon9".to_string(), Coordinates::new("chr22".to_string(), 42139576 - preshift, 42140456 - postshift));
                regions.insert("exon9".to_string(), Coordinates::new("chr22".to_string(), 42140203 - preshift, 42140456 - postshift));

                // spacer region - should start and end with "GGT"
                // currently we do not use this
                // regions.insert("spacer".to_string(), Coordinates::new("chr22".to_string(), 42138115 - preshift, 42139679 - postshift));

                regions
            }
        );

        /*
        These are no longer used, but knowing them may be useful in the future
        // started encountering noise in this region, which was already documented by the paraph-rs config
        pub static ref CYP_NOISY_REGIONS: Vec<Coordinates> = {
            let preshift = 1;
            let postshift = 0;
            vec![
                // in paraphase config
                Coordinates::new("chr22".to_string(), 42132023 - preshift, 42132051 - postshift),
                // discovered via testing
                Coordinates::new("chr22".to_string(), 42127650 - preshift, 42127655 - postshift), // poly-C, often extra C
                Coordinates::new("chr22".to_string(), 42128657 - preshift, 42128662 - postshift), // poly-G, often extra G
            ]
        };
        */

        // these are the coordinates of the deleted region, which looks like REP6 and spans over a gap to near the start of D7
        let star5_start = rep6_start;
        let star5_end = 42135343 - postshift;
        let cyp2d6_star5_del: Coordinates = Coordinates::new("chr22".to_string(), star5_start, star5_end);

        // this is a translator for the hybrids into the actual alleles
        let cyp_translate: BTreeMap<String, String> = [
            // all of the D7::D6 can all safely map to *13, but there are sub-alleles if we eventually want to get specific
            ("CYP2D7::CYP2D6::intron1", "13"),
            ("CYP2D7::CYP2D6::exon2", "13"),
            ("CYP2D7::CYP2D6::intron2", "13"),
            ("CYP2D7::CYP2D6::exon3", "13"),
            ("CYP2D7::CYP2D6::intron3", "13"),
            ("CYP2D7::CYP2D6::exon4", "13"),
            ("CYP2D7::CYP2D6::intron4", "13"),
            ("CYP2D7::CYP2D6::exon5", "13"),
            ("CYP2D7::CYP2D6::intron5", "13"),
            ("CYP2D7::CYP2D6::exon6", "13"),
            ("CYP2D7::CYP2D6::intron6", "13"),
            ("CYP2D7::CYP2D6::exon7", "13"),
            ("CYP2D7::CYP2D6::intron7", "13"),
            ("CYP2D7::CYP2D6::exon8", "13"),
            ("CYP2D7::CYP2D6::intron8", "13"),
            ("CYP2D7::CYP2D6::exon9", "13"),
            // the ones we have seen so far are just *68 and *36; the rest are educated guesses
            ("CYP2D6::CYP2D7::intron1", "68"), // intron 1; category A; found an example where we needed this for *68
            ("CYP2D6::CYP2D7::exon2", "68"), // intron 1; category A; this is our standard *68
            ("CYP2D6::CYP2D7::exon8", "61"), // intron 7; category B
            ("CYP2D6::CYP2D7::intron8", "63") // exon 8; category B
            // ("CYP2D6::CYP2D7::exon9", "") // this one becomes either *4.013, *36, or *83; all of which have full entries for D6
        ].iter()
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect();

        // Any population inferred connections go here, this can include hybrid (e.g. *4 + *68) or duplication (e.g. *2 + *2).
        // Allele pairs that are not on this list will get penalized during the chain assessment.
        let inferred_connections: BTreeSet<(String, String)> = [
                // known dups
                ("*1", "*1"),
                ("*2", "*2"),
                ("*3", "*3"),
                ("*4", "*4"),
                ("*6", "*6"),
                ("*9", "*9"),
                ("*10", "*10"),
                ("*17", "*17"),
                ("*28", "*28"),
                ("*29", "*29"),
                ("*35", "*35"),
                ("*41", "*41"),
                ("*43", "*43"),
                ("*45", "*45"),
                ("*146", "*146"),
                // hybrid connections
                ("*4", "*68"),
                ("*10", "*36")
            ].iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect();

        // This is a set of alleles that we expect to always find with something else; these are basically hybrids that don't fly solo
        let unexpected_singletons: BTreeSet<String> = [
                "*36", "*68"
            ].iter()
            .map(|k| k.to_string())
            .collect();

        /*
        Notes on the complex alleles described here: https://a.storyblok.com/f/70677/x/ecb9681e8d/cyp2d6_structural-variation_v3-0.pdf
        - D6-like and D7-like downstream is based on the absence/presence of the spacer block
        - basic D6->D7 partial conversions are listed as normal alleles; I checked 35.002 and 82
        - ABSENT: *5, the deletion allele; breakpoints supposedly in the REP6/7 regions (homologous) and should have the spacer region still
        - duplications - no specific allele, just represented as *4x2 if known, or xN if unknown; 
            to date, only *1, *2, *4, and *41 have been described with more than 2 alleles
        - D7::D6 hybrids - start as D7 but transition to D6
            - grouped under *13 overall (not in database)
            - must include the frame-shift in exon 1; this is basically the *13 overall definition
            - none of the sub-alleles seem to have a database entry either, we would have to encode the Figure 8 description 
            - seemingly all of them have D6-like 3' ends without the spacer region
        - D6::D7 hybrids - start as D6 but transition to D7, also sometimes back to D6 again
            - Category A - transitions back to D6 in downstream
            - Category B - stays as D7 in downstream
            - apparently both A and B are viable, so we might not even be able to use that information at all other than flagging it
            - some of these have entries already
            - ones without entries: 61, 63, 68
            - Figure 10 seems to have a drawing of some of these
        - Table 5 has a list of commonly co-occuring star alleles
        - reference materials (we probably sequenced some of these) are in Table 6
        - appendix has some information on tie-breaking that we may want to pull in at some point
        */

        Self {
            cyp_coordinates,
            cyp_regions,
            cyp2d6_star5_del,
            cyp_translate,
            inferred_connections,
            unexpected_singletons
        }
    }
}

/// This will join together the reference version of CYP2D6 and CYP2D7 at each exon/intron boundary.
/// Sequences are labeled based on the coding orientation, which is on the rev-comp strand relative to the reference genome.
/// This matches the orientation in the PharmVar PDF describing D6::D7 hybrids.
/// This function will also add other "targets" that we are search for, namely: D6 (full), D7 (full), CYP2D6*5 (deletion), and the upstream regions.
/// # Arguments
/// * `reference_genome` - the reference genome we are encoding
/// # Errors
/// * if we can't decode a reference slice into UTF-8 (we have bigger problems if this happens though)
pub fn generate_cyp_hybrids(reference_genome: &ReferenceGenome, cyp2d6_config: &Cyp2d6Config) -> Result<HashMap<Cyp2d6RegionLabel, String>, Box<dyn std::error::Error>> {
    let mut ret: HashMap<Cyp2d6RegionLabel, String> = Default::default();

    let exon_count: usize = 9;
    let gene1 = "CYP2D6";
    let gene2 = "CYP2D7";

    // these are in reverse 
    let cyp_coordinates = cyp2d6_config.cyp_coordinates();
    let chrom = cyp_coordinates.get(gene1).unwrap().chrom();
    let g1_start = cyp_coordinates.get(gene1).unwrap().start() as usize;
    let g1_end = cyp_coordinates.get(gene1).unwrap().end() as usize;
    let g2_start = cyp_coordinates.get(gene2).unwrap().start() as usize;
    let g2_end = cyp_coordinates.get(gene2).unwrap().end() as usize;

    // insert the full sequences first
    ret.insert(
        Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6, None),
        String::from_utf8(reference_genome.get_slice(chrom, g1_start, g1_end).to_vec())?
    );
    ret.insert(
        Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d7, None),
        String::from_utf8(reference_genome.get_slice(chrom, g2_start, g2_end).to_vec())?
    );

    // now lets insert our *5 target sequence which is the breakpoint +- the buffer on each side
    // the buffer values below were basically the minimum before UCSC BLAT correctly split it, there's a lot of homology downstream of the split
    let cyp2d6_star5_del = cyp2d6_config.cyp2d6_star5_del();
    let star5_signature = String::from_utf8(
        reference_genome.get_slice(
            cyp2d6_star5_del.chrom(),
            cyp2d6_star5_del.start() as usize - STAR5_PRE_BUFFER,
            cyp2d6_star5_del.start() as usize
        ).to_vec())? + std::str::from_utf8(
        reference_genome.get_slice(
            cyp2d6_star5_del.chrom(),
            cyp2d6_star5_del.end() as usize,
            cyp2d6_star5_del.end() as usize + STAR5_POST_BUFFER
        ))?;
    // debug!("star5_signature: {star5_signature}");
    ret.insert(
        Cyp2d6RegionLabel::new(Cyp2d6RegionType::Cyp2d6Deletion, None),
        star5_signature
    );
    
    // now loop through and add in all of the hybrids we can by looking at exon/intron boundaries and splicing them together
    let cyp_regions = cyp2d6_config.cyp_regions();
    for exon_index in 1..(exon_count+1) {
        let current_exon = format!("exon{exon_index}");
        let g1_exon = cyp_regions.get(gene1).unwrap().get(&current_exon).unwrap();
        let g2_exon = cyp_regions.get(gene2).unwrap().get(&current_exon).unwrap();

        // first exon will not have a breakpoint at the start, since that's just the full normal sequence
        if exon_index != 1 {
            // these breakpoints are at the start of the exon, which is "end()" here due to rev-comp
            let breakpoint1 = g1_exon.end() as usize;
            let breakpoint2 = g2_exon.end() as usize;
    
            //D6::D7 splices - remember this means the the early exons are D6 which will be the latter half of our 5'->3' sequence
            let d6_d7_label = format!("{gene1}::{gene2}::exon{exon_index}");
            let d6_d7_seq = String::from_utf8(reference_genome.get_slice(chrom, g2_start, breakpoint2).to_vec())?+
                &String::from_utf8(reference_genome.get_slice(chrom, breakpoint1, g1_end).to_vec())?;
            ret.insert(
                Cyp2d6RegionLabel::new(Cyp2d6RegionType::Hybrid, Some(d6_d7_label)),
                d6_d7_seq
            );

            //D7::D6 splices - remember this means the the early exons are D7 which will be the latter half of our 5'->3' sequence
            let d7_d6_label = format!("{gene2}::{gene1}::exon{exon_index}");
            let d7_d6_seq = String::from_utf8(reference_genome.get_slice(chrom, g1_start, breakpoint1).to_vec())?+
                &String::from_utf8(reference_genome.get_slice(chrom, breakpoint2, g2_end).to_vec())?;
            ret.insert(
                Cyp2d6RegionLabel::new(Cyp2d6RegionType::Hybrid, Some(d7_d6_label)),
                d7_d6_seq
            );
        }

        // last exon will not have a breakpoint at the end, since that's just the full normal sequence
        if exon_index != exon_count {
            // these breakpoints are at the end of the exon (start of intron), which is "start()" here due to rev-comp
            let breakpoint1 = g1_exon.start() as usize;
            let breakpoint2 = g2_exon.start() as usize;

            //D6::D7 splices - remember this means the the early exons are D6 which will be the latter half of our 5'->3' sequence
            let d6_d7_label = format!("{gene1}::{gene2}::intron{exon_index}");
            let d6_d7_seq = String::from_utf8(reference_genome.get_slice(chrom, g2_start, breakpoint2).to_vec())?+
                &String::from_utf8(reference_genome.get_slice(chrom, breakpoint1, g1_end).to_vec())?;
            ret.insert(
                Cyp2d6RegionLabel::new(Cyp2d6RegionType::Hybrid, Some(d6_d7_label)),
                d6_d7_seq
            );

            //D7::D6 splices - remember this means the the early exons are D7 which will be the latter half of our 5'->3' sequence
            let d7_d6_label = format!("{gene2}::{gene1}::intron{exon_index}");
            let d7_d6_seq = String::from_utf8(reference_genome.get_slice(chrom, g1_start, breakpoint1).to_vec())?+
                &String::from_utf8(reference_genome.get_slice(chrom, breakpoint2, g2_end).to_vec())?;
            ret.insert(
                Cyp2d6RegionLabel::new(Cyp2d6RegionType::Hybrid, Some(d7_d6_label)),
                d7_d6_seq
            );
        }
    }

    // finally, add in the surrounding regions; which we need for chaining
    let extras = [
        ("REP6", Cyp2d6RegionType::Rep6),
        ("REP7", Cyp2d6RegionType::Rep7),
        ("spacer", Cyp2d6RegionType::Spacer),
        ("link_region", Cyp2d6RegionType::LinkRegion)
    ];
    for (extra_region, region_type) in extras.into_iter() {
        let extra_start = cyp_coordinates.get(extra_region).unwrap().start() as usize;
        let extra_end = cyp_coordinates.get(extra_region).unwrap().end() as usize;
        let sequence = String::from_utf8(reference_genome.get_slice(chrom, extra_start, extra_end).to_vec())?;
        ret.insert(Cyp2d6RegionLabel::new(region_type, None), sequence);
    }

    Ok(ret)
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::path::PathBuf;

    use crate::util::file_io::load_json;

    #[test]
    fn test_config_full_length() {
        // full file
        let test_fn = PathBuf::from("test_data/CYP2D6_configs/full_length.json");
        let config: Cyp2d6Config = load_json(&test_fn).unwrap();
        assert!(config.validate_config().is_ok());
    }

    #[test]
    fn test_config_missing_regions() {
        // this one is missing a CYP2D6 coordinate
        let test_fn = PathBuf::from("test_data/CYP2D6_configs/missing_regions.json");
        let config: Cyp2d6Config = load_json(&test_fn).unwrap();
        assert!(config.validate_config().is_err());
    }

    #[test]
    fn test_config_missing_exons() {
        // this one is missing a CYP2D6 exon
        let test_fn = PathBuf::from("test_data/CYP2D6_configs/missing_exons.json");
        let config: Cyp2d6Config = load_json(&test_fn).unwrap();
        assert!(config.validate_config().is_err());
    }
}