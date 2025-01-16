
use lazy_static::lazy_static;
use log::debug;
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use serde::{Deserialize, Serialize};
use simple_error::{SimpleError, bail};
use std::collections::btree_map::Keys;
use std::collections::{BTreeMap, BTreeSet};

use crate::cli::core::FULL_VERSION;
use crate::data_types::coordinates::Coordinates;
use crate::data_types::gene_definition::{GeneCollection, GeneDefinition};
use crate::data_types::mapping::MappingStats;
use crate::util::mapping::standard_hifi_aligner;

lazy_static! {
    /// Contains a list of all HLA genes we currently support; BTreeSet keeps them ordered for iterating and searching
    pub static ref SUPPORTED_HLA_GENES: BTreeSet<String> = {
        let supported_vec = [
            // the original ones for PGx
            "HLA-A", "HLA-B",
            // additional Class I: C has the most alleles; optional future additions are EFG
            "HLA-C",
            // Class II: DRB1, DQB1, and DQA1 seems to be relative priorities
            "HLA-DPA1", "HLA-DPB1",
            "HLA-DQA1", "HLA-DQB1",
            "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
            // add any new ones here
        ];
        supported_vec.iter()
            .map(|hla| hla.to_string())
            .collect()
    };

    /// Contains the set of copyable coordinates for genes that do *NOT* appear in RefSeq's GRCh38 coordinates
    pub static ref HLA_COORDINATE_COPIES: BTreeMap<String, String> = {
        // (copy_to, copy_from)
        let copy_keys = [
            // DRB1 was a better match than DRB5, but they're both pretty bad
            ("HLA-DRB3", "HLA-DRB1"),
            ("HLA-DRB4", "HLA-DRB1")
        ];
        copy_keys.iter()
            .map(|(to_key, from_key)| (to_key.to_string(), from_key.to_string()))
            .collect()
    };

    /// Genes we can use for coverage normalization estimates
    pub static ref NORMALIZING_HLA_GENES: BTreeSet<String> = {
        let normalizing_hla_genes = [
            // these three are pretty stable; DP and DQ are also if we want to add more in the future
            // "HLA-A", "HLA-B", "HLA-C"
            // normalizing on this was not successful
            "HLA-DRB1" // I think we're going to have to normalize genes based on their pseudogene, or do something more sophisticated
        ];
        normalizing_hla_genes.iter()
            .map(|hla| hla.to_string())
            .collect()
    };

    /// Contains the set of genes that can be absent on a haplotype; hetero-, homo-, hemi- and no-zygous are all possible
    pub static ref ABSENT_HLA_GENES: BTreeSet<String> = {
        let absent_hla_genes = [
            "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
        ];
        absent_hla_genes.iter()
            .map(|hla| hla.to_string())
            .collect()
    };
}

/// Contains the configuration for the HLA database.
/// These are generally coordinates that we do not expect to change except between reference builds.
/// This is mostly a wrapper with some additional logic around a collection of gene definitions
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct HlaConfig {
    /// The core gene collection data
    #[serde(default="HlaConfig::default_gene_collection")]
    gene_collection: GeneCollection
}

impl HlaConfig {
    /// This function should be called after loading a config to verify that everything required to run the algorithms is present.
    pub fn validate_config(&self) -> Result<(), SimpleError> {
        /*
        TODO: with A and B, we just knew it should be 8; but some of the others are different (e.g. 6)
        do we need to validate this or just trust people to not break it with a custom DB?
        
        // make sure all exon regions are defined
        let expected_num_exons = 8;
        */
        for (gene_name, gene_def) in self.gene_collection.gene_dict().iter() {
            let exons = gene_def.exons();
            
            // we used to verify that we had the right number of exons; instead, we now just verify that it's non-empty
            if exons.is_empty() {
                bail!("Found {} exons for \"{gene_name}\", expected >0.", exons.len());
            }
        }

        Ok(())
    }

    /// Given a collection of HLA sequences, this will extend the RefSeq regions to insure that each one is covered.
    /// # Arguments
    /// * `gene_collection` - the collection of gene definitions that can get extended
    /// * `hla_sequences` - the dict of sequences that we need to handle
    /// * `reference_genome` - the pre-loaded reference genome which we will need to align the alleles to
    pub fn new(mut gene_collection: GeneCollection, hla_sequences: &BTreeMap<String, HlaAlleleDefinition>, reference_genome: &ReferenceGenome) -> Result<Self, Box<dyn std::error::Error>> {
        // shared settings for mappings - we only need cigar and md for debugging
        // other settings for mapping
        let output_cigar: bool = false;
        let output_md: bool = false;
        let max_frag_len: Option<usize> = None;
        let extra_flags = None;

        for (gene_name, gene_def) in gene_collection.gene_dict_mut().iter_mut() {
            // check if this gene can be missing
            if ABSENT_HLA_GENES.contains(gene_name) {
                debug!("Marking {gene_name} as absent capable.");
                gene_def.set_absent_capable();
            }

            // pull the coordinates and add an alignment buffer
            let coordinates = gene_def.coordinates();
            debug!("Initial coordinates for {gene_name}: {coordinates}");
            let buffer_size = 2000;
            let alignment_start = coordinates.start() - buffer_size;
            let alignment_end = coordinates.end() + buffer_size;

            // pull the full buffered region and make a mapper for it
            let reference_sequence = reference_genome.get_slice(
                coordinates.chrom(), 
                alignment_start as usize, alignment_end as usize
            );
            let ref_len = reference_sequence.len();
            let dna_aligner = standard_hifi_aligner()
                .with_seq(reference_sequence)?;

            // start with a perfect match as "worst"
            let mut worst_stats = MappingStats::new(ref_len, 0, 0);

            // iterate over each HLA definition with a matching gene name and DNA sequence
            for (hla_id, hla_def) in hla_sequences.iter()
                .filter(|(_, hla_def)| hla_def.gene_name() == gene_name && hla_def.dna_sequence().is_some()) {
                // map the sequence to the region
                let hla_sequence = hla_def.dna_sequence().unwrap();
                let d_mappings = dna_aligner.map(
                    hla_sequence.as_bytes(),
                    output_cigar, output_md, max_frag_len, extra_flags, None
                )?;

                // find the best mapping in the region
                let mut best_stats = MappingStats::new(ref_len, ref_len, 0);
                let mut best_mapping = None;
                for m in d_mappings.iter() {
                    // scoring is based on the lowest edit distance, including unmapped
                    let nm = m.alignment.as_ref().unwrap().nm as usize;

                    // some sanity checks while we debug
                    let unmapped = hla_sequence.len() - (m.query_end - m.query_start) as usize;
                    let stats = MappingStats::new(hla_sequence.len(), nm, unmapped);
                    if stats.mapping_score().score() < best_stats.mapping_score().score() {
                        best_stats = stats;
                        best_mapping = Some(m);
                    }
                }

                if let Some(bm) = best_mapping {
                    // we have the best mapping of the set
                    let corrected_start = alignment_start + bm.target_start as u64;
                    let corrected_end = alignment_start + bm.target_end as u64;

                    let updated = gene_def.extend_coordinates(corrected_start, corrected_end);
                    if updated {
                        debug!("{gene_name} coordinates updated from {hla_id}: {}", gene_def.coordinates());
                    }

                    // higher scores are worse, and we want to find the worst scoring alignment
                    if best_stats.mapping_score().score() > worst_stats.mapping_score().score() {
                        debug!("{gene_name} worst mapping updated to {hla_id} ({:.4}): {bm:?}", best_stats.mapping_score().score());
                        worst_stats = best_stats;
                    }
                }
            }
        }

        Ok(Self {
            gene_collection
        })
    }

    // getters
    pub fn gene_collection(&self) -> &GeneCollection {
        &self.gene_collection
    }

    pub fn gene_names(&self) -> Keys<String, GeneDefinition> {
        self.gene_collection.gene_dict().keys()
    }

    pub fn gene_definition(&self, gene_name: &str) -> Option<&GeneDefinition> {
        self.gene_collection.gene_dict().get(gene_name)
    }

    // Defaults for strand orientations. Having this as a separate functions allows us to do an update without requiring a DB update.
    fn default_strand() -> BTreeMap<String, bool> {
        let mut hla_is_forward_strand: BTreeMap<String, bool> = Default::default();
        hla_is_forward_strand.insert("HLA-A".to_string(), true);
        hla_is_forward_strand.insert("HLA-B".to_string(), false);
        hla_is_forward_strand
    }

    /// Generates a default gene collection, preserving backwards compatibility with older databases.
    /// While this will be formatted differently, it should match the hard-coded values prior to v1.0 of StarPhase.
    fn default_gene_collection() -> GeneCollection {
        let mut hla_coordinates: BTreeMap<String, Coordinates> = Default::default();
        let preshift = 1; // coordinate below are from UCSC and not 0-base shifted
        let postshift = 0;
        // HLA-A
        // chr6:29,942,532-29,945,870 from UCSC browser
        // 03:01:01:01 blats to chr6:29942254-29945755, partially updated below (old end was higher); TODO: systematically solve this?
        hla_coordinates.insert("HLA-A".to_string(), Coordinates::new("chr6".to_string(), 29942254 - preshift, 29945870 - postshift));
        // HLA-B
        // chr6:31,353,875-31,357,179 from UCSC browser
        // 08:01:01:01 blats to chr6:31353362-31357442, updated below; TODO: systematically solve this?
        hla_coordinates.insert("HLA-B".to_string(), Coordinates::new("chr6".to_string(), 31353362 - preshift, 31357442 - postshift));

        let hla_is_forward_strand: BTreeMap<String, bool> = HlaConfig::default_strand();

        let mut hla_transcript_id: BTreeMap<String, String> = Default::default();
        hla_transcript_id.insert("HLA-A".to_string(), "NM_002116.8".to_string());
        hla_transcript_id.insert("HLA-B".to_string(), "NM_005514.8".to_string());
        /*
        // we were originally using this as a filter, but it's no longer necessary; preserving in case we need this in the future
        // the HLA lengths
        static ref HLA_CDNA_LENGTHS: HashMap<String, usize> = {
            let mut hla_lengths: HashMap<String, usize> = Default::default();
            hla_lengths.insert("HLA-A".to_string(), 1098);
            hla_lengths.insert("HLA-B".to_string(), 1089);
            hla_lengths
        };
        */

        // TODO: figure out how to automate this part in the future?
        let mut hla_exons: BTreeMap<String, Vec<Coordinates>> = Default::default();
        hla_exons.insert("HLA-A".to_string(), 
            vec![
                // 0-based, exclusive, sorted; copied from RefSeq file
                Coordinates::new("chr6".to_string(), 29942532 - preshift, 29942626 - postshift),
                Coordinates::new("chr6".to_string(), 29942757 - preshift, 29943026 - postshift),
                Coordinates::new("chr6".to_string(), 29943268 - preshift, 29943543 - postshift),
                Coordinates::new("chr6".to_string(), 29944122 - preshift, 29944397 - postshift),
                Coordinates::new("chr6".to_string(), 29944500 - preshift, 29944616 - postshift),
                Coordinates::new("chr6".to_string(), 29945059 - preshift, 29945091 - postshift),
                Coordinates::new("chr6".to_string(), 29945234 - preshift, 29945281 - postshift),
                Coordinates::new("chr6".to_string(), 29945451 - preshift, 29945870 - postshift)
            ]
        );
        hla_exons.insert("HLA-B".to_string(), 
            vec![
                // 0-based, exclusive, sorted; copied from RefSeq file
                Coordinates::new("chr6".to_string(), 31353875 - preshift, 31354296 - postshift),
                Coordinates::new("chr6".to_string(), 31354479 - preshift, 31354526 - postshift),
                Coordinates::new("chr6".to_string(), 31354633 - preshift, 31354665 - postshift),
                Coordinates::new("chr6".to_string(), 31355107 - preshift, 31355223 - postshift),
                Coordinates::new("chr6".to_string(), 31355317 - preshift, 31355592 - postshift),
                Coordinates::new("chr6".to_string(), 31356167 - preshift, 31356442 - postshift),
                Coordinates::new("chr6".to_string(), 31356688 - preshift, 31356957 - postshift),
                Coordinates::new("chr6".to_string(), 31357086 - preshift, 31357179 - postshift)
            ]
        );

        // it doesn't like to directly use the reference, go figure
        let full_version = FULL_VERSION.to_string();
        let version = format!("starphase_{full_version}_default");
        
        // copy the definitions into the expected format
        let mut gene_dict: BTreeMap<String, GeneDefinition> = Default::default();
        for (gene, coordinates) in hla_coordinates.into_iter() {
            // pull the aux info out
            let is_forward_strand = *hla_is_forward_strand.get(&gene).unwrap();
            let mod_coords = if is_forward_strand {
                hla_exons.get(&gene).unwrap().clone()
            } else {
                hla_exons.get(&gene).unwrap().iter().rev().cloned().collect()
            };

            // now do the definition
            let mut gene_def = GeneDefinition::new(
                gene.clone(), coordinates, is_forward_strand
            );
            let transcript_id = hla_transcript_id.get(&gene).unwrap().clone();
            gene_def.add_transcript_id(transcript_id, None).unwrap();
            for exon in mod_coords.into_iter() {
                gene_def.add_exon(exon).unwrap();
            }
            gene_dict.insert(gene, gene_def);
        }

        // now we can return
        GeneCollection::new(version, gene_dict)
    }
}

impl Default for HlaConfig {
    // This default is based on our initial hard-coded values and is really just for test backwards compatibility at this point
    fn default() -> Self {
        Self { 
            gene_collection: Self::default_gene_collection()
        }
    }
}

/// Wrapper for all the HLA allele informations
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct HlaAlleleDefinition {
    /// The identifier from IMGT
    hla_id: String,
    /// The gene name this goes to; e.g. HLA-A
    gene_name: String,
    /// The assigned star allele as a Vec; e.g. ["01", "01", "01"] is a three field allele
    star_allele: Vec<String>,
    /// The DNA sequence for the record
    dna_sequence: Option<String>,
    /// The cDNA sequence for the record
    cdna_sequence: String
}

impl HlaAlleleDefinition {
    /// Creates a new HlaAlleleDefinition and performs some checks along the way
    /// # Arguments
    /// * `hla_id` - the identifier, expected to be of form "HLA:HLA00001"
    /// * `description` - basically, the star allele, should be of form "A*01:01:01:01"
    /// * `dna_sequence` - the optional DNA sequence, should be ACGT symbols only
    /// * `cdna_sequence` - the cDNA sequence, should be ACGT symbols only
    pub fn new(hla_id: String, description: &str, dna_sequence: Option<String>, cdna_sequence: String) -> Result<HlaAlleleDefinition, Box<dyn std::error::Error>> {
        let star_split: Vec<&str> = description.split('*').collect();
        if star_split.len() != 2 {
            bail!("Star split length != 2 for allele description: {description}");
        }
        let gene_name: String = format!("HLA-{}", star_split[0]);

        let star_allele: Vec<String> = star_split[1].split(':').map(String::from).collect();
        if star_allele.len() > 4 {
            bail!("Unexpected number of fields for allele description: {description}");
        }

        let allowed_symbols = ['A', 'C', 'G', 'T'];
        if let Some(dna) = dna_sequence.as_ref() {
            if !dna.chars().all(|c| allowed_symbols.contains(&c)) {
                bail!("DNA sequence contains non-ACGT symbols.");
            }
        }
        if !cdna_sequence.chars().all(|c| allowed_symbols.contains(&c)) {
            bail!("cDNA sequence contains non-ACGT symbols.");
        }

        Ok(HlaAlleleDefinition {
            hla_id,
            gene_name,
            star_allele,
            dna_sequence,
            cdna_sequence
        })
    }

    pub fn hla_id(&self) -> &str {
        &self.hla_id
    }

    pub fn gene_name(&self) -> &str {
        &self.gene_name
    }

    pub fn star_allele(&self) -> &[String] {
        &self.star_allele
    }

    pub fn dna_sequence(&self) -> Option<&str> {
        self.dna_sequence.as_deref()
    }

    pub fn cdna_sequence(&self) -> &str {
        &self.cdna_sequence
    }
}

#[cfg(test)]
mod tests {
    use rustc_hash::FxHashMap as HashMap;

    use super::*;

    use std::io::Read;
    use std::path::PathBuf;

    use crate::util::file_io::load_json;

    #[test]
    fn test_config_full_length() {
        // full file
        let test_fn = PathBuf::from("test_data/HLA_configs/full_length.json");
        let config: HlaConfig = load_json(&test_fn).unwrap();
        assert!(config.validate_config().is_ok());
    }

    #[test]
    fn test_config_missing_regions() {
        // this one is missing a CYP2D6 coordinate
        let test_fn = PathBuf::from("test_data/HLA_configs/missing_regions.json");
        let config_result: Result<HlaConfig, Box<dyn std::error::Error>> = load_json(&test_fn);
        assert!(config_result.is_err());
    }

    #[test]
    fn test_config_missing_exons() {
        // this one is missing a CYP2D6 exon
        let test_fn = PathBuf::from("test_data/HLA_configs/missing_exons.json");
        let config: HlaConfig = load_json(&test_fn).unwrap();
        assert!(config.validate_config().is_err());
    }
    
    #[test]
    fn test_good_allele_def() {
        let test_name = "test_name".to_string();
        let test_gene = "A";
        let test_star = "01:01:01:01";
        let test_description = format!("{test_gene}*{test_star}");
        let test_dna = Some("ACGT".to_string());
        let test_cdna = "CG".to_string();
        let test_result = HlaAlleleDefinition::new(
            test_name.clone(),
            &test_description,
            test_dna.clone(),
            test_cdna.clone()
        ).unwrap();
        assert_eq!(test_result, HlaAlleleDefinition {
            hla_id: test_name,
            gene_name: "HLA-A".to_string(),
            star_allele: vec!["01".to_string(); 4],
            dna_sequence: test_dna,
            cdna_sequence: test_cdna
        });
    }

    #[test]
    fn test_bad_fields() {
        // not getting tested
        let test_name = "test_name".to_string();
        let test_dna = Some("ACGT".to_string());
        let test_cdna = "CG".to_string();

        // too many fields
        let test_gene = "A";
        let test_star = "01:01:01:01:01";
        let test_description = format!("{test_gene}*{test_star}");
        let test_result = HlaAlleleDefinition::new(
            test_name,
            &test_description,
            test_dna,
            test_cdna
        );
        assert!(test_result.is_err());
    }

    #[test]
    fn test_bad_alleles() {
        let test_name = "test_name".to_string();
        let test_gene = "A";
        let test_star = "01:01:01:01";
        let test_description = format!("{test_gene}*{test_star}");
        
        // bad dna
        let test_bad = "BOB".to_string();
        let test_good = "CG".to_string();
        let test_result = HlaAlleleDefinition::new(
            test_name.clone(),
            &test_description,
            Some(test_bad.clone()),
            test_good.clone()
        );
        assert!(test_result.is_err());

        // bad cnda
        let test_result = HlaAlleleDefinition::new(
            test_name.clone(),
            &test_description,
            Some(test_good),
            test_bad
        );
        assert!(test_result.is_err());
    }

    #[test]
    fn test_hlaconfig_new() {
        // load a faux gene collection with just A and B
        let filename = PathBuf::from("./test_data/refseq_faux/refseq_small.gff.gz");
        let mut targets: BTreeSet<String> = Default::default();
        targets.insert("HLA-A".to_string());
        targets.insert("HLA-B".to_string());
        let gene_collection = GeneCollection::load_refseq_file(&filename, Some(&targets)).unwrap();

        // read the files into memory
        let dna_filename = "./test_data/HLA-faux/hla_gen.fa";
        let mut dna_file = std::fs::File::open(&dna_filename).unwrap();
        let mut dna_result: String = Default::default();
        dna_file.read_to_string(&mut dna_result).unwrap();

        let cdna_filename = "./test_data/HLA-faux/hla_nuc.fa";
        let mut cdna_file = std::fs::File::open(&cdna_filename).unwrap();
        let mut cdna_result: String = Default::default();
        cdna_file.read_to_string(&mut cdna_result).unwrap();
        
        // parse them and collapse them
        let dna_data: HashMap<String, (String, String)> = crate::build_database::convert_fasta_str_to_map(&dna_result, false).unwrap();
        let cdna_data: HashMap<String, (String, String)> = crate::build_database::convert_fasta_str_to_map(&cdna_result, false).unwrap();
        
        // now load the HLA sequences that extend the region outwards
        let hla_sequences = crate::build_database::collapse_hla_lookup(dna_data, cdna_data).unwrap();
        println!("Loaded {} HLA sequences.", hla_sequences.len());

        // load the config
        let fasta_fn = PathBuf::from("./test_data/refseq_faux/hg38_chr6_masked.fa.gz");
        let reference_genome = ReferenceGenome::from_fasta(&fasta_fn).unwrap();
        let config = HlaConfig::new(gene_collection, &hla_sequences, &reference_genome).unwrap();

        // test the results
        let default_config = HlaConfig::default();
        assert_eq!(config.gene_collection().version(), "NCBI RefSeq GCF_000001405.40-RS_2024_08");
        assert_eq!(config.gene_collection().gene_dict(), default_config.gene_collection().gene_dict());

        /*
        // quick way to get the compressible chr6 we want, preserved for future reference
        let buffer = 2000;
        let mut ref_seq = reference_genome.get_full_chromosome("chr6").to_vec();
        let coords = default_config.gene_definition("HLA-A").unwrap().coordinates();
        ref_seq[0..(coords.start() as usize - buffer)].fill(b'N');
        let coords_b = default_config.gene_definition("HLA-B").unwrap().coordinates();
        ref_seq[(coords.end() as usize + buffer)..(coords_b.start() as usize - buffer)].fill(b'N');
        ref_seq[(coords_b.end() as usize + buffer)..].fill(b'N');
        let out_filename = PathBuf::from("./test_data/refseq_faux/hg38_chr6_masked.fa"); // manually gzip afterwards
        let mut out_file = std::fs::File::create(&out_filename).unwrap();
        out_file.write_all(b">chr6\n").unwrap();
        out_file.write_all(&ref_seq).unwrap();
        */
    }
}