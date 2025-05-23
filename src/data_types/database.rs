
use log::{debug, warn};
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use rustc_hash::FxHashMap as HashMap;
use serde::{Deserialize, Serialize};
use simple_error::{SimpleError, bail};
use std::collections::{BTreeMap, BTreeSet};
use std::collections::btree_map::Entry::{Occupied, Vacant};

use crate::cyp2d6::definitions::Cyp2d6Config;
use crate::data_types::alleles::AlleleDefinition;
use crate::data_types::db_const::{CPIC_IGNORED_GENES, CPIC_FULL_DELETIONS, CPIC_PARTIAL_DELETIONS};
use crate::data_types::cpic_api_results::CpicAlleleDefinition;
use crate::data_types::gene_definition::GeneCollection;
use crate::data_types::hgvs::ParsedHgvs;
use crate::data_types::pgx_structural_variants::{PgxStructuralVariants, FullDeletion};
use crate::hla::alleles::{HlaAlleleDefinition, HlaConfig, SUPPORTED_HLA_GENES};

use super::pgx_structural_variants::PartialDeletion;
use super::pharmvar_api_results::PharmvarAlleleDefinition;

#[derive(Clone, Copy, Deserialize, Debug, Default, Eq, PartialEq, Serialize)]
pub enum PgxDataSource {
    #[default]
    Unknown,
    Cpic,
    PharmVar
}

/// This is the full set of PGx information that we have available
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct PgxDatabase {
    /// Metadata for the database
    database_metadata: PgxMetadata,
    /// RefSeq gene entries
    #[serde(default)] // old database versions did not have gene_collection, set to default if missing
    gene_collection: GeneCollection,
    /// This is a map from gene name to all of the relevant information for that gene
    gene_entries: BTreeMap<String, PgxGene>,
    /// The configuration for HLA genes
    #[serde(default)] // will populate with our default() if not in the file
    hla_config: HlaConfig,
    /// The sequences for the HLA alleles
    hla_sequences: BTreeMap<String, HlaAlleleDefinition>,
    /// The configuration for the CYP2D6 gene
    #[serde(default)] // will populate with our default() if not in the file
    cyp2d6_config: Cyp2d6Config,
    /// The sequences for CYP2D6
    cyp2d6_gene_def: BTreeMap<String, AlleleDefinition>
}

impl PgxDatabase {
    /// Creates a new database from the CPIC allele definitions
    /// # Arguments
    /// * `gene_to_chrom` - a simple hashmap from gene name to chromosome
    /// * `cpic_allele_definitions` - the full set of allele definitions from CPIC
    /// * `pharmvar_allele_definitions` - the full set of allele definitions from PharmVar
    /// * `hla_version` - the HLA version from GitHub
    /// * `hla_sequences` - the hashmap from HLA identifier to the definitions (e.g., sequences)
    /// * `pharmvar_version` - the PharmVar version
    /// * `cyp2d6_gene_def` - the CYP2D6 gene definition
    /// * `reference_genome` - the pre-loaded reference genome data
    /// * `opt_refseq_fn` - a path to a local RefSeq file; if not specified, the latest will be downloaded
    /// # Errors
    /// * if CPIC API inconsistencies are detected
    /// * if there is an error while adding a new CPIC allele definition to our database
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        gene_to_chrom: &HashMap<String, String>, cpic_allele_definitions: &[CpicAlleleDefinition],
        pharmvar_allele_definitions: &[PharmvarAlleleDefinition],
        hla_version: String, hla_sequences: BTreeMap<String, HlaAlleleDefinition>,
        pharmvar_version: String, cyp2d6_gene_def: BTreeMap<String, AlleleDefinition>,
        reference_genome: &ReferenceGenome, opt_refseq_fn: Option<&std::path::Path>
    ) -> Result<PgxDatabase, Box<dyn std::error::Error>> {
        // initialize all the gene entries
        let mut gene_entries: BTreeMap<String, PgxGene> = Default::default();
        for (gene_name, chrom) in gene_to_chrom.iter() {
            if CPIC_IGNORED_GENES.contains(gene_name.as_str()) {
                warn!("Gene {gene_name} is on the CPIC ignored genes list, skipping it and all allele definitions for it.");
                continue;
            }
            gene_entries.insert(
                gene_name.clone(),
                PgxGene::new(gene_name, chrom, PgxDataSource::Cpic)
            );
        }
        
        // now add the allele definitions
        for allele_def in cpic_allele_definitions.iter() {
            // make sure are not ignoring this gene
            let gene: &str = &allele_def.gene_symbol;
            if CPIC_IGNORED_GENES.contains(gene) {
                continue;
            }

            // make sure this isn't an SV, we ignore those currently
            if allele_def.is_sv {
                warn!("SV allele detected, ignoring: {gene}, {}", allele_def.allele_name);
                continue;
            }

            // add it now
            let gene_entry: &mut PgxGene = gene_entries.get_mut(gene)
                .ok_or(format!("An allele definition was provided for {gene}, but it was not found in the gene to chromosome list."))?;
            gene_entry.add_cpic_allele(allele_def)?;
        }

        // go through all the default SV events and add them
        debug!("Adding CPIC structural variants...");
        let mut cpic_sv_genes: BTreeSet<String> = Default::default();
        for ((gene, allele), event) in CPIC_FULL_DELETIONS.iter() {
            // make sure are not ignoring this gene
            if CPIC_IGNORED_GENES.contains(gene.as_str()) {
                continue;
            }

            // check if the gene has an entry; this should really only fail in our unit testing
            if let Some(gene_entry) = gene_entries.get_mut(gene) {
                // add it now
                debug!("\tAdding: {gene}{allele} => {event:?}");
                gene_entry.add_full_deletion(allele.clone(), event.clone())?;

                // add any relevent genes to our RefSeq lookup
                for gene in event.full_genes_deleted().iter() {
                    cpic_sv_genes.insert(gene.clone());
                }
            } else {
                warn!("An SV allele definition was provided for {gene}, but it was not found in the gene to chromosome list; skipping.");
            }
        }

        for ((gene, allele), event) in CPIC_PARTIAL_DELETIONS.iter() {
            // make sure are not ignoring this gene
            if CPIC_IGNORED_GENES.contains(gene.as_str()) {
                continue;
            }

            // check if the gene has an entry; this should really only fail in our unit testing
            if let Some(gene_entry) = gene_entries.get_mut(gene) {
                debug!("\tAdding: {gene}{allele} => {event:?}");
                gene_entry.add_partial_deletion(allele.clone(), event.clone())?;

                // add any relevent genes to our RefSeq lookup
                for gene in event.exons_deleted().keys() {
                    cpic_sv_genes.insert(gene.clone());
                }
            } else {
                warn!("An SV allele definition was provided for {gene}, but it was not found in the gene to chromosome list; skipping.");
            }
        }

        // build the metadata for saving to the database output
        let build_time = chrono::Utc::now();
        let cpic_version = format!("API-{}", build_time.to_rfc3339_opts(chrono::SecondsFormat::Nanos, true));
        let database_metadata = PgxMetadata { 
            pbstarphase_version: crate::cli::core::FULL_VERSION.to_string(),
            cpic_version,
            hla_version,
            pharmvar_version,
            build_time
        };

        // add the HLA genes as well
        cpic_sv_genes.extend(SUPPORTED_HLA_GENES.iter().cloned());

        // lastly add in all the PharmVar genes, which need chromosome lookups
        cpic_sv_genes.extend(pharmvar_allele_definitions.iter().map(|a| a.gene_symbol.clone()));

        // generate the default configs here
        // the CPIC config requires coordinates for interpreting any SVs into known forms
        // TODO: do we want to allow users to ultimately specify a file and/or URL? maybe if requested
        let full_gene_collection = if let Some(filename) = opt_refseq_fn {
            // this is really just for unit testing locally
            GeneCollection::load_refseq_file(filename, Some(&cpic_sv_genes))?
        } else {
            // standard path pulls the latest refseq
            GeneCollection::load_refseq_url(None, Some(&cpic_sv_genes))?
        };

        // now we should know which chromosome each gene is on, so we need to add in all the PharmVar allele definitions
        for allele_def in pharmvar_allele_definitions.iter() {
            let gene = &allele_def.gene_symbol;
            if !gene_entries.contains_key(gene) {
                let gene_dict = full_gene_collection.gene_dict().get(gene)
                    .ok_or(format!("{gene} was not found in the RefSeq gene definitions."))?;
                let chrom = gene_dict.coordinates().chrom();
                gene_entries.insert(gene.clone(), PgxGene::new(gene, chrom, PgxDataSource::PharmVar));
            }

            // add the pharmvar allele
            let gene_entry = gene_entries.get_mut(gene).unwrap();
            let chrom = gene_entry.chromosome();
            let reference = reference_genome.get_full_chromosome(chrom);
            gene_entry.add_pharmvar_allele(allele_def, reference)?;
        }

        // all the PharmVar entries are missing a reference allele, so just add one in
        for gene_entry in gene_entries.values_mut() {
            if gene_entry.reference_allele().is_none() {
                gene_entry.add_pharmvar_reference_allele()?;
            }
        }

        // build the HLA and CYP2D6 configs
        let hla_config = HlaConfig::new(&full_gene_collection, &hla_sequences, reference_genome)?;
        let cyp2d6_config = Cyp2d6Config::default();

        // go through the list and remove any genes with no defined alleles
        let filtered_entries: BTreeMap<String, PgxGene> = BTreeMap::from_iter(
            gene_entries.into_iter()
            .filter(|(k, v)| {
                let is_empty: bool = v.defined_haplotypes.is_empty();
                let no_variants = v.variants().is_empty();
                if is_empty {
                    debug!("No defined haplotypes detected for {k}, removing from gene list.");
                } else if no_variants {
                    debug!("No defined variants detected for {k}, removing from gene list.");
                } else {
                    debug!("{k} stats: {} variants, {} haplotypes", v.variants().len(), v.defined_haplotypes.len())
                }
                !is_empty && !no_variants
            })
        );

        // if we made it here, we're all good yo
        Ok(PgxDatabase {
            database_metadata,
            gene_collection: full_gene_collection,
            gene_entries: filtered_entries,
            hla_config,
            hla_sequences,
            cyp2d6_config,
            cyp2d6_gene_def
        })
    }

    /// Validates the loaded database where possible.
    /// This does not prevent data errors, but it will prevent crashes due to missing information.
    /// # Errors
    /// * if the HLA configuration is missing information
    /// * if the CYP2D6 configuration is missing information
    pub fn validate(&self) -> Result<(), SimpleError> {
        self.hla_config.validate_config()?;
        self.cyp2d6_config.validate_config()
    }

    pub fn database_metadata(&self) -> &PgxMetadata {
        &self.database_metadata
    }

    pub fn gene_collection(&self) -> &GeneCollection {
        &self.gene_collection
    }

    pub fn gene_entries(&self) -> &BTreeMap<String, PgxGene> {
        &self.gene_entries
    }

    pub fn hla_config(&self) -> &HlaConfig {
        &self.hla_config
    }

    pub fn hla_sequences(&self) -> &BTreeMap<String, HlaAlleleDefinition> {
        &self.hla_sequences
    }

    pub fn cyp2d6_config(&self) -> &Cyp2d6Config {
        &self.cyp2d6_config
    }

    pub fn cyp2d6_gene_def(&self) -> &BTreeMap<String, AlleleDefinition> {
        &self.cyp2d6_gene_def
    }
}

/// Contains metadata about the construction of the database
#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct PgxMetadata {
    /// The version of pbstarphase we're running
    pbstarphase_version: String,
    /// The version of the CPIC database
    cpic_version: String,
    /// The version of the HLA database
    hla_version: String,
    /// The version of the PharmVar database
    pharmvar_version: String,
    /// The time the database was constructed
    build_time: chrono::DateTime<chrono::Utc>
}

/// A PGx gene has defined variants as well as alleles that are composites of the defined variants.
/// There is also always one "reference" allele
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct PgxGene {
    /// The name of the gene
    gene_name: String,
    /// Source of the definitions
    #[serde(default)] // enabled Unknown for backwards compatibility
    data_source: PgxDataSource,
    /// The chromosome the gene is located on
    chromosome: String,
    /// Variants by their provided ID
    variants: BTreeMap<u64, PgxVariant>,
    /// Any defined structural variants associated with this gene
    #[serde(default)] // will populate with our default() if not in the file; which is None
    structural_variants: Option<PgxStructuralVariants>,
    /// Each name points to an allele that contains variants in the same order as "ordered_variants"
    defined_haplotypes: BTreeMap<String, PgxHaplotype>,
    /// The allele corresponding to reference, it is the only one that will typically be "full" on the PgxHaplotype definition
    reference_allele: Option<String>,
}

impl PgxGene {
    /// Creates a new, blank PgxGene with minimum info. Alleles need to be added afterwards.
    /// # Arguments
    /// * `gene_name` - the name of the gene
    /// * `chromosome` - the chromosome this gene and all of it's variants are on
    /// * `data_source` - origin of this definition
    fn new(gene_name: &str, chromosome: &str, data_source: PgxDataSource) -> PgxGene {
        PgxGene { 
            gene_name: gene_name.to_string(), 
            chromosome: chromosome.to_string(), 
            variants: Default::default(),
            structural_variants: None,
            defined_haplotypes: Default::default(),
            reference_allele: None,
            data_source
        }
    }

    /// This will add a new allele definition for this gene and perform any sanity checks along the way.
    /// # Arguments
    /// * `allele_definition` - the CPIC allele definition
    /// # Errors
    /// * if this allele was already provided
    /// * if this is flagged as the reference allele, but we already have a reference allele
    /// * if this is flagged as containing an SV, we are not handling those currently
    /// * if there is a variant definition conflict
    /// # Panics
    /// * if the gene names do not match
    fn add_cpic_allele(&mut self, allele_definition: &CpicAlleleDefinition) -> Result<(), Box<dyn std::error::Error>> {
        // make sure gene symbol matches
        assert_eq!(self.gene_name, allele_definition.gene_symbol);

        // make sure this allele is not already defined
        if self.defined_haplotypes.contains_key(&allele_definition.allele_name) {
            bail!("Duplicate allele definition found for {}: {}", self.gene_name, allele_definition.allele_name);
        }

        // make sure we did not get a second reference allele, we can't handle that pressure
        if self.reference_allele.is_some() && allele_definition.is_reference {
            bail!("Multiple reference alleles provided for {}: {} and {}", self.gene_name, self.reference_allele.as_ref().unwrap(), allele_definition.allele_name);
        }

        // make sure this is not an SV
        if allele_definition.is_sv {
            bail!("SV allele detected for {}, these are not handled: {}", self.gene_name, allele_definition.allele_name);
        }

        if self.data_source != PgxDataSource::Cpic {
            bail!("Adding a CPIC allele to a gene without a CPIC data source");
        }

        // check if we should mark this as reference allele
        let is_reference = allele_definition.is_reference;
        if is_reference {
            self.reference_allele = Some(allele_definition.allele_name.clone());
        }

        // initialize the defined alleles as "None" for all known variants (this may change in length and assignment)
        let mut haplotype: BTreeMap<u64, String> = Default::default();

        // okay, now we can check the variants and add them as necessary
        for ad_variant in allele_definition.variants.iter() {
            let variant_name: String = ad_variant.sequence_location.name.clone();
            let dbsnp_id: Option<String> = ad_variant.sequence_location.dbsnp_id.clone();
            let variant_id: u64 = ad_variant.sequence_location.id;
            let position: usize = ad_variant.sequence_location.position;
            let variant_sequence: &str = &ad_variant.variant_allele;

            match self.variants.entry(variant_id) {
                Occupied(mut entry) => {
                    // we have this variant loaded, do sanity checks
                    let variant: &mut PgxVariant = entry.get_mut();
                    if variant.position != position {
                        bail!("Encountered variant with id {} but different positions: {} != {}", variant_id, variant.position, position);
                    }
                    if variant.dbsnp_id != dbsnp_id {
                        bail!("Encountered variants with id {} but different dbsnp IDs: {:?} != {:?}", variant_id, dbsnp_id, variant.dbsnp_id);
                    }
                    if is_reference {
                        match &variant.alleles[0] {
                            Some(ra) => if ra != variant_sequence { 
                                bail!("Encountered variant with id {} but different reference alleles: {} != {}", variant_id, ra, variant_sequence); 
                            },
                            None => variant.alleles[0] = Some(variant_sequence.to_string())
                        };
                    } else {
                        // check if any alleles already match the sequence
                        let match_index: Option<usize> = variant.alleles.iter()
                            .position(|a| a.as_ref().unwrap_or(&"".to_string()) == variant_sequence);
                        match match_index {
                            Some(i) => {
                                // we already have this sequence in our list
                                // this should never be the REF allele
                                assert_ne!(i, 0);
                            },
                            None => {
                                // this is a new sequence to our list, append it
                                variant.alleles.push(Some(variant_sequence.to_string()));
                            }
                        };
                    }
                },
                Vacant(entry) => {
                    // this is a new variant for us
                    let var_alleles: Vec<Option<String>> = if is_reference {
                        vec![Some(variant_sequence.to_string())]
                    } else {
                        vec![None, Some(variant_sequence.to_string())]
                    };
                    let new_variant = PgxVariant {
                        name: variant_name,
                        dbsnp_id,
                        position,
                        alleles: var_alleles
                    };

                    // store it and update our alleles
                    entry.insert(new_variant);
                }
            };

            // now we need to save the allele for this particular haplotype
            match haplotype.entry(variant_id) {
                Vacant(entry) => {
                    entry.insert(variant_sequence.to_string());
                },
                Occupied(_entry) => {
                    bail!("Detected CPIC allele with same variant assigned multiple times: {} {}", allele_definition.allele_name, variant_id);
                }
            }
        }

        // finally, add the allele set we parsed
        self.defined_haplotypes.insert(
            allele_definition.allele_name.clone(),
            PgxHaplotype { haplotype }
        );

        Ok(())
    }

    /// This will add a new PharmVar allele definition for this gene and perform any sanity checks along the way.
    /// # Arguments
    /// * `allele_definition` - the PharmVar allele definition
    /// * `reference` - the reference sequence (e.g. chromosome) that this allele resides on
    /// # Errors
    /// * if this allele was already provided
    /// * if this is flagged as the reference allele, but we already have a reference allele
    /// * if this is flagged as containing an SV, we are not handling those currently
    /// * if there is a variant definition conflict
    /// # Panics
    /// * if the gene names do not match
    pub fn add_pharmvar_allele(&mut self, allele_definition: &PharmvarAlleleDefinition, reference: &[u8]) -> Result<(), Box<dyn std::error::Error>> {
        // make sure gene symbol matches
        assert_eq!(self.gene_name, allele_definition.gene_symbol);

        // make sure this allele is not already defined
        if self.defined_haplotypes.contains_key(allele_definition.star_allele()) {
            bail!("Duplicate allele definition found for {}: {} => {}", self.gene_name, allele_definition.allele_name, allele_definition.star_allele());
        }

        // PharmVar does not directly provide reference alleles, so make sure we have variants
        if allele_definition.variants.is_empty() {
            bail!("Non-reference allele definition provided with no variants.");
        }

        if self.data_source != PgxDataSource::PharmVar {
            bail!("Adding a PharmVar allele to a gene without a PharmVar data source");
        }

        // initialize the defined alleles as "None" for all known variants (this may change in length and assignment)
        let mut haplotype: BTreeMap<u64, String> = Default::default();

        // okay, now we can check the variants and add them as necessary
        for ad_variant in allele_definition.variants.iter() {
            let variant_name = ad_variant.variant_id.clone();
            let dbsnp_id = ad_variant.rs_id.clone();
            let variant_id: u64 = ad_variant.variant_id.parse()?;

            // originally we were using ad_variant.hgvs, but it has some very annoying things in it
            let parsed_hgvs = ParsedHgvs::new(&ad_variant.position)?;
            let (position, ref_seq, alt_seq) = parsed_hgvs.generate_ref_alt(reference)?;

            debug!("{} => {parsed_hgvs:?} => ({position}, {ref_seq:?}, {alt_seq:?})", ad_variant.position);

            match self.variants.entry(variant_id) {
                Occupied(mut entry) => {
                    // we have this variant loaded, do sanity checks
                    let variant: &mut PgxVariant = entry.get_mut();
                    if variant.position != position {
                        bail!("Encountered variant with id {} but different positions: {} != {}", variant_id, variant.position, position);
                    }
                    if variant.dbsnp_id != dbsnp_id {
                        bail!("Encountered variants with id {} but different dbsnp IDs: {:?} != {:?}", variant_id, dbsnp_id, variant.dbsnp_id);
                    }

                    // check if any alleles already match the sequence
                    let match_index: Option<usize> = variant.alleles.iter()
                        .position(|a| a.as_ref().unwrap_or(&"".to_string()) == &alt_seq);
                    match match_index {
                        Some(i) => {
                            // we already have this sequence in our list
                            // this should never be the REF allele
                            assert_ne!(i, 0);
                        },
                        None => {
                            // this is a new sequence to our list, append it
                            // variant.alleles.push(Some(alt_seq.clone()));

                            // Note for future Matt: if this pops up, we can probably just append like above
                            //       but we don't have a test for this since we have no examples
                            bail!("Unhandled conflicting alleles in PharmVar");
                        }
                    };
                },
                Vacant(entry) => {
                    // this is a new variant for us
                    let var_alleles = vec![
                        Some(ref_seq), Some(alt_seq.clone())
                    ];
                    let new_variant = PgxVariant {
                        name: variant_name,
                        dbsnp_id,
                        position,
                        alleles: var_alleles
                    };

                    // store it and update our alleles
                    entry.insert(new_variant);
                }
            };

            // now we need to save the allele for this particular haplotype
            match haplotype.entry(variant_id) {
                Vacant(entry) => {
                    entry.insert(alt_seq);
                },
                Occupied(_entry) => {
                    bail!("Detected CPIC allele with same variant assigned multiple times: {} {}", allele_definition.allele_name, variant_id);
                }
            };
        }

        // finally, add the allele set we parsed
        self.defined_haplotypes.insert(
            allele_definition.star_allele().to_string(),
            PgxHaplotype { haplotype }
        );

        Ok(())
    }

    /// Creates a PharmVar reference allele as *1
    pub fn add_pharmvar_reference_allele(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let allele_name = "*1";

        // make sure this allele is not already defined
        if self.defined_haplotypes.contains_key(allele_name) {
            bail!("Duplicate allele definition found for {}: {}", self.gene_name, allele_name);
        }

        if self.reference_allele.is_some() {
            bail!("Reference allele has already been set.");
        }

        // finally, add the allele and set it as the reference
        self.defined_haplotypes.insert(
            allele_name.to_string(),
            PgxHaplotype::default() // empty set of variants
        );
        self.reference_allele = Some(allele_name.to_string());

        Ok(())
    }

    /// Adds a full deletion event to the alleles for this gene
    /// * `label` - the label (i.e. genotype) for this event
    /// * `event` - the full deletion to add
    pub fn add_full_deletion(&mut self, label: String, event: FullDeletion) -> Result<(), Box<dyn std::error::Error>> {
        if self.structural_variants.is_none() {
            self.structural_variants = Some(Default::default());
        }

        let sv_db = self.structural_variants.as_mut().unwrap();
        sv_db.add_full_deletion(label, event)?;
        Ok(())
    }

    /// Adds a partial deletion event to the alleles for this gene
    /// * `label` - the label (i.e. genotype) for this event
    /// * `event` - the full deletion to add
    pub fn add_partial_deletion(&mut self, label: String, event: PartialDeletion) -> Result<(), Box<dyn std::error::Error>> {
        if self.structural_variants.is_none() {
            self.structural_variants = Some(Default::default());
        }

        let sv_db = self.structural_variants.as_mut().unwrap();
        sv_db.add_partial_deletion(label, event)?;
        Ok(())
    }

    pub fn gene_name(&self) -> &str {
        &self.gene_name
    }

    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }

    pub fn reference_allele(&self) -> Option<&str> {
        self.reference_allele.as_deref()
    }

    pub fn variants(&self) -> &BTreeMap<u64, PgxVariant> {
        &self.variants
    }

    pub fn structural_variants(&self) -> Option<&PgxStructuralVariants> {
        self.structural_variants.as_ref()
    }

    pub fn defined_haplotypes(&self) -> &BTreeMap<String, PgxHaplotype> {
        &self.defined_haplotypes
    }
}

/// This is the core information needed to identify a variant in PGx land
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct PgxVariant {
    /// the name of the variant
    name: String,
    /// DBSNP ID if available
    dbsnp_id: Option<String>,
    /// The 1-based coordinate of the variant
    position: usize,
    /// All of the alleles, index 0 is *always* reference allele
    alleles: Vec<Option<String>>
}

impl PgxVariant {
    pub fn new(name: String, dbsnp_id: Option<String>, position: usize, alleles: Vec<Option<String>>) -> PgxVariant {
        PgxVariant {
            name, dbsnp_id, position, alleles
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn dbsnp_id(&self) -> &Option<String> {
        &self.dbsnp_id
    }

    pub fn position(&self) -> usize {
        self.position
    }

    pub fn alleles(&self) -> &[Option<String>] {
        &self.alleles
    }
}

/// For the database representation, it makes more sense to have a sparser format that is just ID -> allele value
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct PgxHaplotype {
    /// The PGx variant ID to the allele for this haplotype, only defined ones are stored
    haplotype: BTreeMap<u64, String>
}

impl PgxHaplotype {
    pub fn haplotype(&self) -> &BTreeMap<u64, String> {
        &self.haplotype
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::path::PathBuf;

    use crate::data_types::db_const::CYP2D6;
    use crate::util::file_io::load_json;

    fn create_masked_reference() -> ReferenceGenome {
        let mut ref_gen = ReferenceGenome::empty_reference();
        ref_gen.add_contig("chr6".to_string(), &"N".repeat(200000000)).unwrap();
        ref_gen.add_contig("chr8".to_string(), &"N".repeat(200000000)).unwrap();
        ref_gen
    }

    #[test]
    fn test_simple_cacna1s() {
        // consts we can tweak
        let gene_name = "CACNA1S";
        let chrom = "chr1";
        let cacna1s_fn = PathBuf::from("test_data/CACNA1S/CPIC_API.json");
        let nat2_fn = PathBuf::from("test_data/NAT2/PharmVar_API.json");
        
        // set up a simple single gene hashmap
        let mut gene_to_chrom = HashMap::default();
        gene_to_chrom.insert(gene_name.to_string(), chrom.to_string());

        // load allele definitions
        let cacna1s_allele_defs: Vec<CpicAlleleDefinition> = load_json(&cacna1s_fn).unwrap();
        let pharmvar_allele_defs: Vec<PharmvarAlleleDefinition> = load_json(&nat2_fn).unwrap();

        let mut simple_hla: BTreeMap<String, HlaAlleleDefinition> = Default::default();
        let allele_name: String = "HLA00001".to_string();
        simple_hla.insert(
            allele_name.clone(), 
            HlaAlleleDefinition::new(allele_name.clone(), "A*01:01:01:01", Some("ACGT".to_string()), "TGCA".to_string()).unwrap()
        );

        let mut simple_cyp: BTreeMap<String, AlleleDefinition> = Default::default();
        let allele_name: String = "PV00124".to_string();
        simple_cyp.insert(
            allele_name.clone(), 
            AlleleDefinition::new(Some(allele_name.clone()), "CYP2D6*1", vec![]).unwrap()
        );

        // pre-load our small reference genome
        let reference_genome = create_masked_reference();
        
        // build the database
        let hla_version: String = "hla_v1".to_string();
        let pharmvar_version: String = "pharmvar_v1".to_string();
        let refseq_fn = std::path::PathBuf::from("./test_data/refseq_faux/refseq_small.gff.gz");
        let pgx_database = PgxDatabase::new(
            &gene_to_chrom,
            &cacna1s_allele_defs,
            &pharmvar_allele_defs,
            hla_version.clone(),
            simple_hla.clone(),
            pharmvar_version.clone(),
            simple_cyp.clone(),
            &reference_genome,
            Some(&refseq_fn)
        ).unwrap();

        // check that one gene inside
        assert_eq!(pgx_database.gene_entries.len(), 2);

        // check the cacna1s data
        let cacna1s_entry = pgx_database.gene_entries.get(gene_name).unwrap();
        assert_eq!(cacna1s_entry.gene_name, gene_name);
        assert_eq!(cacna1s_entry.chromosome, chrom);
        assert_eq!(cacna1s_entry.reference_allele.as_ref().unwrap(), "Reference");
        
        // check one of the variants for cacna1s
        assert_eq!(cacna1s_entry.variants.len(), 2);
        let variant = cacna1s_entry.variants.get(&777260).unwrap();
        assert_eq!(variant.dbsnp_id.as_ref().unwrap(), "rs772226819");
        assert_eq!(variant.position, 201091993);
        assert_eq!(variant.alleles, vec![Some("G".to_string()), Some("A".to_string())]);

        // check the alleles as well
        assert_eq!(cacna1s_entry.defined_haplotypes.len(), 3);

        let reference = cacna1s_entry.defined_haplotypes.get("Reference").unwrap();
        assert_eq!(reference.haplotype.len(), 2);
        assert_eq!(reference.haplotype.get(&777260).unwrap(), "G");
        assert_eq!(reference.haplotype.get(&777261).unwrap(), "C");

        let alt1 = cacna1s_entry.defined_haplotypes.get("c.520C>T").unwrap();
        assert_eq!(alt1.haplotype.len(), 1);
        assert_eq!(alt1.haplotype.get(&777260).unwrap(), "A");

        let alt2 = cacna1s_entry.defined_haplotypes.get("c.3257G>A").unwrap();
        assert_eq!(alt2.haplotype.len(), 1);
        assert_eq!(alt2.haplotype.get(&777261).unwrap(), "T");

        // check NAT2 stuff from our partial PharmVar
        let nat2_entry = pgx_database.gene_entries().get("NAT2").unwrap();
        assert_eq!(nat2_entry.gene_name.as_str(), "NAT2");
        assert_eq!(nat2_entry.chromosome, "chr8");
        assert_eq!(nat2_entry.reference_allele.as_ref().unwrap(), "*1");

        // check one of the variants for cacna1s
        assert_eq!(nat2_entry.variants.len(), 3);
        let variant = nat2_entry.variants.get(&2704).unwrap();
        assert_eq!(variant.dbsnp_id.as_ref().unwrap(), "rs1208");
        assert_eq!(variant.position, 18400806);
        assert_eq!(variant.alleles, vec![Some("N".to_string()), Some("A".to_string())]);

        // check that we auto-created a reference *1 allele
        assert_eq!(nat2_entry.defined_haplotypes.len(), 2);
        let reference = nat2_entry.defined_haplotypes.get("*1").unwrap();
        assert!(reference.haplotype.is_empty());

        // and check the allele from our mock JSON
        let alt1 = nat2_entry.defined_haplotypes.get("*36").unwrap();
        assert_eq!(alt1.haplotype.len(), 3);
        assert_eq!(alt1.haplotype.get(&2704).unwrap(), "A");

        // check the HLA stuff
        assert_eq!(pgx_database.database_metadata().hla_version, hla_version);
        assert_eq!(pgx_database.hla_sequences(), &simple_hla);

        // check the CYP2D6 stuff
        assert_eq!(pgx_database.database_metadata().pharmvar_version, pharmvar_version);
        let base_entry = pgx_database.cyp2d6_gene_def().get("PV00124").unwrap();
        assert_eq!(base_entry.gene_name(), CYP2D6);
        assert_eq!(base_entry.star_allele(), "1");
        assert_eq!(base_entry.variants(), &[]);
    }

    #[test]
    #[should_panic]
    fn test_error_sv() {
        // test that SVs cause a panic
        let gene_name = "CACNA1S";
        let chrom = "chr1";
        let mut pgx_gene = PgxGene::new(gene_name, chrom, PgxDataSource::Cpic);

        // load allele definitions
        let cacna1s_fn = PathBuf::from("test_data/CACNA1S/CPIC_API.json");
        let mut cacna1s_allele_defs: Vec<CpicAlleleDefinition> = load_json(&cacna1s_fn).unwrap();
        
        // mark as an SV and add it
        cacna1s_allele_defs[0].is_sv = true;
        pgx_gene.add_cpic_allele(&cacna1s_allele_defs[0]).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_error_duplicate() {
        // test that SVs cause a panic
        let gene_name = "CACNA1S";
        let chrom = "chr1";
        let mut pgx_gene = PgxGene::new(gene_name, chrom, PgxDataSource::Cpic);

        // load allele definitions
        let cacna1s_fn = PathBuf::from("test_data/CACNA1S/CPIC_API.json");
        let cacna1s_allele_defs: Vec<CpicAlleleDefinition> = load_json(&cacna1s_fn).unwrap();
        
        // add the same ID twice
        pgx_gene.add_cpic_allele(&cacna1s_allele_defs[0]).unwrap();
        pgx_gene.add_cpic_allele(&cacna1s_allele_defs[0]).unwrap();
    }

    #[test]
    #[should_panic]
    fn test_error_double_reference() {
        // test that SVs cause a panic
        let gene_name = "CACNA1S";
        let chrom = "chr1";
        let mut pgx_gene = PgxGene::new(gene_name, chrom, PgxDataSource::Cpic);

        // load allele definitions
        let cacna1s_fn = PathBuf::from("test_data/CACNA1S/CPIC_API.json");
        let mut cacna1s_allele_defs: Vec<CpicAlleleDefinition> = load_json(&cacna1s_fn).unwrap();
        
        // add two different alleles both flagged as reference
        cacna1s_allele_defs[0].is_reference = true;
        cacna1s_allele_defs[1].is_reference = true;
        pgx_gene.add_cpic_allele(&cacna1s_allele_defs[0]).unwrap();
        pgx_gene.add_cpic_allele(&cacna1s_allele_defs[1]).unwrap();
    }
}