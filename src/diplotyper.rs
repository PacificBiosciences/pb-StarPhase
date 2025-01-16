
use log::{debug, error, info, trace, warn};
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use simple_error::bail;
use std::collections::hash_map::Entry::{Occupied, Vacant};
use std::path::{Path, PathBuf};

use crate::cli::diplotype::DiplotypeSettings;
use crate::cyp2d6::caller::diplotype_cyp2d6;
use crate::data_types::database::{PgxDatabase, PgxGene, PgxVariant};
use crate::data_types::normalized_variant::{Genotype, NormalizedGenotype, NormalizedVariant, NormalizedPgxHaplotype};
use crate::data_types::pgx_diplotypes::{PgxDiplotypes, PgxGeneDetails, PgxVariantDetails, Diplotype};
use crate::hla::caller::{diplotype_hla_batch, diplotype_hla};
use crate::util::file_io::load_file_lines;
use crate::visualization::debug_bam_writer::DebugBamWriter;
use crate::visualization::igv_session_writer::IgvSessionWriter;

/// This is the main function to call all of the diplotypes.
/// It handles all the VCF parsing and variant normalization.
/// # Arguments
/// * `database` - the pre-loaded database of PGx data
/// * `opt_vcf_fn` - the optional VCF that we will scan for variants from the database, we assume it is indexed
/// * `reference_genome` - the pre-loaded reference genome; if None, some normalization steps will not happen
/// * `bam_fns` - any BAM files that will be used for HLA calling
/// * `cli_settings` - the full settings for diplotyping
/// # Errors
/// * if there are errors normalizing database entries
/// * if there are errors loading the VCF file
/// * if there are errors calling the diplotypes
pub fn call_diplotypes(
    database: &PgxDatabase, opt_vcf_fn: Option<&Path>, reference_genome: Option<&ReferenceGenome>,
    bam_fns: &[PathBuf], cli_settings: &DiplotypeSettings
) -> Result<PgxDiplotypes, Box<dyn std::error::Error>> {
    let mut diplotypes: PgxDiplotypes = PgxDiplotypes::new(database.database_metadata().clone());

    // figure out the set of genes to include / exclude
    let opt_exclude_set: Option<HashSet<String>> = if let Some(efn) = cli_settings.exclude_fn.as_deref() {
        Some(load_file_lines(efn)?)
    } else {
        None
    };
    let opt_include_set: Option<HashSet<String>> = if let Some(ifn) = cli_settings.include_fn.as_deref() {
        Some(load_file_lines(ifn)?)
    } else {
        None
    };

    if let Some(vcf_fn) = opt_vcf_fn {
        for (gene_name, gene_entry) in database.gene_entries().iter() {
            // assume include set has everything UNLESS it is specified
            let included = opt_include_set.as_ref().map(|include_set| include_set.contains(gene_name)).unwrap_or(true);
            if !included {
                debug!("Skipping {gene_name}, not in include set");
                continue;
            }
            
            // assume exclude set has nothing UNLESS it is specified
            let excluded = opt_exclude_set.as_ref().map(|exclude_set| exclude_set.contains(gene_name)).unwrap_or(false);
            if excluded {
                debug!("Skipping {gene_name}, part of exclude set");
                continue;
            }

            info!("Solving {gene_name}...");

            // we need to normalize all of the variants in our database and load haplotypes as well
            let (variant_hash, normalized_haplotypes) = load_database_haplotypes(gene_entry, reference_genome)?;
            debug!("Loaded {} normalized variants.", variant_hash.len());
            debug!("Loaded {} normalized haplotypes.", normalized_haplotypes.len());

            // make sure we have variants, otherwise we can't do anything
            if variant_hash.is_empty() {
                warn!("No variants found for {gene_name}, returning default reference allele.");
                let reference_name: &str = gene_entry.reference_allele().unwrap_or("NO_REFERENCE_ALLELE");
                let all_ref_diplotype = Diplotype::new(reference_name, reference_name);

                debug!("\t{:?}", all_ref_diplotype.diplotype());
                let gene_details = PgxGeneDetails::new(
                    vec![all_ref_diplotype],
                    None,
                    vec![]
                )?;
                diplotypes.insert(gene_name.clone(), gene_details)?;
                continue;
            }

            // load the normalize variants from the VCF
            let vcf_variants: HashMap<NormalizedVariant, NormalizedGenotype> = load_vcf_variants(vcf_fn, &variant_hash, reference_genome)?;
            debug!("Loaded {} normalized genotypes.", vcf_variants.len());
            
            // finally, call the diplotype from the loaded variants
            let diplotype = solve_diplotype(&normalized_haplotypes, &vcf_variants)?;
            debug!("Diplotype for {gene_name} => {:?}", diplotype.iter().map(|d| d.diplotype()).collect::<Vec<&str>>());
            let variant_details: Vec<PgxVariantDetails> = vcf_variants.into_iter()
                .map(|(nv, ng)| {
                    let variant_meta = variant_hash.get(&nv).unwrap();
                    PgxVariantDetails::new(
                        variant_meta.variant_id,
                        variant_meta.name.clone(),
                        variant_meta.dbsnp_id.clone(),
                        nv,
                        ng
                    )
                })
                .collect();
            let gene_details = PgxGeneDetails::new(
                diplotype,
                None,
                variant_details
            )?;
            diplotypes.insert(gene_name.clone(), gene_details)?;
        }
    } else {
        info!("No VCF file provided, all variant based diplotyping was skipped.");
    }

    if !bam_fns.is_empty() {
        // this check is not _really_ necessary, but we have an Option so lets keep it safe
        if reference_genome.is_none() {
            bail!("Reference genome is required for reading alignment files");
        }

        let mut debug_bam_writer: Option<DebugBamWriter> = match cli_settings.debug_folder.as_ref() {
            Some(debug_folder) => {
                let extension = "debug_consensus.bam";
                let consensus_fn = debug_folder.join(extension);
                Some(DebugBamWriter::new(consensus_fn, reference_genome.unwrap())?)
            },
            None => None
        };

        let mut debug_custom_writer: Option<IgvSessionWriter> = cli_settings.debug_folder.as_ref().map(|debug_folder| {
            let extension = "hla_igv_custom";
            let session_folder = debug_folder.join(extension);
            IgvSessionWriter::new(session_folder, true)
        });

        if !cli_settings.debug_skip_hla {
            // the HLA list now comes from the database config, no longer hard coded
            let initial_hla_list: Vec<String> = database.hla_config().gene_names()
                .cloned().collect();

            // filter out anything that is excluded due to CLI parameters
            let final_hla_list: Vec<String> = initial_hla_list.into_iter()
                .filter(|gene_name| {
                    // assume include set has everything UNLESS it is specified
                    let included = opt_include_set.as_ref().map(|include_set| include_set.contains(gene_name)).unwrap_or(true);
                    if !included {
                        debug!("Skipping {gene_name}, not in include set");
                    }
                    
                    // assume exclude set has nothing UNLESS it is specified
                    let excluded = opt_exclude_set.as_ref().map(|exclude_set| exclude_set.contains(gene_name)).unwrap_or(false);
                    if excluded {
                        debug!("Skipping {gene_name}, part of exclude set");
                    }

                    // this one is only kept if it is both included AND NOT excluded
                    included && !excluded
                })
                .collect();

            // only call this if something is in the list to diplotype
            if !final_hla_list.is_empty() {
                // user gave us BAM files, so lets add in the HLA genes we have ready
                let batch_diplotyping = !cli_settings.hla_revert_method;
                let hla_calls = if batch_diplotyping {
                    diplotype_hla_batch(
                        &final_hla_list,
                        database,
                        bam_fns,
                        reference_genome.unwrap(),
                        debug_bam_writer.as_mut(),
                        debug_custom_writer.as_mut(),
                        cli_settings
                    )?
                } else {
                    // old method that has only been tested on A and B
                    diplotype_hla(
                        &final_hla_list,
                        database,
                        bam_fns,
                        reference_genome.unwrap(),
                        debug_bam_writer.as_mut(),
                        cli_settings
                    )?
                };

                // add in each result, insert will make sure we do not duplicate
                for (gene_name, gene_details) in hla_calls.into_iter() {
                    diplotypes.insert(gene_name, gene_details)?;
                }
            }
        }

        // assume include set has everything UNLESS it is specified
        let gene_name = "CYP2D6";
        let included = opt_include_set.as_ref().map(|include_set| include_set.contains(gene_name)).unwrap_or(true);
        if !included {
            debug!("Skipping {gene_name}, not in include set");
        }
        
        // assume exclude set has nothing UNLESS it is specified
        let excluded = opt_exclude_set.as_ref().map(|exclude_set| exclude_set.contains(gene_name)).unwrap_or(false);
        if excluded {
            debug!("Skipping {gene_name}, part of exclude set");
        }

        if included && !excluded {
            // CYP2D6 also requires a BAM file
            match diplotype_cyp2d6(
                database,
                bam_fns,
                reference_genome.unwrap(),
                debug_bam_writer.as_mut(),
                cli_settings
            ) {
                // happy path, we produced a D6 call
                Ok(cyp2d6_call) => diplotypes.insert("CYP2D6".to_string(), cyp2d6_call)?,
                // unhappy path - this could be an "expected" error due to something like low coverage
                //                OR it could be something where we want to propagate something unexpected so we get a user report
                Err(e) => {
                    // first, check if this is an error we "expect" to happen
                    use crate::cyp2d6::errors::CallerError;
                    if let Some(caller_error) = e.downcast_ref::<CallerError>() {
                        error!("Received error while calling CYP2D6: {caller_error}");
                        error!("Setting result to NO_MATCH state");
                        diplotypes.insert("CYP2D6".to_string(), PgxGeneDetails::no_match())?;
                    } else {
                        // not an expected error, propagate and maybe a user will message us for debugging
                        return Err(e);
                    }
                }
            };
        }

        // we finished all processing, finalize the debug BAM if it exists
        if let Some(dbw) = debug_bam_writer.as_mut() {
            // write all the records we have saved from the HLA/CYP2D6 processes
            match dbw.write_all_records() {
                Ok(()) => {},
                Err(e) => {
                    error!("Error while writing debug BAM: {e}");
                    error!("Continuing processes...");
                }
            }
        }

        if let Some(dcw) = debug_custom_writer.as_mut() {
            match dcw.write_session() {
                Ok(()) => {},
                Err(e) => {
                    error!("Error while writing HLA custom session: {e}");
                    error!("Continuing processess...");
                }
            }
        }

    } else {
        info!("No BAM files were provided, all alignment based diplotyping was skipped.");
    }

    Ok(diplotypes)
}

/// This is just a basic wrapper for variant-level metadata that we can tag on NormalizedVariants.
#[derive(Clone, Debug, Default, PartialEq)]
struct VariantMeta {
    /// CPIC variant ID
    pub variant_id: u64,
    /// CPIC variant name
    pub name: String,
    /// DBSNP ID, if available
    pub dbsnp_id: Option<String>
}

/// This will load all the haplotypes from the database and normalize them.
/// It returns the set of loaded variants AND the haplotypes.
/// Note that if a variant fails to normalize, the whole haplotype will be ignored.
/// 
/// Returns a tuple (`variant_hash`, `normalized_haplotypes`):
/// * `variant_hash` - a HashMap from a NormalizedVariant to the original variant metadata
/// * `normalized_haplotypes` - a Vec of NormalizedPgxHaplotypes loaded from the database
/// # Arguments
/// * `gene_entry` - the single gene entry from our database
/// * `reference_genome` - the pre-loaded reference genome; if None, some normalization steps will not happen
/// # Errors
/// * if a variant_id is present in a haplotype definition but not in the variant set (aka, undefined)
/// * if a variant has fewer than two alleles
/// * if a variant in the database is incomplete
#[allow(clippy::type_complexity)]
fn load_database_haplotypes(gene_entry: &PgxGene, reference_genome: Option<&ReferenceGenome>) 
    -> Result<(HashMap<NormalizedVariant, VariantMeta>, Vec<NormalizedPgxHaplotype>), Box<dyn std::error::Error>> {
    let mut normalized_haplotypes: Vec<NormalizedPgxHaplotype> = vec![];
    let mut normalized_variants: HashMap<NormalizedVariant, VariantMeta> = Default::default();

    let pgx_variants = gene_entry.variants();
    for (haplotype_name, pgx_haplotype) in gene_entry.defined_haplotypes() {
        // initialize a haplotype
        let mut normalized_haplotype = NormalizedPgxHaplotype::new(haplotype_name.clone());
        let mut normalized_variant_meta: Vec<VariantMeta> = vec![];
        let mut normalized: bool = true;
        for (variant_id, variant_allele) in pgx_haplotype.haplotype().iter() {
            let variant: &PgxVariant = pgx_variants.get(variant_id).ok_or(format!("variant {variant_id} is referenced but not defined"))?;
            let dbsnp: &Option<String> = variant.dbsnp_id();
            let variant_name: String = variant.name().to_string();
            
            let alleles = variant.alleles();
            if alleles.len() < 2 {
                bail!("Encountered variant {variant_id} with fewer than two alleles.");
            }
            for allele in alleles.iter() {
                if allele.is_none() {
                    bail!("Encountered variant {variant_id} with undefined alleles.");
                }
            }

            // ref is allele 0
            let ref_allele = alleles[0].as_deref().unwrap();
            let alt_allele = variant_allele.as_str();

            // check if this is a reference allele
            if ref_allele != alt_allele {
                match NormalizedVariant::multi_new(
                    gene_entry.chromosome().to_string(),
                    // 1-based -> 0-based
                    variant.position() - 1,
                    ref_allele,
                    alt_allele,
                    reference_genome
                ) {
                    Ok(nv) => {
                        // the allele was successfully constructed, add it to this haplotype
                        normalized_haplotype.add_variant(nv);
                        normalized_variant_meta.push(VariantMeta{
                            variant_id: *variant_id, 
                            name: variant_name,
                            dbsnp_id: dbsnp.clone()
                        });
                    },
                    Err(e) => {
                        warn!("Error while normalizing database variant {variant_id}: {e}, {variant:?}");
                        warn!("Ignoring {haplotype_name:?} due to variant incompatibility.");
                        normalized = false;
                        break;
                    }
                }
            } else {
                // this is a reference allele, we can just ignore for now
            }
        }

        if normalized {
            // we need to add every variant we found to the hash set
            assert_eq!(normalized_haplotype.variants().len(), normalized_variant_meta.len());
            for (or_variants, variant_meta) in normalized_haplotype.variants().iter().zip(normalized_variant_meta.into_iter()) {
                // this is a list of logically OR'ed variants, they should all have the same metadata though (except None)
                for opt_nv in or_variants.iter() {
                    match opt_nv {
                        Some(nv) => {
                            // add this variant if it is not already added
                            match normalized_variants.entry(nv.clone()) {
                                Occupied(entry) => {
                                    // we already saved this one, make sure there are no surprises
                                    assert_eq!(entry.get(), &variant_meta);
                                },
                                Vacant(entry) => {
                                    // new one to add
                                    entry.insert(variant_meta.clone());
                                }
                            };
                        },
                        None => {
                            // do nothing, it's just indiciating optional
                        }
                    };
                }
            }

            // now we need to save this haplotype
            normalized_haplotypes.push(normalized_haplotype);
        }
    }

    Ok(
        (normalized_variants, normalized_haplotypes)
    )
}

/// This will load all the identified variants from a VCF file.
/// # Arguments
/// * `vcf_fn` - the path to the VCF file, which must have an index
/// * `variant_hash` - the set of normalized variants we are looking for
/// * `reference_genome` - the pre-loaded reference genome; if None, some normalization steps will not happen
/// # Errors
/// * if there are issues opening the indexed VCF
/// * if the chromosome for the variants cannot be found in the VCF
/// # Panics
/// * if the variant_hash is empty, the chrom.unwrap() will panic
fn load_vcf_variants(vcf_fn: &Path, variant_hash: &HashMap<NormalizedVariant, VariantMeta>, reference_genome: Option<&ReferenceGenome>) 
    -> Result<HashMap<NormalizedVariant, NormalizedGenotype>, Box<dyn std::error::Error>> {
    // first we need to open up the vcf and pull out the header
    let mut vcf_reader: bcf::IndexedReader = bcf::IndexedReader::from_path(vcf_fn)?;
    let vcf_header: bcf::header::HeaderView = vcf_reader.header().clone();

    // prep our return struct
    let mut ret: HashMap<NormalizedVariant, NormalizedGenotype> = Default::default();

    // iterate over each variant, search for the corresponding entry in the VCF
    for (variant, _v_meta) in variant_hash.iter() {
        trace!("Searching for variant: {variant:?}");
        let chrom: &str = variant.chrom();
        let chrom_index: u32 = vcf_header.name2rid(chrom.as_bytes())?;
        let position: usize = variant.position();
        const BUFFER: usize = 50;
        let min_search: usize = position.saturating_sub(BUFFER);
        let max_search: usize = position.saturating_add(BUFFER);

        // if we find anything, this gets changed
        let mut search_genotype: Option<NormalizedGenotype> = None;
        match vcf_reader.fetch(chrom_index, min_search as u64, Some(max_search as u64)) {
            Ok(()) => {
                // we can iterate as normal
                for record_result in vcf_reader.records() {
                    let record: rust_htslib::bcf::Record = record_result?;
                    let alleles: Vec<&str> = record.alleles().iter()
                        .map(|a| std::str::from_utf8(a).unwrap_or("UTF8_ERROR"))
                        .collect();
                    let ref_allele = alleles[0];
                    
                    // TODO: if we want to allow for multi-sample VCFs, we need to adjust this
                    let all_genotypes = record.genotypes()?;
                    let genotype = all_genotypes.get(0);
                    if genotype.len() != 2 {
                        warn!("Error while parsing genotype.len() != 2, ignoring: {} {} {:?} => {:?}", chrom, record.pos(), alleles, genotype);
                        continue;
                    }
                    
                    // figure out what the genotype is
                    let mut is_phased = false;
                    let gt1 = match genotype[0] {
                        GenotypeAllele::Unphased(at) => Some(at),
                        GenotypeAllele::Phased(at) => {
                            is_phased = true;
                            Some(at)
                        },
                        //TODO: ignore these for now, not sure how to handle it?
                        GenotypeAllele::UnphasedMissing |
                        GenotypeAllele::PhasedMissing=> None
                    };
                    let gt2 = match genotype[1] {
                        GenotypeAllele::Unphased(at) => Some(at),
                        GenotypeAllele::Phased(at) => {
                            is_phased = true;
                            Some(at)
                        },
                        //TODO: ignore these for now, not sure how to handle it?
                        GenotypeAllele::UnphasedMissing |
                        GenotypeAllele::PhasedMissing=> None
                    };

                    // if we encounter empty genotypes, we will ignore them for now
                    if gt1.is_none() || gt2.is_none() {
                        warn!("Error while parsing incomplete genotype, ignoring: {} {} {:?} => {:?}", chrom, record.pos(), alleles, genotype);
                        continue;
                    }
                    let gt1 = gt1.unwrap() as usize;
                    let gt2 = gt2.unwrap() as usize;

                    let phase_set: Option<usize> = if is_phased {
                        // check for a phase set ID if we are phased
                        match record.format(b"PS").integer() {
                            Ok(all_ps_tag) => {
                                // phase set parsing was fine
                                let ps_tag = all_ps_tag[0];
                                assert_eq!(ps_tag.len(), 1);
                                Some(ps_tag[0] as usize)
                            },
                            Err(e) => {
                                // phase set parsing failed, which is weird
                                warn!("Failed to parse \"PS\" tag for variant, setting unphased: {} {} {:?} => {}", chrom, record.pos(), alleles, e);
                                is_phased = false;
                                None
                            }
                        }
                    } else {
                        // marked as unphased, so no ID
                        None
                    };
                    
                    for (alt_index, &alt_allele) in alleles.iter().enumerate().skip(1) {
                        match NormalizedVariant::new(
                            chrom.to_string(),
                            record.pos() as usize,
                            ref_allele,
                            alt_allele,
                            reference_genome
                        ) {
                            Ok(nv) => {
                                if nv == *variant {
                                    if alt_index == gt1 && alt_index == gt2 {
                                        // homozygous call - we know this can never be reference because we did .skip(1)
                                        if phase_set.is_some() {
                                            // we do not expect homozygous records to have a phase set
                                            bail!("Homozygous record detected with a phase set ID (PS): {}", record.desc());
                                        }
                                        assert!(search_genotype.is_none());
                                        search_genotype = Some(NormalizedGenotype::new(
                                            Genotype::HomozygousAlternate,
                                            phase_set
                                        ));
                                    } else if alt_index == gt1 && is_phased {
                                        // heterozygous call like 1|0
                                        if phase_set.is_none() {
                                            // heterozygous and phased, we need a PS tag
                                            bail!("Phased record detected without a phase set ID (PS): {}", record.desc());
                                        }
                                        assert!(search_genotype.is_none());
                                        search_genotype = Some(NormalizedGenotype::new(
                                            Genotype::HeterozygousPhasedFlip,
                                            phase_set
                                        ));
                                    } else if alt_index == gt2 && is_phased {
                                        // heterozygous call like 0|1
                                        if phase_set.is_none() {
                                            // heterozygous and phased, we need a PS tag
                                            bail!("Phased record detected without a phase set ID (PS): {}", record.desc());
                                        }
                                        assert!(search_genotype.is_none());
                                        search_genotype = Some(NormalizedGenotype::new(
                                            Genotype::HeterozygousPhased,
                                            phase_set
                                        ));
                                    } else if (alt_index == gt1 || alt_index == gt2) && !is_phased {
                                        // heterozygous call like 0/1 (this can handle situations like 1/2 also though)
                                        if phase_set.is_some() {
                                            // heterozygous and marked as unphased, so if we have a PS tag, something is weird
                                            bail!("Unphased heterozygous record detected with a phase set ID (PS): {}", record.desc());
                                        }
                                        assert!(search_genotype.is_none());
                                        search_genotype = Some(NormalizedGenotype::new(
                                            Genotype::HeterozygousUnphased,
                                            phase_set
                                        ));
                                    } else {
                                        // neither, possibly hom-reference or a different form of the allele
                                    }
                                } else {
                                    // this is a variant, just not one that matches what we are looking for
                                }
                            },
                            Err(e) => {
                                warn!("Error parsing VCF variant {} {} {:?}: {e}", chrom, record.pos(), alleles);
                            }
                        };
                    }
                }
            },
            Err(e) => {
                // this usually happens when there are no entries for the chromosome
                // error in the search, we will handle these later if they pop up
                warn!("Received \'{}\', while seeking to {}:{}-{}, assuming no variants present", e, chrom, min_search, max_search);
            }
        };

        // check if the search found anything
        if let Some(ng) = search_genotype {
            debug!("Genotype found for {variant:?}: {ng:?}");
            ret.insert(variant.clone(), ng);
        } else {
            trace!("Genotype not found.");
        }
    }

    Ok(ret)
}

/// This is the workhorse function for computing the diplotype.
/// This initial implementation requires an *exact* match of the two haplotypes.
/// There can be ambiguity though, even with exact matches, so a Vec is returned.
/// # Arguments
/// * `normalized_haplotypes` - all possible haplotypes that are defined
/// * `variant_calls` - the set of variants that were identified as well as their genotype
/// # Errors
/// * None so far
fn solve_diplotype(normalized_haplotypes: &[NormalizedPgxHaplotype], variant_calls: &HashMap<NormalizedVariant, NormalizedGenotype>) -> Result<Vec<Diplotype>, Box<dyn std::error::Error>> {
    // build up our base homozygous haplotype and also order the het variants at the same time
    let mut base_haplotype: Vec<NormalizedVariant> = vec![];
    let mut het_variants: Vec<NormalizedVariant> = vec![];
    let mut null_haplogroups: usize = 0;
    let mut identified_haplogroups: HashSet<usize> = Default::default();
    for (variant, genotype) in variant_calls.iter() {
        match genotype.genotype() {
            Genotype::HomozygousReference => panic!("we do not actually do anything with this"),
            Genotype::HomozygousAlternate => base_haplotype.push(variant.clone()),
            Genotype::HeterozygousUnphased | 
            Genotype::HeterozygousPhased | 
            Genotype::HeterozygousPhasedFlip => {
                // always save the variant
                het_variants.push(variant.clone());

                // figure out the haplogroup
                match genotype.phase_set() {
                    Some(ps) => { identified_haplogroups.insert(*ps); },
                    None => null_haplogroups += 1
                };
            }
        };
    }

    // now the magic
    let diplotype: Vec<Diplotype> = if het_variants.is_empty() {
        // there are no heterozygous variants, so we are returning a homozygous allele
        let mut matched: Option<String> = None;
        for haplotype in normalized_haplotypes.iter() {
            if haplotype.matches(&base_haplotype) {
                assert!(matched.is_none());
                matched = Some(haplotype.haplotype_name().to_string());
            }
        }

        // since we require exact matching, it's possible we fail to find something
        let match_name: String = matched.unwrap_or("NO_MATCH".to_string());
        vec![
            Diplotype::new(&match_name, &match_name)
        ]
    } else {
        // complicated path - we have heterozygous variants to resolve
        let total_haplogroups: usize = null_haplogroups + identified_haplogroups.len();

        // if we have 'x' hets, there are 2^(x-1) combinations due to symmetry; iterate over that combination count
        let max_combinations = 2_usize.pow(total_haplogroups as u32 - 1);
        let mut valid_diplotypes: Vec<Diplotype> = vec![];
        for combination in 0..max_combinations {
            // first resolve this combination into the two haplotypes
            let mut h1: Vec<NormalizedVariant> = base_haplotype.clone();
            let mut h2: Vec<NormalizedVariant> = base_haplotype.clone();

            let mut combo_index: usize = 0;
            let mut ps_lookup: HashMap<usize, bool> = Default::default();
            for hv in het_variants.iter() {
                // let is_h1: bool = ((combination >> i) & 0x1) != 0;
                let genotype = variant_calls.get(hv).unwrap();

                // this basically controls the permutation we are on
                // if we have 3 iterations of unphased variants, the first one is h1=000, h2=111; i.e. all REF on h1, all ALT on h2
                let is_h1: bool = match genotype.phase_set() {
                    Some(ps) => {
                        match ps_lookup.entry(*ps) {
                            Occupied(entry) => {
                                // we already assigned this PS
                                *entry.get()
                            },
                            Vacant(entry) => {
                                // we have not assigned this PS, so do so now
                                let r = ((combination >> combo_index) & 0x1) != 0;
                                entry.insert(r);
                                combo_index += 1;
                                r
                            }
                        }
                    },
                    None => {
                        // unphased, so we always combo bump
                        let r = ((combination >> combo_index) & 0x1) != 0;
                        combo_index += 1;
                        r
                    }
                };

                // controls orientation
                let orientation01: bool = match genotype.genotype() {
                    // unphased 0/1 and phased 0|1 orientations are treated the same
                    Genotype::HeterozygousUnphased |
                    Genotype::HeterozygousPhased => true,
                    // flipped 1|0 orientation
                    Genotype::HeterozygousPhasedFlip => false,
                    _ => panic!("we should not have homs in the het list")
                };

                if is_h1 == orientation01 {
                    // if (is_h1 AND normal orientation) OR (is_h2 AND flipped orientation)
                    h1.push(hv.clone());
                } else {
                    // if (is_h1 AND flipped orientatio) OR (is_h2 AND normal orientation)
                    h2.push(hv.clone());
                }
            }

            assert_eq!(combo_index, total_haplogroups);
            debug!("\t combination {} = {:?}, {:?}", combination, h1, h2);
            
            // now compare them to see if they both match an allele
            let mut h1_matched: Option<String> = None;
            let mut h2_matched: Option<String> = None;
            for haplotype in normalized_haplotypes.iter() {
                if haplotype.matches(&h1) {
                    assert!(h1_matched.is_none());
                    h1_matched = Some(haplotype.haplotype_name().to_string());
                }
                if haplotype.matches(&h2) {
                    assert!(h2_matched.is_none());
                    h2_matched = Some(haplotype.haplotype_name().to_string());
                }
            }

            if let (Some(h1_name), Some(h2_name)) = (h1_matched, h2_matched) {
                // both matched, save this combination
                valid_diplotypes.push(Diplotype::new(&h1_name, &h2_name));
            }
        }

        if valid_diplotypes.is_empty() {
            // no exact matching combinations were found
            vec![
                Diplotype::new("NO_MATCH", "NO_MATCH")
            ]
        } else {
            // we found one or more
            valid_diplotypes
        }
    };

    Ok(diplotype)
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::path::PathBuf;

    use crate::util::file_io::load_json;

    fn load_test_reference() -> ReferenceGenome {
        let ref_fn = PathBuf::from("test_data/test_reference.fa");
        ReferenceGenome::from_fasta(&ref_fn).unwrap()
    }

    fn create_dummy_cli() -> DiplotypeSettings {
        Default::default()
    }

    /// mainly test that the database values loaded matches expectations, CACNA1S is a relatively simple starter
    #[test]
    fn test_load_database_haplotypes() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/CACNA1S/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        let gene_entry: &PgxGene = database.gene_entries().get("CACNA1S").unwrap();

        // get the haplotypes
        let (normalized_variants, normalized_haplotypes) = load_database_haplotypes(gene_entry, None).unwrap();

        // check the variants
        let v1 = NormalizedVariant::new("chr1".to_string(), 201091992, "G", "A", None).unwrap();
        let v1_meta = VariantMeta { variant_id: 777260, name: "faux".to_string(), dbsnp_id: Some("rs772226819".to_string()) };
        let v2 = NormalizedVariant::new("chr1".to_string(), 201060814, "C", "T", None).unwrap();
        let v2_meta = VariantMeta { variant_id: 777261, name: "faux".to_string(), dbsnp_id: Some("rs1800559".to_string()) };
        let expected_variants = HashMap::from_iter(vec![
            (v1.clone(), v1_meta),
            (v2.clone(), v2_meta)
        ].into_iter());
        assert_eq!(normalized_variants, expected_variants);

        // check the haplotypes
        assert_eq!(normalized_haplotypes.len(), 3);
        let h1 = NormalizedPgxHaplotype::new("Reference".to_string());
        let mut h2 = NormalizedPgxHaplotype::new("c.3257G>A".to_string());
        h2.add_variant(vec![Some(v2)]);
        let mut h3 = NormalizedPgxHaplotype::new("c.520C>T".to_string());
        h3.add_variant(vec![Some(v1)]);
        let expected_haplotypes = vec![h1, h2, h3];
        assert_eq!(normalized_haplotypes, expected_haplotypes);
    }

    /// tests that we load variants from a homozygous only VCF correctly
    #[test]
    fn test_load_vcf_variants() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/CACNA1S/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        let gene_entry: &PgxGene = database.gene_entries().get("CACNA1S").unwrap();

        // get the haplotypes
        let (normalized_variants, _normalized_haplotypes) = load_database_haplotypes(gene_entry, None).unwrap();

        // now we need to load from VCF
        let vcf_fn = PathBuf::from("./test_data/CACNA1S/hom.vcf.gz");
        let vcf_variants = load_vcf_variants(&vcf_fn, &normalized_variants, None).unwrap();

        // make sure we have this variant as homozygous
        let expected_variant = NormalizedVariant::new("chr1".to_string(), 201060814, "C", "T", None).unwrap();
        let expected_genotype = NormalizedGenotype::new(Genotype::HomozygousAlternate, None);

        // we just expect the one hom call
        assert_eq!(vcf_variants.len(), 1);
        assert_eq!(*vcf_variants.get(&expected_variant).unwrap(), expected_genotype);
    }

    /// tests that bad variants through errors
    #[test]
    fn test_invalid_ps_vcf() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/CACNA1S/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        let gene_entry: &PgxGene = database.gene_entries().get("CACNA1S").unwrap();

        // get the haplotypes
        let (normalized_variants, _normalized_haplotypes) = load_database_haplotypes(gene_entry, None).unwrap();

        // now we need test the bad VCF
        let vcf_fn = PathBuf::from("./test_data/CACNA1S/bad_hom_ps.vcf.gz");
        let result = load_vcf_variants(&vcf_fn, &normalized_variants, None);
        assert!(result.is_err());
    }

    /// This is basically a full test of the single CACNA1S homozygous call
    #[test]
    fn test_solve_diplotype_hom() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/CACNA1S/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();

        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/CACNA1S/hom.vcf.gz");

        let diplotype = call_diplotypes(&database, Some(&vcf_fn), None, &[], &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("CACNA1S").unwrap().diplotypes(), vec![Diplotype::new("c.3257G>A", "c.3257G>A")]);
    }

    /// This is basically a full test of the single CACNA1S heterozygous call
    #[test]
    fn test_solve_diplotype_het() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/CACNA1S/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/CACNA1S/het.vcf.gz");

        let diplotype = call_diplotypes(&database, Some(&vcf_fn), None, &[], &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("CACNA1S").unwrap().diplotypes(), vec![Diplotype::new("Reference", "c.3257G>A")]);
    }

    /// This is basically a full test of the CACNA1S compound heterozygous call
    #[test]
    fn test_solve_diplotype_compound_het() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/CACNA1S/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/CACNA1S/compound_het.vcf.gz");

        let diplotype = call_diplotypes(&database, Some(&vcf_fn), None, &[], &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("CACNA1S").unwrap().diplotypes(), vec![Diplotype::new("c.520C>T", "c.3257G>A")]);
    }

    /// This is basically a full test of the RNR1 compound het for overlapping variants
    #[test]
    fn test_solve_diplotype_overlapping_compound_het() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/RNR1-faux/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        let reference_genome = load_test_reference();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/RNR1-faux/compound_het.vcf.gz");

        let diplotype = call_diplotypes(&database, Some(&vcf_fn), Some(&reference_genome), &[], &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("MT-RNR1").unwrap().diplotypes(), vec![Diplotype::new("961T>del", "961T>del+Cn")]);
    }

    /// This is basically a full test of the RNR1 homozygous for overlapping variants
    #[test]
    fn test_solve_diplotype_overlapping_hom() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/RNR1-faux/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        let reference_genome = load_test_reference();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/RNR1-faux/hom.vcf.gz");

        let diplotype = call_diplotypes(&database, Some(&vcf_fn), Some(&reference_genome), &[], &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("MT-RNR1").unwrap().diplotypes(), vec![Diplotype::new("961T>del+Cn", "961T>del+Cn")]);
    }

    /// This is basically a full test of the UGT1A1 variant *1/*80+*28
    #[test]
    fn test_solve_same_phase_001() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/UGT1A1-faux/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        let reference_genome = load_test_reference();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/UGT1A1-faux/same_phase_001.vcf.gz");
        let diplotype = call_diplotypes(&database, Some(&vcf_fn), Some(&reference_genome), &[], &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("UGT1A1").unwrap().diplotypes(), vec![Diplotype::new("*1", "*80+*28")]);
    }

    /// Same as above, but reversed; for now it keeps the alleles ordered by the phasing
    #[test]
    fn test_solve_same_phase_002() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/UGT1A1-faux/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        let reference_genome = load_test_reference();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/UGT1A1-faux/same_phase_002.vcf.gz");
        let diplotype = call_diplotypes(&database, Some(&vcf_fn), Some(&reference_genome), &[], &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("UGT1A1").unwrap().diplotypes(), vec![Diplotype::new("*80+*28", "*1")]);
    }

    /// Opposite phased alleles
    #[test]
    fn test_solve_opposite_phase_001() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/UGT1A1-faux/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        let reference_genome = load_test_reference();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/UGT1A1-faux/opposite_phase_001.vcf.gz");
        let diplotype = call_diplotypes(&database, Some(&vcf_fn), Some(&reference_genome), &[], &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("UGT1A1").unwrap().diplotypes(), vec![Diplotype::new("*28", "*80")]);
    }
    
    /// Opposite phased alleles, and swap to *37
    #[test]
    fn test_solve_opposite_phase_002() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/UGT1A1-faux/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        let reference_genome = load_test_reference();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/UGT1A1-faux/opposite_phase_002.vcf.gz");
        let diplotype = call_diplotypes(&database, Some(&vcf_fn), Some(&reference_genome), &[], &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("UGT1A1").unwrap().diplotypes(), vec![Diplotype::new("*80", "*37")]);
    }

    /// Homozygous allele + phased allele
    #[test]
    fn test_solve_hethom_phase_001() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/UGT1A1-faux/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        let reference_genome = load_test_reference();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/UGT1A1-faux/hethom_phase_001.vcf.gz");
        let diplotype = call_diplotypes(&database, Some(&vcf_fn), Some(&reference_genome), &[], &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("UGT1A1").unwrap().diplotypes(), vec![Diplotype::new("*80+*28", "*80+*37")]);
    }

    /// Different phase sets, so effectively unphased; should return two options
    #[test]
    fn test_solve_different_phaseset_001() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/UGT1A1-faux/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        let reference_genome = load_test_reference();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/UGT1A1-faux/different_phaseset_001.vcf.gz");
        let diplotype = call_diplotypes(&database, Some(&vcf_fn), Some(&reference_genome), &[], &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("UGT1A1").unwrap().diplotypes(), vec![
            Diplotype::new("*1", "*80+*28"),
            Diplotype::new("*28", "*80")
        ]);
    }

    /// Different phase sets, so effectively unphased; should return two options that shifts the *80 around
    #[test]
    fn test_solve_different_phaseset_002() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/UGT1A1-faux/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        let reference_genome = load_test_reference();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/UGT1A1-faux/different_phaseset_002.vcf.gz");
        let diplotype = call_diplotypes(&database, Some(&vcf_fn), Some(&reference_genome), &[], &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("UGT1A1").unwrap().diplotypes(), vec![
            Diplotype::new("*28", "*80+*37"),
            Diplotype::new("*37", "*80+*28")
        ]);
    }

    // test when we don't provide a VCF, we should get no calls
    // we prevent this in the CLI, but worth checking that we don't crash
    #[test]
    fn test_no_files() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("./data/v0.9.0/cpic_20240404.json.gz");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        
        // and here's the VCF
        let no_ref_fn = None;
        let no_vcf_fn = None;
        let no_bams = [];
        let diplotype = call_diplotypes(&database, no_vcf_fn, no_ref_fn, &no_bams, &create_dummy_cli()).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        // make sure we didn't call anything
        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 0);
    }

    #[test]
    fn test_include_set() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/CACNA1S/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/CACNA1S/compound_het.vcf.gz");

        // test a CLI that just includes CACNA1S
        let mut cli = create_dummy_cli();
        cli.include_fn = Some(PathBuf::from("test_data/CACNA1S/CACNA1S_gene_list.txt"));

        let diplotype = call_diplotypes(&database, Some(&vcf_fn), None, &[], &cli).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("CACNA1S").unwrap().diplotypes(), vec![Diplotype::new("c.520C>T", "c.3257G>A")]);

        // now do the same test, but with an empty include list
        let mut cli = create_dummy_cli();
        cli.include_fn = Some(PathBuf::from("test_data/empty_gene_list.txt"));

        let diplotype = call_diplotypes(&database, Some(&vcf_fn), None, &[], &cli).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 0);
    }

    #[test]
    fn test_exclude_set() {
        // load the database
        let database_fn: PathBuf = PathBuf::from("test_data/CACNA1S/database.json");
        let database: PgxDatabase = load_json(&database_fn).unwrap();
        
        // and here's the VCF
        let vcf_fn = PathBuf::from("./test_data/CACNA1S/compound_het.vcf.gz");

        // test a CLI that excludes nothing from the list
        let mut cli = create_dummy_cli();
        cli.exclude_fn = Some(PathBuf::from("test_data/empty_gene_list.txt"));

        let diplotype = call_diplotypes(&database, Some(&vcf_fn), None, &[], &cli).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 1);
        assert_eq!(*diplotypes.get("CACNA1S").unwrap().diplotypes(), vec![Diplotype::new("c.520C>T", "c.3257G>A")]);

        // now do the same test, but with our gene in the exclude list
        let mut cli = create_dummy_cli();
        cli.exclude_fn = Some(PathBuf::from("test_data/CACNA1S/CACNA1S_gene_list.txt"));
        
        let diplotype = call_diplotypes(&database, Some(&vcf_fn), None, &[], &cli).unwrap();

        // make sure the metadata matches exactly
        assert_eq!(database.database_metadata(), diplotype.database_metadata());

        let diplotypes = diplotype.gene_details();
        assert_eq!(diplotypes.len(), 0);
    }
}