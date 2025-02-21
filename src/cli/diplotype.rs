

use clap::Args;
use log::{debug, info, warn};
use simple_error::bail;
use std::path::PathBuf;

use crate::cli::core::{AFTER_HELP, check_optional_filename, check_required_filename};

#[derive(Args, Clone, Default)]
#[clap(author, about, 
    after_help = &**AFTER_HELP)]
pub struct DiplotypeSettings {
    /// Input database file (JSON)
    #[clap(required = true)]
    #[clap(short = 'd')]
    #[clap(long = "database")]
    #[clap(value_name = "JSON")]
    #[clap(help_heading = Some("Input/Output"))]
    pub input_database: PathBuf,

    /// Reference FASTA file
    #[clap(required = true)]
    #[clap(short = 'r')]
    #[clap(long = "reference")]
    #[clap(value_name = "FASTA")]
    #[clap(help_heading = Some("Input/Output"))]
    pub reference_filename: PathBuf,

    /// Input variant file in VCF format
    #[clap(short = 'c')]
    #[clap(long = "vcf")]
    #[clap(value_name = "VCF")]
    #[clap(help_heading = Some("Input/Output"))]
    pub vcf_filename: Option<PathBuf>,

    /// Input structural variant file in VCF format
    #[clap(short = 's')]
    #[clap(long = "sv-vcf")]
    #[clap(value_name = "VCF")]
    #[clap(help_heading = Some("Input/Output"))]
    pub sv_vcf_filename: Option<PathBuf>,

    /// Input alignment file in BAM format, can be specified multiple times; required for HLA diplotyping
    #[clap(short = 'b')]
    #[clap(long = "bam")]
    #[clap(value_name = "BAM")]
    #[clap(help_heading = Some("Input/Output"))]
    pub bam_filenames: Vec<PathBuf>,

    /// Output diplotype call file (JSON)
    #[clap(required = true)]
    #[clap(short = 'o')]
    #[clap(long = "output-calls")]
    #[clap(value_name = "JSON")]
    #[clap(help_heading = Some("Input/Output"))]
    pub diplotype_filename: PathBuf,

    /// Output file that can be provided to PharmCAT for further call interpretation
    #[clap(long = "pharmcat-tsv")]
    #[clap(value_name = "TSV")]
    #[clap(help_heading = Some("Input/Output"))]
    pub pharmcat_tsv: Option<PathBuf>,

    /// Optional file indicating the list of genes to include in diplotyping, one per line
    #[clap(long = "include-set")]
    #[clap(value_name = "TXT")]
    #[clap(help_heading = Some("Input/Output"))]
    pub include_fn: Option<PathBuf>,

    /// Optional file indicating the list of genes to exclude from diplotyping, one per line
    #[clap(long = "exclude-set")]
    #[clap(value_name = "TXT")]
    #[clap(help_heading = Some("Input/Output"))]
    pub exclude_fn: Option<PathBuf>,

    /// Optional output debug folder
    #[clap(long = "output-debug")]
    #[clap(value_name = "DIR")]
    #[clap(help_heading = Some("Input/Output"))]
    pub debug_folder: Option<PathBuf>,

    /// Enables scoring by cDNA and tie-breaking with DNA
    #[clap(hide = true)]
    #[clap(long = "disable-cdna-scoring")]
    #[clap(help_heading = Some("HLA calling"))]
    pub disable_cdna_scoring: bool,

    /// Requires HLA alleles to have a DNA sequence definition
    #[clap(long = "hla-require-dna")]
    #[clap(help_heading = Some("HLA calling"))]
    pub hla_require_dna: bool,

    /// The maximum error rate for a read to the HLA reference allele
    #[clap(long = "max-error-rate")]
    #[clap(value_name = "FLOAT")]
    #[clap(default_value = "0.07")]
    #[clap(help_heading = Some("HLA calling"))]
    pub max_error_rate: f64,

    /// The minimum cumulative distribution function probability for a heterozygous call
    #[clap(long = "min-cdf-prob")]
    #[clap(value_name = "FLOAT")]
    #[clap(default_value = "0.001")]
    #[clap(help_heading = Some("HLA calling"))]
    pub min_cdf: f64,

    /// Expected minor allele frequency; reduce to account for skew from sequencing bias
    #[clap(long = "expected-maf")]
    #[clap(value_name = "FLOAT")]
    #[clap(default_value = "0.45")]
    #[clap(help_heading = Some("HLA calling"))]
    pub expected_maf: f64,

    /// Reverts to the old HLA calling algorithm for just HLA-A and HLA-B
    #[clap(hide = true)]
    #[clap(long = "hla-revert-method")]
    #[clap(help_heading = Some("HLA calling"))]
    pub hla_revert_method: bool,

    /// Additional HLA targets for the debug BAM file
    #[clap(hide = true)]
    #[clap(long = "debug-hla-target")]
    #[clap(value_name = "HLA_ID")]
    #[clap(help_heading = Some("HLA debug"))]
    pub debug_hla_targets: Vec<String>,

    /// Allows us to skip HLA for the purpose of debugging quickly
    #[clap(hide = true)]
    #[clap(long = "debug-skip-hla")]
    #[clap(help_heading = Some("HLA debug"))]
    pub debug_skip_hla: bool,

    /// (Deprecated) Optional output realignment file for CYP2D6 consensus sequences
    #[clap(hide = true)]
    #[clap(long = "output-cyp2d6-bam")]
    #[clap(value_name = "BAM")]
    #[clap(help_heading = Some("CYP2D6 calling"))]
    pub cyp2d6_bam_filename: Option<PathBuf>,

    /// Enables inferrence of connected alleles based on population observations
    #[clap(long = "infer-connections")]
    #[clap(help_heading = Some("CYP2D6 calling"))]
    pub infer_connections: bool,

    /// Disables normalizing coverage with D7 and hybrid alleles
    #[clap(long = "normalize-d6-only")]
    #[clap(help_heading = Some("CYP2D6 calling"))]
    pub normalize_d6_only: bool,

    /// The minimum fraction of sequences required to split into multiple consensuses (e.g. MAF)
    #[clap(long = "min-consensus-fraction")]
    #[clap(value_name = "FLOAT")]
    #[clap(default_value = "0.10")]
    #[clap(help_heading = Some("Consensus (HLA and CYP2D6)"))]
    pub min_consensus_fraction: f64,

    /// The minimum counts of sequences required to split into multiple consensuses
    #[clap(long = "min-consensus-count")]
    #[clap(value_name = "COUNT")]
    #[clap(default_value = "3")]
    #[clap(help_heading = Some("Consensus (HLA and CYP2D6)"))]
    pub min_consensus_count: u64,

    /// The edit distance delta threshold to stop tracking divergent sequences (efficiency heuristic)
    #[clap(long = "dual-max-ed-delta")]
    #[clap(value_name = "COUNT")]
    #[clap(default_value = "100")]
    #[clap(help_heading = Some("Consensus (HLA and CYP2D6)"))]
    pub dual_max_ed_delta: usize,

    /// Number of threads to use for phasing.
    #[clap(hide = true)]
    #[clap(short = 't')]
    #[clap(long = "threads")]
    #[clap(value_name = "THREADS")]
    #[clap(default_value = "1")]
    pub threads: usize,

    /// Enable verbose output.
    #[clap(short = 'v')]
    #[clap(long = "verbose")]
    #[clap(action = clap::ArgAction::Count)]
    pub verbosity: u8,
}

pub fn check_diplotype_settings(mut settings: DiplotypeSettings) -> Result<DiplotypeSettings, Box<dyn std::error::Error>> {
    info!("Inputs:");

    // check for all the required input files
    check_required_filename(&settings.input_database, "Database JSON");
    check_required_filename(&settings.reference_filename, "Reference FASTA");
    check_optional_filename(settings.vcf_filename.as_deref(), "VCF file");
    check_optional_filename(settings.sv_vcf_filename.as_deref(), "SV VCF file");

    // these are optional, but make sure that any specified exist
    for bam_fn in settings.bam_filenames.iter() {
        check_required_filename(bam_fn, "Alignment file");
    }

    // dump stuff to the logger
    info!("\tDatabase: {:?}", settings.input_database);
    info!("\tReference: {:?}", &settings.reference_filename);
    if let Some(vcf_fn) = settings.vcf_filename.as_ref() {
        info!("\tVCF: {:?}", vcf_fn);
        if let Some(sv_fn) = settings.sv_vcf_filename.as_ref() {
            info!("\tSV VCF: {:?}", sv_fn);
        } else {
            info!("\tSV VCF: None");
        }
    } else {
        warn!("\tVCF: No variant call files provided, all variant-based diplotyping is disabled")
    }
    if settings.bam_filenames.is_empty() {
        warn!("\tBAM: No alignment files provided, HLA and CYP2D6 diplotyping is disabled");
    } else {
        // these are optional, but make sure that any specified exist
        for bam_fn in settings.bam_filenames.iter() {
            info!("\tBAM: {:?}", bam_fn);
        }
    }

    if settings.vcf_filename.is_none() && settings.bam_filenames.is_empty() {
        // user didn't provide any data, bail out
        bail!("Must provide a VCF file and/or aligned BAM file to perform diplotyping.");
    }

    if settings.include_fn.is_some() && settings.exclude_fn.is_some() {
        bail!("Only one of --exclude-set and --include-set can be specified.");
    }
    if let Some(ifn) = settings.include_fn.as_ref() {
        check_required_filename(ifn, "Include set");
        info!("\tInclude set file: {ifn:?}");
    }
    if let Some(efn) = settings.exclude_fn.as_ref() {
        check_required_filename(efn, "Exclude set");
        info!("\tExclude set file: {efn:?}");
    }

    // outputs
    info!("Outputs:");
    info!("\tDiplotype calls: {:?}", settings.diplotype_filename);
    if let Some(filename) = settings.pharmcat_tsv.as_ref() {
        info!("\tPharmCAT TSV: {:?}", filename);
    }
    if let Some(debug_folder) = settings.debug_folder.as_ref() {
        debug!("\tDebug folder: {debug_folder:?}");
    }

    // miscellaneous settings
    if !settings.bam_filenames.is_empty() {
        info!("HLA settings:");
        if settings.disable_cdna_scoring {
            info!("\tDisable cDNA scoring: {}", settings.disable_cdna_scoring);
            if !settings.hla_require_dna {
                settings.hla_require_dna = true;
                info!("\tAutomatically enabling HLA DNA requirement")
            }
        }
        
        if settings.hla_require_dna {
            info!("\tRequire DNA for HLA calls: {}", settings.hla_require_dna);
        }

        if !(0.0..=1.0).contains(&settings.max_error_rate) {
            bail!("--max-error-rate must be between 0.0 and 1.0");
        }
        info!("\tMax read error rate: {}", settings.max_error_rate);
        
        if !(0.0..=1.0).contains(&settings.min_cdf) {
            bail!("--min-cdf-prob must be between 0.0 and 1.0");
        }
        info!("\tMinimum CDF probability: {}", settings.min_cdf);

        if !(0.01..=0.5).contains(&settings.expected_maf) {
            bail!("--expected-maf must be between 0.01 and 0.5");
        }
        info!("\tExpected MAF: [{}, 0.5]", settings.expected_maf);

        if settings.debug_folder.is_some() {
            debug!("\tHLA debug targets: {:?}", settings.debug_hla_targets);
        }

        if settings.hla_revert_method {
            debug!("\tHLA revert method: {}", settings.hla_revert_method);
        }

        if settings.cyp2d6_bam_filename.is_some() {
            warn!("The --output-cyp2d6-bam option is deprecated, use --output-debug instead.");
        }

        info!("CYP2D6 settings:");
        info!("\tConnection inferrence: {}", if settings.infer_connections { "ENABLED" } else { "DISABLED" });
        info!("\tNormalize D6 only: {}", if settings.normalize_d6_only { "ENABLED" } else { "DISABLED" });

        info!("Consensus settings:");
        if !(0.0..=1.0).contains(&settings.min_consensus_fraction) {
            bail!("--min-consensus-fraction must be between 0.0 and 1.0");
        }
        info!("\tMinimum consensus fraction: {}", settings.min_consensus_fraction);
        info!("\tMinimum consensus count: {}", settings.min_consensus_count);
        info!("\tDual max edit distance delta: {}", settings.dual_max_ed_delta);
        
        if settings.threads == 0 {
            settings.threads = 1;
        }
        if settings.threads != 1 {
            warn!("Threads (deprecated): {}", settings.threads);
        }
    }

    Ok(settings)
}
