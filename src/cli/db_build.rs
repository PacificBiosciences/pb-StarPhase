

use clap::Args;
use log::info;
use std::path::PathBuf;

use crate::cli::core::{check_optional_filename, check_required_filename, AFTER_HELP};

#[derive(Clone, Args)]
#[clap(author, about, 
    after_help = &**AFTER_HELP)]
pub struct BuildSettings {
    /// Reference FASTA file
    #[clap(short = 'r')]
    #[clap(long = "reference")]
    #[clap(value_name = "FASTA")]
    #[clap(help_heading = Some("Input/Output"))]
    pub reference_filename: PathBuf,

    /// Output database location (JSON)
    #[clap(short = 'o')]
    #[clap(long = "output-db")]
    #[clap(value_name = "JSON")]
    #[clap(help_heading = Some("Input/Output"))]
    pub output_database: PathBuf,

    /// User-provided build options (JSON)
    #[clap(short = 'b')]
    #[clap(long = "build-options")]
    #[clap(value_name = "JSON")]
    #[clap(help_heading = Some("Input/Output"))]
    pub build_options: Option<PathBuf>,

    /// Enable verbose output.
    #[clap(short = 'v')]
    #[clap(long = "verbose")]
    #[clap(action = clap::ArgAction::Count)]
    pub verbosity: u8,
}

pub fn check_build_settings(settings: BuildSettings) -> BuildSettings {
    // dump stuff to the logger
    check_required_filename(&settings.reference_filename, "Reference FASTA");
    
    info!("Reference: {:?}", &settings.reference_filename);
    info!("Output database: {:?}", settings.output_database);

    // check the build options
    check_optional_filename(settings.build_options.as_deref(), "Build options");
    if let Some(build_options) = settings.build_options.as_ref() {
        info!("Build options: {build_options:?}");
    } else {
        info!("Build options: default");
    }

    settings
}
