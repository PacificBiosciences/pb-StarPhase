

use clap::Args;
use log::info;
use std::path::PathBuf;

use crate::cli::core::AFTER_HELP;

#[derive(Clone, Args)]
#[clap(author, about, 
    after_help = &**AFTER_HELP)]
pub struct BuildSettings {
    /// Output database location (JSON)
    #[clap(required = true)]
    #[clap(short = 'o')]
    #[clap(long = "output-db")]
    #[clap(value_name = "JSON")]
    #[clap(help_heading = Some("Input/Output"))]
    pub output_database: PathBuf,

    /// Enable verbose output.
    #[clap(short = 'v')]
    #[clap(long = "verbose")]
    #[clap(action = clap::ArgAction::Count)]
    pub verbosity: u8,
}

pub fn check_build_settings(settings: BuildSettings) -> BuildSettings {
    // dump stuff to the logger
    info!("Output database: {:?}", settings.output_database);
    settings
}
