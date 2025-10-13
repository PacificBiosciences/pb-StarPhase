use clap::Args;
use log::info;
use std::path::PathBuf;

use crate::cli::core::{check_required_filename, AFTER_HELP};

#[derive(Clone, Args)]
#[clap(author, about, 
    after_help = &**AFTER_HELP)]
pub struct DbStatSettings {
    /// Input database file (JSON)
    #[clap(required = true)]
    #[clap(short = 'd')]
    #[clap(long = "database")]
    #[clap(value_name = "JSON")]
    #[clap(help_heading = Some("Input/Output"))]
    pub input_database: PathBuf,

    /// Enable verbose output.
    #[clap(short = 'v')]
    #[clap(long = "verbose")]
    #[clap(action = clap::ArgAction::Count)]
    pub verbosity: u8,
}

pub fn check_db_stat_settings(settings: DbStatSettings) -> DbStatSettings {
    // dump stuff to the logger
    check_required_filename(&settings.input_database, "Database JSON");
    
    info!("Input database: {:?}", &settings.input_database);
    
    settings
}
