
use clap::{Parser, Subcommand};
use chrono::Datelike;
use lazy_static::lazy_static;
use log::error;
use std::path::Path;

use crate::cli::diplotype::DiplotypeSettings;
use crate::cli::db_build::BuildSettings;
use crate::cli::db_stat::DbStatSettings;

lazy_static! {
    /// Stores the full version string we plan to use, which is generated in build.rs
    /// # Examples
    /// * `0.11.0-6bb9635-dirty` - while on a dirty branch
    /// * `0.11.0-6bb9635` - with a fresh commit
    pub static ref FULL_VERSION: String = format!("{}-{}", env!("CARGO_PKG_VERSION"), env!("VERGEN_GIT_DESCRIBE"));

    /// Shared after help string containing the legalese.
    pub static ref AFTER_HELP: String = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year());
}

#[derive(Parser)]
#[clap(author, 
    version = &**FULL_VERSION, 
    about, 
    after_help = &**AFTER_HELP)]
#[command(propagate_version = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands
}

// Here lies PharmGOAT, you were truly the greatest of all time.
//               ,,~~--___---,
//              /            .~,
//        /  _,~             )
//       (_-(~)   ~, ),,,(  /'
//        Z6  .~`' ||     \ |
//        /_,/     ||      ||
//  ~~~~~~~~~~~~~~~W`~~~~~~W`~~~~~~~~~
// PharmGOAT, a tool for diplotyping PGx genes from HiFi data.

/// pb-StarPhase, a tool for diplotyping PGx genes from HiFi data.
/// Select a subcommand to see more usage information:
#[derive(Subcommand)]
pub enum Commands {
    /// Download and build the database for StarPhase
    Build(Box<BuildSettings>),
    /// Generate statistics about a database file
    DbStat(Box<DbStatSettings>),
    /// Run the diplotyper on a dataset
    Diplotype(Box<DiplotypeSettings>),
}

pub fn get_cli() -> Cli {
    Cli::parse()
}

/// Checks if a file exists and will otherwise exit
/// # Arguments
/// * `filename` - the file path to check for
/// * `label` - the label to use for error messages
pub fn check_required_filename(filename: &Path, label: &str) {
    if !filename.exists() {
        error!("{} does not exist: \"{}\"", label, filename.display());
        std::process::exit(exitcode::NOINPUT);
    } else {
        // file exists, we're good
    }
}

/// Checks if a file exists and will otherwise exit
/// # Arguments
/// * `filename` - the file path to check for
/// * `label` - the label to use for error messages
pub fn check_optional_filename(opt_filename: Option<&Path>, label: &str) {
    if let Some(filename) = opt_filename {
        if !filename.exists() {
            error!("{} does not exist: \"{}\"", label, filename.display());
            std::process::exit(exitcode::NOINPUT);
        } else {
            // file exists, we're good
        }
    }
}