
use log::{LevelFilter, error, info};
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use serde::Serialize;
use std::fs::File;
use std::path::Path;

use pbstarphase::cli::diplotype::{DiplotypeSettings, check_diplotype_settings};
use pbstarphase::cli::core::{Commands, get_cli};
use pbstarphase::cli::db_build::{BuildSettings, check_build_settings};
use pbstarphase::cli::db_stat::{DbStatSettings, check_db_stat_settings};
use pbstarphase::data_types::pgx_diplotype::Diplotype;
use pbstarphase::data_types::starphase_json::StarphaseJson;
use pbstarphase::database::pgx_database::PgxDatabase;
use pbstarphase::database::db_config::{DatabaseBuildOptions};
use pbstarphase::util::file_io::{load_json, save_json};

/// This will run the "build" mode of the tool
/// # Arguments
/// * `settings` - the BuildSettings object
fn run_build(settings: BuildSettings) {
    // get the settings
    // let settings: Settings = get_raw_settings();
    let filter_level: LevelFilter = match settings.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace
    };

    // immediately setup logging first
    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        .init();

    // okay, now we can check all the other settings
    let cli_settings: BuildSettings = check_build_settings(settings);

    // now create our build options, this should also be fast and likely never fail
    let build_options: DatabaseBuildOptions = if let Some(bo_filename) = cli_settings.build_options.as_ref() {
        match load_json(bo_filename) {
            Ok(bo) => bo,
            Err(e) => {
                error!("Error while loading database build options: {e}");
                std::process::exit(exitcode::IOERR);
            }
        }
    } else {
        DatabaseBuildOptions::default()
    };
    info!("Build options: {build_options:#?}");

    // pre-load the reference genome also
    info!("Loading reference genome from {:?}...", cli_settings.reference_filename);
    let reference_genome: ReferenceGenome = match ReferenceGenome::from_fasta(&cli_settings.reference_filename) {
        Ok(rg) => rg,
        Err(e) => {
            error!("Error while loading reference genome file: {e}");
            std::process::exit(exitcode::IOERR);
        }
    };

    // all the work
    let pgx_db: PgxDatabase = match pbstarphase::build_database::build_database_via_api(
        &build_options,
        &reference_genome
    ) {
        Ok(pdb) => pdb,
        Err(e) => {
            error!("Error while building StarPhase database: {e}");
            std::process::exit(exitcode::IOERR);
        }
    };

    // debug!("Full database:\n{pgx_db:#?}");
    // save the database to the defined file
    info!("Saving database to {:?}", cli_settings.output_database);
    match save_json(&pgx_db, &cli_settings.output_database) {
        Ok(()) => {},
        Err(e) => {
            error!("Error while writing database to file: {e}");
            std::process::exit(exitcode::IOERR);
        }
    };
}

/// This will run the "build" mode of the tool
/// # Arguments
/// * `settings` - the BuildSettings object
fn run_diplotype(settings: DiplotypeSettings) {
    // get the settings
    // let settings: Settings = get_raw_settings();
    let filter_level: LevelFilter = match settings.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace
    };

    // immediately setup logging first
    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        .init();

    // okay, now we can check all the other settings
    let cli_settings: DiplotypeSettings = match check_diplotype_settings(settings) {
        Ok(s) => s,
        Err(e) => {
            error!("Error while processing CLI settings: {e}");
            std::process::exit(exitcode::USAGE);
        }
    };

    // create a debug folder if specified
    if let Some(debug_folder) = cli_settings.debug_folder.as_ref() {
        info!("Creating debug folder at {debug_folder:?}...");
        match std::fs::create_dir_all(debug_folder) {
            Ok(()) => {},
            Err(e) => {
                error!("Error while creating debug folder: {e}");
                std::process::exit(exitcode::IOERR);
            }
        }
    }

    // first load the database
    info!("Loading PGx database from {:?}...", cli_settings.input_database);
    let database: PgxDatabase = match load_json(&cli_settings.input_database) {
        Ok(db) => db,
        Err(e) => {
            error!("Error while loading PGx database file: {e}");
            std::process::exit(exitcode::IOERR);
        }
    };

    // we also need to validate that the database is complete enough to run
    if let Err(e) = database.validate() {
        error!("Error while validating PGx database file: {e}");
        std::process::exit(exitcode::IOERR);
    }

    // pre-load the reference genome also
    info!("Loading reference genome from {:?}...", cli_settings.reference_filename);
    let reference_genome: ReferenceGenome = match ReferenceGenome::from_fasta(&cli_settings.reference_filename) {
        Ok(rg) => rg,
        Err(e) => {
            error!("Error while loading reference genome file: {e}");
            std::process::exit(exitcode::IOERR);
        }
    };

    // now hand it to the diplotype caller
    let diplotypes: StarphaseJson = match pbstarphase::diplotyper::call_diplotypes(
        &database,
        cli_settings.vcf_filename.as_deref(),
        Some(&reference_genome),
        &cli_settings.bam_filenames,
        &cli_settings
    ) {
        Ok(dc) => dc,
        Err(e) => {
            error!("Error while calling diplotypes: {e}");
            std::process::exit(exitcode::DATAERR);
        }
    };

    // debug!("Full diplotypes:\n{diplotypes:#?}");
    // save the diplotypes to the defined file
    info!("Saving diplotypes to {:?}", cli_settings.diplotype_filename);
    match save_json(&diplotypes, &cli_settings.diplotype_filename) {
        Ok(()) => {},
        Err(e) => {
            error!("Error while writing diplotypes to file: {e}");
            std::process::exit(exitcode::IOERR);
        }
    };

    if let Some(filename) = cli_settings.pharmcat_tsv.as_ref() {
        info!("Saving PharmCAT diplotypes to {:?}", filename);
        match save_pharmcat_tsv(&diplotypes, filename) {
            Ok(()) => {},
            Err(e) => {
                error!("Error while writing PharmCAT diplotypes to file: {e}");
                std::process::exit(exitcode::IOERR);
            }
        };
    }
}

/// Wrapper for the pharmcat output
#[derive(Serialize)]
struct PharmCatRow {
    #[serde(rename = "#gene")]
    gene: String,
    diplotype: String
}

/// Helper function to save the basic TSV file for feeding into PharmCAT
/// # Arguments
/// * `diplotypes` - our reported diplotypes
/// * `filename` - the output filename, TSV
/// # Errors
/// * if we have any errors opening or writing to the file
fn save_pharmcat_tsv(diplotypes: &StarphaseJson, filename: &Path) -> Result<(), Box<dyn std::error::Error>> {
    // modify the delimiter to "," if it ends with .csv
    let delimiter: u8 = b'\t';
    let mut csv_writer: csv::Writer<File> = csv::WriterBuilder::new()
        .delimiter(delimiter)
        .from_path(filename)?;

    // make sure we go through the blocks in order
    for (gene, details) in diplotypes.gene_details().iter() {
        // check if we have one result or multiple
        let diplotypes = details.dedup_simple_diplotypes();
        let diplotype_dt: Diplotype = if diplotypes.len() > 1 {
            Diplotype::new("Multiple", "Multiple")
        } else {
            diplotypes[0].clone()
        };

        // check if it's a haplotype or diplotype gene
        let diplotype: String = if gene == "MT-RNR1" {
            // PharmCAT only accepts a single haplotype for MT (makes sense)
            diplotype_dt.homozygous_haplotype().unwrap_or("Unknown").to_string()
        } else {
            // all others are a diplotype
            diplotype_dt.pharmcat_diplotype()
        };

        // write the row out
        let block_row = PharmCatRow {
            gene: gene.clone(),
            diplotype
        };
        csv_writer.serialize(&block_row)?;
    }
    csv_writer.flush()?;
    Ok(())
}

/// This will run the "db_stat" mode of the tool
/// # Arguments
/// * `settings` - the DbStatSettings object
fn run_db_stat(settings: DbStatSettings) {
    // get the settings
    let filter_level: LevelFilter = match settings.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace
    };

    // immediately setup logging first
    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        .init();

    // okay, now we can check all the other settings
    let cli_settings: DbStatSettings = check_db_stat_settings(settings);

    // first load the database
    info!("Loading PGx database from {:?}...", cli_settings.input_database);
    let database: PgxDatabase = match load_json(&cli_settings.input_database) {
        Ok(db) => db,
        Err(e) => {
            error!("Error while loading PGx database file: {e}");
            std::process::exit(exitcode::IOERR);
        }
    };

    // we also need to validate that the database is complete enough to run
    if let Err(e) = database.validate() {
        error!("Error while validating PGx database file: {e}");
        std::process::exit(exitcode::IOERR);
    }
    info!("Database loaded successfully.");

    // display the database statistics
    pbstarphase::db_stat::print_stats(&database);
}

fn main() {
    let cli = get_cli();
    match cli.command {
        Commands::Build(settings) => {
            run_build(*settings);
        },
        Commands::Diplotype(settings) => {
            run_diplotype(*settings);
        },
        Commands::DbStat(settings) => {
            run_db_stat(*settings);
        }
    }

    info!("Process finished successfully.");
}