
use bio::io::fasta;
use rustc_hash::FxHashSet as HashSet;
use simple_error::bail;
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::fs::File;
use std::path::Path;

/// Helper function that loads a file into some type, helpful generic
/// # Arguments
/// * `filename` - the file path to open and parse
/// # Errors
/// * if the file does not open properly
/// * if the deserialization throws errors
pub fn load_json<T: serde::de::DeserializeOwned>(filename: &Path) -> Result<T, Box<dyn std::error::Error>> {
    let fp: Box<dyn std::io::Read> = if filename.extension().unwrap_or_default() == "gz" {
        Box::new(
            flate2::read::MultiGzDecoder::new(
                File::open(filename)?
            )
        )
    } else {
        Box::new(File::open(filename)?)
    };
    let result: T = serde_json::from_reader(fp)?;
    Ok(result)
}

/// This will save a generic serializable struct to JSON.
/// # Arguments
/// * `data` - the data in memory
/// * `out_filename` - user provided path to write to 
/// # Errors
/// * if opening or writing to the file throw errors
/// * if JSON serialization throws errors
pub fn save_json<T: serde::Serialize>(data: &T, out_filename: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let file: Box<dyn std::io::Write> = if out_filename.extension().unwrap_or_default() == "gz" {
        Box::new(
            flate2::write::GzEncoder::new(
                File::create(out_filename)?,
                flate2::Compression::best()
            )
        )
    } else {
        Box::new(File::create(out_filename)?)
    };
    let mut writer = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut writer, data)?;
    writer.flush()?;
    Ok(())
}

/// Helper function that will just read a file line-by-line and return the list as a HashSet.
/// # Arguments
/// * `filename` - The file to load into the hash set
/// # Errors
/// * if a file is provided but cannot be opened or read
pub fn load_file_lines(filename: &Path) -> Result<HashSet<String>, Box<dyn std::error::Error>> {
    // open the file and throw into a buffered reader
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    // now add each line
    let mut hashset: HashSet<String> = Default::default();
    for line in reader.lines() {
        hashset.insert(line?);
    }
    Ok(hashset)
}

/// Given a map of keys and values, this will write out them out as a FASTA file
/// # Arguments
/// * `data` - each key is the entry label, each value is the DNA sequence
/// * `filename` - location to save fasta file to
pub fn save_fasta(data: &BTreeMap<String, String>, filename: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let mut fasta_writer = fasta::Writer::to_file(filename)?;
    for (k, v) in data.iter() {
        fasta_writer.write(k, None, v.as_bytes())?;
    }
    Ok(())
}

/// Fasta indexer, mirrored from https://docs.rs/rust-htslib/latest/src/rust_htslib/faidx/mod.rs.html#37-48
/// Ideally, we would call that directly, but it requires an htslib update and that seems to be going poorly...
pub fn index_fasta(filename: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let os_path = std::ffi::CString::new(filename.display().to_string())?;
    let rc = unsafe { rust_htslib::htslib::fai_build(os_path.as_ptr()) };
    if rc < 0 {
        bail!("Error {rc} while building index for {filename:?}");
    }
    Ok(())
}