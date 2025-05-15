
use bio::io::fasta;
use itertools::Itertools;
use log::{debug, info, trace, warn};
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use rustc_hash::FxHashMap as HashMap;
use serde::Deserialize;
use std::collections::BTreeMap;
use std::collections::hash_map::Entry::{Occupied, Vacant};
use std::io::Read;
use simple_error::bail;

use crate::data_types::alleles::{AlleleDefinition, VariantDefinition};
use crate::data_types::cpic_api_results::CpicAlleleDefinition;
use crate::data_types::database::PgxDatabase;
use crate::data_types::db_const::PHARMVAR_IGNORED_GENES;
use crate::data_types::pharmvar_api_results::{PharmvarGeneDefinition, PharmvarAlleleDefinition};
use crate::hla::alleles::{HlaAlleleDefinition, SUPPORTED_HLA_GENES};

// CPIC API quickstart: https://github.com/cpicpgx/cpic-data/wiki
// CPI API full book: https://documenter.getpostman.com/view/1446428/Szt78VUJ?version=latest
// Useful postgrest reference: https://postgrest.org/en/v7.0.0/api.html#horizontal-filtering-rows

/// Base API addresses
const CPIC_API_URL: &str = "https://api.cpicpgx.org/v1";

// fortunately, HLA is version controlled on GitHub, makes life a little better for versioning
const HLA_REPO_LOOKUP: &str = "https://api.github.com/repos/ANHIG/IMGTHLA/releases/latest";
const HLA_GITHUB_PREFIX: &str = "https://raw.githubusercontent.com/ANHIG/IMGTHLA";
const HLA_GENOME_FASTA: &str = "fasta/hla_gen.fasta";
const HLA_GENOME_FASTA_ZIP: &str = "fasta/hla_gen.fasta.zip"; // started with v3.57.0-alpha
const HLA_CDNA_FASTA: &str = "fasta/hla_nuc.fasta";

// PharmVar API: https://www.pharmvar.org/documentation
// PharmVar gene information with useful download: https://www.pharmvar.org/gene/CYP2D6
//   from there, you can get this zip file: https://www.pharmvar.org/get-download-file?name={gene}&refSeq=ALL&fileType=zip&version={version}
//   version can be "current" or numbered like "6.0.8"
// PharmVar API link: https://www.pharmvar.org/api-service/alleles?exclude-sub-alleles=false&include-reference-variants=false&include-retired-alleles=false&include-retired-reference-sequences=false&reference-sequence=NC_000022.11
const PHARMVAR_API_URL: &str = "https://www.pharmvar.org/api-service";

/// This is the primary call to build out our database locally via multiple API queries.
/// # Arguments
/// * `reference_genome` - required now for checking the HLA alignments automatically during DB construction
/// # Errors
/// * if there are errors retrieving the CPIC gene list
/// * if there are errors retrieving allele definitions for a gene
pub fn build_database_via_api(reference_genome: &ReferenceGenome) -> Result<PgxDatabase, Box<dyn std::error::Error>> {
    // first get all the CPI genes
    info!("Starting CPIC API queries...");
    let all_genes: HashMap<String, String> = get_all_cpic_genes()?;
    info!("\tCPIC gene list: {:?}", all_genes.keys().sorted().collect::<Vec<_>>());

    // If testing, you can limit this to a particular gene
    let query_limit = None; //Some("CACNA1S");

    // get the alleles for the gene
    let cpic_alleles: Vec<CpicAlleleDefinition> = query_gene_cpic_api(query_limit)?;
    info!("CPIC API queries complete.");

    // now handle the PharmVar genes
    info!("Starting PharmVar gene queries...");
    let pharmvar_genes = get_all_pharmvar_genes()?;
    debug!("\tFull PharmVar gene list: {pharmvar_genes:?}");
    let filtered_pharmvar_genes: Vec<String> = pharmvar_genes.into_iter()
        .filter(|g| !all_genes.contains_key(g) && !PHARMVAR_IGNORED_GENES.contains(g.as_str())) // remove anything from CPIC and CYP2D6 (handled separate)
        .sorted()
        .collect();
    info!("\tFiltered PharmVar gene list: {filtered_pharmvar_genes:?}");

    // now get all the PharmVar alleles, which will be missing the all REF alleles
    let pharmvar_alleles = query_gene_pharmvar_api(&filtered_pharmvar_genes)?;
    info!("Found {} PharmVar alleles via API.", pharmvar_alleles.len());

    // now we need to pull down HLA data as well
    info!("Starting HLA queries...");
    let latest_hla_version: String = get_latest_hla_tag()?;
    info!("Found latest HLA version: {latest_hla_version}");
    let hla_data: BTreeMap<String, HlaAlleleDefinition> = get_hla_sequences(&latest_hla_version)?;

    // finally, get the PharmVar CYP2D6 data
    let (pharmvar_version, cyp2d6_data) = get_pharmvar_variants("CYP2D6", "current")?;
    info!("Found latest PharmVar version: {pharmvar_version}");

    // now build our database and ship it back
    // let opt_refseq_fn = Some(std::path::PathBuf::from("./GRCh38_latest_genomic.gff.gz"));
    let full_database: PgxDatabase = PgxDatabase::new(
        &all_genes, 
        &cpic_alleles,
        &pharmvar_alleles,
        latest_hla_version,
        hla_data,
        pharmvar_version,
        cyp2d6_data,
        reference_genome,
        None
    )?;

    // todo!("remove query filters and remove opt_refseq_fn");

    Ok(full_database)
}

/// This pulls the list of genes that are available from CPIC and stores useful metadata like chromosome
/// # Errors
/// * if the URL request has issues connecting or converting to JSON
/// * if duplicate gene names are detected while parsing
/// * if required entries are missing or fail to parse
fn get_all_cpic_genes() -> Result<HashMap<String, String>, Box<dyn std::error::Error>> {
    // this endpoint gets the list of genes, ordered by symbol, where the URL field is not empty
    //   this tends to correlate with genes that have allele definitions
    //   if we ever find that does not hold, we can remove the url= filter component and just accept extra queries downstream
    let gene_url: String = format!("{CPIC_API_URL}/gene?url=not.eq.null&order=symbol");
    info!("\tQuerying CPIC gene list via {gene_url}");

    // hit the end point so we can parse it
    let result: String = reqwest::blocking::get(gene_url)?.text()?;
    debug!("Response received.");
    
    // now parse it via serde
    // we are using a generic Value here because we really just need one field right now
    let parsed: Vec<serde_json::Value> = serde_json::from_str(&result)?;
    debug!("Parsing complete.");
    
    // now pull out the chromosome for each gene we care about
    let mut ret: HashMap<String, String> = Default::default();
    for gene_entry in parsed.iter() {
        // make sure we get a gene name and a chromosome
        let gene_name: String = match gene_entry["symbol"].as_str() {
            Some(s) => s.to_string(),
            None => bail!("Error while parsing field \"symbol\" as a string")
        };
        let chromosome: String = match gene_entry["chr"].as_str() {
            Some(s) => s.to_string(),
            None => {
                warn!("Error while parsing field \"chr\" for {gene_name}, ignoring: {:?}", gene_entry["chr"]);
                continue;
            }
        };

        // the clippy warning here is far less readable IMO, disabling it
        #[allow(clippy::map_entry)]
        if ret.contains_key(&gene_name) {
            bail!("Detected duplicate gene name during parsing: {gene_name}");
        } else {
            debug!("\t\t{gene_name} -> {chromosome}");
            ret.insert(gene_name, chromosome);
        }
    }
    
    Ok(ret)
}

/// This will pull all the CPIC allele definitions via a single API query.
/// # Arguments
/// * `gene` - the gene name to query; if None, then all definitions are pulled
/// # Errors
/// * if the URL fails to get
/// * if the response fails to parse into JSON or our allele definition
fn query_gene_cpic_api(gene: Option<&str>) -> Result<Vec<CpicAlleleDefinition>, Box<dyn std::error::Error>> {
    // this will pull allele definitions as well as the variants that go with them
    // let definition_url: String = format!("{CPIC_API_URL}/allele_definition?genesymbol=eq.{gene}&select=*,%20allele_location_value(*,%20sequence_location(*))&order=name");
    let definition_url: String = match gene {
        Some(g) => format!("{CPIC_API_URL}/allele_definition?genesymbol=eq.{g}&select=*,%20allele_location_value(*,%20sequence_location(*))&order=name"),
        None => format!("{CPIC_API_URL}/allele_definition?select=*,%20allele_location_value(*,%20sequence_location(*))&order=name")
    };
    let label = gene.unwrap_or("all_genes");
    info!("\tQuerying \"{label}\" via {definition_url}");

    // hit the end point so we can parse it
    let result: String = reqwest::blocking::get(definition_url)?.text()?;
    debug!("Response received.");
    
    // now parse it via serde
    let parsed: Vec<CpicAlleleDefinition> = serde_json::from_str(&result)?;
    debug!("Parsing complete.");
    Ok(parsed)
}

/// Gets the latest version tag of the HLA sequences
/// # Errors
/// * if the URL request fails
/// * if "tag_name" is not present in the response or cannot be converted into a string
fn get_latest_hla_tag() -> Result<String, Box<dyn std::error::Error>> {
    // we need the User Agent specified for GitHub queries, set it to our tool name
    let client = reqwest::blocking::Client::builder()
        .user_agent(env!("CARGO_PKG_NAME"))
        .build()?;

    // hit the end point so we can parse it
    info!("\tQuerying lastest HLA tag via {HLA_REPO_LOOKUP}");
    let result: String = client.get(HLA_REPO_LOOKUP)
        .send()?
        .error_for_status()?
        .text()?;
    debug!("Response received.");
    
    // now parse it via serde
    let parsed: serde_json::Value = serde_json::from_str(&result)?;
    match parsed.get("tag_name") {
        Some(v) => {
            match v.as_str() {
                Some(str_form) => {
                    Ok(str_form.to_string())
                },
                None => {
                    bail!("Key \"tag_name\" could not be converted to a String.");
                }
            }
        },
        None => {
            bail!("Key \"tag_name\" was not found in GitHub latest response for HLA repository.");
        }
    }
}

/// This will pull the files for the corresponding version string and get them ready
/// # Arguments
/// * `version` - the version we are pulling out of GitHub, must be a valid tag
/// # Errors
/// if the URLs requested have an error
/// if we cannot convert the response into FASTA sequences
/// if we get an error while collapsing the DNA and cDNA entries together
fn get_hla_sequences(version: &str) -> Result<BTreeMap<String, HlaAlleleDefinition>, Box<dyn std::error::Error>> {
    // hit the DNA end point
    let dna_url: String = format!("{HLA_GITHUB_PREFIX}/{version}/{HLA_GENOME_FASTA_ZIP}");
    info!("\tQuerying HLA DNA sequences via {dna_url}");
    let client = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(300))
        .build()?;
    let dna_result: String = match client.get(dna_url).send()?.error_for_status() {
        Ok(r) => {
            let result = r.bytes()?;
            
            // convert into a cursor for the archive and open it
            let cursor_result = std::io::Cursor::new(result.to_vec());
            let mut archive = zip::ZipArchive::new(cursor_result)?;

            // get the exact file we want and read into the string
            let mut zip_file = archive.by_name("hla_gen.fasta")?;
            let mut text_form: String = Default::default();
            zip_file.read_to_string(&mut text_form)?;
            text_form
        },
        Err(e) => {
            debug!("\tFailed to find zipped HLA fasta, error: {e}");
            let dna_url_unzip = format!("{HLA_GITHUB_PREFIX}/{version}/{HLA_GENOME_FASTA}");
            info!("\tQuerying HLA DNA sequences via backup URL: {dna_url_unzip}");
            let unzip_result = client.get(dna_url_unzip).send()?.error_for_status()?;
            unzip_result.text()?
        }
    };
    
    debug!("Response received.");
    let dna_data: HashMap<String, (String, String)> = convert_fasta_str_to_map(&dna_result, false)?;
    debug!("Parsing complete.");

    // now hit the cDNA (e.g., exon) end point
    let cdna_url: String = format!("{HLA_GITHUB_PREFIX}/{version}/{HLA_CDNA_FASTA}");
    info!("\tQuerying HLA cDNA sequences via {cdna_url}");
    let response = client.get(cdna_url).send()?.error_for_status()?;
    let cdna_result: String = response.text()?;
    debug!("Response received.");
    let cdna_data: HashMap<String, (String, String)> = convert_fasta_str_to_map(&cdna_result, false)?;
    debug!("Parsing complete.");

    let collapsed_lookup = collapse_hla_lookup(dna_data, cdna_data)?;
    Ok(collapsed_lookup)
}

/// This will take a loaded FASTA file as a string and convert it into a HashMap, performing some minor checks.
/// The key of this HashMap is an HLA identifer and the values is the (star-allele, sequence).
/// # Arguments
/// * `raw_fasta` - the raw FASTA sequence to convert
/// * `reversed_ids` - if True, then we need to swap the "id" and "star_allele" in the lookups
/// # Errors
/// * if the FASTA is invalid format
/// * if the sequence cannot be converted from UTF-8
/// * if a duplicate entry is detected in the FASTA
pub fn convert_fasta_str_to_map(raw_fasta: &str, reversed_ids: bool) -> Result<HashMap<String, (String, String)>, Box<dyn std::error::Error>> {
    let mut ret: HashMap<String, (String, String)> = Default::default();
    let reader = fasta::Reader::new(raw_fasta.as_bytes());
    for result in reader.records() {
        let record = result?;
        let id: String = record.id().to_string();
        let seq: String = String::from_utf8(record.seq().to_vec())?;

        // this converts from desc() = Some("A*01:01:01:01 3503 bp") -> "A*01:01:01:01"
        let star_allele: String = record.desc().unwrap_or_default()
            .split_whitespace().next().unwrap_or_default()
            .to_string();
        
        // handle a swap if necessary
        let (id, star_allele) = if reversed_ids {
            (star_allele, id)
        } else {
            (id, star_allele)
        };

        trace!("Found record {}({}) with length {}", id, star_allele, seq.len());
        
        // apparently there are some duplicates in the FASTA file for some reason
        match ret.entry(id.clone()) {
            Occupied(entry) => {
                // make sure the entry is just a duplicate
                if entry.get() != &(star_allele, seq) {
                    bail!("FASTA record with multiple IDs/sequences detected: {id}");
                }
            },
            Vacant(entry) => {
                // normal path, we insert the new entry
                entry.insert((star_allele, seq));
            }
        };
    }
    Ok(ret)
}

/// Helper function that collapse the DNA and cDNA entries into a single HlaAlleleDefinition for each allele.
/// # Arguments
/// * `dna_data` - a map from HLA ID to (star allele ID, DNA sequence)
/// * `cdna_data` - a map from HLA ID to (star allele ID, cDNA sequence)
/// # Errors
/// * if the key sets for each HashMap do not match
/// * if the star alleles for an HLA ID are different in the two maps
/// * if there is an error parsing an HLA allele definition
pub fn collapse_hla_lookup(dna_data: HashMap<String, (String, String)>, cdna_data: HashMap<String, (String, String)>) 
    -> Result<BTreeMap<String, HlaAlleleDefinition>, Box<dyn std::error::Error>> {
    // every DNA has a cDNA, but not all cDNAs have a DNA; count them before we drain the `cdna_data`
    let mut missed_dna: usize = 0;
    for hla_id in dna_data.keys() {
        if !cdna_data.contains_key(hla_id) {
            missed_dna += 1;
        }
    }
    if missed_dna > 0 {
        warn!("Detected {missed_dna} DNA entries that do not have a cDNA, ignoring them.");
    }

    // now we can build up all the usable entries
    let mut ret: BTreeMap<String, HlaAlleleDefinition> = Default::default();
    let mut ignored_alleles = 0;
    for (hla_id, (cdna_desc, cdna_seq)) in cdna_data.into_iter() {
        // future note: this is similar to a .map() function but we do it this way so we can propagate the bail! correctly
        let opt_dna_seq = match dna_data.get(&hla_id) {
            Some((dna_desc, dna_seq)) => {
                if dna_desc != &cdna_desc {
                    bail!("{hla_id} has description \"{dna_desc}\" for DNA and \"{cdna_desc}\" for cDNA.");
                };
                Some(dna_seq.clone())
            },
            None => None
        };

        // checks out so far, make the allele and insert it
        let new_allele = HlaAlleleDefinition::new(
            hla_id.clone(), &cdna_desc, opt_dna_seq, cdna_seq
        )?;

        // restrict our database to the alleles we plan to match, we can relax this if we ever do an update to HLA
        // if new_allele.gene_name() == "HLA-A" || new_allele.gene_name() == "HLA-B" {
        if SUPPORTED_HLA_GENES.contains(new_allele.gene_name()) {
            ret.insert(hla_id, new_allele);
        } else {
            ignored_alleles += 1;
        }
    }
    let supported_vec: Vec<&String> = SUPPORTED_HLA_GENES.iter().collect(); // purely for print statement
    debug!("Removed {ignored_alleles} alleles that are not in supported HLA gene set: {supported_vec:?}");

    Ok(ret)
}

/// This pulls the list of genes that are available from PharmVar.
/// # Errors
/// * if the URL request has issues connecting or converting to JSON
fn get_all_pharmvar_genes() -> Result<Vec<String>, Box<dyn std::error::Error>> {
    // this endpoint gets the list of genes, ordered by symbol, where the URL field is not empty
    //   this tends to correlate with genes that have allele definitions
    //   if we ever find that does not hold, we can remove the url= filter component and just accept extra queries downstream
    let gene_url: String = format!("{PHARMVAR_API_URL}/genes/list");
    info!("\tQuerying PharmVar gene list via {gene_url}");

    // hit the end point so we can parse it
    let result: String = reqwest::blocking::get(gene_url)?.text()?;
    debug!("Response received.");

    // now parse it via serde
    // we are using a generic Value here because we really just need one field right now
    let parsed: Vec<String> = serde_json::from_str(&result)?;
    debug!("Parsing complete.");

    Ok(parsed)
}

/// This will pull all the PharmVar allele definitions for the genes in the list.
/// # Arguments
/// * `gene_list` - the list of genes to query
/// # Errors
/// * if the URL fails to get
/// * if the response fails to parse into JSON or our allele definition
fn query_gene_pharmvar_api(gene_list: &[String]) -> Result<Vec<PharmvarAlleleDefinition>, Box<dyn std::error::Error>> {
    let mut ret = vec![];
    for gene in gene_list.iter() {
        // example URL: https://www.pharmvar.org/api-service/genes/NAT2?exclude-sub-alleles=false&include-reference-variants=false&include-retired-alleles=false&include-retired-reference-sequences=false&reference-collection=GRCh38
        // TODO: if we do not specify the gene, it pulls the whole DB in one query; is that better?
        // TODO: we have exclude-sub-alleles=true, which excludes sub-alleles
        let definition_url = format!("{PHARMVAR_API_URL}/genes/{gene}?exclude-sub-alleles=true&include-reference-variants=false&include-retired-alleles=false&include-retired-reference-sequences=false&reference-collection=GRCh38");
        info!("\tQuerying \"{gene}\" via {definition_url}");

        // hit the end point so we can parse it
        let result: String = reqwest::blocking::get(definition_url)?.text()?;
        debug!("Response received.");

        // now parse it via serde
        let parsed: PharmvarGeneDefinition = serde_json::from_str(&result)?;
        debug!("Parsing complete.");

        ret.extend(parsed.alleles.iter().cloned());
    }
    Ok(ret)
}

/// Gets a PharmVar zip file from the website and returns a tuple with the version and the PgxGene definition
/// # Arguments
/// * `gene` - the gene to retrieve
/// * `version` - the version to fetch; can be "current" to get the latest
/// # Errors
/// * if the URL fetch fails
/// * if parsing the ZIP archive fails
/// * if the file containing fasta sequence cannot be found
fn get_pharmvar_variants(gene: &str, version: &str) -> Result<(String, BTreeMap<String, AlleleDefinition>), Box<dyn std::error::Error>> {
    // the URL is pretty standard
    let gene_url = format!("https://www.pharmvar.org/get-download-file?name={gene}&refSeq=ALL&fileType=zip&version={version}");
    info!("Querying PharmVar({gene}, {version}) via {gene_url}");

    // hit the end point so we can parse it; returns as bytes::Bytes
    let client = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(300))
        .build()?;
    let response = client.get(gene_url).send()?.error_for_status()?;
    let result = response.bytes()?;
    debug!("Response received.");

    // convert into a cursor for the archive
    let cursor_result = std::io::Cursor::new(result.to_vec());
    // now open up the archive
    let mut archive = zip::ZipArchive::new(cursor_result)?;
    let mut version: Option<String> = None;
    let coordinate_version = "GRCh38";

    // time to populate the alleles
    let mut ret: BTreeMap<String, AlleleDefinition> = Default::default();
    for i in 0..archive.len() {
        let mut file = archive.by_index(i)?;
        let outpath = match file.enclosed_name() {
            Some(path) => path.to_owned(),
            None => continue,
        };
        if outpath.is_dir() {
            // ignore directories
            continue;
        }

        // make sure version is consistent
        let parent_split: Vec<&str> = outpath.components()
            .map(|c| c.as_os_str().to_str().unwrap())
            .collect();

        // we only care about files that are two folders deep
        if parent_split.len() != 3 {
            continue;
        }

        // root folder should be "{gene}-{version}"
        let root_folder = parent_split[0];
        assert!(root_folder.starts_with(gene));
        let file_version = root_folder[gene.len()+1..].to_string();
        match version.as_ref() {
            Some(v) => {
                if v != &file_version {
                    bail!("Found mismatched versions in ZIP file: {v}, {file_version}");
                }
            },
            None => {
                trace!("file_version={file_version}");
                version = Some(file_version);
            }
        };

        let filename = outpath.file_name().unwrap().to_str().unwrap();
        
        // structure is "{gene}-{version}/{coordinate_version}/{gene}_{allele}.vcf"
        if parent_split[1] == coordinate_version && filename.ends_with(".vcf") {
            // pull out the actual allele name from the filename
            assert_eq!(&filename[..gene.len()], gene);
            let allele = filename[gene.len()+1..filename.len()-4].to_string();

            // now we need to actually parse the contained VCF file
            let mut vcf_content: String = Default::default();
            let _bytes = file.read_to_string(&mut vcf_content)?;
            let start_index = vcf_content.find("#CHROM").unwrap();
            
            // now put this into a VCF reader for ease of use
            let variants = load_vcf_from_bytes(vcf_content[start_index..].as_bytes())?;
            
            // create an allele definition using the full star ID; we do not get the special identifiers with this approach
            let full_star = format!("{gene}*{allele}");
            let allele_def = AlleleDefinition::new(None, &full_star, variants)?;

            // insert it using the ID; make sure the previous value was empty
            assert!(ret.insert(allele_def.id().to_string(), allele_def).is_none());
        }
    }

    // we need to add in the reference key here at the end
    let ref_allele = format!("{gene}*1.001");
    let allele_def = AlleleDefinition::new(None, &ref_allele, vec![])?;
    assert!(ret.insert(allele_def.id().to_string(), allele_def).is_none());

    if let Some(v) = version {
        Ok((v, ret))
    } else {
        bail!("No files or version identified in ZIP file");
    }
}

/// This assists in parsing the CSV reader of the VCF
#[derive(Debug, Deserialize)]
struct VcfRow {
    #[serde(alias = "#CHROM")]
    pub chrom: String,
    #[serde(alias = "POS")]
    pub pos: usize,
    #[serde(alias = "ID")]
    pub id: String,
    #[serde(alias = "REF")]
    pub reference: String,
    #[serde(alias = "ALT")]
    pub alternate: String,
    // #[serde(alias = "QUAL")]
    // quality: String,
    // #[serde(alias = "FILTER")]
    // filter: String,
    #[serde(alias = "INFO")]
    pub info: String
}

/// This is a wrapper function to the load the VCF files for CYP2D6 from memory (e.g., from the ZIP in memory).
/// In theory, you wouldn't normally need this, but rust_htslib does not have a function to parse VCFs from bytes in memory.
/// # Arguments
/// * `vcf_content` - the bytes that make up the VCF starting at "#CHROM" (i.e., the header is skipped)
/// # Errors
/// * if the deserializing throws errors
/// # Panics
/// * we baked in some specific assumptions about VCF content and those currently panic if they fail
fn load_vcf_from_bytes(vcf_content: &[u8]) -> Result<Vec<VariantDefinition>, Box<dyn std::error::Error>> {
    /*
    rust_htslib does not have a way to read from memory, I guess we can do one of these:
    1. write the files out and read them as normal - I'd like to avoid creating temp files if possible
    2. read them via noodles_bcf which can load from anything implementing Read - turns out noodles does not like these file formats
    3. parse it out ourselves via csv - "easiest" for now, could have scaling issues if VCFs get complicated
    */
    // open the tab-delimited reader
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(vcf_content);
    let mut ret = vec![];

    // go through each row and build a variant entry from it
    for row in reader.deserialize() {
        let record: VcfRow = row?;
        let mut extras: BTreeMap<String, String> = Default::default();
        for key_value in record.info.split(';') {
            if !key_value.is_empty() && key_value != "." {
                // get the key value extra here
                let kv_split: Vec<&str> = key_value.split('=').collect();
                assert_eq!(kv_split.len(), 2);

                // insert and make sure we do not have that key already
                let k = kv_split[0].to_string();
                let v = kv_split[1].to_string();
                assert!(extras.insert(k, v).is_none());
            }
        }

        // set the record ID if it's not empty
        let record_id = if record.id != "." {
            Some(record.id.clone())
        } else {
            None
        };

        // now make the variant and save it
        let var_def = VariantDefinition::new(
            record_id, record.chrom.clone(), record.pos - 1, //POS is 1-based inside the VCF file, so shift it to 0-based
            record.reference.clone(), record.alternate.clone(),
            extras
        )?;
        ret.push(var_def);
    }

    Ok(ret)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_get_all_cpic_genes() {
        let all_genes: HashMap<String, String> = get_all_cpic_genes().unwrap();

        // this list can change, but presumably things will not get removed
        assert!(all_genes.len() >= 24);

        // check a gene
        assert_eq!(all_genes.get("CACNA1S").unwrap(), "chr1");
    }

    #[test]
    fn test_get_latest_hla_tag() {
        let latest_tag = get_latest_hla_tag();
        assert!(latest_tag.is_ok());
    }

    #[test]
    fn test_get_hla_sequences() {
        // note to future self: this particular version appear to be malformed; most others would be v3.54.0-alpha
        let fixed_version = "v3.54.0-alpha";
        let hla_db = get_hla_sequences(fixed_version).unwrap();
        // assert_eq!(hla_db.len(), 38408); // this was before we restricted to just A and B
        // assert_eq!(hla_db.len(), 17585); // this was with just A and B
        assert_eq!(hla_db.len(), 36287); // v1.1.0 gene list
        
        let first_entry = hla_db.get("HLA:HLA00001").unwrap();
        assert_eq!(first_entry.hla_id(), "HLA:HLA00001");
        assert_eq!(first_entry.gene_name(), "HLA-A");
        assert_eq!(first_entry.star_allele(), vec!["01".to_string(); 4]);
        assert_eq!(first_entry.dna_sequence().unwrap().len(), 3503);
        assert_eq!(first_entry.cdna_sequence().len(), 1098);
    }

    #[test]
    fn test_zip_hla_sequences() {
        // note to future self: this particular version appear to be malformed; most others would be v3.54.0-alpha
        let fixed_version = "v3.57.0-alpha";
        let hla_db = get_hla_sequences(fixed_version).unwrap();
        // assert_eq!(hla_db.len(), 38408); // this was before we restricted to just A and B
        // assert_eq!(hla_db.len(), 18461); // this was with just A and B
        assert_eq!(hla_db.len(), 38316); // v1.1.0 gene list
        
        let first_entry = hla_db.get("HLA:HLA00001").unwrap();
        assert_eq!(first_entry.hla_id(), "HLA:HLA00001");
        assert_eq!(first_entry.gene_name(), "HLA-A");
        assert_eq!(first_entry.star_allele(), vec!["01".to_string(); 4]);
        assert_eq!(first_entry.dna_sequence().unwrap().len(), 3503);
        assert_eq!(first_entry.cdna_sequence().len(), 1098);
    }

    #[test]
    fn test_get_pharmvar_sequences() {
        let fixed_version = "6.0.8";
        let (version, cyp2d6_db) = get_pharmvar_variants("CYP2D6", fixed_version).unwrap();
        assert_eq!(&version, fixed_version);
        assert_eq!(cyp2d6_db.len(), 510); // we got 511 when we did the FASTA based, may need to resolve that in the future

        // make sure we have the entry that is missing a VCF due to no variants
        let first_entry = cyp2d6_db.get("CYP2D6*1.001").unwrap();
        assert_eq!(first_entry.id(), "CYP2D6*1.001");
        assert_eq!(first_entry.gene_name(), "CYP2D6");
        assert_eq!(first_entry.star_allele(), "1.001");
        assert_eq!(first_entry.variants().len(), 0);

        // check *2 also
        let second_entry = cyp2d6_db.get("CYP2D6*2").unwrap();
        assert_eq!(second_entry.id(), "CYP2D6*2");
        assert_eq!(second_entry.gene_name(), "CYP2D6");
        assert_eq!(second_entry.star_allele(), "2");
        assert_eq!(second_entry.variants().len(), 2);
    }
}