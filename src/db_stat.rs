
use std::collections::BTreeMap;

use crate::data_types::alleles::VariantDefinition;
use crate::database::pgx_database::PgxDatabase;

/// Prints the statistics for a given database
/// # Arguments
/// * `database` - the database to print the statistics for
pub fn print_stats(database: &PgxDatabase) {
    // display the database metadata
    let db_metadata = database.database_metadata();
    println!("Database metadata:");
    println!("\tVersion: {}", db_metadata.pbstarphase_version);
    println!("\tCPIC version: {}", db_metadata.cpic_version);
    println!("\tHLA version: {}", db_metadata.hla_version);
    println!("\tPharmVar version: {}", db_metadata.pharmvar_version);
    println!("\tBuild time: {}", db_metadata.build_time);

    // display the database gene statistics
    let gene_entries = database.gene_entries();
    let hla_config = database.hla_config();
    let hla_sequences = database.hla_sequences();
    let cyp2d6_gene_def = database.cyp2d6_gene_def();

    println!("Database gene statistics:");
    let total_genes = gene_entries.len() + hla_config.gene_collection().gene_dict().len() + 1;
    println!("\tTotal genes: {} = {} (variant-based genes) + {} (HLA) + 1 (CYP2D6)", total_genes, gene_entries.len(), hla_config.gene_collection().gene_dict().len());
    println!("\tGene entries:");
    println!("\t\tTotal alleles: {}", gene_entries.values().map(|g| g.defined_haplotypes().len()).sum::<usize>());
    println!("\t\tTotal variants: {}", gene_entries.values().map(|g| g.variants().len()).sum::<usize>());
    println!("\tHLA:");
    println!("\t\tTotal HLA alleles: {}", hla_sequences.len());
    println!("\tCYP2D6:");
    println!("\t\tTotal alleles: {}", cyp2d6_gene_def.len());
    let d6_variants: std::collections::HashSet<&VariantDefinition> = cyp2d6_gene_def.values()
        .flat_map(|gd| gd.variants())
        .collect();
    println!("\t\tTotal variants: {}", d6_variants.len());

    // now do per-gene statistics, but these are just if we have elevated verbosity
    if log::log_enabled!(log::Level::Debug) {
        println!();
        println!("Gene entry statistics:");
        println!("gene\tcore_alleles\tsub_alleles\tcore_variants\tsub_variants\tdata_source");
        for (gene, gene_entry) in gene_entries.iter() {
            let core_alleles = gene_entry.defined_haplotypes().iter()
                .filter(|(_, h)| h.is_core_haplotype())
                .count();
            let sub_alleles = gene_entry.defined_haplotypes().iter()
                .filter(|(_, h)| !h.is_core_haplotype()).count();
            let core_variants = gene_entry.variants().iter()
                .filter(|(_, v)| v.is_core_variant())
                .count();
            let sub_variants = gene_entry.variants().iter()
                .filter(|(_, v)| !v.is_core_variant())
                .count();
            let data_source = gene_entry.data_source();
            println!("{gene}\t{core_alleles}\t{sub_alleles}\t{core_variants}\t{sub_variants}\t{data_source}");
        }
        println!();

        println!("HLA gene statistics:");
        println!("gene\talleles\tdna_sequences");
        let mut allele_counts: BTreeMap<String, usize> = BTreeMap::new();
        let mut dna_counts: BTreeMap<String, usize> = BTreeMap::new();
        for (_allele_name, gene_entry) in hla_sequences.iter() {
            let gene_name = gene_entry.gene_name();
            let dna_sequence = gene_entry.dna_sequence();
            *allele_counts.entry(gene_name.to_string()).or_insert(0) += 1;
            if dna_sequence.is_some() {
                *dna_counts.entry(gene_name.to_string()).or_insert(0) += 1;
            }
        }
        for (gene_name, count) in allele_counts.iter() {
            println!("{}\t{}\t{}", gene_name, count, dna_counts.get(gene_name).unwrap_or(&0));
        }
        println!();
    }
}