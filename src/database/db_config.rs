
use std::collections::BTreeMap;
use serde::{Deserialize, Serialize};

use crate::database::db_const::DPYD;

/// Source for a gene definition
#[derive(Clone, Copy, Deserialize, Debug, Default, Eq, PartialEq, Serialize, strum_macros::Display)]
pub enum PgxDataSource {
    #[default]
    Unknown,
    #[strum(to_string = "CPIC")]
    #[serde(alias = "CPIC")]
    Cpic,
    PharmVar
}

/// Options for building a database
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct DatabaseBuildOptions {
    /// Whether to prefer PharmVar over CPIC
    pub default_gene_source: PgxDataSource,
    /// Override the gene source for a specific gene
    pub gene_source_overrides: BTreeMap<String, PgxDataSource>
}

impl Default for DatabaseBuildOptions {
    fn default() -> Self {
        let mut default_overrides: BTreeMap<String, PgxDataSource> = BTreeMap::new();
        default_overrides.insert(DPYD.to_string(), PgxDataSource::Cpic);
        Self {
            default_gene_source: PgxDataSource::PharmVar,
            gene_source_overrides: default_overrides
        }
    }
}
