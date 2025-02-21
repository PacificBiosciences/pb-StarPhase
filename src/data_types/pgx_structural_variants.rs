
use std::collections::{BTreeMap, BTreeSet};
use std::collections::btree_map::Entry;

use serde::{Deserialize, Serialize};
use simple_error::{SimpleError, bail};

/// Structure containing all known PGx SVs for a particular gene
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct PgxStructuralVariants {
    /// Map from a full gene deletion ID to the definition; at most one "generic" full deletion in the set
    full_gene_deletions: BTreeMap<String, FullDeletion>,
    /// Map from a full gene deletion ID to the definition; at most one "generic" partial deletion in the set
    partial_gene_deletions: BTreeMap<String, PartialDeletion>
}

impl PgxStructuralVariants {
    /// Adds a full deletion event to the alleles for this gene
    /// * `label` - the label (i.e. genotype) for this event
    /// * `event` - the full deletion to add
    /// # Errors
    /// * if a duplicate label is provided
    pub fn add_full_deletion(&mut self, label: String, event: FullDeletion) -> Result<(), SimpleError> {
        // make sure this definition is not in the partial deletions
        if self.partial_gene_deletions.contains_key(&label) {
            bail!("Full deletion label is already in partial gene deletion definitions: {label}");
        }

        // now check the full deletions
        match self.full_gene_deletions.entry(label) {
            Entry::Vacant(vacant_entry) => {
                vacant_entry.insert(event);
            },
            Entry::Occupied(occupied_entry) => {
                bail!("Duplicate full deletion detected: {}", occupied_entry.key());
            },
        };
        Ok(())
    }

    /// Adds a partial deletion event to the alleles for this gene
    /// * `label` - the label (i.e. genotype) for this event
    /// * `event` - the full deletion to add
    /// # Errors
    /// * if a duplicate label is provided
    pub fn add_partial_deletion(&mut self, label: String, event: PartialDeletion) -> Result<(), SimpleError> {
        // make sure this definition is not in the full deletions
        if self.full_gene_deletions.contains_key(&label) {
            bail!("Partial deletion label is already in full gene deletion definitions: {label}");
        }

        // now check the full deletions
        match self.partial_gene_deletions.entry(label) {
            Entry::Vacant(vacant_entry) => {
                vacant_entry.insert(event);
            },
            Entry::Occupied(occupied_entry) => {
                bail!("Duplicate partial deletion detected: {}", occupied_entry.key());
            },
        };
        Ok(())
    }

    /// This will return a list of all genes that are potentially impacted by the defined structural variants
    pub fn impacted_gene_set(&self) -> BTreeSet<String> {
        let mut ret: BTreeSet<String> = Default::default();
        for fd in self.full_gene_deletions.values() {
            ret.extend(fd.full_genes_deleted().iter().cloned());
        }

        for pd in self.partial_gene_deletions.values() {
            ret.extend(pd.exons_deleted().keys().cloned());
        }

        ret
    }

    // getters
    pub fn full_gene_deletions(&self) -> &BTreeMap<String, FullDeletion> {
        &self.full_gene_deletions
    }

    pub fn partial_gene_deletions(&self) -> &BTreeMap<String, PartialDeletion> {
        &self.partial_gene_deletions
    }
}

/// Describes a full deletion structural variant
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct FullDeletion {
    /// If True, this is a generic full deletion, which is generally a core allele in PGx land
    is_generic: bool,
    /// The set of genes that must be deleted to match; if `is_generic`, then this should just be the primary gene
    full_genes_deleted: BTreeSet<String>
}

impl FullDeletion {
    /// Constructor
    pub fn new(is_generic: bool, full_genes_deleted: BTreeSet<String>) -> Self {
        Self {
            is_generic,
            full_genes_deleted
        }
    }

    // getters
    pub fn is_generic(&self) -> bool {
        self.is_generic
    }

    pub fn full_genes_deleted(&self) -> &BTreeSet<String> {
        &self.full_genes_deleted
    }
}

/// Describes a partial deletion structural variant
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct PartialDeletion {
    /// If True, this is a generic partial deletion, which is generally a core allele in PGx land
    is_generic: bool,
    /// A list of genes and the exons that are deleted; this exon list is always relative to the reference orientation
    exons_deleted: BTreeMap<String, std::ops::Range<usize>>
}

impl PartialDeletion {
    /// Constructor
    pub fn new(is_generic: bool, exons_deleted: BTreeMap<String, std::ops::Range<usize>>) -> Self {
        Self {
            is_generic,
            exons_deleted
        }
    }

    // getters
    pub fn is_generic(&self) -> bool {
        self.is_generic
    }

    pub fn exons_deleted(&self) -> &BTreeMap<String, std::ops::Range<usize>> {
        &self.exons_deleted
    }
}