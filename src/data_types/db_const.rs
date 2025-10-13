
use lazy_static::lazy_static;
use rustc_hash::FxHashSet as HashSet;
use std::collections::BTreeMap;

use crate::data_types::pgx_structural_variants::{FullDeletion, PartialDeletion};

// gene names to prevent dev typos
pub const CYP2A6: &str = "CYP2A6";
pub const CYP2B6: &str = "CYP2B6";
pub const CYP2B7P: &str = "CYP2B7P";
pub const CYP2C18: &str = "CYP2C18";
pub const CYP2C19: &str = "CYP2C19";
pub const CYP2D6: &str = "CYP2D6";
pub const CYP4F2: &str = "CYP4F2";
pub const HELLS: &str = "HELLS";
pub const HLA_A: &str = "HLA-A";
pub const HLA_B: &str = "HLA-B";
pub const NAT2: &str = "NAT2";
pub const SLCO1B1: &str = "SLCO1B1";
pub const TBC1D12: &str = "TBC1D12";

lazy_static!{
    /// List of genes that are ignored from CPIC
    pub static ref CPIC_IGNORED_LIST: Vec<&'static str> = vec![
        CYP2D6, HLA_A, HLA_B, // handled on separate path
        NAT2 // added to CPIC, but not the genes list; we want to use PharmVar here
    ];

    /// Same list, but in set format
    pub static ref CPIC_IGNORED_GENES: HashSet<&'static str> = {
        CPIC_IGNORED_LIST.iter().cloned().collect()
    };

    // hard-coded set of known deletions
    // CYP2B6 reference: https://a.storyblok.com/f/70677/x/97da7caf47/cyp2b6_variation_v1-2.pdf
    //      TODO: *30 should resemble a tandem duplication, which is not currently coded; we may need advanced features to detect that correctly
    //      TODO: duplications with little information mentioned for CYP2B6, maybe a future addition
    // CYP2C19 reference: https://a.storyblok.com/f/70677/x/9225fa2a03/variation_cyp2c19.pdf
    // CYP4F2 reference: https://a.storyblok.com/f/70677/x/bca90e3cbf/cyp4f2_structural-variation_v1-0.pdf
    //      TODO: there are duplications with little information mentioned for CYP4F2, maybe a future addition
    // SLCO1B1 reference: https://a.storyblok.com/f/70677/x/33548d2d14/slco1b1_variation_v1-1.pdf

    /// Contains the list of hard-coded full-length deletions
    pub static ref CPIC_FULL_DELETIONS: BTreeMap<(String, String), FullDeletion> = {
        // array contains (gene name, haplotype label, is_generic, list of deleted genes)
        let full_deletions = [
            // CYP2C19 generic full deletions
            (CYP2C19, "*36", true, vec![CYP2C19]),
            // CYP2C19 specific full deletions
            (CYP2C19, "*36.001", false, vec![CYP2C19, CYP2C18, HELLS]),
            (CYP2C19, "*36.002", false, vec![CYP2C19, CYP2C18, HELLS, TBC1D12]),
            // CYP4F2 just has a full generic deletion
            (CYP4F2, "*16", true, vec![CYP4F2]),
            // technically CYP4F2*16.001 describes a specific deletion event, but we do not have coding for that currently
            // SLCO1B1 just has generic full and partial deletions
            (SLCO1B1, "*48", true, vec![SLCO1B1])
        ];
        full_deletions.into_iter()
            .map(|(gene, hap, is_generic, deleted_genes)| {
                let event = FullDeletion::new(is_generic, deleted_genes.iter().map(|s| s.to_string()).collect());
                ((gene.to_string(), hap.to_string()), event)
            })
            .collect()
    };

    // now the partial deletions
    pub static ref CPIC_PARTIAL_DELETIONS: BTreeMap<(String, String), PartialDeletion> = {
        let partial_deletions = [
            // this is a hybrid between CYP2B6 and CYP2B7 caused by homology in intron 4
            (CYP2B6, "*29", false, vec![
                (CYP2B7P, 4..9),
                (CYP2B6, 0..4)
            ]),
            // CYP2C19 generic partial deletion
            (CYP2C19, "*37", true, vec![
                (CYP2C19, 0..9)
            ]),
            // CYP2C19 more precise partial deletions
            (CYP2C19, "*37.001", false, vec![
                // nothing deleted in C18
                (CYP2C19, 0..5)
            ]),
            (CYP2C19, "*37.002", false, vec![
                (CYP2C18, 7..9),
                (CYP2C19, 0..4)
            ]),
            (CYP2C19, "*37.003", false, vec![
                (CYP2C18, 0..9),
                (CYP2C19, 0..1)
            ]),
            (CYP2C19, "*37.004", false, vec![
                (CYP2C18, 4..9),
                (CYP2C19, 0..7)
            ]),
            (CYP2C19, "*37.005", false, vec![
                (CYP2C18, 1..9),
                (CYP2C19, 0..7)
            ]),
            // SLCO1B1 just has generic full and partial deletions
            (SLCO1B1, "*49", true, vec![
                (SLCO1B1, 0..15)
            ])
        ];
        partial_deletions.into_iter()
            .map(|(gene, hap, is_generic, deleted_exons)| {
                let event = PartialDeletion::new(is_generic, deleted_exons.into_iter().map(|(s, r)| (s.to_string(), r)).collect());
                ((gene.to_string(), hap.to_string()), event)
            })
            .collect()
    };

    /// List of genes that are ignored from PharmVar
    pub static ref PHARMVAR_IGNORED_LIST: Vec<&'static str> = vec![
        CYP2A6, // this one has some D6-like behaviors with SVs; lets ignore it for now
        CYP2D6,
    ];

    /// Same list, but in set format
    pub static ref PHARMVAR_IGNORED_GENES: HashSet<&'static str> = {
        PHARMVAR_IGNORED_LIST.iter().cloned().collect()
    };
}
