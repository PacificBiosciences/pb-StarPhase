
use lazy_static::lazy_static;
use regex::Regex;
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use rustc_hash::FxHashSet as HashSet;
use serde::Serialize;
use simple_error::{SimpleError, bail};

lazy_static! {
    /// This matches tandem-repeat like patterns: AC(8) OR ACGTAGT(3).
    /// "seq" matches the bases and "count" matches the contained repeat value.
    pub static ref TR_REGEX: Regex = Regex::new(r"^(?<seq>[A-Z]+)\((?<count>[0-9]+)\)$").unwrap();
}

/// A normalized variant is unambiguously defined and left-aligned to the reference genome
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct NormalizedVariant {
    /// chromosome of the variant
    chrom: String,
    /// 0-based position of the variant
    position: usize,
    /// ref allele
    reference: String,
    /// alt allele
    alternate: String,
    /// contains relevant SV statistics if this is an SV variant
    sv_stats: Option<StructuralVariantStats>
}

impl NormalizedVariant {
    /// This will take an allele definition and normalize it using the reference genome
    /// # Arguments
    /// * `chrom` - the chromosome the variant is positioned on
    /// * `position` - the 0-based coordinate of the variant
    /// * `ref_allele` - the given reference sequence
    /// * `alt_allele` - the given alternate sequence
    /// * `reference_genome` - the pre-loaded reference genome; if None, some normalization steps will not happen and it may fail to normalize
    /// # Errors
    /// * if the reference allele is empty
    /// * if a variant is not yet supported
    /// * if the ref_allele does not match the given reference
    /// * if there is unexpected sequence
    pub fn new(chrom: String, position: usize, ref_allele: &str, alt_allele: &str, reference_genome: Option<&ReferenceGenome>) -> Result<NormalizedVariant, Box<dyn std::error::Error>> {
        // I think CPIC is smart enough to always put "del", but best to double check
        if ref_allele.is_empty() {
            bail!("ref_allele cannot be empty");
        }

        // this shouldn't happen from what I've seen
        if ref_allele == "del" && !alt_allele.starts_with("ins") {
            bail!("Unexpected non-ins alt sequence with a del reference");
        }
        
        // these may change as we normalize the alleles, put into byte vecs for manipulation
        let mut position: usize = position;
        let mut ref_allele: Vec<u8> = parse_sequence(ref_allele);
        let mut alt_allele: Vec<u8> = parse_sequence(alt_allele);
        
        // check if we were given a reference genome to normalize with
        let chrom_seq = match reference_genome {
            Some(rg) => {
                // TODO: this will panic if the chromosome is not in the reference genome
                if rg.contig_keys().contains(&chrom) {
                    // we have the contig, make sure our reference allele matches
                    let cs = rg.get_full_chromosome(&chrom);
                    let rg_seq = &cs[position..position+ref_allele.len()];
                    if ref_allele != rg_seq {
                        let ref_allele = std::str::from_utf8(&ref_allele).unwrap_or("utf8-error");
                        let rg_seq = std::str::from_utf8(rg_seq).unwrap_or("utf8-error");
                        bail!("At {chrom}:{position}, provided reference allele has {ref_allele:?} but reference genome has {rg_seq:?}");
                    }
                    Some(cs)
                } else {
                    bail!("Reference genome does not contain contig {chrom:?}");
                }
            },
            None => None
        };

        // make sure that we have some sequences
        if ref_allele.is_empty() && alt_allele.is_empty() {
            bail!("ref_allele and alt_allele cannot both be empty");
        } else if ref_allele.is_empty() {
            // we are inserting _after_ this position (see https://www.ncbi.nlm.nih.gov/snp/rs777311140 for example)
            // position does not change, we just need to pre-pend this position
            if let Some(cs) = chrom_seq {
                ref_allele.insert(0, cs[position]);
                alt_allele.insert(0, cs[position]);
            }
        } else if alt_allele.is_empty() {
            if position == 0 {
                bail!("alt_allele is empty at position 0");
            }

            // we need to pre-pend the anchor base because one allele is empty
            if let Some(cs) = chrom_seq {
                position -= 1;
                ref_allele.insert(0, cs[position]);
                alt_allele.insert(0, cs[position]);
            }
        }

        // at this point, we should only have ACGTs
        // first, chop off any duplicate sequence at the end
        while ref_allele.len() > 1 && alt_allele.len() > 1 && ref_allele[ref_allele.len() - 1] == alt_allele[alt_allele.len() - 1] {
            // position does not change here
            ref_allele.pop();
            alt_allele.pop();
        }

        // similarly, chop off any duplicate sequence at the start
        while ref_allele.len() > 1 && alt_allele.len() > 1 && ref_allele[0] == alt_allele[0] {
            // we need to shift position as we chop
            position += 1;
            ref_allele.remove(0);
            alt_allele.remove(0);
        }

        // now we need to see if there is any shifting that needs to happen
        while ref_allele[ref_allele.len() - 1] == alt_allele[alt_allele.len() - 1] {
            // TODO: I think this assertion is true, lets leave it in for ourselves
            assert!(ref_allele.len() == 1 || alt_allele.len() == 1);

            if position == 0 {
                // we cannot prepend, just break out?
                break;
            } else if let Some(cs) = chrom_seq {
                // this is pre-pending the reference base
                position -= 1;
                ref_allele.insert(0, cs[position]);
                alt_allele.insert(0, cs[position]);
            } else {
                // we cannot pre-pend because no reference
                break;
            }

            // remove the trailing character that matched
            ref_allele.pop();
            alt_allele.pop();
        }

        // we have finished the byte normalizing, convert everything back into a string
        let reference: String = String::from_utf8(ref_allele)?;
        let alternate: String = String::from_utf8(alt_allele)?;

        // one last sanity check on the allowed bases
        let allowed_bases: HashSet<char> = HashSet::from_iter(['A', 'C', 'G', 'T']);
        if reference.chars().all(|c| allowed_bases.contains(&c)) && 
            alternate.chars().all(|c| allowed_bases.contains(&c)) {
            // all sequences are composed of ACGT
            Ok(NormalizedVariant {
                chrom,
                position,
                reference,
                alternate,
                sv_stats: None
            })
        } else {
            bail!("ACGT alleles only");
        }
    }

    /// This will take an allele definition and normalize it using the reference genome, potentially producing multiple variants that can match.
    /// If one of the results is "None", it indicates that one of the variants is the reference.
    /// For example, if you have A->R, this is the same as (A->A OR A->G) where A->A is just the reference allele.
    /// # Arguments
    /// * `chrom` - the chromosome the variant is positioned on
    /// * `position` - the 0-based coordinate of the variant
    /// * `ref_allele` - the given reference sequence
    /// * `alt_allele` - the given alternate sequence
    /// * `reference_genome` - the pre-loaded reference genome; if None, some normalization steps will not happen and it may fail to normalize
    /// # Errors
    /// * if there are errors creating the sub-alleles
    pub fn multi_new(chrom: String, position: usize, ref_allele: &str, alt_allele: &str, reference_genome: Option<&ReferenceGenome>) -> Result<Vec<Option<NormalizedVariant>>, Box<dyn std::error::Error>> {
        // do any conversions
        let multi_alt = match alt_allele {
            // IUPAC code: https://www.bioinformatics.org/sms/iupac.html
            "K" => vec!["G", "T"],
            "M" => vec!["A", "C"],
            "R" => vec!["A", "G"],
            "S" => vec!["C", "G"],
            "W" => vec!["A", "T"],
            "Y" => vec!["C", "T"],
            "B" => vec!["C", "G", "T"],
            "D" => vec!["A", "G", "T"],
            "H" => vec!["A", "C", "T"],
            "V" => vec!["A", "C", "G"],
            // there is one allele with this pattern: "delinsCC; delinsCCC; delinsCCCC; delinsCCCCC; delinsCCCCCC; delinsCCCCCCC"
            _ => alt_allele.split("; ").collect()
        };

        // now iterate over each possible allele and save it
        let mut ret: Vec<Option<NormalizedVariant>> = vec![];
        for aa in multi_alt {
            if ref_allele == aa {
                // ref and alt match, so just store None; this allows us to match later
                ret.push(None);
            } else {
                // they're different, so convert it
                ret.push(
                    Some(
                        Self::new(
                            chrom.clone(),
                            position,
                            ref_allele,
                            aa,
                            reference_genome
                        )?
                    )
                );
            }
        }
        Ok(ret)
    }

    /// Constructor for a structural variant.
    /// No normalization actually happens here, but it collects stats for translating into haplotypes later.
    /// # Arguments
    /// * `sv_type` - the type of SV detected
    /// * `chrom` - chromosome
    /// * `position` - the start of the event
    /// * `end` - the end coordinate of the event
    /// * `haplotype_label` - the PGx haplotype that this event matches
    pub fn new_sv(sv_type: SvType, chrom: String, position: usize, end: usize, haplotype_label: String) -> Result<Self, SimpleError> {
        let sv_stats = Some(StructuralVariantStats::new(sv_type, position, end, haplotype_label)?);
        Ok(Self {
            chrom,
            position,
            reference: Default::default(),
            alternate: Default::default(),
            sv_stats
        })
    }

    /// Returns true if this variant has SV statistics associated with it
    pub fn is_sv(&self) -> bool {
        self.sv_stats.is_some()
    }

    // getters
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    pub fn position(&self) -> usize {
        self.position
    }

    pub fn sv_stats(&self) -> Option<&StructuralVariantStats> {
        self.sv_stats.as_ref()
    }

    pub fn variant_name(&self) -> String {
        format!("{}:{}{}>{}", self.chrom, self.position, self.reference, self.alternate)
    }
}

/// Converts a sequence from CPIC into a Vec representation, while parsing any regex along the way.
/// STR regex is expanded and indel events are parsed into the corresponding base-only sequence.
/// # Arguments
/// * `sequence` - the sequence we are converting / parsing
fn parse_sequence(sequence: &str) -> Vec<u8> {
    if let Some((_, [seq, count])) = TR_REGEX.captures(sequence).map(|c| c.extract()) {
        // this is a pattern like ACGT(8), so we make ACGT * 8 and return that
        let count = count.parse::<usize>().unwrap();
        seq.repeat(count).into_bytes()
    } else if sequence.starts_with("delins") {
        // just skip the "delins" sequence
        sequence.as_bytes()[6..].to_vec()
    } else if sequence.starts_with("ins") {
        // just skip the "ins" sequence
        sequence.as_bytes()[3..].to_vec()
    } else if sequence.starts_with("del") {
        // this is basically saying to delete the whole thing, so return empty vec
        vec![]
    } else {
        sequence.as_bytes().to_vec()
    }
}

/// Enum for SV types
#[derive(Clone, Copy, Debug, strum_macros::Display, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum SvType {
    /// Anything we haven't set an enum for yet
    Unknown,
    /// Deletion is the only one we currently care about
    Deletion
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct StructuralVariantStats {
    /// The associated SV type
    sv_type: SvType,
    /// The 0-based, inclusive start of the SV
    start: usize,
    /// The 0-based, exclusive end of the SV
    end: usize,
    /// Label that matches a known PGx haplotype
    haplotype_label: String
}

impl StructuralVariantStats {
    /// Constructor
    pub fn new(sv_type: SvType, start: usize, end: usize, haplotype_label: String) -> Result<Self, SimpleError> {
        if start >= end {
            bail!("SV definition requires that start < end");
        }
        Ok(Self {
            sv_type,
            start,
            end,
            haplotype_label
        })
    }

    // getters
    pub fn sv_type(&self) -> SvType {
        self.sv_type
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn end(&self) -> usize {
        self.end
    }

    pub fn haplotype_label(&self) -> &str {
        &self.haplotype_label
    }
}

/// The possible genotypes we allow
#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub enum Genotype {
    /// 0/0
    #[serde(rename = "0/0")]
    HomozygousReference,
    /// 0/1
    #[serde(rename = "0/1")]
    HeterozygousUnphased,
    /// 0|1
    #[serde(rename = "0|1")]
    HeterozygousPhased,
    /// 1|0
    #[serde(rename = "1|0")]
    HeterozygousPhasedFlip,
    /// 1/1
    #[serde(rename = "1/1")]
    HomozygousAlternate
}

/// A normalized genotype is composed of both the genotype and the phase set 
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct NormalizedGenotype {
    /// The genotype
    genotype: Genotype,
    /// The phase set ID, only for phased alleles
    phase_set: Option<usize>
}

impl NormalizedGenotype {
    /// Basic constructor
    pub fn new(genotype: Genotype, phase_set: Option<usize>) -> NormalizedGenotype {
        NormalizedGenotype {
            genotype, phase_set
        }
    }

    pub fn genotype(&self) -> Genotype {
        self.genotype
    }

    pub fn phase_set(&self) -> &Option<usize> {
        &self.phase_set
    }
}

/// Normalized PGx haplotype has both a haplotype name and a list of required variants.
/// This is complicated by the fact that multiple variants can "match" an allele, for example a "C" -> "R" is ["C" -> "A" OR "C" -> "G"].
/// This struct handles the logic around comparing these mixtures of AND and OR composite haplotype definitions.
#[derive(Debug, Eq, PartialEq)]
pub struct NormalizedPgxHaplotype {
    // the name of this haplotype
    haplotype_name: String,
    /// if true, this is an SV haplotype
    is_sv: bool,
    // variants defining the haplotype - the outer Vec is AND, the inner Vec is OR; None indicates that it doesn't have to match this
    variants: Vec<Vec<Option<NormalizedVariant>>>
}

impl NormalizedPgxHaplotype {
    /// Basic constructor
    pub fn new(haplotype_name: String) -> NormalizedPgxHaplotype {
        NormalizedPgxHaplotype {
            haplotype_name,
            is_sv: false,
            variants: vec![]
        }
    }

    /// Adds a single required variant that matches one of the provided NormalizedVariants
    /// # Arguments
    /// * `variant` - a Vec of normalized variants, any of them would be considered a match for this variant
    pub fn add_variant(&mut self, variant: Vec<Option<NormalizedVariant>>) {
        // flatten will skip any None values
        for nv in variant.iter().flatten() {
            if nv.is_sv() {
                self.is_sv = true;
            }
        }
        self.variants.push(variant);
    }

    /// Returns true if a collection of NormalizedVariants matches this haplotype.
    /// A match requires that both all of the provided variants match a variant in this haplotype AND that all variants in the haplotype have a match.
    /// If double matches are detected (e.g. same variant twice), then this will return false.
    /// # Arguments
    /// * `other_variants` - the list of other variants we want to match
    pub fn matches(&self, other_variants: &[NormalizedVariant]) -> bool {
        // initialize that all variants are unmatched so far
        let mut match_vec: Vec<bool> = vec![false; self.variants.len()];
        for ov in other_variants.iter() {
            // find the first variant that matches this one
            // if we ever get multiple matches, this logic will need to change
            let match_index = self.variants.iter()
                .position(|v| {
                    v.iter().any(|sub_v| sub_v.as_ref() == Some(ov))
                });
            match match_index {
                Some(mi) => {
                    if match_vec[mi] {
                        // something has already matched this loci, we cannot double match
                        return false;
                    } else {
                        // mark this one as matched
                        match_vec[mi] = true;
                    }
                },
                None => {
                    // there is no match for this allele
                    return false;
                }
            }
        }

        // all the other_variants had a match, make sure everything in this haplotype also had a match
        match_vec.iter()
            .enumerate()
            .all(|(i, &matched)| {
                // allowed if we either had a match OR None is allowed
                matched || self.variants[i].contains(&None)
            })
    }

    /// Quantifies the match between this haplotype and a list of other variants.
    /// Returns a tuple of (missing_variants, extra_variants)
    ///     missing_variants - the variants that are in this haplotype but not in the other list
    ///     extra_variants - the variants that are in the other list, but not in this haplotype
    /// # Arguments
    /// * `other_variants` - the list of other variants we want to match
    pub fn quant_match(&self, other_variants: &[NormalizedVariant]) -> (Vec<NormalizedVariant>, Vec<NormalizedVariant>) {
        // make sure this is not an SV haplotype
        assert!(!self.is_sv, "SV haplotypes should not be quantified");
        
        let mut missing_variants = vec![];
        let mut extra_variants = vec![];
        // initialize that all variants are unmatched so far
        let mut match_vec: Vec<bool> = vec![false; self.variants.len()];
        for ov in other_variants.iter() {
            // find the first variant that matches this one
            // if we ever get multiple matches, this logic will need to change
            let match_index = self.variants.iter()
                .position(|v| {
                    v.iter().any(|sub_v| sub_v.as_ref() == Some(ov))
                });
            match match_index {
                Some(mi) => {
                    if match_vec[mi] {
                        // something has already matched this loci, we cannot double match
                        extra_variants.push(ov.clone());
                    } else {
                        // mark this one as matched
                        match_vec[mi] = true;
                    }
                },
                None => {
                    // there is no match for this allele
                    extra_variants.push(ov.clone());
                }
            }
        }

        // all the other_variants had a match, make sure everything in this haplotype also had a match
        for (matched, variant) in match_vec.into_iter().zip(self.variants.iter()) {
            if !(matched || variant.contains(&None)) {
                // get the first variant that is not None
                let first_variant = variant.iter().find(|v| v.is_some())
                    .unwrap().clone().unwrap();
                missing_variants.push(first_variant);
            }
        }
        (missing_variants, extra_variants)
    }

    pub fn haplotype_name(&self) -> &str {
        &self.haplotype_name
    }

    pub fn is_sv(&self) -> bool {
        self.is_sv
    }

    pub fn variants(&self) -> &[Vec<Option<NormalizedVariant>>] {
        &self.variants
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    /// Utility that loads our tiny reference for us
    fn load_test_reference() -> ReferenceGenome {
        let ref_fn = PathBuf::from("test_data/test_reference.fa");
        ReferenceGenome::from_fasta(&ref_fn).unwrap()
    }

    /// Checks a basic SNP, no reference
    #[test]
    fn test_normalize_snp() {
        let chrom = "chr1".to_string();
        let position = 10;
        let reference = "A".to_string();
        let alternate = "C".to_string();
        let na = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, None).unwrap();
        assert_eq!(na, NormalizedVariant {
            chrom,
            position,
            reference,
            alternate,
            sv_stats: None
        });
    }

    /// Checks an A->R (non-multi) to verify we get an error.
    #[test]
    fn test_normalize_multisnp() {
        let chrom = "chr1".to_string();
        let position = 10;
        let reference = "A".to_string();
        let alternate = "R".to_string();
        let na = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, None);
        assert!(na.is_err());
    }

    /// Check basic indel normalization (no reference).
    #[test]
    fn test_normalize_indel() {
        let chrom = "chr1".to_string();
        let position = 10;
        let reference = "AC".to_string();
        let alternate = "ACC".to_string();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, None).unwrap();
        let expected = NormalizedVariant {
            chrom,
            position,
            reference: "A".to_string(),
            alternate: "AC".to_string(),
            sv_stats: None
        };
        assert_eq!(nv, expected);
    }

    /// Same as above, but with reference
    #[test]
    fn test_normalize_ins_ref_001() {
        // this case should get normalize to chr1:10 A -> AC
        let chrom = "chr1".to_string();
        let position = 10;
        let reference = "AC".to_string();
        let alternate = "ACC".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = NormalizedVariant {
            chrom,
            position,
            reference: "A".to_string(),
            alternate: "AC".to_string(),
            sv_stats: None
        };
        assert_eq!(nv, expected);
    }

    /// tests for removing redundant bases in prefix
    #[test]
    fn test_normalize_ins_ref_002() {
        // this case should get normalize to chr1:12 A -> AC
        let chrom = "chr1".to_string();
        let position = 10;
        let reference = "ACAC".to_string();
        let alternate = "ACACC".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = NormalizedVariant {
            chrom,
            position: 12,
            reference: "A".to_string(),
            alternate: "AC".to_string(),
            sv_stats: None
        };
        assert_eq!(nv, expected);
    }

    /// tests for left shifting a large insertion in a TR
    #[test]
    fn test_normalize_ins_ref_003() {
        // this case should get normalize to chr1:9 A -> AACAC
        let chrom = "chr1".to_string();
        let position = 10;
        let reference = "ACACACACAC".to_string();
        let alternate = "ACACACACACACAC".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = NormalizedVariant {
            chrom,
            position: 9,
            reference: "A".to_string(),
            alternate: "AACAC".to_string(),
            sv_stats: None
        };
        assert_eq!(nv, expected);
    }

    /// test for normalizing a deletion in a TR
    #[test]
    fn test_normalize_del_ref_001() {
        // this case should get normalize to chr1:9 AAC -> A
        let chrom = "chr1".to_string();
        let position = 16;
        let reference = "ACAC".to_string();
        let alternate = "AC".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = NormalizedVariant {
            chrom,
            position: 9,
            reference: "AAC".to_string(),
            alternate: "A".to_string(),
            sv_stats: None
        };
        assert_eq!(nv, expected);
    }

    /// test for normalizing a del in an TR
    #[test]
    fn test_normalize_del_ref_002() {
        // this case should get normalize to chr1:9 AAC -> A
        let chrom = "chr1".to_string();
        let position = 16;
        let reference = "ACAC".to_string();
        let alternate = "AC".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = NormalizedVariant {
            chrom,
            position: 9,
            reference: "AAC".to_string(),
            alternate: "A".to_string(),
            sv_stats: None
        };
        assert_eq!(nv, expected);
    }

    /// make sure that empty ref & alt throws an error
    #[test]
    fn test_empty_refalt() {
        // this case should get normalize to chr1:9 CAGT -> C
        let chrom = "chr2".to_string();
        let position = 13;
        let reference = "".to_string();
        let alternate = "".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome));
        assert!(nv.is_err());
    }

    /// make sure that an empty alt will add the anchor base and left-shift
    #[test]
    fn test_empty_alt() {
        // this case should get normalize to chr1:9 CAGT -> C
        let chrom = "chr2".to_string();
        let position = 13;
        let reference = "AGT".to_string();
        let alternate = "".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = NormalizedVariant {
            chrom,
            position: 9,
            reference: "CAGT".to_string(),
            alternate: "C".to_string(),
            sv_stats: None
        };
        assert_eq!(nv, expected);
    }

    /// make sure the reference always matches
    #[test]
    fn test_ref_mismatch() {
        // this case should get fail to normalize due to reference mismatch
        let chrom = "chr2".to_string();
        let position = 13;
        let reference = "MISS".to_string();
        let alternate = "A".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome));
        assert!(nv.is_err());
    }

    /// CPIC has an annoying syntax for insertions, this checks that we appropriately anchor it.
    #[test]
    fn test_cpic_ins() {
        // this case should get normalize to chr1:9 C -> CAGT
        let chrom = "chr2".to_string();
        let position = 12;
        let reference = "del".to_string();
        let alternate = "insAGT".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = NormalizedVariant {
            chrom,
            position: 9,
            reference: "C".to_string(),
            alternate: "CAGT".to_string(),
            sv_stats: None
        };
        assert_eq!(nv, expected);
    }

    /// tests a CPIC-style deletion
    #[test]
    fn test_cpic_del() {
        // this case should get normalize to chr1:9 C -> CAGT
        let chrom = "chr2".to_string();
        let position = 13;
        let reference = "AGT".to_string();
        let alternate = "delAGT".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = NormalizedVariant {
            chrom,
            position: 9,
            reference: "CAGT".to_string(),
            alternate: "C".to_string(),
            sv_stats: None
        };
        assert_eq!(nv, expected);
    }

    /// tests a CPIC-style indel
    #[test]
    fn test_cpic_delins() {
        // this case should get normalize to chr1:10 A -> CGG
        let chrom = "chr2".to_string();
        let position = 10;
        let reference = "A".to_string();
        let alternate = "delinsCGG".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = NormalizedVariant {
            chrom,
            position,
            reference: "A".to_string(),
            alternate: "CGG".to_string(),
            sv_stats: None
        };
        assert_eq!(nv, expected);
    }

    /// Tests a CPIC deletion in a TR
    #[test]
    fn test_cpic_tr_del() {
        // this case should get normalize to chr1:10 A -> CGG
        let chrom = "chr2".to_string();
        let position = 10;
        let reference = "AGT(3)".to_string();
        let alternate = "AGT(2)".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = NormalizedVariant {
            chrom,
            position: 9,
            reference: "CAGT".to_string(),
            alternate: "C".to_string(),
            sv_stats: None
        };
        assert_eq!(nv, expected);
    }

    /// Tests a CPIC inseriton in a TR
    #[test]
    fn test_cpic_tr_ins() {
        // this case should get normalize to chr1:10 A -> CGG
        let chrom = "chr2".to_string();
        let position = 10;
        let reference = "AGT(3)".to_string();
        let alternate = "AGT(4)".to_string();
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = NormalizedVariant {
            chrom,
            position: 9,
            reference: "C".to_string(),
            alternate: "CAGT".to_string(),
            sv_stats: None
        };
        assert_eq!(nv, expected);
    }

    /// Tests a CPIC IUPAC multi-match with the reference allele
    #[test]
    fn test_multinew_iupac_ref_included() {
        // this case should get normalize to chr1:10 A -> [None, G]
        let chrom = "chr1".to_string();
        let position = 10;
        let reference = "A".to_string();
        let alternate = "R".to_string(); // A or G
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::multi_new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = vec![
            None,
            Some(NormalizedVariant {
                chrom,
                position: 10,
                reference: "A".to_string(),
                alternate: "G".to_string(),
                sv_stats: None
            })
        ];
        assert_eq!(nv, expected);
    }

    /// Tests a CPIC IUPAC multi-match where both are non-reference
    #[test]
    fn test_multinew_iupac_double_alt() {
        // this case should get normalize to chr1:10 A -> [C, T]
        let chrom = "chr1".to_string();
        let position = 10;
        let reference = "A".to_string();
        let alternate = "Y".to_string(); // C or T
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::multi_new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = vec![
            Some(NormalizedVariant {
                chrom: chrom.clone(),
                position: 10,
                reference: "A".to_string(),
                alternate: "C".to_string(),
                sv_stats: None
            }),
            Some(NormalizedVariant {
                chrom,
                position: 10,
                reference: "A".to_string(),
                alternate: "T".to_string(),
                sv_stats: None
            })
        ];
        assert_eq!(nv, expected);
    }

    /// Tests a CPIC multi-match using the "; " delimiter. There is only one in CPIC thus far.
    #[test]
    fn test_multinew_semicolon() {
        // this case should get normalize to chr1:10 A -> [C, CC, CCC]
        let chrom = "chr1".to_string();
        let position = 10;
        let reference = "A".to_string();
        let alternate = "delinsC; delinsCC; delinsCCC".to_string(); // C or CC or CCC
        let reference_genome = load_test_reference();
        let nv = NormalizedVariant::multi_new(chrom.clone(), position, &reference, &alternate, Some(&reference_genome)).unwrap();
        let expected = vec![
            Some(NormalizedVariant {
                chrom: chrom.clone(),
                position: 10,
                reference: "A".to_string(),
                alternate: "C".to_string(),
                sv_stats: None
            }),
            Some(NormalizedVariant {
                chrom: chrom.clone(),
                position: 10,
                reference: "A".to_string(),
                alternate: "CC".to_string(),
                sv_stats: None
            }),
            Some(NormalizedVariant {
                chrom,
                position: 10,
                reference: "A".to_string(),
                alternate: "CCC".to_string(),
                sv_stats: None
            })
        ];
        assert_eq!(nv, expected);
    }

    /// make sure an empty haplotype only matches empty list
    #[test]
    fn test_ref_matches() {
        // can only match reference allele
        let pgx_hap = NormalizedPgxHaplotype::new("test".to_string());
        let test_variant = NormalizedVariant::new("chr1".to_string(), 10, "A", "C", None).unwrap();
        assert!(pgx_hap.matches(&[]));
        assert!(!pgx_hap.matches(&[
            test_variant.clone()
        ]));

        // now do a quant test
        assert_eq!(pgx_hap.quant_match(&[]), (vec![], vec![]));
        assert_eq!(pgx_hap.quant_match(&[
            test_variant.clone()
        ]), (vec![], vec![test_variant]));
    }

    /// test a simple single-variant haplotype match
    #[test]
    fn test_alt_matches() {
        // can only match alternate allele
        let mut pgx_hap = NormalizedPgxHaplotype::new("test".to_string());
        let test_variant = NormalizedVariant::new("chr1".to_string(), 10, "A", "C", None).unwrap();
        pgx_hap.add_variant(vec![Some(test_variant.clone())]);
        assert!(!pgx_hap.matches(&[]));
        assert!(pgx_hap.matches(&[
            test_variant.clone()
        ]));

        // now do a quant test
        assert_eq!(pgx_hap.quant_match(&[]), (vec![test_variant.clone()], vec![]));
        assert_eq!(pgx_hap.quant_match(&[
            test_variant.clone()
        ]), (vec![], vec![]));
    }

    /// tests that optional matches work both with and without the optional variant
    #[test]
    fn test_optional_matches() {
        // optional variant, so both with and without will match
        let mut pgx_hap = NormalizedPgxHaplotype::new("test".to_string());
        let test_variant = NormalizedVariant::new("chr1".to_string(), 10, "A", "C", None).unwrap();
        pgx_hap.add_variant(vec![None, Some(test_variant.clone())]);
        assert!(pgx_hap.matches(&[]));
        assert!(pgx_hap.matches(&[
            test_variant.clone()
        ]));

        // now do a quant test - here, presence and absence should both be empty since the variant is optional
        assert_eq!(pgx_hap.quant_match(&[]), (vec![], vec![]));
        assert_eq!(pgx_hap.quant_match(&[
            test_variant.clone()
        ]), (vec![], vec![]));
    }

    /// tests a multi-variant list, make sure that each individually will match, but both together should not
    #[test]
    fn test_multivariant_matches() {
        // optional variant, so both with and without will match
        let mut pgx_hap = NormalizedPgxHaplotype::new("test".to_string());
        let test_variant_form1 = NormalizedVariant::new("chr1".to_string(), 10, "A", "C", None).unwrap();
        let test_variant_form2 = NormalizedVariant::new("chr1".to_string(), 10, "A", "T", None).unwrap();
        pgx_hap.add_variant(vec![Some(test_variant_form1.clone()), Some(test_variant_form2.clone())]);

        // either one by itself should match, both together should not
        assert!(pgx_hap.matches(&[
            test_variant_form1.clone()
        ]));
        assert!(pgx_hap.matches(&[
            test_variant_form2.clone()
        ]));
        assert!(!pgx_hap.matches(&[
            test_variant_form1.clone(),
            test_variant_form2.clone()
        ]));

        // now do a quant test
        assert_eq!(pgx_hap.quant_match(&[]), (vec![test_variant_form1.clone()], vec![])); // only the first is marked as missing
        assert_eq!(pgx_hap.quant_match(&[
            test_variant_form1.clone()
        ]), (vec![], vec![]));
        assert_eq!(pgx_hap.quant_match(&[
            test_variant_form2.clone()
        ]), (vec![], vec![]));
        assert_eq!(pgx_hap.quant_match(&[
            test_variant_form1.clone(),
            test_variant_form2.clone()
        ]), (vec![], vec![test_variant_form2.clone()]));
    }
}
