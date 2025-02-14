
use std::collections::BTreeMap;
use std::fmt::Display;

/// Core region types that are associated with our CYP2D6 locus
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, strum_macros::Display)]
pub enum Cyp2d6RegionType {
    /// Generic unknown type
    #[strum(to_string="UNKNOWN")]
    Unknown,
    /// REP6 region, before CYP2D6 typically
    #[strum(to_string="REP6")]
    Rep6,
    /// Main CYP2D6 allele type
    #[strum(to_string="CYP2D6")]
    Cyp2d6,
    /// Link region, follows CYP2D6
    #[strum(to_string="link_region")]
    LinkRegion,
    /// REP7, follow link region typically
    #[strum(to_string="REP7")]
    Rep7,
    /// Spacer region, between REP7 and CYP2D7
    #[strum(to_string="spacer")]
    Spacer,
    /// Main CYP2D7 type
    #[strum(to_string="CYP2D7")]
    Cyp2d7,
    /// CYP2D6*5
    #[strum(to_string="CYP2D6*5")]
    Cyp2d6Deletion,
    /// Currently capturing both D6::D7 and D7::D6 hybrids
    Hybrid,
    /// Sentinel for when we call an allele but it turns out to be a false positive
    FalseAllele
}

impl Cyp2d6RegionType {
    /// Returns true if this is a "classic" CYP2D region
    pub fn is_cyp2d(&self) -> bool {
        match self {
            // these do not count
            Cyp2d6RegionType::Unknown | 
            Cyp2d6RegionType::Rep6 |
            Cyp2d6RegionType::LinkRegion |
            Cyp2d6RegionType::Rep7 |
            Cyp2d6RegionType::Spacer |
            Cyp2d6RegionType::FalseAllele => false,
            // D6, D7, *5, and hybrids count
            Cyp2d6RegionType::Cyp2d6 |
            Cyp2d6RegionType::Cyp2d7 |
            Cyp2d6RegionType::Cyp2d6Deletion |
            Cyp2d6RegionType::Hybrid => true
        }
    }

    /// Returns true if this is a "REP" region
    pub fn is_rep(&self) -> bool {
        matches!(self, Cyp2d6RegionType::Rep6 | Cyp2d6RegionType::Rep7)
    }

    /// Returns true if this type would show up in the final diplotype
    pub fn is_reported_allele(&self) -> bool {
        matches!(self,
            Cyp2d6RegionType::Cyp2d6 |
            Cyp2d6RegionType::Cyp2d6Deletion |
            Cyp2d6RegionType::Hybrid
        )
    }
}

/// Fully describes a region we identify, first by the region type (D6 + nearby components) and then optionally with a more descriptive label.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Cyp2d6RegionLabel {
    /// The high level type of region this represents
    region_type: Cyp2d6RegionType,
    /// If Some, this will generally contain a star allele; e.g. "*4.001".
    subtype_label: Option<String>
}

impl Cyp2d6RegionLabel {
    /// Constructor
    pub fn new(region_type: Cyp2d6RegionType, subtype_label: Option<String>) -> Cyp2d6RegionLabel {
        Cyp2d6RegionLabel {
            region_type,
            subtype_label
        }
    }

    /// Converts this into a FalseAllele, preserving anything else that may be present.
    pub fn mark_false_allele(&mut self) {
        self.region_type = Cyp2d6RegionType::FalseAllele;
    }

    /// This constructs a simplified version of a particular allele.
    /// E.g. CYP2D6*4.001 -> *4.001 (detailed) OR *4 (!detailed).
    /// CYP2D7 and other non-descript hybrid should retain the full allele def.
    /// # Arguments
    /// * `detailed` - if False, then this will reduce any D6 subunits into their integer form; e.g., *4.001 -> *4
    /// * `cyp_translate` - a translation hashmap from internal representation to user-friendly
    pub fn simplify_allele(&self, detailed: bool, cyp_translate: &BTreeMap<String, String>) -> String {
        match self.region_type {
            Cyp2d6RegionType::Cyp2d6 |
            Cyp2d6RegionType::Hybrid => {
                if let Some(subtype_label) = self.subtype_label.as_deref() {
                    // we should not have to strip out the '*' here
                    if let Some(translation) = cyp_translate.get(subtype_label) {
                        // we have a direct translation already, so do that
                        format!("*{translation}")
                    } else if detailed {
                        // they want detailed, so no changes to the output
                        format!("*{subtype_label}")
                    } else {
                        // non-detailed, simplify any floats into ints
                        match subtype_label.parse::<f64>() {
                            Ok(float_value) => {
                                let int_value = float_value.floor() as i64;
                                format!("*{int_value}")
                            },
                            Err(_e) => {
                                // we failed to parse this one into a float, so just strip the prefix and replace with an asterisk
                                format!("*{subtype_label}")
                            }
                        }
                    }
                } else {
                    // we don't expect this to happen normally
                    self.full_allele()
                }
            },
            // this will always be *5
            Cyp2d6RegionType::Cyp2d6Deletion => "*5".to_string(),
            // all the other should not get reported, so just spit out the internal string format
            _ => self.full_allele()
        }
    }

    // Assembles a string version of the full allele name
    pub fn full_allele(&self) -> String {
        match self.region_type {
            // pass through the strum format
            Cyp2d6RegionType::Unknown |
            Cyp2d6RegionType::Rep6 |
            Cyp2d6RegionType::LinkRegion |
            Cyp2d6RegionType::Rep7 |
            Cyp2d6RegionType::Spacer |
            Cyp2d6RegionType::Cyp2d7 |
            Cyp2d6RegionType::Cyp2d6Deletion => format!("{}", self.region_type),
            // e.g. CYP2D6*4.001
            Cyp2d6RegionType::Cyp2d6 => if let Some(stl) = self.subtype_label.as_ref() {
                format!("{}*{}", self.region_type, stl)
            } else {
                self.region_type.to_string()
            },
            // Hybrids will usually get translated later
            Cyp2d6RegionType::Hybrid => if let Some(stl) = self.subtype_label.as_ref() {
                stl.clone()
            } else {
                self.region_type.to_string()
            }
            // Hybrids will usually get translated later
            Cyp2d6RegionType::FalseAllele => if let Some(stl) = self.subtype_label.as_ref() {
                format!("{}_{}", self.region_type, stl)
            } else {
                self.region_type.to_string()
            }
        }
    }

    /// Returns true if this label is allowed to be a part of a chain
    pub fn is_allowed_label(&self) -> bool {
        !matches!(self.region_type, Cyp2d6RegionType::Unknown | Cyp2d6RegionType::FalseAllele)
    }

    /// Checks a link candidate is allowed to be linked from this label.
    /// # Arguments
    /// * `link_candidate` - the label we want to link to
    pub fn is_allowed_label_pair(&self, link_candidate: &Cyp2d6RegionLabel) -> bool {
        use Cyp2d6RegionType::*;

        // we can't have two *5s in one allele; this sometimes is "observed" when we have *5 but no D7 alleles identified
        let type1 = self.region_type();
        let type2 = link_candidate.region_type();
        let double_star5 = type1 == Cyp2d6Deletion && type2 == Cyp2d6Deletion;

        let unexpected_order = // there may be a bunch of these
            // in normal land, we expect REP6 to be the start
            // Note: its very similar to REP7, so this restriction may need to change to a penalty if this becomes an issue
            type2 == Rep6 ||
            // if we have a CYP2D allele (except *5), we should always expect a link region to follow
            (type1.is_cyp2d() && type1 != Cyp2d6Deletion && type2 != LinkRegion) ||
            // if we are joining a link region, we should always have a full allele before it
            (type2 == LinkRegion && !type1.is_cyp2d()) ||
            // if we have a link_region, we should always expect a REP to follow
            (type1 == LinkRegion && !type2.is_rep()) ||
            // REP should always be preceeded by link_region
            (type2.is_rep() && type1 != LinkRegion) ||
            // REP should always be followed by either the spacer or a CYP2D
            (type1.is_rep() && !(type2 == Spacer || type2.is_cyp2d())) ||
            // if the chained-to is "spacer", then we expect a REP before it OR special case deletion: *5
            (type2 == Spacer && !(type1.is_rep() || type1 == Cyp2d6Deletion)) ||
            // if the chained-from is "spacer", then we should have a CYP2D alleles that follows
            (type1 == Spacer && !type2.is_cyp2d()) ||
            // if the second allele is a D7 allele, then we should have a spacer before it
            (type2 == Cyp2d7 && type1 != Spacer) ||
            // we do not expect any extensions after D7; it's the final link in a chain
            type1 == Cyp2d7
            // TODO: we add the above because of a problem in NA12877 where a region gets skipped in a read, leading to problematic mapping
            //       this mapping has a higher "error" because it's just the wrong allele entirely, and it dominates the compute
            //       what we really want is a mechanism to recognize and ignore when that happens, likely with a penalty still, but not like what it's currently doing
            //       ideas: put a cap on error contribution per read? ignore anything > X, 
            //              and add a static penalty because ignored?
            //              detect and allow for a "skip" of the bad allele, again with a penalty?
            //              add back in the explained reads threshold as primary?
            // TODO: any others?
        ;

        // only allowed if not a double *5 AND
        !double_star5 && 
            // not in an unexpected order
            !unexpected_order
    }

    /// This will return true if this allele is a valid start point for a chain.
    /// # Arguments
    /// * `normalize_all_alleles` - if True, then all CYP2D alleles are used for coverage normalization
    pub fn is_candidate_chain_head(&self, normalize_all_alleles: bool) -> bool {
        use Cyp2d6RegionType::*;
        match self.region_type {
            // REP6 is normally the start, unless we have a deletion allele; both are explicitly allowed
            Rep6 |
            Cyp2d6Deletion => true,

            // Sometimes we don't get a REP6, in which case we should allow these alleles if they're used for normalizing
            Cyp2d6 |
            Hybrid => self.is_normalizing_allele(normalize_all_alleles),
            
            // these should NEVER start a chain
            Unknown |
            LinkRegion |
            Rep7 |
            Spacer |
            Cyp2d7 |
            FalseAllele => false
        }
    }

    /// Checks if a label matches one of those for normalizing coverage under targeted conditions where some alleles are under-captured.
    /// In general, we expect all D6 alleles to get captured fine.
    /// However, D7 is often an "off-target" allele, and sometimes the hybrids fall into that off-target bucket as well.
    /// # Arguments
    /// * `normalize_all_alleles` - if true, this will include all full length alleles in normalization
    pub fn is_normalizing_allele(&self, normalize_all_alleles: bool) -> bool {
        if normalize_all_alleles {
            // normalize all valid CYP2D alleles; included D6, D7, hybrids, and deletion
            self.region_type.is_cyp2d()
        } else {
            // only allow the D6 alleles; hybrids, deletions, and D7 have wonky coverage
            self.region_type == Cyp2d6RegionType::Cyp2d6
        }
    }

    /// Returns true if this type would show up in the final diplotype
    pub fn is_reported_allele(&self) -> bool {
        self.region_type.is_reported_allele()
    }

    // wrappers
    /// Returns true if this allele is considered a D6 or D7 allele
    pub fn is_cyp2d(&self) -> bool {
        self.region_type.is_cyp2d()
    }

    /// Creates a new unknown label.
    pub fn new_unknown() -> Cyp2d6RegionLabel {
        Cyp2d6RegionLabel::new(Cyp2d6RegionType::Unknown, None)
    }

    // getters
    pub fn region_type(&self) -> Cyp2d6RegionType {
        self.region_type
    }

    pub fn subtype_label(&self) -> Option<&str> {
        self.subtype_label.as_deref()
    }
}

impl Display for Cyp2d6RegionLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.full_allele())?;
        Ok(())
    }
}
