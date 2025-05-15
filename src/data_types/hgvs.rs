
use regex::Regex;
use simple_error::bail;

/// Represent HGVS nomenclature but parsed into something usable. We are only considering g. for now.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParsedHgvs {
    /// Chromosome string, but usually in annoying form
    chrom: String,
    /// 1-based position due to HGVS nomenclature
    position: u64,
    /// Type of variant and any relevant sub-fields
    variant: HgvsVariant
}

impl ParsedHgvs {
    /// Constructor, which just parses a string
    pub fn new(hgvs_str: &str) -> Result<Self, Box<dyn std::error::Error>> {
        // we have a lot of Regexs to test here, this is probably sub-optimal but we do this once in a blue moon so who cares
        let snv_regex = Regex::new("^(?<chrom>.+):g.(?<pos>\\d+)(?<ref>[ACGNT])>(?<alt>[ACGNT])$")?;
        let del_regex = Regex::new("^(?<chrom>.+):g.(?<pos>\\d+)del(?<alt>[ACGNT]+)$")?;
        let ins_regex = Regex::new("^(?<chrom>.+):g.(?<start>\\d+)_(?<end>\\d+)ins(?<alt>[ACGNT]+)$")?;

        let result = if let Some(captures) = snv_regex.captures(hgvs_str) {
            let chrom = captures.name("chrom").unwrap().as_str();
            let pos: u64 = captures.name("pos").unwrap().as_str().parse()?;
            let ref_allele= captures.name("ref").unwrap().as_str();
            let alt_allele = captures.name("alt").unwrap().as_str();

            assert_eq!(ref_allele.len(), 1);
            assert_eq!(alt_allele.len(), 1);

            ParsedHgvs {
                chrom: chrom.to_string(),
                position: pos,
                variant: HgvsVariant::Snv {
                    ref_allele: ref_allele.to_string(),
                    alt_allele: alt_allele.to_string()
                }
            }
        } else if let Some(captures) = del_regex.captures(hgvs_str) {
            let chrom = captures.name("chrom").unwrap().as_str();
            let pos: u64 = captures.name("pos").unwrap().as_str().parse()?;
            let deleted = captures.name("alt").unwrap().as_str().to_string();

            ParsedHgvs {
                chrom: chrom.to_string(),
                position: pos,
                variant: HgvsVariant::Deletion { deleted }
            }
        } else if let Some(captures) = ins_regex.captures(hgvs_str) {
            let chrom = captures.name("chrom").unwrap().as_str();
            let start: u64 = captures.name("start").unwrap().as_str().parse()?;
            let end: u64 = captures.name("end").unwrap().as_str().parse()?;
            let inserted = captures.name("alt").unwrap().as_str().to_string();

            ParsedHgvs {
                chrom: chrom.to_string(),
                position: start,
                variant: HgvsVariant::Insertion { end, inserted }
            }
        } else {
            bail!("failed to parse {hgvs_str}");
        };

        Ok(result)
    }

    /// Converts this into REF and ALT sequences using the provided reference sequence.
    /// It also returns the 1-based position for those REF/ALT sequences.
    /// # Arguments
    /// * `reference` - the reference chromosome sequence
    /// # Errors
    /// * if any REF sequences we know about do not match what is provided in `reference`
    pub fn generate_ref_alt(&self, reference: &[u8]) -> Result<(usize, String, String), Box<dyn std::error::Error>> {
        // zero-based position
        let zpos = self.position as usize - 1;
        let (pos, r, a) = match &self.variant {
            HgvsVariant::Snv { ref_allele, alt_allele } => {
                // verify our reference base matches
                if &reference[zpos..(zpos+1)] != ref_allele.as_bytes() {
                    bail!("REF allele base does not match provided reference");
                }
                (zpos+1, ref_allele.clone(), alt_allele.clone())
            },
            HgvsVariant::Insertion { end, inserted } => {
                if *end as usize != zpos+2 {
                    bail!("Unexpected end on insertion");
                }
                let ref_seq = std::str::from_utf8(&reference[zpos..(zpos+1)])?.to_string();
                let alt_seq = format!("{ref_seq}{inserted}");
                (zpos+1, ref_seq, alt_seq)
            },
            HgvsVariant::Deletion { deleted } => {
                if &reference[zpos..(zpos+deleted.len())] != deleted.as_bytes() {
                    bail!("Deleted based do not match provided reference");
                }
                let ref_seq = std::str::from_utf8(&reference[(zpos-1)..(zpos+deleted.len())])?.to_string();
                let alt_seq = ref_seq[0..1].to_string();
                (zpos, ref_seq, alt_seq)
            },
        };

        Ok((pos, r, a))
    }

    // getters
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    pub fn position(&self) -> u64 {
        self.position
    }

    pub fn variant(&self) -> &HgvsVariant {
        &self.variant
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum HgvsVariant {
    Snv {
        ref_allele: String,
        alt_allele: String
    },
    Insertion {
        end: u64,
        inserted: String
    },
    Deletion {
        deleted: String
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hgvs_variant() {
        /*
        original HGVS
        todo!("NC_000015.10:g.74749863C>G");
        todo!("NC_000015.10:g.74749841del");
        todo!("NC_000019.10:g.40843713_40843714del");
        todo!("NC_000022.11:g.42126658_42126666dup");
        todo!("NC_000019.10:g.40843720_40843721insG");
        NC_000019.10:g.40848264GC[1] <- this one is super annoying
        */
        let chrom = "mock";
        let reference = "ACGTACGTACGT";

        // test the SNV format
        let snv_hgvs = format!("{chrom}:g.6C>G");
        let snv_var = ParsedHgvs::new(&snv_hgvs).unwrap();
        assert_eq!(snv_var, ParsedHgvs {
            chrom: chrom.to_string(),
            position: 6,
            variant: HgvsVariant::Snv { ref_allele: "C".to_string(), alt_allele: "G".to_string() }
        });
        assert_eq!(snv_var.generate_ref_alt(reference.as_bytes()).unwrap(), (6, "C".to_string(), "G".to_string()));

        // test the INS format
        let ins_hgvs = format!("{chrom}:g.8_9insTTG");
        let ins_var = ParsedHgvs::new(&ins_hgvs).unwrap();
        assert_eq!(ins_var, ParsedHgvs {
            chrom: chrom.to_string(),
            position: 8,
            variant: HgvsVariant::Insertion { end: 9, inserted: "TTG".to_string() }
        });
        assert_eq!(ins_var.generate_ref_alt(reference.as_bytes()).unwrap(), (8, "T".to_string(), "TTTG".to_string()));

        // test the DEL format
        let del_hgvs = format!("{chrom}:g.6delCGT");
        let del_var = ParsedHgvs::new(&del_hgvs).unwrap();
        assert_eq!(del_var, ParsedHgvs {
            chrom: chrom.to_string(),
            position: 6,
            variant: HgvsVariant::Deletion { deleted: "CGT".to_string() }
        });
        assert_eq!(del_var.generate_ref_alt(reference.as_bytes()).unwrap(), (5, "ACGT".to_string(), "A".to_string()));
    }
}