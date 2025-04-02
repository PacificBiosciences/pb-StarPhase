# Debug Outputs
This file contains details on some of the outputs from the `--output-debug` option.

## CYP2D6 alleles
Filename: `cyp2d6_alleles.json`

Contains the detailed variants information for identified _CYP2D6_ alleles.
This is primarily for users looking to discover differences relative to the database alleles.

Fields:
* `hap1` and `hap2` - The assigned haplotypes from StarPhase
  * `deep_form` - A deep haplotype string containing all alleles on the haplotype. Each allele will be of the form `({index}_{label} {variant_delta})`. The `{index}` is a unique ID for when multiple `{label}`s are identical. Variant deltas will be included if the observed haplotype does not exactly match the assigned database haplotype. In the example below, we see `"(4_CYP2D6*41.004 +rs28735595)"` indicating that the allele most closely matches *41.004, but it also includes an additional unexpected variant rs28735595.
  * `suballele_form` - The sub-allele form of the haplotype string.
  * `core_form` - The core allele form of the haplotype string.
* `alleles` - A dictionary where each key corresponds to a uniquely discovered allele in the dataset, with the form `{index}_{label}`. Each value is a list of relevant variants that were identified and their state information.
  * `label` - The variant label, which is typically an rsID. If no rsID is included in the database, then this will match the following format: `{chrom}:{position}{ref}>{alt}`.
  * `is_vi` - Indicates that the variant is flagged with the `VI` field, indicating that it impacts the cDNA and protein. VI variants determine the core allele while non-VI variant influence the sub-allele.
  * `variant_state` - Indicates the state of the variant in the sample. It will have one of the following entries:
    * `Match` - The variant was found in the same state as expected by the allele. Typically, this means both are ALT state.
    * `Unexpected` - The variant was identified in the ALT state, but the database allele expected the REF state.
    * `Missing` - The variant was identified in the REF state, but the database allele expected the ALT state.
    * `AmbiguousUnexpected` - The variant was found in an ambiguous state, but the database allele expected the REF state.
    * `AmbiguousMissing` - The variant was found in an ambiguous state, but the database allele expected the ALT state.
    * `UnknownUnexpected` - The variant was found in an unknown state, but the database allele expected the REF state. This is common when variants overlap and one of them is found in an ALT state.
    * `UnknownMissing` - The variant was found in an unknown state, but the database allele expected the ALT state.

Example:
```
{
  "hap1": {
    "deep_form": "(4_CYP2D6*41.004 +rs28735595)",
    "suballele_form": "*41.004",
    "core_form": "*41"
  },
  "hap2": {
    "deep_form": "(0_CYP2D6::CYP2D7::exon2) + (2_CYP2D6*4.001 +rs28735595)",
    "suballele_form": "*68 + *4.001",
    "core_form": "*68 + *4"
  },
  "alleles": {
    "2_CYP2D6*4.001": [
      {
        "label": "rs28371738",
        "is_vi": false,
        "variant_state": "Match"
      },
      ...
    ],
    "4_CYP2D6*41.004": [
      ...
    ]
  }
}
```

## HLA alleles
Filename: `hla_debug.json`

Contains the alignment information for HLA database sequences against the StarPhase consensus sequences.
This is primarily for users looking to discover new HLA alleles or debug HLA allele assignment in StarPhase.

Fields:
* `read_mapping_stats` - Contains the mapping for each database sequence against each consensus sequence. All sub-keys of this dictionary are gene names, `{gene}`, with the following fields:
  * `consensus1` - Statistics for the mapping database sequences against consensus sequence 1
    * `best_match_id` - The IMGT/HLA database ID of the best matching HLA sequence. When determining the best match, alleles are scored first by the cDNA mapping, and then second by the DNA mapping.
    * `best_match_star` - The star allele of the best matching HLA sequence
    * `mapping_stats` - A dictionary of all mapping comparisons between the consensus sequence and the database entries for the `{gene}`. All sub-keys of this dictionary are star alleles, `{star_allele}`, with the following fields:
      * `cdna_mapping` - The mapping statistics for the database cDNA sequence against the consensus cDNA sequence. The consensus sequences are typically longer than the database sequences to allow for a buffer during mapping.
        * `query_len` - The length of the database sequence
        * `target_len` - The length of the consensus sequence
        * `match_len` - The number of matching bases from the alignment
        * `nm` - NM tag from minimap2 mapping. This represents the edit distance of the aligned sequence.
        * `query_unmapped` - The number of bases from the database sequence that were left unmapped. If `nm + query_unmapped == 0`, then the consensus sequence is an **exact match** to the database sequence (see example below).
        * `target_unmapped` - The number of bases from the consensus sequence that were unused by the mapping
        * `cigar` - The CIGAR string of the mapping. This is useful for determining where differences between the consensus and database sequences are.
        * `md` - The MD string of the mapping
      * `dna_mapping` - The same statistics as `cdna_mapping`, but for the full-length DNA sequence. These may be unavailable if the database entry is missing a DNA sequence.
  * `consensus2` - Statistics for the mapping database sequences against consensus sequence 2. All sub-fields are identical to those of `consensus1`.
* `dual_passing_stats` - Contains the statistics for determining if a diplotype is reported as homozygous or heterozygous. All sub-keys of this dictionary are gene names, `{gene}`, with the following fields:
  * `is_passing` - if True, then the diplotype is a passing heterozygous call
  * `is_dual` - if True, then a dual consensus was identified (it may not be passing); if False, the following values will be empty because only a single consensus was identified
  * `counts1` - the number of reads assigned to consensus 1
  * `counts2` - the number of reads assigned to consensus 2
  * `maf` - Minor allele frequency = `min(counts1, counts2) / (counts1 + counts2)`
  * `cdf` - Cumulative distribution function value for the provided expected MAF

Example:
```
{
  "read_mapping_stats": {
    "HLA-A": {
      "consensus1": {
        "best_match_id": "HLA:HLA00090",
        "best_match_star": "30:02:01:01",
        "mapping_stats": {
          ...
          "30:02:01:01": {
            "cdna_mapping": {
              "query_len": 1098,
              "target_len": 1535,
              "match_len": 1098,
              "nm": 0,
              "query_unmapped": 0,
              "target_unmapped": 437,
              "cigar": "1098=",
              "md": "1098"
            },
            "dna_mapping": {
              "query_len": 3503,
              "target_len": 3818,
              "match_len": 3503,
              "nm": 0,
              "query_unmapped": 0,
              "target_unmapped": 315,
              "cigar": "3503=",
              "md": "3503"
            }
          },
          ...
        }
      },
      "consensus2": {
        ...
      }
    },
    ...
  },
  "dual_passing_stats": {
    "HLA-A": {
      "is_passing": true,
      "is_dual": true,
      "counts1": 27,
      "counts2": 10,
      "maf": 0.2702702702702703,
      "cdf": 0.019406414321609413
    },
    ...
  }
}
```
