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