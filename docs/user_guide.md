# User Guide
Table of contents:

* [Quickstart](#quickstart)
* [Supported upstream processes](#supported-upstream-processes)
* [Output files](#output-files)
* [FAQ](#faq)

# Quickstart
The following command will generate a single output JSON file with diplotype calls:

```bash
pharmgoat diplotype \
    --database {DATABASE} \
    --reference {REFERENCE} \
    --vcf {IN_VCF} \
    --output-calls {OUT_RESULTS}
```

Parameters:
* `--database {DATABASE}` - the path to the PGx allele database; these can be found in our provided [Data](../data) or [generated with a separate command](#how-do-i-update-the-database-for-pharmgoat)
* `--reference {REFERENCE}` - a FASTA file containing the reference genome, gzip allowed
* `--vcf {IN_VCF}` - path to a VCF file containing the variants for PharmGOAT to convert into diplotypes; any absent variants are assumed to be homozygous for the reference allele; the VCF must be indexed (e.g., `{path}.vcf.gz` with `{path}.vcf.gz.tbi`)
* `--output-calls {OUT_RESULTS}` - path to where the summarize output will be written

## Quickstart Example
This example also includes some auxiliary outputs that may be useful for tracking phase result statistics.
For details on each, refer to the [Output files](#output-files) section.

```
pharmgoat diplotype \
    --database ./data/v0.6.1/cpic_20230914.json \
    --reference ./reference/human_GRCh38_no_alt_analysis_set.fasta \
    --vcf ./HG001/hiphase/HG001.GRCh38.hiphase.vcf.gz \
    --output-calls ./output/HG001.pharmgoat.json
        
[2023-09-15T13:25:28.530Z INFO  pharmgoat::cli::diplotype] Input database: "./data/v0.6.1/cpic_20230914.json"
[2023-09-15T13:25:28.530Z INFO  pharmgoat::cli::diplotype] Input reference: "./reference/human_GRCh38_no_alt_analysis_set.fasta"
[2023-09-15T13:25:28.530Z INFO  pharmgoat::cli::diplotype] Input VCF: "./HG001/hiphase/HG001.GRCh38.hiphase.vcf.gz"
[2023-09-15T13:25:28.530Z INFO  pharmgoat::cli::diplotype] Output diplotype calls: "./HG001/pharmgoat/HG001.pharmgoat.json"
[2023-09-15T13:25:28.530Z INFO  pharmgoat] Loading PGx database from "/data/v0.6.1/cpic_20230914.json"...
[2023-09-15T13:25:28.716Z INFO  pharmgoat] Loading reference genome from "./reference/human_GRCh38_no_alt_analysis_set.fasta"...
[2023-09-15T13:25:58.847Z INFO  pharmgoat::diplotyper] Solving ABCG2...
[2023-09-15T13:25:59.241Z INFO  pharmgoat::diplotyper] Solving CACNA1S...
[2023-09-15T13:25:59.338Z INFO  pharmgoat::diplotyper] Solving CFTR...
[2023-09-15T13:25:59.527Z INFO  pharmgoat::diplotyper] Solving CYP2B6...
[2023-09-15T13:25:59.734Z INFO  pharmgoat::diplotyper] Solving CYP2C19...
[2023-09-15T13:25:59.910Z INFO  pharmgoat::diplotyper] Solving CYP2C9...
[2023-09-15T13:26:00.152Z INFO  pharmgoat::diplotyper] Solving CYP3A5...
[2023-09-15T13:26:00.280Z INFO  pharmgoat::diplotyper] Solving CYP4F2...
[2023-09-15T13:26:00.399Z INFO  pharmgoat::diplotyper] Solving DPYD...
[2023-09-15T13:26:00.614Z INFO  pharmgoat::diplotyper] Solving G6PD...
[2023-09-15T13:26:00.614Z WARN  pharmgoat::diplotyper] Error while normalizing database variant 778641: At chrX:154532438, provided reference allele has "G" but reference genome has "A", PgxVariant { name: "c.1311C>T", dbsnp_id: Some("rs2230037"), position: 154532439, alleles: [Some("G"), Some("A")] }
[2023-09-15T13:26:00.614Z WARN  pharmgoat::diplotyper] Ignoring "Mediterranean Haplotype" due to variant incompatibility.
[2023-09-15T13:26:01.034Z INFO  pharmgoat::diplotyper] Solving IFNL3...
[2023-09-15T13:26:01.131Z INFO  pharmgoat::diplotyper] Solving MT-RNR1...
[2023-09-15T13:26:01.131Z WARN  pharmgoat::diplotyper] Error while normalizing database variant 828080: At chrM:929, provided reference allele has "A" but reference genome has "G", PgxVariant { name: "930A>G", dbsnp_id: None, position: 930, alleles: [Some("A"), Some("G")] }
[2023-09-15T13:26:01.131Z WARN  pharmgoat::diplotyper] Ignoring "930A>G" due to variant incompatibility.
[2023-09-15T13:26:01.235Z INFO  pharmgoat::diplotyper] Solving NUDT15...
[2023-09-15T13:26:01.412Z INFO  pharmgoat::diplotyper] Solving RYR1...
[2023-09-15T13:26:01.571Z INFO  pharmgoat::diplotyper] Solving SLCO1B1...
[2023-09-15T13:26:01.768Z INFO  pharmgoat::diplotyper] Solving TPMT...
[2023-09-15T13:26:01.994Z INFO  pharmgoat::diplotyper] Solving UGT1A1...
[2023-09-15T13:26:02.103Z INFO  pharmgoat::diplotyper] Solving VKORC1...
[2023-09-15T13:26:02.204Z INFO  pharmgoat] Saving diplotypes to "./HG001/pharmgoat/HG001.pharmgoat.json"
[2023-09-15T13:26:02.450Z INFO  pharmgoat] Process finished successfully.
```

# Supported upstream processes
The following upstream processes are supported as inputs to PharmGOAT:

* Data types
  * PacBio whole genome sequencing
  * PacBio PGx panel
* Variant callers
  * [DeepVariant](https://github.com/google/deepvariant) - for SNV/indel
* Phasers
  * [HiPhase](https://github.com/PacificBiosciences/HiPhase) - for phased inputs; phased inputs are generally recommended as this can reduce or remove ambiguity in the diplotype assignment

Other upstream data types and processes may work with PharmGOAT, but there is no official support for them at this time.

# Output files
## Output call file
The output call file is a JSON file containing both the diplotype call for each gene as well as the supported information.
Fields are described below, with a partial example further down:

* `pharmgoat_version` - the version of PharmGOAT that generated the output file; this will match the version from `pharmgoat -V`
* `database_metadata` - copy of the database metadata describing how the database was generated
  * `pharmgoat_version` - the version of PharmGOAT that generated the database
  * `cpic_version` - a tag indicating the CPIC version; the API does not currently provide a global version tag, so it is instead labeled as `API-{build_time}`
  * `build_time` - the time this database was built; UTC
* `gene_details` - the core output; each gene will have key in this dictionary with the following information:
  * `diplotypes` - a list of all exact matching diplotype combinations; ambiguous combinations (i.e., lack of phase information) may lead to more than one possible diplotype combination
    * `hap1` - name of the first haplotype in the diplotype
    * `hap2` - name of the second haplotype in the diplotype
    * `diplotype` - name of the combined diplotype; `{hap1}/{hap2}`
  * `variant_details` - contains the list of all identified variants from the VCF file that match a variant definition
    * `cpic_variant_id` - the CPIC assigned variant identifier (unsigned integer)
    * `cpic_name` - the CPIC assigned name for this variant; this is typically a human-readable identifier that has been historically used in pharmacogenomic literature
    * `dbsnp` - a DBSNP identifier, if available
    * `normalized_variant` - describes the identified variant after variant normalization; this will almost always match VCF specification for describing a variant (left-shifted, no redundant bases); this may not match the exact CPIC definition due oddities in how CPIC describes alleles
      * `chrom` - the chromosome from the reference
      * `position` - the **0-based** position along the chromosome
      * `reference` - the reference sequence at this position
      * `alternate` - the alternate (variant) sequence at this position
    * `normalized_genotype` - describes the associated genotype after genotype normalization
      * `genotype` - this will almost always match the genotype (GT) from the VCF file; multi-allelic sites will get converted into one of the following: `[0/0, 0/1, 0|1, 1|0, 1/1]`
      * `phase_set` - a phase set ID (PS) from the VCF file; this will be `null` for unphased or homozygous variants; variants on different phase sets are _not_ considered phased with each other (they are effectively, unphased)

Partial example:
```
{
  "pharmgoat_version": "0.6.1-c952b5e",
  "database_metadata": {
    "pharmgoat_version": "0.6.1-c0d01e4",
    "cpic_version": "API-2023-09-14T20:48:01.821179769Z",
    "build_time": "2023-09-14T20:48:01.821179769Z"
  },
  "gene_details": {
    ...
    "CYP4F2": {
      "diplotypes": [
        {
          "hap1": "*5",
          "hap2": "*4",
          "diplotype": "*5/*4"
        }
      ],
      "variant_details": [
        {
          "cpic_variant_id": 778236,
          "cpic_name": "34T>G",
          "dbsnp": "rs3093105",
          "normalized_variant": {
            "chrom": "chr19",
            "position": 15897577,
            "reference": "A",
            "alternate": "C"
          },
          "normalized_genotype": {
            "genotype": "0|1",
            "phase_set": 15871890
          }
        },
        ...
      ]
    }
    ...
  }
}
```

## PharmCAT ingestible TSV
PharmGOAT will generate diplotype calls, but it does not assign phenotypes, activity scores, etc. to those calls.
PharmCAT is a tool that can [accept outside calls](https://pharmcat.org/using/Outside-Call-Format/) and add them to a report with further interpretation.
The `--pharmcat-tsv {filename}` option can be used in PharmGOAT to generate this file.
PharmGOAT will generate a basic two-column TSV file with one column for the gene and one column for the diplotype call.
If ambiguity is present, it will report "Multiple/Multiple" in place of the diplotype.
If no matching diplotype is generated, it will report "NO_MATCH/NO_MATCH" in place of the diplotype.
An example partial/truncated output is below:

```
#gene	diplotype
CYP2B6	Multiple/Multiple
CYP4F2	*3/*4
DPYD	NO_MATCH/NO_MATCH
MT-RNR1	Reference
...
```

# FAQ
## How do I update the database for PharmGOAT?
A new database can be generated automatically using the available data with the following command:

```bash
pharmgoat build \
  --output-json {path_to_new_database}.json
```

This requires an internet connection that can query the CPIC API for the latest genes, allele definitions, and variants.

## Why are some of the haplotypes ignored?
Some genes, like _G6PD_, include reference variants that are alternate sequences relative to the GRCh38 reference genome.
This can cause some haplotypes, like _G6PD_'s "Mediterranean Haplotype", to require variants that would be a REF allele in a standard VCF file.
Unfortunately, this can generate ambiguity between an implied homozygous reference call (i.e., "0/0") and a lack of coverage of that variant site.
If these haplotypes are not ignored, there can be systematic errors in the downstream reporting due to the assumption of REF allele when no VCF entry is identified.
Any haplotype that is ignored this way is reported as a warning in the PharmGOAT stderr.
As of v0.6.1, there are only two ignored haplotypes in the PharmGOAT database.
Future work may try to resolve this ambiguity.
