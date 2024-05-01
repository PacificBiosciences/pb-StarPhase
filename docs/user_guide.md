# User Guide
Table of contents:

* [Quickstart](#quickstart)
* [Common uses cases](#common-use-cases)
* [Supported upstream processes](#supported-upstream-processes)
* [Output files](#output-files)
* [FAQ](#faq)

# Quickstart
The following command will generate a single output JSON file with diplotype calls:

```bash
pbstarphase diplotype \
    --database {DATABASE} \
    --reference {REFERENCE} \
    --vcf {IN_VCF} \
    --output-calls {OUT_RESULTS}
```

Parameters:
* `--database {DATABASE}` - The path to the PGx allele database; these can be found in our provided [Data](../data) or [generated with a separate command](#how-do-i-update-the-database-for-pb-starphase).
* `--reference {REFERENCE}` - A FASTA file containing the reference genome; gzip allowed.
* `--vcf {IN_VCF}` - Path to a VCF file containing the variants for pb-StarPhase to convert into diplotypes; any absent variants are assumed to be homozygous for the reference allele; the VCF must be indexed (e.g., `{path}.vcf.gz` with `{path}.vcf.gz.tbi`).
* `--output-calls {OUT_RESULTS}` - Path to where the summarize output will be written.

## Quickstart Example
This example also includes some auxiliary outputs that may be useful for tracking phase result statistics.
For details on each, refer to the [Output files](#output-files) section.

```
pbstarphase diplotype \
    --database ./data/v0.6.1/cpic_20230914.json \
    --reference ./reference/human_GRCh38_no_alt_analysis_set.fasta \
    --vcf ./HG001/hiphase/HG001.GRCh38.hiphase.vcf.gz \
    --output-calls ./output/HG001.pbstarphase.json
        
[2023-09-15T13:25:28.530Z INFO  pbstarphase::cli::diplotype] Input database: "./data/v0.6.1/cpic_20230914.json"
[2023-09-15T13:25:28.530Z INFO  pbstarphase::cli::diplotype] Input reference: "./reference/human_GRCh38_no_alt_analysis_set.fasta"
[2023-09-15T13:25:28.530Z INFO  pbstarphase::cli::diplotype] Input VCF: "./HG001/hiphase/HG001.GRCh38.hiphase.vcf.gz"
[2023-09-15T13:25:28.530Z INFO  pbstarphase::cli::diplotype] Output diplotype calls: "./HG001/pbstarphase/HG001.pbstarphase.json"
[2023-09-15T13:25:28.530Z INFO  pbstarphase] Loading PGx database from "/data/v0.6.1/cpic_20230914.json"...
[2023-09-15T13:25:28.716Z INFO  pbstarphase] Loading reference genome from "./reference/human_GRCh38_no_alt_analysis_set.fasta"...
[2023-09-15T13:25:58.847Z INFO  pbstarphase::diplotyper] Solving ABCG2...
[2023-09-15T13:25:59.241Z INFO  pbstarphase::diplotyper] Solving CACNA1S...
[2023-09-15T13:25:59.338Z INFO  pbstarphase::diplotyper] Solving CFTR...
[2023-09-15T13:25:59.527Z INFO  pbstarphase::diplotyper] Solving CYP2B6...
[2023-09-15T13:25:59.734Z INFO  pbstarphase::diplotyper] Solving CYP2C19...
[2023-09-15T13:25:59.910Z INFO  pbstarphase::diplotyper] Solving CYP2C9...
[2023-09-15T13:26:00.152Z INFO  pbstarphase::diplotyper] Solving CYP3A5...
[2023-09-15T13:26:00.280Z INFO  pbstarphase::diplotyper] Solving CYP4F2...
[2023-09-15T13:26:00.399Z INFO  pbstarphase::diplotyper] Solving DPYD...
[2023-09-15T13:26:00.614Z INFO  pbstarphase::diplotyper] Solving G6PD...
[2023-09-15T13:26:00.614Z WARN  pbstarphase::diplotyper] Error while normalizing database variant 778641: At chrX:154532438, provided reference allele has "G" but reference genome has "A", PgxVariant { name: "c.1311C>T", dbsnp_id: Some("rs2230037"), position: 154532439, alleles: [Some("G"), Some("A")] }
[2023-09-15T13:26:00.614Z WARN  pbstarphase::diplotyper] Ignoring "Mediterranean Haplotype" due to variant incompatibility.
[2023-09-15T13:26:01.034Z INFO  pbstarphase::diplotyper] Solving IFNL3...
[2023-09-15T13:26:01.131Z INFO  pbstarphase::diplotyper] Solving MT-RNR1...
[2023-09-15T13:26:01.131Z WARN  pbstarphase::diplotyper] Error while normalizing database variant 828080: At chrM:929, provided reference allele has "A" but reference genome has "G", PgxVariant { name: "930A>G", dbsnp_id: None, position: 930, alleles: [Some("A"), Some("G")] }
[2023-09-15T13:26:01.131Z WARN  pbstarphase::diplotyper] Ignoring "930A>G" due to variant incompatibility.
[2023-09-15T13:26:01.235Z INFO  pbstarphase::diplotyper] Solving NUDT15...
[2023-09-15T13:26:01.412Z INFO  pbstarphase::diplotyper] Solving RYR1...
[2023-09-15T13:26:01.571Z INFO  pbstarphase::diplotyper] Solving SLCO1B1...
[2023-09-15T13:26:01.768Z INFO  pbstarphase::diplotyper] Solving TPMT...
[2023-09-15T13:26:01.994Z INFO  pbstarphase::diplotyper] Solving UGT1A1...
[2023-09-15T13:26:02.103Z INFO  pbstarphase::diplotyper] Solving VKORC1...
[2023-09-15T13:26:02.204Z INFO  pbstarphase] Saving diplotypes to "./HG001/pbstarphase/HG001.pbstarphase.json"
[2023-09-15T13:26:02.450Z INFO  pbstarphase] Process finished successfully.
```

# Common use cases
## HLA and _CYP2D6_ diplotyping
With v0.10.0, pb-StarPhase supports diplotyping of _HLA-A_, _HLA-B_, and _CYP2D6_ from an aligned BAM file.
If using targeted sequencing datasets, see [our recommended parameters](#can-i-diplotype-using-targeted-sequencing-data).
To enable HLA and _CYP2D6_ diplotyping, simply provide the BAM file(s) in addition to the normal parameters.
Both HLA and _CYP2D6_ diplotyping is more computationally expensive than the CPIC genes.
If run-time is an issue, we recommend using the `--threads` option to provide additional cores to StarPhase, which will improve the HLA diplotyping components.

```bash
pbstarphase diplotype \
    --bam ${BAM} \
    --threads ${THREADS} \
    ...
```

# Supported upstream processes
The following upstream processes are supported as inputs to pb-StarPhase:

* Data types
  * PacBio whole genome sequencing
  * PacBio PGx panel - _CYP2D6_ currently unsupported
* Reference genomes
  * GRCh38 - we recommend the [human_GRCh38_no_alt_analysis_set](https://github.com/PacificBiosciences/reference_genomes/tree/main/reference_genomes/human_GRCh38_no_alt_analysis_set)
* Aligners (BAM files):
  * [pbmm2](https://github.com/PacificBiosciences/pbmm2) (recommended)
  * [minimap2](https://github.com/lh3/minimap2)
* Variant callers
  * [DeepVariant](https://github.com/google/deepvariant) - For SNV/indel.
* Phasers
  * [HiPhase](https://github.com/PacificBiosciences/HiPhase) - For phased inputs; phased inputs are generally recommended as this can reduce or remove ambiguity in the diplotype assignment.

Other upstream data types and processes may work with pb-StarPhase, but there is no official support for them at this time.

# Output files
## Output call file
The output call file is a JSON file containing both the diplotype call for each gene as well as the supported information.
Fields are described below, with a partial example further down:

* `pbstarphase_version` - The version of pb-StarPhase that generated the output file; this will match the version from `pbstarphase -V`.
* `database_metadata` - Copy of the database metadata describing how the database was generated:
  * `pbstarphase_version` - The version of pb-StarPhase that generated the database.
  * `cpic_version` - A tag indicating the CPIC version; the API does not currently provide a global version tag, so it is instead labeled as `API-{build_time}`.
  * `hla_version` - A tag indicating the IMGTHLA version; these are identical to GitHub tags on [https://github.com/ANHIG/IMGTHLA](https://github.com/ANHIG/IMGTHLA).
  * `build_time` - The time this database was built; in UTC format.
* `gene_details` - The core output; each gene will have a key in this dictionary with the following information:
  * `diplotypes` - A list of all exact matching diplotype combinations; ambiguous combinations (i.e., lack of phase information) may lead to more than one possible diplotype combination.
    * `hap1` - Name of the first haplotype in the diplotype.
    * `hap2` - Name of the second haplotype in the diplotype.
    * `diplotype` - Name of the combined diplotype; `{hap1}/{hap2}`.
  * `variant_details` - Contains the list of all identified variants from the VCF file that match a variant definition. Will be `null` for HLA genes.
    * `cpic_variant_id` - The CPIC assigned variant identifier (unsigned integer).
    * `cpic_name` - The CPIC assigned name for this variant; this is typically a human-readable identifier that has been historically used in pharmacogenomic literature.
    * `dbsnp` - A DBSNP identifier, if available.
    * `normalized_variant` - Describes the identified variant after variant normalization; this will almost always match VCF specification for describing a variant (left-shifted, no redundant bases); this may not match the exact CPIC definition due to oddities in how CPIC describes alleles.
      * `chrom` - The chromosome from the reference.
      * `position` - The **0-based** position along the chromosome.
      * `reference` - The reference sequence at this position.
      * `alternate` - The alternate (variant) sequence at this position.
    * `normalized_genotype` - Describes the associated genotype after genotype normalization.
      * `genotype` - This will almost always match the genotype (GT) from the VCF file; multi-allelic sites will get converted into one of the following: `[0/0, 0/1, 0|1, 1|0, 1/1]`.
      * `phase_set` - A phase set ID (PS) from the VCF file; this will be `null` for unphased or homozygous variants. Note that variants on different phase sets are _not_ considered phased with each other (they are effectively, unphased).
  * `mapping_details` - Contains the list of identified HLA read mappings from the BAM files. Will be `null` for non-HLA genes.
    * `read_qname` - The read query name.
    * `best_hla_id` - The best matching HLA ID, does not have to be one of the haplotypes reported in the best diplotype.
    * `best_star_allele` - The HLA star allele corresponding to the `best_hla_id`.
    * `best_mapping_stats` - The length, number of mismatches (nm), and number of unmapped bases for the `best_hla_id` sequences for cDNA and DNA mapped against the read. Note that cDNA/DNA will only be present if used as part of the scoring.
    * `is_ignored` - True if the read was ignored when solving the diplotype due to a high error rate.
  * `multi_mapping_details` - Contains the list of identified _CYP2D6_ read mapping segments from the BAM files. Will be `null` for all other genes.
    * `read_qname` - The read query name.
    * `read_position` - The portion (or slice) of the read corresponding to the allele.
    * `consensus_id` - The assigned consensus ID for this read segment.
    * `consensus_star_allele` - The assigned star allele for this read segment.

Partial example:
```
{
  "pbstarphase_version": "0.8.0-fa50d82",
  "database_metadata": {
    "pbstarphase_version": "0.8.0-79d4679",
    "cpic_version": "API-2023-12-19T16:11:50.938951041Z",
    "hla_version": "v.354.0-alpha",
    "build_time": "2023-12-19T16:11:50.938951041Z"
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
    },
    ...
    "CYP2D6": {
      "diplotypes": [
        {
          "hap1": "*4.001+*68",
          "hap2": "*4.001",
          "diplotype": "*4.001+*68/*4.001"
        }
      ],
      "variant_details": null,
      "mapping_details": null,
      "multi_mapping_details": [
        {
          "read_qname": "mfake/1234/ccs",
          "read_position": {
            "start": 0,
            "end": 2701
          },
          "consensus_id": 1,
          "consensus_star_allele": "CYP2D7"
        },
        ...
      ]
    },
    ...
    "HLA-A": {
      "diplotypes": [
        {
          "hap1": "01:01:01:01",
          "hap2": "26:01:01:01",
          "diplotype": "01:01:01:01/26:01:01:01"
        }
      ],
      "variant_details": null,
      "mapping_details": [
        {
          "read_qname": "mfake/5678/ccs",
          "best_hla_id": "HLA:HLA00073",
          "best_star_allele": "26:01:01:01",
          "best_mapping_stats": {
            "cdna_len": null,
            "cdna_nm": null,
            "cdna_unmapped": null,
            "dna_len": 3517,
            "dna_nm": 4,
            "dna_unmapped": 0
          },
          "is_ignored": false
        },
        ...
      ]
    }
    ...
  }
}
```

## PharmCAT ingestible TSV
pb-StarPhase will generate diplotype calls, but it does not assign phenotypes, activity scores, and so on, to those calls.
PharmCAT is a tool that can [accept outside calls](https://pharmcat.org/using/Outside-Call-Format/) and add them to a report with further interpretation.
The `--pharmcat-tsv {filename}` option can be used in pb-StarPhase to generate this file.
pb-StarPhase will generate a basic two-column TSV file with one column for the gene and one column for the diplotype call.
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
## How do I update the database for pb-StarPhase?
A new database can be generated automatically using the available data with the following command:

```bash
pbstarphase build \
  --output-json {path_to_new_database}.json
```

This requires an internet connection that can query the CPIC API and IMGTHLA GitHub for the latest genes, allele definitions, and variants.

## Why are some of the haplotypes ignored?
Some CPIC genes, like _G6PD_, include reference variants that are alternate sequences relative to the GRCh38 reference genome.
This can cause some haplotypes, like _G6PD_'s "Mediterranean Haplotype", to require variants that would be a REF allele in a standard VCF file.
Unfortunately, this can generate ambiguity between an implied homozygous reference call (i.e., "0/0") and a lack of coverage of that variant site.
If these haplotypes are not ignored, there can be systematic errors in the downstream reporting due to the assumption of REF allele when no VCF entry is identified.
Any haplotype that is ignored this way is reported as a warning in the pb-StarPhase stderr.
As of v0.6.1, there are only two ignored CPIC haplotypes in the pb-StarPhase database.
Future work may try to resolve this ambiguity.

IMGTHLA also includes several incomplete HLA definitions that may be missing DNA sequences, cDNA sequences, or be heavily truncated (e.g., several only include exons 2 and 3).
By default, pb-StarPhase requires a DNA sequence to be present for each HLA haplotype.
Any haplotypes without a DNA sequence are ignored.

## Why am I getting a weird hybrid allele call for _CYP2D6_?
Internally, pb-StarPhase searches for hybrid D6/D7 alleles by generating representative sequences at the intron/exon boundaries and labeling them accordingly.
For example, the allele "CYP2D6::CYP2D7::intron1" indicates an allele that looks like _CYP2D6_ until the start of intron1 (i.e., end of exon1), and then it looks more like _CYP2D7_ (note that these are not precise breakpoints for the hybrid).
Prior to output to JSON, we attempt to resolve these alleles to known star alleles.
For example, both "CYP2D6::CYP2D7::intron1" and "CYP2D6::CYP2D7::exon2" are re-mapped to CYP2D6*68.
While all "CYP2D7::CYP2D6" alleles are currently mapped to *13, most "CYP2D6::CYP2D7" alleles do _not_ have a known re-mapping.
Those without a known re-mapping are left in the pb-StarPhase internal format.
If you encounter an allele that you think should be re-mapped, please open an issue on our GitHub.

## Can I diplotype using targeted sequencing data?
In general yes: in our internal tests the CPIC and HLA genes behave similar to their WGS counterparts.
However, _CYP2D6_ tends to be more difficult to accurately call with targeted sequencing.
This is typically due to shorter read lengths, increased coverage variation across alleles, and full-allele drop out due to the capture.
For _CYP2D6_, this is particular problematic due to the presence of deletion, duplication, and hybrid alleles that may influence the final diplotype.
For targeted sequencing, we recommend using the following _CYP2D6_-specific additional parameters, which attempt to account for these complicating factors: `--infer-connections --normalize-d6-only`.
