# StarPhase database guide
Table of contents:

* [Overview](#overview)
* [Building the database](#building-the-database)
* [Database statistics](#database-statistics)

## Overview
The StarPhase database file is a JSON file (optionally, gzipped) that includes all the necessary information required for StarPhase to generate the diplotype calls.
The required information varies by gene, but it can be broadly categorized as variant-based genes, HLA genes, and CYP2D6.
This information is gathered from multiple databases (CPIC, PharmVar, IMGT/HLA) and consolidated into a single portable file.
We periodically provide public versions that have been built and tested in our [Data folder](../data/).

## Building the database
StarPhase includes the ability to build the database as part of the toolkit.
In general, this will query multiple APIs available by the upstream databases and consolidate the information into a single file.
Unfortunately, not all information is currently available in APIs and is instead hard-coded with fixed values by StarPhase during database creation.
This primarily impacts structural variation (multiple genes) and the hybrid allele definitions (CYP2D6).

A new database can be generated automatically using the available data with the following command:
```bash
pbstarphase build \
  --output-json {path_to_new_database}.json
```

This requires an internet connection that can query the various APIs.
Additionally, this command relies on upstream databases maintaining a known structure.
If that structure changes, this command may fail and require an update to the software to resolve it.
If you encounter and issue with building the database, please open an issue on GitHub so we can investigate it.

## Database configuration
For some genes, there are entries across multiple databases with the greatest overlap between CPIC and PharmVar.
By default, StarPhase will prefer sourcing genes from PharmVar because the database tracks deeper sub-allele information.
The one exception to this default rule is _DPYD_, where StarPhase will prefer the CPIC entries for downstream compatibility with tools like PharmCAT.
These preferences can be overridden by providing a database configuration with the `--build-options` (or `-b`) flag when running `pbstarphase build`.
We provide the current default database configuration in [default_db_config.json](../data/default_db_config.json).

Options for configuration:
* `default_gene_source` (string): Sets the default data source for genes present in both CPIC and PharmVar. Can be `"CPIC"` or `"PharmVar"`. Default is `"PharmVar"`.
* `gene_source_overrides` (object): A map of gene names to their preferred data source (`"CPIC"` or `"PharmVar"`). These overrides take precedence over the `default_gene_source` setting for the specified genes. For example, `{"DPYD": "CPIC"}` forces _DPYD_ to use CPIC data regardless of the `default_gene_source` setting.

## Database statistics
StarPhase includes a `db-stat` function, which will provide high-level summary information about the contents of a provided database JSON file.
Additional per-gene statistics can be output with the `-v` flag.

Here is an example `db-stat` run:
```
Database metadata:
	Version: 1.5.0-b25ef87
	CPIC version: API-2025-10-02T16:35:06.217945084Z
	HLA version: v3.61.0-alpha
	PharmVar version: 6.2.15
	Build time: 2025-10-02 16:35:06.217945084 UTC
Database gene statistics:
	Total genes: 35 = 23 (gene entries) + 11 (HLA) + 1 (CYP2D6)
	Gene entries:
		Total alleles: 1251
		Total variants: 1145
	HLA:
		Total HLA alleles: 40673
	CYP2D6:
		Total alleles: 552
		Total variants: 403
```
