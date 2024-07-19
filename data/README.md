# Auxiliary data files
## Database files
Databases were pre-generated for ease-of-use, sharing, and backwards compatibility when using older versions of pb-StarPhase.
Though we will do our best to avoid breaking changes, database files are **not** guaranteed to work on versions that do not match.
Additionally, these databases represent a snapshot in time of upstream data sources (e.g., CPIC, IMGTHLA, and PharmVar).
Running the same command at a later date may produce a different database with updated annotations.

As of v0.12.0, each file is labeled as `{version}/pbstarphase_{YYYYMMDD}.json.gz` and represents a run of the following command using the specified `{version}` of pb-StarPhase on the corresponding date (`{YYYYMMDD}`):

```bash
pbstarphase build \
    --output-db {version}/pbstarphase_{YYYYMMDD}.json.gz
```

# Data sources and citations
## CPIC data citations
All CPIC variants are sourced via the [CPIC API](https://cpicpgx.org/api-and-database/).
As of this writing, versioning is *not* available via the API, so we instead tag the creation date in the `database_metadata/cpic_version` for all database files.
CPIC data is released under [CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/) license.
If you use the CPIC data within our database files, the authors and maintainers of CPIC ask that you follow their instructions for citing the relevant publications and resources: [https://cpicpgx.org/license/](https://cpicpgx.org/license/).

## HLA data citations
All HLA sequences are sourced and unmodified from [https://github.com/ANHIG/IMGTHLA](https://github.com/ANHIG/IMGTHLA).
The exact git tag version is stored in the `database_metadata/hla_version` for all database versions after v0.8.0.
If you use the HLA sequences within our database files, the authors and maintainers of IMGTHLA ask that you cite the following:

1. Robinson J, Barker DJ, Georgiou X, Cooper MA, Flicek P, Marsh SGE: IPD-IMGT/HLA Database. Nucleic Acids Research (2020), 48:D948-55
2. Robinson J, Malik A, Parham P, Bodmer JG, Marsh SGE: IMGT/HLA - a sequence database for the human major histocompatibility complex Tissue Antigens (2000), 55:280-287

For more details, see the [IMGTHLA license file](https://github.com/ANHIG/IMGTHLA/blob/Latest/LICENCE.md).

## PharmVar data citations
All CYP2D6 sequences are sourced and unmodified from PharmVar's "Download Gene Data" ZIP file on [https://www.pharmvar.org/gene/CYP2D6](https://www.pharmvar.org/gene/CYP2D6).
The version is stored in the `database_metadata/pharmvar_version` for all database versions after v0.9.0.
PharmVar data is released under a [Creative Commons Attribution-ShareAlike 4.0 International license](https://creativecommons.org/licenses/by-sa/4.0/) license, see [https://www.pharmvar.org/terms-and-conditions](https://www.pharmvar.org/terms-and-conditions) for more details.
If you use the CYP2D6 sequences within our database files, the authors and maintainers of PharmVar ask that you cite the relevant publications from their website: [https://www.pharmvar.org/publications](https://www.pharmvar.org/publications).
