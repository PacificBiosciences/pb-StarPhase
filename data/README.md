# Auxiliary data files
## Database files
Databases were pre-generated for ease-of-use, sharing, and backwards compatibility when using older versions of pb-StarPhase.
Though we will do our best to avoid breaking changes, database files are **not** guaranteed to work on versions that do not match.
Additionally, these databases represent a snapshot in time of upstream data sources (e.g., CPIC).
Running the same command at a later date may produce a different database with updated annotations.

Each file is labeled as `{version}/cpic_{YYYYMMDD}.json` and represents a run of the following command using the specified `{version}` of pb-StarPhase on the corresponding date (`{YYYYMMDD}`):

```bash
pbstarphase build \
    --output-db {version}/cpic_{YYYYMMDD}.json
```
