<h1 align="center"><img width="300px" src="images/logo_pb-StarPhase.svg"/></h1>

<h1 align="center">pb-StarPhase</h1>

<p align="center">A phase-aware pharmacogenomic diplotyper for PacBio sequencing data</p>

***

The pb-StarPhase tool will diplotype pharmacogenomic (PGx) genes from [PacBio](https://www.pacb.com/technology/) datasets.
Key features include:

* Ability to create a database from latest CPIC and IMGTHLA information
* Ability to diplotype most genes from CPIC and as well as _HLA-A_, _HLA-B_, and _CYP2D6_
* Works on PacBio datasets from targeted and whole genome sequencing

Authors: [Matt Holt](https://github.com/holtjma), [John Harting](https://github.com/jrharting), [Zev Kronenberg](https://github.com/zeeev)

## Early release warning
Please note that pb-StarPhase is in early development. 
We are still tweaking the input and output file formats and making changes that can affect the behavior of the program.

## Availability
* [Latest release with binary](https://github.com/PacificBiosciences/pb-StarPhase/releases/latest)
* [Auxiliary data files](./data)

## Documentation
* [Supported genes](docs/supported_genes.md)
* [Installation instructions](docs/install.md)
* [User guide with quickstart](docs/user_guide.md)
* [Output files](docs/user_guide.md#output-files)
* [Methods](docs/methods.md)
* [Performance](docs/performance.md)

## Citation
If you use StarPhase, please cite our bioRxiv pre-print:

[Holt, J. M. et al. (2024). StarPhase: Comprehensive Phase-Aware Pharmacogenomic Diplotyper for Long-Read Sequencing Data. _bioRxiv_, 2024-12.](https://doi.org/10.1101/2024.12.10.627527)

## Need help?
If you notice any missing features, bugs, or need assistance with analyzing the output of pb-StarPhase, 
please don't hesitate to open a GitHub issue.

## Support information
pb-StarPhase is a pre-release software intended for research use **only** and not for use in diagnostic procedures. 
While efforts were made to ensure that pb-StarPhase lives up to the quality that PacBio strives for, we make no warranty regarding this software.

As pb-StarPhase is **not** covered by any service level agreement or the like, please do not contact a PacBio Field Applications Scientists or PacBio Customer Service for assistance with any pb-StarPhase release. 
Please report all issues through GitHub instead. 
We make no warranty that any such issue will be addressed, to any extent or within any time frame.

### DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
