# Performance
Table of Contents:
* [Summary results](#results)
* [Tool definitions](#tool-definitions)

## Results
### PharmCAT comparison
We tested 64 targeted HiFi datasets that were sequenced on Revio (40) and Sequel IIe (24) systems.
For each dataset, we ran both PharmCAT and pb-StarPhase to calculate the diplotypes.
Additionally, we ran each dataset on an unphased VCF (DeepVariant) and one that was phased (DeepVariant + pbsv + HiPhase).
The following table shows the comparison for all CPIC genes.
_RNR1_ was excluded because it is not called by PharmCAT.
_HLA-A_, _HLA-B_, and _CYP2D6_ were excluded because they are not called by either tool.

| Gene | Identical (Unphased) | Identical (phased) |
| -- | -- | -- |
| _ABCG2_ | 100% | 100% |
| _CACNA1S_ | 100% | 100% |
| _CFTR_ | 100% | 100% |
| _CYP2B6_ | 100% | 100.0% |
| _CYP2C19_ | 100% | 100.0% |
| _CYP2C9_ | 100% | 98.4% |
| _CYP3A5_ | 100% | 100% |
| _CYP4F2_ | 100% | 100% |
| _DPYD_ | 82.8% | 68.8% |
| _G6PD_ | 100% | 100% |
| _IFNL3_ | 100% | 100% |
| _NUDT15_ | 100% | 100% |
| _RNR1 (MT)_ | 100% | 100% |
| _SLCO1B1_ | 100% | 100% |
| _TPMT_ | 100% | 100% |
| _UGT1A1_ | 100% | 100% |
| _VKORC1_ | 100% | 100% |
| Overall | 99.0% | 98.1% |

We manually curated the remaining differences and determined they are due to either differences in reporting or mishandling of phased alleles in PharmCAT:

* _CYP2C9_ - There was one total discrepancies in the phased solutions. In this case, PharmCAT reported "Unknown/Unknown" while pb-StarPhase reported a single diplotype. Interestingly, PharmCAT reported the same answer as pb-StarPhase when unphased variants were provided. When we inspected the variants, it appeared that PharmCAT was erroneously treating variants in different phase blocks as being in phase with each other. This created haplotypes that have not been described before and led to the "Unknown/Unknown" reports.
* _DPYD_ - This gene is a little different from the other genes in that traditional haplotypes have not been built into CPIC yet. Instead, each "haplotype" is just a single variant. 
  * Unphased - If more than two variants are identified, PharmCAT handles these by reporting all variants as "{variant}/None". In contrast, pb-StarPhase will report "NO_MATCH/NO_MATCH" as there is no special treatment for _DPYD_. All reported diplotypes with <=2 variants were identical between PharmCAT and pb-StarPhase.
  * Phased - The same phase block issue with _CYP2C9_ impacts the _DPYD_ results. In addition to the mismatches from 3+ variants, phase blocks can get misinterpreted in PharmCAT leading to reduced identity scores in this gene.

In summary, we expect most unphased results to be identical between the two tools, with _DPYD_ reporting being the notable exception.
With phased HiFi data, pb-StarPhase is more likely to generate correct results as it properly handles variants that are on different phase blocks.
However, we note from our experiment that these situations are relatively rare.

## Tool definitions
### pb-StarPhase
This is the tool released with this repository.

* Version - v0.7.2
* Database version - v0.6.1-20230914
* Extra parameters - None

### PharmCAT
The PharmCAT pipeline (`/pharmcat/pharmcat_pipeline`) was run from the provided docker image:

* Version - v2.7.1
* Image - `docker://pgkb/pharmcat:2.7.1`
* Extra parameters - `--matcher-all-results --missing-to-ref --reporter-save-json`
