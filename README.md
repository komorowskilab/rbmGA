# Rule-Based Modeling for Gene Expression Data Analysis

## Description

A robust systems biology approach to represent biomolecular interactions,
including gene regulation, as a series of logical "if-then" rules instead
of relying on rigid, fully quantified kinetic equations. This method is
particularly effective at addressing the combinatorial complexity found in
immune system and gene expression networks.

## Prerequisites

The method is implemented and tested in R, so it requires **R version 4.4.1
or later**.

The rule-based modeling using R.ROSETTA is compatible with both UNIX and
Windows operating systems. However, UNIX systems (such as macOS or Linux)
require a 32-bit version of [Wine](https://www.winehq.org/), a free and
open-source compatibility layer. Please note that macOS Catalina and later
versions no longer support 32-bit applications.

## Development Setup and Installation

To install the package, run the following in R:

```r
install.packages("devtools")
devtools::install_github("GirishPulinkala/rbmGA")
```

To reproduce the results and check the analysis:

```bash
git clone git@github.com:GirishPulinkala/rbmGA.git
cd analysis
```

### Example Script

Make sure you are in the `analysis` folder:

```bash
pwd
cd R.ROSETTA_analysis
Rscript external_validation.R
```

## Citation

If you find this package useful, please cite us:

```bibtex
@article{,
  title={},
  author={},
  journal={},
  year={},
  pages={},
  doi={},
  url={}
}
```




