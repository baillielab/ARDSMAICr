
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ARDSMAICR <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/JonathanEMillar/ARDSMAICr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JonathanEMillar/ARDSMAICr/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/JonathanEMillar/ARDSMAICr/branch/main/graph/badge.svg?token=5ILIXWJ542)](https://codecov.io/gh/JonathanEMillar/ARDSMAICr)

<!-- badges: end -->

The ARDSMAICR package contains the results of a meta-analysis by
information content (MAIC) of genome-wide studies of the host response
to acute respiratory distress syndrome (ARDS). These data are
accompanied by a range of helper functions useful for analysis.

<br />

## Meta-Analysis by Information Content (MAIC)

<br />

[MAIC](https://github.com/baillielab/maic) was developed in the [Baillie
Lab](https://baillielab.net), Roslin Institute, University of Edinburgh
as a method for combining lists of genes arising from diverse
experimental sources.

It has been used to study the host response to
[Influenza](https://doi.org/10.1038/s41467-019-13965-x) and
[SARS-CoV-2](https://doi.org/10.1038/s41586-020-03065-y).

MAIC consistently out performs similar algorithms in the case of ranked
and unranked data sources and in the presence of heterogeneity in the
quality of studies. Further details can be found
[here](https://doi.org/10.1093/bioinformatics/btac621).

The source code for MAIC is hosted
[here](https://github.com/baillielab/maic).

<br />

## ARDSMAICR

<br />

### Installation

You can install the latest version of ARDSMAICR from
[GitHub](https://github.com/) with:

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("baillielab/ARDSMAICR")
```

<br />

### Getting started

``` r
library(ARDSMAICR)
```

<br />

### ARDS MAIC results

The systematic review, data extraction, and MAIC are described in detail
on the study [website](https://ardsmaic.site44.com).

The package contains the output of the MAIC covering the period
**1<sup>st</sup> January 1967 to 1<sup>st</sup> December 2022**. Future
releases will match regular updates of the systematic review.

These data are contained in `data_genes`. It has the following format:

| gene   | study_1_id | study_2_id | study_n_id | maic_score | contributors                               |
|--------|------------|------------|------------|------------|--------------------------------------------|
| GENE_A | 0.000      | 1.234      | 0.000      | 1.234      | METHOD_1: study_2_id                       |
| GENE_B | 1.345      | 1.234      | 0.000      | 1.456      | METHOD_1: study_2_id, METHOD_2: study_id_1 |

Additional data are:

- `data_study` - A summary of the studies identified by the systematic
  review, including their methods.
- `data_contributions` - Calculated study contributions to MAIC.
- `data_biolitmine` - The results of a [BioLitMine]() search for the
  MeSH “Respiratory Distress Syndrome, Adult”.
- `data_covidmaic` - A ranked list of genes from a MAIC of [SARS-CoV-2
  studies]().

<br />

### Helper functions

Several groups of functions useful in the analysis of the MAIC results
are included. They fall into the following broad families:

- Summary tables
- Helper functions - counts, genes, lists, methods and categories
- Gene prioritisation
- Information content and contribution

The majority of these functions can be applied to the standard
[MAIC](https://github.com/baillielab/maic) output of any analysis.
