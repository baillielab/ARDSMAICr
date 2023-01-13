
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ARDSMAICR <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/JonathanEMillar/ARDSMAICr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JonathanEMillar/ARDSMAICr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/JonathanEMillar/ARDSMAICr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/JonathanEMillar/ARDSMAICr?branch=main)
<!-- badges: end -->

ARDSMAICR is a package containing helper functions and data useful to
analyse the results of a Meta-analysis by Information Content of ARDS
whole-genome studies.

[MAIC](https://github.com/baillielab/maic) was developed in the [Baillie
Lab](https://baillielab.net), Roslin Institute, University of Edinburgh
as a method for combining lists of genes arising from diverse
experimental sources.

It has been used to study the host response to
[Influenza](https://doi.org/10.1038/s41467-019-13965-x) and
[SARS-CoV-2](https://doi.org/10.1038/s41586-020-03065-y).

MAIC consistently out performs similar algorithms in the case of ranked
and unranked data sources and in the presence of high quality
heterogeneity. Further details can be found
[here](https://doi.org/10.1093/bioinformatics/btac621).

The source code for MAIC can be found
[here](https://github.com/baillielab/maic).

<br />

## ARDS Genomics Systematic Review

<br />

The ARDS Genomics Systematic Review is designed to identify original
studies of human host genomics in Acute Respiratory Distress Syndrome.

We include whole-genome studies, regardless of methodology or the
specific research question. The primary intention is to recover lists of
genes implicated in the host response.

Details of the systematic review and the review protocol can be found on
the study [website](https://ardsmaic.site44.com).

<br />

## ARDS Genomics MAIC

<br />

The results of our systematic review are passed to MAIC.

The output is a ranked list of genes, each with an associated MAIC
score. The list summarises the total information supporting the
association of any given gene with the host response to ARDS.

Our current systematic review and MAIC covers the period
**1<sup>st</sup> January 1967 to 1<sup>st</sup> December 2022**.

<br />

## ARDSMAICR

<br />

### Installation

You can install the development version of ARDSMAICR from
[GitHub](https://github.com/) with:

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("JonathanEMillar/ARDSMAICr")
```

<br />

### Getting started

``` r
library(ARDSMAICR)
```

<br />

ARDSMAICR includes:

- “Data” which provides the results of the most up to date ARDS Genomics
  systematic review and MAIC.

- “Helper functions” which assist in analysing these data.

<br />

#### Data

1.  `data_genes`

The standard output of the [MAIC
algorithm](https://github.com/baillielab/maic).

It has the following format:

| gene   | study_1\_id | study_2\_id | study_n\_id | maic_score | contributors                                |
|--------|-------------|-------------|-------------|------------|---------------------------------------------|
| GENE_A | 0.000       | 1.234       | 0.000       | 1.234      | METHOD_1: study_2\_id                       |
| GENE_B | 1.345       | 1.234       | 0.000       | 1.456      | METHOD_1: study_2\_id, METHOD_2: study_id_1 |

2.  `data_study`

Summary data for studies included by the ARDS Genomics systematic
review.

3.  `data_contributions`

Study contributions to MAIC derived from the
`contributions_calculation()` function.

4.  `data_biolitmine`

Results of a [BioLitMine](https://www.flyrnai.org/tools/biolitmine/web/)
search for the MeSH “Acute Respiratory Distress Syndrome, Adult”.

5.  `data_covidmaic`

A data frame of the genes identified by a MAIC of
[SARS-CoV-2](https://doi.org/10.1038/s41586-020-03065-y).

<br />

#### Helper functions

These functions are provided to aid the analysis of MAIC. Most can be
used for any output in the MAIC format.
