---
author: Kyoko Yamaguchi
date: 23 November 2026
title: R package "methylationhelper"
---

# Overview
This code repository holds code for installing and using the custom R package "methylationhelper".
This package contains functions that are useful for DNA methylation array analysis.

Install the latest version of the package with:

```r
# install.packages("remotes")
remotes::install_github("yamachan18/methylationhelper")
```


# Structure
* code: Code used for development/testing (scratch.Rmd) and code used for illustrating examples (TBD) are found here.
* reports: The corresponding knitted HTML/PDF documents if applicable.
* R: R scripts holding each custom function of the methylationhelper package.
* man: manual pages for each R function auto-rendered by roxygen2.
* data: This directory holds data used to illustrate the workflow with the custom functions. GSE171424.rda and GSE174666.rda were taken from The Gene Expression Ombibus (GEO) are are used for the query and reference datasets, respectively.
<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171424">GSE171424</a>
<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174666">GSE174666</a>

