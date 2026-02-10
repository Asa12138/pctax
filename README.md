

<!-- README.md is generated from README.qmd. Please edit that file -->
<!-- badges: start -->

[![](https://img.shields.io/badge/blog-@asa-blue.svg)](https://asa-blog.netlify.app/)
[![](http://cranlogs.r-pkg.org/badges/grand-total/pctax)](https://cran.r-project.org/package=pctax)
[![](http://cranlogs.r-pkg.org/badges/last-month/pctax)](https://cran.r-project.org/package=pctax)
[![](https://www.r-pkg.org/badges/version/pctax?color=green)](https://cran.r-project.org/package=pctax)
[![](https://img.shields.io/badge/devel%20version-0.1.8-green.svg)](https://github.com/Asa12138/pctax)

<!-- badges: end -->

# pctax

`pctax` provides a comprehensive suite of tools for analyzing omics
data.

## Install

``` r
install.packages("devtools")
devtools::install_github("Asa12138/pcutils")
devtools::install_github("Asa12138/pctax")
```

## 🚀 NEWS 🚀

Recently I added a function to plot element cycling because element
cycling genes are important in the microbiome (especially the
environmental microbiome). Supports simple cycle diagram drawing of C,
N, P, S, Fe (manual arrangement, there must be some missing parts, will
be continuously added in the future):

``` r
plot_element_cycle(cycle = "Nitrogen cycle")
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the pctax package.
#>   Please report the issue at <https://github.com/Asa12138/pctax/issues>.
#> recommend ggsave(width = 12,height = 10)
```

![Nitrogen cycle](README_files/figure-commonmark/unnamed-chunk-3-1.png)

## Usage

For the full vignette, please visit [pctax: Analyzing Omics Data with
R](https://bookdown.org/Asa12138/pctax_book/).

**Some Functionalities of `pctax`:**

``` mermaid
flowchart LR
  B(pctax)--> C{Functionalities}
  C --> D[Visualization]
  C --> E[Diversity Analysis]
  C --> F[Differential Abundance Analysis]
  C --> G[Community Assembly Analysis]
  C --> H[Functional Enrichment Analysis]
  C --> I[Network Analysis]
  C --> J[Elemental Cycling Analysis]
```

## Citation

Please cite:

Chen Peng, Chao Jiang (2023). *pctax: Professional Comprehensive
Microbiome Data Analysis Pipeline*. R package,
<https://github.com/Asa12138/pctax>.
