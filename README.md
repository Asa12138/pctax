
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pctax <img src="man/figures/pctax.png" align="right" width="120" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/Asa12138/pctax/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Asa12138/pctax/actions/workflows/R-CMD-check.yaml)
[![](https://img.shields.io/badge/doi-waiting-yellow.svg)](https://doi.org/waiting)
[![](https://img.shields.io/badge/blog-@asa-blue.svg)](https://asa-blog.netlify.app/)
[![](http://cranlogs.r-pkg.org/badges/grand-total/pctax)](https://cran.r-project.org/package=pctax)
[![](http://cranlogs.r-pkg.org/badges/last-month/pctax)](https://cran.r-project.org/package=pctax)
[![](https://www.r-pkg.org/badges/version/pctax?color=green)](https://cran.r-project.org/package=pctax)
[![](https://img.shields.io/badge/devel%20version-0.1.2-green.svg)](https://github.com/Asa12138/pctax)
<!-- badges: end -->

`pctax` provides a comprehensive suite of tools for analyzing omics
data.

The HTML documentation of the latest version is available at [Github
page](https://asa12138.github.io/pctax/).

Please go to <https://bookdown.org/Asa12138/pctax_book/> for the full
vignette.

## Installation

The stable version is available on CRAN:

    install.packages("pctax")

Or you can install the development version of pctax from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Asa12138/pctax")
```

## 🚀 NEWS 🚀

Recently, I added a function to draw two trees and their relationships:

``` r
data(otutab, package = "pcutils")
df2tree(taxonomy[1:50, ]) -> tax_tree
df2tree(taxonomy[51:100, ]) -> tax_tree2
link <- data.frame(from = sample(tax_tree$tip.label, 20), to = sample(tax_tree2$tip.label, 20))
plot_two_tree(tax_tree, tax_tree2, link,
  tree1_tip = T, tree2_tip = T,
  tip1_param = list(size = 2), tip2_param = list(size = 2)
)
```

<div class="figure">

<img src="man/figures/README-unnamed-chunk-2-1.png" alt="Two trees and their relationships" width="100%" />
<p class="caption">
Two trees and their relationships
</p>

</div>

Recently I added a function to plot element cycling because element
cycling genes are important in the microbiome (especially the
environmental microbiome). Supports simple cycle diagram drawing of C,
N, P, S, Fe (manual arrangement, there must be some missing parts, will
be continuously added in the future):

``` r
plot_element_cycle(cycle = "Nitrogen cycle")
#> recommend ggsave(width = 12,height = 10)
```

<div class="figure">

<img src="man/figures/README-unnamed-chunk-3-1.png" alt="Nitrogen cycle" width="100%" />
<p class="caption">
Nitrogen cycle
</p>

</div>

## Usage

For the full vignette, please visit [pctax: Analyzing Omics Data with
R](https://bookdown.org/Asa12138/pctax_book/).

**Some Functionalities of `pctax`:**

![](man/figures/pctax1.png)

## Citation

Please cite:

Chen Peng, Chao Jiang (2023). *pctax: Professional Comprehensive
Microbiome Data Analysis Pipeline*. R package,
<https://github.com/Asa12138/pctax>.
