---
title: "Quick start: pctax"
output: 
  prettydoc::html_pretty:
    keep_md: true
    toc: true
    toc_depth: 3
    theme: architect
vignette: >
  %\VignetteIndexEntry{pctax}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---





# pctax

`pctax` provides a comprehensive suite of tools for analyzing microbiome data.

## Install


```r
install.packages("devtools")
devtools::install_github("Asa12138/pcutils", dependencies = T)
devtools::install_github("Asa12138/pctax", dependencies = T)
```

## 🚀 NEWS 🚀

Recently I added a function to plot element cycling because element cycling genes are important in the microbiome (especially the environmental microbiome). Supports simple cycle diagram drawing of C, N, P, S, Fe (manual arrangement, there must be some missing parts, will be continuously added in the future):


```r
plot_element_cycle(cycle = "Nitrogen cycle")
#> Loading required namespace: ggchicklet
#> recommend ggsave(width = 12,height = 10)
```

![Nitrogen cycle](README_files/figure-html/unnamed-chunk-3-1.png)

## Usage
It includes functionalities for α-diversity analysis, β-diversity analysis, differential analysis, community assembly, visualization of phylogenetic tree and functional enrichment analysis... 

Look at the test data:

```r
library(pctax)
library(pcutils)
data(otutab, package = "pcutils")
# help(otutab)

head(otutab)
#>                               NS1  NS2  NS3  NS4  NS5  NS6  WS1  WS2  WS3  WS4
#> s__un_f__Thermomonosporaceae 1092 1920  810 1354 1064 1070 1252 1597 1330  941
#> s__Pelomonas_puraquae        1962 1234 2362 2236 2903 1829  644  495 1230 1284
#> s__Rhizobacter_bergeniae      588  458  889  901 1226  853  604  470 1070 1028
#> s__Flavobacterium_terrae      244  234 1810  673 1445  491  318 1926 1493  995
#> s__un_g__Rhizobacter         1432  412  533  759 1289  506  503  590  445  620
#>                               WS5  WS6  CS1  CS2  CS3  CS4  CS5  CS6
#> s__un_f__Thermomonosporaceae 1233 1011 2313 2518 1709 1975 1431 1527
#> s__Pelomonas_puraquae         953  635 1305 1516  844 1128 1483 1174
#> s__Rhizobacter_bergeniae      846  670 1029 1802 1002 1200 1194  762
#> s__Flavobacterium_terrae      577  359 1080 1218  754  423 1032 1412
#> s__un_g__Rhizobacter          657  429 1132 1447  550  583 1105  903
#>  [ reached 'max' / getOption("max.print") -- omitted 1 rows ]
head(metadata)
#>      Id Group     env1     env2     env3      env4     env5        env6
#> NS1 NS1    NS 3.057248 10.23571 5.554576  8.084997 25.00795 -1.15456682
#> NS2 NS2    NS 4.830219 11.13453 5.613455  8.556829 16.67690  0.81168745
#> NS3 NS3    NS 3.753133 10.06232 5.582916 10.226572 21.68926  1.40733211
#> NS4 NS4    NS 4.262264 10.84401 5.258419  9.002256 24.81046  1.47805320
#> NS5 NS5    NS 2.476135  7.52584 6.255314  9.357587 19.70553  0.05813095
#> NS6 NS6    NS 5.131004 10.82761 5.180966  8.141506 18.39021 -1.70032569
#>          lat     long
#> NS1 26.94422 103.4767
#> NS2 29.08733 109.6196
#> NS3 28.25164 104.0361
#> NS4 33.82415 106.8651
#> NS5 33.51011 105.4571
#> NS6 31.86864 102.7832
head(taxonomy)
#>                                  Kingdom            Phylum
#> s__un_f__Thermomonosporaceae k__Bacteria p__Actinobacteria
#> s__Pelomonas_puraquae        k__Bacteria p__Proteobacteria
#> s__Rhizobacter_bergeniae     k__Bacteria p__Proteobacteria
#> s__Flavobacterium_terrae     k__Bacteria  p__Bacteroidetes
#> s__un_g__Rhizobacter         k__Bacteria p__Proteobacteria
#> s__un_o__Burkholderiales     k__Bacteria p__Proteobacteria
#>                                               Class               Order
#> s__un_f__Thermomonosporaceae      c__Actinobacteria  o__Actinomycetales
#> s__Pelomonas_puraquae         c__Betaproteobacteria  o__Burkholderiales
#> s__Rhizobacter_bergeniae     c__Gammaproteobacteria  o__Pseudomonadales
#> s__Flavobacterium_terrae          c__Flavobacteriia o__Flavobacteriales
#> s__un_g__Rhizobacter         c__Gammaproteobacteria  o__Pseudomonadales
#> s__un_o__Burkholderiales      c__Betaproteobacteria  o__Burkholderiales
#>                                                Family
#> s__un_f__Thermomonosporaceae   f__Thermomonosporaceae
#> s__Pelomonas_puraquae               f__Comamonadaceae
#> s__Rhizobacter_bergeniae          f__Pseudomonadaceae
#> s__Flavobacterium_terrae         f__Flavobacteriaceae
#> s__un_g__Rhizobacter              f__Pseudomonadaceae
#> s__un_o__Burkholderiales     f__un_o__Burkholderiales
#>                                                     Genus
#> s__un_f__Thermomonosporaceae g__un_f__Thermomonosporaceae
#> s__Pelomonas_puraquae                        g__Pelomonas
#> s__Rhizobacter_bergeniae                   g__Rhizobacter
#> s__Flavobacterium_terrae                g__Flavobacterium
#> s__un_g__Rhizobacter                       g__Rhizobacter
#> s__un_o__Burkholderiales         g__un_o__Burkholderiales
#>                                                   Species
#> s__un_f__Thermomonosporaceae s__un_f__Thermomonosporaceae
#> s__Pelomonas_puraquae               s__Pelomonas_puraquae
#> s__Rhizobacter_bergeniae         s__Rhizobacter_bergeniae
#> s__Flavobacterium_terrae         s__Flavobacterium_terrae
#> s__un_g__Rhizobacter                 s__un_g__Rhizobacter
#> s__un_o__Burkholderiales         s__un_o__Burkholderiales
```

### α-diversity analysis

Calculate a_diversity of otutab then link to experiment group or environment variable.

```r
a_diversity(otutab) -> a_res
#> Loading required namespace: vegan
plot(a_res, "Group", metadata)
```

![α-diversity](README_files/figure-html/a-diversity-1.png)

```r
plot(a_res, "env1", metadata)
#> Loading required namespace: ggpmisc
```

![α-diversity](README_files/figure-html/a-diversity-2.png)

### β-diversity analysis

There are a range of dimensionality reduction methods available for analysis, including Constrained and non-Constrained.

Like PCA, PCoA, NMDS, RDA, CCA... For example:

PCA:

```r
b_analyse(otutab, method = "pca") -> b_res
#> Loading required namespace: ade4
#> four dataframes in a list, 1 is eig, 2 is sample_site, 3 is var, 4 is var contribution
plot(b_res, "Group", metadata, bi = T, rate = 0.5)
#> Loading required namespace: ggnewscale
#> Loading required namespace: ggrepel
#> Loading required namespace: RColorBrewer
```

![PCA for β-diversity](README_files/figure-html/b-diversity-1.png)

```r
plot(b_res, "Group", metadata, mode = 3)
```

![PCA for β-diversity](README_files/figure-html/b-diversity-2.png)

RDA:

```r
env <- metadata[, 6:10]
# RDA
myRDA(otutab, env) -> phy.rda
#> 
#> Call:
#> vegan::decorana(veg = dat.h) 
#> 
#> Detrended correspondence analysis with 26 segments.
#> Rescaling of axes with 4 iterations.
#> Total inertia (scaled Chi-square): 0.3207 
#> 
#>                         DCA1    DCA2    DCA3    DCA4
#> Eigenvalues          0.03127 0.02265 0.01916 0.01729
#> Additive Eigenvalues 0.03127 0.02265 0.01917 0.01727
#> Decorana values      0.03150 0.02146 0.01701 0.01035
#> Axis lengths         0.74268 0.74498 0.57253 0.52361
#> 
#> DCA analysis, select the sorting analysis model according to the first value of the Axis lengths row
#>    Axis Lengths >4.0-CCA (based on unimodal model, canonical correspondence analysis);
#>    If it is between 3.0-4.0 - both RDA/CCA;
#>    If less than 3.0-RDA (based on linear model, redundancy analysis)
#> [1] "===============Initial Model================"
#> [1] "Initial cca, vif>20 indicates serious collinearity:"
#>     env4     env5     env6      lat     long 
#> 2.574997 2.674671 1.252002 1.381839 1.211392 
#> Initial Model R-square: 0.04817859 
#> [1] "=============Statistics==========="
#> 0.3281261 Constrained indicates the degree to which environmental factors explain differences in community structure
#> 0.6718739 unconstrained means that the environmental factors cannot explain the part of the community structure
RDA_plot(phy.rda, "Group", metadata)
```

![RDA for β-diversity associated environmental variables](README_files/figure-html/rda-1.png)


### Differential analysis

There are also lots of statistic methods for differential analysis:
ALDEX, ANCOM2, randomForest, t.test, wilcox.test... or deseq2, limma...(Commonly used in transcriptome)

```r
diff_da(otutab, metadata["Group"]) -> res
#> Loading required namespace: DESeq2
#> converting counts to integer mode
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
volcano_p(res)
```

![Volcano plot of differential analysis](README_files/figure-html/diff-1.png)

```r
volcano_p(res, mode = 2)
#> Using compare as id variables
#> Using compare as id variables
```

![Volcano plot of differential analysis](README_files/figure-html/diff-2.png)

### Community assembly

Community assembly in microbiome refers to the processes that shape the composition, diversity, and structure of microbial communities in a particular environment or host. 
Microbiome consist of diverse microbial populations that interact with each other and their surroundings, and understanding how these communities assemble is crucial for comprehending their ecological dynamics and functional implications.


```r
ncm(otutab) -> ncm_res
#> Loading required namespace: Hmisc
#> Loading required namespace: minpack.lm
plot(ncm_res)
#> Loading required namespace: ggpubr
#> Loading required namespace: patchwork
#> Scale for x is already present.
#> Adding another scale for x, which will replace the existing scale.
```

![NCM model](README_files/figure-html/ncm-1.png)
### Rtaxonkit

[Taxonkit](https://github.com/shenwei356/taxonkit) is a Practical and Efficient NCBI Taxonomy Toolkit.[(Shen, W. et al.  Journal of Genetics and Genomics 2021.)](https://www.sciencedirect.com/science/article/abs/pii/S1673852721000837)

We recommend you download this excellent software to help next analysis. Or you can use Taxonkit in R by pctax interface as followed:


```r
# 1. This function help you install suitable version taxonkit
install_taxonkit()
# taxonkit has been successfully installed!

# 2. Then download the NCBI taxonomy database.
download_taxonkit_dataset()
# Taxonkit files downloaded and copied successfully.

# 3. Check whether taxonkit is ready
check_taxonkit()
# ==============Taxonkit is available if there is help message above==============
# =========================Taxonkit dataset is available!=========================
```

Then you can use taxonkit in R just like in terminal.

```r
?taxonkit_lineage

# taxonkit_list
# taxonkit_reformat
# taxonkit_name2taxid
# taxonkit_filter
# taxonkit_lca
```


### Phylogenetic tree

```r
ann_tree(taxonomy, otutab) -> tree
#> Loading required namespace: ggtree
#> Joining with `by = join_by(label)`
easy_tree(tree, add_abundance = FALSE)
#> Loading required namespace: ggtreeExtra
```

![Phylogenetic tree](README_files/figure-html/unnamed-chunk-6-1.png)

## Cite
Please cite:

Chen Peng, Chao Jiang (2023). _pctax: Professional Comprehensive Microbiome Data Analysis Pipeline_. R package, <https://github.com/Asa12138/pctax>.
