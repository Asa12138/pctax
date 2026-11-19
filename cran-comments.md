## Test environments

* local macOS, R 4.5.1
* R-devel failure of 0.1.7: https://www.r-project.org/nosvn/R.check/r-devel-linux-x86_64-debian-gcc/pctax-00check.html

## R CMD check results

There were no ERRORs or WARNINGs.

Local `--as-cran` notes come from this machine, not the package:

* unable to verify the current time
* HTML Tidy is not recent enough to validate the manual

## Comments

This is a resubmission. Version 0.1.7 failed on R-devel in the `plot_two_tree` example because `treelayout` 0.1.0 wrote `NaN` y-coordinates for trees with singleton nodes. `plot_two_tree()` now recomputes those coordinates with `ape::node.depth(method = 2)`.

Braces in the taxonkit help pages are escaped.
