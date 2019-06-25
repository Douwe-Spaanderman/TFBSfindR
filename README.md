# TFBSfindR
#### Documentation
TFBSfindR investigates changes between Transcription Factor (TFs) binding sites in single nucleotide polymorphisms. TFBSfindR uses several methods to calculate changes in TF motifs, first calculating Position Weight Matrixes scores for both reference and variant. Secondly, two methods can be assessed for the identification of important variants, suggested by Kumasaka et al. 2019 and Touzet et al. 2007, TFBSfindR also comes with some plot functions for single variant score analysis.

For an indepth overview on how to use `TFBSfindR`, check out our vignette [TFBSfindR vignette](https://github.com/Douwe-Spaanderman/TFBSfindR/tree/master/vignettes/TFBS.findR_vignette.pdf).

#### Install

```{r}
install.packages("devtools")
devtools::install_github("Douwe-Spaanderman/TFBSfindR")
```
