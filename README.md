# Genetic shifts in viral populations: how sudden vs gradual changes in host species composition affect the stability and dynamics of genetic variability in viral populations

### Description
R scripts to analyze data of evolving populations of Sindbis virus under different rates of host replacement.

### Overview
To run the analysis, processed data (approxwf results and allele frequency tables in .xlsx format) must be obtained from the [Zenodo repository](https://zenodo.org/) and the contents of the tar.gz file must be extracted here.

#### Contents
##### source_me.R
Loading of the data and main functions.

##### approxwf_analysis.R
Plot and analyze population parameters estimated with approxwf. Look for files in approxwfout/.

##### allele_freq_diff.R
Analysis of mean allele frequency differences and richness through passages.

##### fluctuation_scaling.R
Analysis of allele fluctuations across passages. Fits to Taylor's law and related analyses.

##### hurst.R
Analysis of persistent behavior on the placement of mutations.

The analysis was performed on the following R environment:

```
> sessionInfo()
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin22.6.0
Running under: macOS Ventura 13.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /opt/homebrew/Cellar/r/4.4.1/lib/R/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Madrid
tzcode source: internal

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] factoextra_1.0.7       cluster_2.1.7          Hmisc_5.2-1            marginaleffects_0.24.0 ggbreak_0.1.2          gridExtra_2.3         
 [7] statcomp_0.1.0         pracma_2.4.4           minpack.lm_1.2-4       igraph_2.1.2           ggrepel_0.9.6          viridis_0.6.5         
[13] viridisLite_0.4.2      ggpubr_0.6.0           segmented_2.1-3        nlme_3.1-166           MASS_7.3-61            readxl_1.4.3          
[19] writexl_1.5.1          ggplot2_3.5.1          tidyr_1.3.1            dplyr_1.1.4           

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1   farver_2.1.2       fastmap_1.2.0      coop_0.6-3         bayestestR_0.15.0  digest_0.6.37      rpart_4.1.23       estimability_1.5.1
 [9] lifecycle_1.0.4    magrittr_2.0.3     compiler_4.4.1     rlang_1.1.4        tools_4.4.1        utf8_1.2.4         data.table_1.16.4  knitr_1.49        
[17] ggsignif_0.6.4     labeling_0.4.3     htmlwidgets_1.6.4  aplot_0.2.3        abind_1.4-8        withr_3.0.2        foreign_0.8-87     purrr_1.0.2       
[25] datawizard_0.13.0  nnet_7.3-19        fansi_1.0.6        xtable_1.8-4       colorspace_2.1-1   emmeans_1.10.5     scales_1.3.0       insight_1.0.0     
[33] mvtnorm_1.3-2      cli_3.6.3          rmarkdown_2.29     generics_0.1.3     rstudioapi_0.17.1  parameters_0.24.0  stringr_1.5.1      splines_4.4.1     
[41] effectsize_1.0.0   ggplotify_0.1.2    cellranger_1.1.0   base64enc_0.1-3    yulab.utils_0.1.8  vctrs_0.6.5        Matrix_1.7-1       carData_3.0-5     
[49] car_3.1-3          gridGraphics_0.5-1 patchwork_1.3.0    rstatix_0.7.2      Formula_1.2-5      htmlTable_2.4.3    glue_1.8.0         cowplot_1.1.3     
[57] stringi_1.8.4      gtable_0.3.6       munsell_0.5.1      tibble_3.2.1       pillar_1.9.0       htmltools_0.5.8.1  R6_2.5.1           evaluate_1.0.1    
[65] lattice_0.22-6     backports_1.5.0    broom_1.0.7        ggfun_0.1.8        Rcpp_1.0.13-1      coda_0.19-4.1      checkmate_2.3.2    mgcv_1.9-1        
[73] xfun_0.49          fs_1.6.5           pkgconfig_2.0.3      
```

