# pollinatordistanceTSF_SRMA

This repository contains the data and R scripts for our systematic review and meta-analysis on the relationship between distance to natural habitat and pollination services in tropical smallholder farming landscapes.

## Author

Ennia Bosshard

## Repository Contents

The repository is organised in the following structure:

- **`raw_data/`** - Original datasets for the abundance, richness, and fruitset meta-analysis used in the analysis
- **`code/`** - R scripts for data processing and analysis for each of the three meta-analysis
- **`outputs/`** - Processed results, figures, and summary tables, structured into three subfolders:
  - `abundance/` - Contains all outputs related to the pollinator abundance meta-analysis
  - `richness/` - Contains all outputs related to the pollinator richness meta-analysis
  - `fruitset/` - Contains all outputs related to the fruitset meta-analysis

In the main repository, you will also find:  
- **`README.md`** (this file)  

## Dependencies

This analysis was conducted in **R** (version 4.2.2). The following packages are required:  

```r
install.packages(c("ggplot2", "dplyr", "MASS", "lme4", "metafor", "ggeffects", "glmmTMB", "here"))
```

## Variable Explanations

- **study**: Unique study abbreviation/ID for each study
- **report**: Unique abbreviation/ID for each individual report, i.e. published or unpublished article
- **authors**: Author and year for each study (e.g., `author et al. year`)
- **crop**: Common name of the study crop
- **p_dependency**: Pollinator-dependency level of each study crop, following the information from Siopa et al. (2024)
- **farm_size**: Total size of farm in hectares (ha) or nearest classification of farm size available (e.g. 'small <2ha' or 'medium 2-15ha')
- **habitat_description**: Description of semi-natural habitat (e.g., 'forest', 'savannah', 'grassland')
- **habitat_type**: Categorised type of natural habitat based on original description (binary; either 'natural forest' or 'other)
- **agr_intensity**: Agricultural intensity (e.g., 'high' = intense use of pesticides, monocultures; 'low' = minimal pesticide use and no monocultures)
- **distance_m**: Euclidean distance (in metres) to the nearest natural habitat
- **distance_measure**: Whether distance measurement was reported in the original study ('reported') or estimated using remote sensing ('estimated')
- **abundance_all**: Observed number of pollinators visiting crop flowers for all observed insects
- **abundance_wild**: Observed number of wild pollinators excluding managed honeybees (A. mellifera, and in some cases A. cerana)
- **richness_all**: Observed number of species for all pollinators
- **richness_wild**: Observed number of species for wild pollinators excluding managed honeybees (A. mellifera, and in some cases A. cerana)
- **fruitset**: Observed proportion (0-1) of flowers that set fruit
- **treatment**: Fruit set treatment (only used data from 'open' treatments for unbagged/untreated flowers in our synthesis)
- **pollinator**: Focal pollinator taxa of the study to the most precise taxonomic level (e.g., 'Insecta', 'Hymenoptera', 'Diptera', 'bees', 'stingless bees', etc.)
- **sites**: Total number of study sites/farms in the study
- **location**: Name or code for each study site/farm
- **study_design**: Description of the study design for each dataset; either 'single distance per site', 'nested design' or 'paired sites'
- **repeat_measures**: Number of visits/replicate observations for active pollinator observations, number of pan traps for passive sampling, number of plants, branches or flowers observed for fruit set
- **sampling_method**: Method for insect observations (e.g., 'active' for focal observations/sweep netting, 'passive' for pan traps, 'combined' for both methods)
- **sampling_time_min**: Time (in minutes) over which pollinators were sampled per sampling event
- **total_sampling_time_min**: Total time (in minutes) over which pollinators were sampled at each study site/data point (accummulative across total sampling effort)
- **source**: Where the data was obtained from (e.g., authors or publication directly, open access data, or from figures)

## References

Siopa, C., L. G. Carvalheiro, H. Castro, J. Loureiro, and S. Castro. 2024. Animal-pollinated crops and cultivars—A quantitative assessment of pollinator dependence values and evaluation of methodological approaches. *Journal of Applied Ecology*, 61:1279-1288.

## R Session Info

R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.6.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/London
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] glmmTMB_1.1.12      ggeffects_2.3.0     here_1.0.1          metafor_4.6-0       numDeriv_2016.8-1.1 metadat_1.2-0       lme4_1.1-35.5       Matrix_1.7-1       
 [9] MASS_7.3-61         dplyr_1.1.4         ggplot2_3.5.1      

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       TMB_1.9.16         insight_1.2.0      lattice_0.22-6     mathjaxr_1.6-0     vctrs_0.6.5        tools_4.4.2        Rdpack_2.6.2       generics_0.1.4    
[10] datawizard_1.0.2   sandwich_3.1-1     tibble_3.3.0       pkgconfig_2.0.3    RColorBrewer_1.1-3 lifecycle_1.0.4    compiler_4.4.2     farver_2.1.2       textshaping_1.0.0 
[19] codetools_0.2-20   pillar_1.11.0      nloptr_2.1.1       crayon_1.5.3       reformulas_0.4.1   boot_1.3-31        multcomp_1.4-26    nlme_3.1-166       tidyselect_1.2.1  
[28] mvtnorm_1.3-2      labeling_0.4.3     forcats_1.0.0      splines_4.4.2      rprojroot_2.0.4    grid_4.4.2         cli_3.6.5          magrittr_2.0.3     survival_3.7-0    
[37] utf8_1.2.6         TH.data_1.1-2      withr_3.0.2        scales_1.4.0       estimability_1.5.1 emmeans_1.10.7     ragg_1.3.3         zoo_1.8-12         hms_1.1.3         
[46] coda_0.19-4.1      haven_2.5.5        rbibutils_2.3      mgcv_1.9-1         rlang_1.1.6        Rcpp_1.1.0         xtable_1.8-4       glue_1.8.0         rstudioapi_0.17.1 
[55] minqa_1.2.8        R6_2.6.1           systemfonts_1.2.1 
> 
