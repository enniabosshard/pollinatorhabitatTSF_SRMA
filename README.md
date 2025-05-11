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
install.packages(c("metafor", "MASS", "ggplot2", "dplyr", "here"))
```

## Variable Explanations

- **study**: Unique study abbreviation/ID for each study
- **report**: Unique abbreviation/ID for each individual report, i.e. published or unpublished article
- **authors**: Author and year for each study (e.g., `author et al. year`)
- **crop**: Common name of the study crop
- **p_dependency**: Pollinator-dependency level of each study crop, following the information from Siopa et al. (2024)
- **farm_size**: Total size of farm in hectares (ha)
- **habitat**: Definition/description of semi-natural habitat (e.g., 'forest', 'savannah', 'grassland')
- **agr_intensity**: Agricultural intensity (e.g., 'high' = intense use of pesticides, monocultures; 'low' = minimal pesticide use and no monocultures)
- **distance_m**: Euclidean distance (in metres) to the nearest natural habitat
- **distance_measure**: Whether distance measurement was reported in the original study ('reported') or estimated using remote sensing ('estimated')
- **abundance_all**: Observed number of pollinators visiting crop flowers for all observed insects
- **abundance_wild**: Observed number of wild pollinators excluding managed honeybees (A. mellifera, and in some cases A. cerana)
- **richness_all**: Observed number of species for all pollinators
- **richness_wild**: Observed number of species for wild pollinators excluding managed honeybees
- **fruitset**: Observed proportion (0-1) of flowers that set fruit
- **treatment**: Fruit set treatment ('open' for unbagged/untreated flowers)
- **pollinator**: Focal pollinator taxa of the study to the most precise taxonomic level (e.g., 'Insecta', 'Hymenoptera', 'Diptera', 'bees', 'stingless bees', etc.)
- **sites**: Total number of study sites/farms in the study
- **location**: Name or code for each study site/farm
- **sampling_effort**: Sample size (e.g., number of visits/replicate observations for active pollinator observations, number of pan traps for passive sampling, number of flowers counted for fruit set)
- **sampling_method**: Method for insect observations (e.g., 'active' for focal observations/sweep netting, 'passive' for pan traps, 'combined' for both methods)
- **sampling_time_min**: Time (in minutes) over which pollinators were sampled per sampling event
- **total_sampling_time_min**: Total time (in minutes) over which pollinators were sampled at each study site/data point (accummulative across total sampling effort)
- **source**: Where the data was obtained from (e.g., authors or publication directly, open access data, or from figures)

## References

Siopa, C., L. G. Carvalheiro, H. Castro, J. Loureiro, and S. Castro. 2024. Animal-pollinated crops and cultivars—A quantitative assessment of pollinator dependence values and evaluation of methodological approaches. *Journal of Applied Ecology*, 61:1279-1288.

