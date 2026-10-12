# Multi-trait selection defines a plant-architecture ideotype for soybean in the Russian Far East

[![Paper](https://img.shields.io/badge/Paper-10.3389%2Ffpls.2026.1922888-blue.svg)](https://doi.org/10.3389/fpls.2026.1922888)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21792792.svg)](https://doi.org/10.5281/zenodo.21792792)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Analysis code and data for the study of plant-architecture traits in 63 soybean
(*Glycine max* (L.) Merr.) accessions evaluated over three growing seasons (2018–2020) at the
Far Eastern Research Institute of Agriculture, Vostochnoye, Khabarovsk Krai, Russia.

Six traits were recorded — plant height, first pod height, number of branches, productive
nodes on the main stem, productive nodes on lateral branches, and days to maturity. The
analysis partitions their phenotypic variance, estimates genetic parameters and genotypic
values, describes the multivariate structure of the panel, and selects genotypes towards a
regional architecture ideotype using the multi-trait genotype-ideotype distance index (MGIDI)
and the multi-trait stability index (MTSI).

Every table and figure in the paper is produced by the scripts in this repository from the
data file included here.

## Contents

```
data/
  Architecture_63genotypes_2018_2020_clean.csv   plot-level data, 492 observations
  rp5_daily.csv                                  daily weather record, 2018-2020
Scripts/
  01_data_preparation.R                          load and verify the dataset
  02_weather_figure.R                            Figure 1
  03_anova_and_descriptive_stats.R               Tables 1-2, Figure 2, model diagnostics
  04_genetic_parameters_and_blups.R              Table 3, Figures 3-4, cross-validation
  05_correlations_pca_clustering.R               Figures 5-6
  06_mgidi_mtsi_selection.R                      Tables 4-5, Figures 8-9
  07_factor_structure_chord.R                    Figure 7
output/
  figures/                                       generated, not tracked
  tables/                                        generated, not tracked
```

## Data

`Architecture_63genotypes_2018_2020_clean.csv` holds 492 plot-level observations: 63
accessions × 3 years, with two or three replications per accession-year cell (75 cells with
two, 114 with three; harmonic mean 2.503). Columns are `Year`, `Replication`, `Genotype` and
the six traits. Heights are in centimetres, days to maturity in days from sowing, and the
branch and node traits are counts.

The accessions were retained from a larger multi-year trial by restricting to 2018–2020 and
keeping only those recorded in at least two replications in every one of the three years.
Accession names, VIR catalogue numbers and countries of origin are given in Supplementary
Table S3 of the paper.

`rp5_daily.csv` holds the daily weather record for the Khabarovsk station from the Reliable
Prognosis archive (rp5.ru), covering 1 May to 30 September of each year: date, number of
observations per day, mean, minimum and maximum air temperature (°C), and precipitation (mm).

## Running the analysis

```bash
git clone https://github.com/yawarhabib/soybean-architecture-ideotype.git
```

Set the working directory to the root of the repository — the folder containing `data/`,
`Scripts/` and `output/` — and run the scripts in order:

```r
source("Scripts/01_data_preparation.R")
source("Scripts/02_weather_figure.R")
source("Scripts/03_anova_and_descriptive_stats.R")
source("Scripts/04_genetic_parameters_and_blups.R")
source("Scripts/05_correlations_pca_clustering.R")
source("Scripts/06_mgidi_mtsi_selection.R")
source("Scripts/07_factor_structure_chord.R")
```

Scripts 01 to 06 are self-contained: each reads the data file directly and can be run on its
own. Script 07 draws the factor loadings written by script 06, so it must be run after it.
Figures are written to `output/figures/` as 600-dpi LZW-compressed TIFF and tables to
`output/tables/` as CSV; both folders are created on the first run.

A full run takes about four minutes, most of it in script 06, where the two multi-trait
indices each fit a mixed model for all six traits.

## Requirements

R 4.4.0 or later and the following packages:

| | | |
|---|---|---|
| car 3.1.3 | emmeans 1.10.6 | psych 2.4.12 |
| lme4 1.1.36 | lmerTest 3.1.3 | metan 1.19.0 |
| dplyr 1.1.4 | tidyr 1.3.1 | lubridate 1.9.4 |
| ggplot2 4.0.3 | ggpubr 0.6.0 | ggrepel 0.9.6 |
| patchwork 1.3.0 | GGally 2.2.1 | circlize 0.4.16 |
| factoextra 1.0.7 | cluster 2.1.8 | zoo 1.8.12 |

The versions listed are those used to produce the published results. To install everything:

```r
install.packages(c("car", "circlize", "cluster", "dplyr", "emmeans", "factoextra", "GGally",
                   "ggplot2", "ggpubr", "ggrepel", "lme4", "lmerTest", "lubridate", "metan",
                   "patchwork", "psych", "tidyr", "zoo"))
```

## Notes on the analysis

The combined analysis of variance uses Type III sums of squares under the sum-to-zero
constraint, because the replication is unbalanced and the genotype × year interaction is
significant for every trait.

Genetic parameters come from a mixed model with genotype, year and their interaction as
random effects. Equation (1) of the paper also carries a replication-within-year term; that
variance is estimated at or near zero for every trait and the fit is singular, so the reported
parameters come from the reduced model. Script 04 fits both and prints the comparison.

The gap statistic in script 05 is resampled, so a seed is set to make that figure
reproducible. Everything else in the analysis is deterministic.

## Citation

If you use the paper, the code or the data, please cite the paper:

> Habib Y, Aseeva T, Shepel O and Potokina E (2026) Multi-trait selection defines a
> plant-architecture ideotype for soybean in the Russian Far East. *Front. Plant Sci.*
> 17:1922888. doi: [10.3389/fpls.2026.1922888](https://doi.org/10.3389/fpls.2026.1922888)

## Contact

Yawar Habib — Center for Bio and Medical Technologies, Skolkovo Institute of Science and
Technology, Moscow, Russia. y.habib@skoltech.ru — ORCID
[0009-0008-1812-3910](https://orcid.org/0009-0008-1812-3910)

## License

Released under the MIT License. See [LICENSE](LICENSE).
