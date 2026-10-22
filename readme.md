# Holocene climate reconstructions for the Western part of the Arid Central Asia (WCA)

## Overview
This repository contains the R code and data for the studies "*Holocene hydroclimate trend in Western Arid Central Asia: evidence of coherent precipitation patterns from the Iranian Plateau to the western Tian Shan piedmonts*" under reviews in *Earth-Science Reviews* in 2027 (Dugerdil et al., in reviews).

**Author**: **Lucas Dugerdil**<sup>1,2,3</sup>

**Affiliations**: 
1. Key Laboratory of Western China's Environmental Systems (Ministry of Education), College of Earth and Environmental Sciences, Lanzhou University, Lanzhou, China
2. Univ. Lyon, ENS de Lyon, Université Lyon 1, CNRS, UMR 5276 LGL-TPE, F-69364, Lyon, France
3. Université de Montpellier, CNRS, IRD, EPHE, UMR 5554 ISEM, Montpellier, France
**ORCID:** [0000-0003-0266-564X](https://orcid.org/0000-0003-0266-564X)

## Description
This R script permits the reproducibility of all data analyses, quantitative reconstructions, and climate simulation extractions. It produce the stacked curves for the five studied regions (Caucasus Moutains, Iranian Plateau, Caspian/Aral Basin, western and eastern Tian Shan areas) for both data and models outputs, assessed their inter-correlations, their linkage with seasonal insolation and the spatial sampling biases over the study area. All the figures provides in the manuscript are produced by this repository, except Figs. 1 and 14 realised in `QGIS`, and Fig.~ 13B processed in `Python`.

## Structure
- `Import/Laskar_2004/` – seasonal insolation at 35°N extracted from Laskar et al. (2004)
- `Import/SIG/` – the spatial boundaries of the WCA in shapefile
- `Import/test/` – the climate reconstructions for each individual site
- `Import/TraCE21K/` – the TraCE-21K simulations for the different regions extracted from He et al. (2010)

## Install/run R
1. Install [R](https://larmarange.github.io/analyse-R/installation-de-R-et-RStudio.html)
2. It is easier to use [Rstudio](https://posit.co/downloads/)
3. Download this GitHub repository from ZIP file (by clicing on the green button `<> Code` beyond.

## Package management and reproducibility
This project uses [`renv`](https://rstudio.github.io/renv/) to ensure that the R package environment required to run the analyses can be reproduced. The file `renv.lock` records the package versions used for the project.

### Default: use `renv`
The recommended way to initialize the project is simply to run:

```r
source("Import/Script/setup_packages.R")
```

By default, `setup_packages.R` uses `renv` to load the project's package environment. If the `renv` environment has not yet been restored, the required packages will automatically be restored from `renv.lock`.

Therefore, for a fresh download or clone of the repository, no manual package installation should be necessary.

```r
source("Import/Script/setup_packages.R")
```

will automatically:

1. Check whether `renv` is available.
2. Load the project's `renv` environment.
3. Restore the environment from `renv.lock` when necessary.
4. Load the packages required by the project.

### Alternative: work outside `renv`
If you prefer not to use `renv`, the project also provides an alternative package-installation mechanism.

In:

```text
Import/Script/setup_packages.R
```

change:

```r
Install.PKG = F
```

to:

```r
Install.PKG = T
```

With `Install.PKG = T`, the script will install the mandatory packages listed in the project's package configuration using the standard R package installation mechanism rather than relying on the `renv` environment.

You can then run the project normally:

```r
source("Import/Script/setup_packages.R")
```

> **Recommendation:** Using `renv` is recommended for reproducibility, as it restores the package versions recorded in `renv.lock`. The alternative installation mode is provided for users who prefer to work with their regular R library.

