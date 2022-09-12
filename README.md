# September Workshop
This repository stores the code, presentations, and material used in the R-HTA in LMICs September Workshop. The following sections provide a breakdown of the primary documents and guidance on how to use them for your own personal training. The tutorial is based on the open-source DARTH group materials, and the majority of the code used for modelling in this tutorial is explained in the following manuscript:

-   Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM, Pechlivanoglou P, Jalal H. [A Tutorial on Time-Dependent Cohort State-Transition Models in R using a Cost-Effectiveness Analysis Example](https://arxiv.org/abs/2108.13552). arXiv:2108.13552v2. 2022:1-37.

The [`analysis`](https://github.com/R-HTA-in-LMICs/September-Workshop/tree/main/analysis) folder includes the scripts with all the code, description and comments to reproduce the CEA, probabilistic sensitivity analysis (PSA) and generation of epidemiological measures of the manuscript:

-   [`cSTM_time_dep_simulation.R`](https://github.com/R-HTA-in-LMICs/September-Workshop/blob/main/analysis/cSTM_time_dep_simulation.R): Code to replicate all simulation-time dependent cSTMs analyses of the manuscript.
-   [`cSTM_time_dep_state_residence.R`](https://github.com/R-HTA-in-LMICs/September-Workshop/blob/main/analysis/Extra%20material%20-%20cSTM_time_dep_state_residence/cSTM_time_dep_state_residence.R): Code to replicate all state-residence time dependent cSTMs analyses of the manuscript.
-   [`cSTM_torn_diag.R`](https://github.com/R-HTA-in-LMICs/September-Workshop/blob/main/analysis/cSTM_torn_diag.R): Code to implement a tornado diagram analysis for Deterministic Sensitivity Analysis.

The R scripts require loading functions that synthesize cSTMs outputs and conduct several sensitivity analyses included in the [`R`](https://github.com/R-HTA-in-LMICs/September-Workshop/tree/main/R) folder:

-   [`Funtions.R`](https://github.com/R-HTA-in-LMICs/September-Workshop/blob/main/R/Functions.R): Functions to generate epidemiological measures from time-dependent cSTMs.
-   [`Functions_cSTM_time_dep_simulation.R`](https://github.com/R-HTA-in-LMICs/September-Workshop/blob/main/R/Functions_cSTM_time_dep_simulation.R): These functions wrap the simulation-time dependent cSTMs, compute CEA and epidemiological measures, and generate probabilistic sensitivity analysis (PSA) input datasets.
-   [`Functions_cSTM_time_dep_state_residence.R`](https://github.com/R-HTA-in-LMICs/September-Workshop/blob/main/R/Functions_cSTM_time_dep_state_residence.R): These functions wrap the state-residence time dependent cSTMs, compute CEA and epidemiological measures, and generate probabilistic sensitivity analysis (PSA) input datasets.

## Preliminaries

-   Install [RStudio](https://www.rstudio.com/products/rstudio/download/)
-   Install [`dampack`](https://cran.r-project.org/web/packages/dampack/index.html) R package from CRAN

```{r, eval=FALSE}
# Install release version from CRAN
install.packages("dampack")

# Or install development version from GitHub
# devtools::install_github("DARTH-git/dampack")
```

-   Install `devtools` to install [`darthtools`](https://github.com/DARTH-git/darthtools) R package from [DARTH's GitHub](https://github.com/DARTH-git)

```{r, eval=FALSE}
# Install release version from CRAN
install.packages("devtools")

# Or install development version from GitHub
# devtools::install_github("r-lib/devtools")
```

-   Install `darthtools` using `devtools`

```{r, eval=FALSE}
# Install development version from GitHub
devtools::install_github("DARTH-git/darthtools")
```

To run the CEA, you require [`dampack`: Decision-Analytic Modeling Package](https://cran.r-project.org/web/packages/dampack/index.html), an R package for analysing and visualizing the health economic outputs of decision models.

We strongly recommend reading the DARTH group's introductory tutorial on time-independent cSTMs in R:

-   Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM, Pechlivanoglou P, Jalal H. [An Introductory Tutorial on Cohort State-Transition Models in R Using a Cost-Effectiveness Analysis Example](https://journals.sagepub.com/doi/full/10.1177/0272989X221103163). [Medical Decision Making](https://journals.sagepub.com/home/mdm), 2022 (Online First):1-18. <https://doi.org/10.1177/0272989X221103163>

Although more technical, try to also understand the use of multidimensional arrays to represent cSTM dynamics in R described in:

-   Krijkamp EM, Alarid-Escudero F, Enns EA, Pechlivanoglou P, Hunink MGM, Yang A, Jalal HJ. [A multidimensional array representation of state-transition model dynamics](https://journals.sagepub.com/doi/full/10.1177/0272989X19893973). [Medical Decision Making](https://journals.sagepub.com/home/mdm), 2020;40(2):242-248. <https://doi.org/10.1177/0272989X19893973>,

Lastly, we recommend familiarising with the useful [DARTH](http://darthworkgroup.com) coding framework described in:

-   Alarid-Escudero F, Krijkamp EM, Pechlivanoglou P, Jalal HJ, Kao SYZ, Yang A, Enns EA. [A Need for Change! A Coding Framework for Improving Transparency in Decision Modeling](https://link.springer.com/article/10.1007/s40273-019-00837-x). [PharmacoEconomics](https://www.springer.com/journal/40273), 2190;37(11):1329--1339. <https://doi.org/10.1007/s40273-019-00837-x>

## Citation

This tutorial is based on the R code from Alarid-Escudero F et al. (2022)":

> Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM, Pechlivanoglou P, Jalal H. A Tutorial on Time-Dependent Cohort State-Transition Models in R using a Cost-Effectiveness Analysis Example (<https://arxiv.org/abs/2108.13552>). arXiv:2108.13552v2. 2022:1-37.

> Alarid-Escudero F, Krijkamp EM, Enns EA, Yang A, Hunink MGM, Pechlivanoglou P, Jalal H (2022). R Code for A Tutorial on Time-Dependent Cohort State-Transition Models in R using a Cost-Effectiveness Analysis Example (Version v0.2.0). Zenodo. [10.5281/zenodo.6620902](https://www.doi.org/10.5281/zenodo.6620902). Last accessed 7 June 2022.

# Additional Information
Visit our [webpage](https://r-hta-in-lmics.github.io/) and follow the links to our social media to keep up-to-date on our latest tutorials. Alternatively, follow us on [EventBrite](https://www.eventbrite.co.uk/o/r-hta-in-lmics-46016978693) to receive notifications for when new events go live!