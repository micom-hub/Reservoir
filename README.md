© 2025 The Regents of the University of Michigan  
Carson Dudley, Reiden Magdaleno -- Michigan Public Health Integrated Center for Outbreak Analytics and Modeling

# Reservoir: A Large-Scale Simulated Dataset for Training and Evaluating Epidemiological Models

![License: PolyForm Noncommercial 1.0.0](https://img.shields.io/badge/license-PolyForm--Noncommercial%201.0.0-blue)
![R](https://img.shields.io/badge/r-4.4%2B-276DC3)

---

## Try It


---

## Paper

- **Reservoir Paper:** [Reservoir: A Large-Scale Simulated Dataset for Training and Evaluating Epidemiological Models]() <br>

---

## Installation

Reservoir contains C++ that has to be compiled on your machine, so you need a compiler:

| OS      | what to run                                       |
|---------|---------------------------------------------------|
| macOS   | `xcode-select --install` in Terminal              |
| Windows | install Rtools from CRAN, matching your R version |
| Linux   | `sudo apt install r-base-dev`                     |


### Downloading Reservoir through Github

There are 2 ways to download Reservoir.

In R:
```{r}
install.packages("remotes")
remotes::install_github("micom-hub/Reservoir")
```

You can also directly download `Reservoir-main` as a ZIP from GitHub (green **Code** button to **Download ZIP**) and unzip it. That folder **is** the package.

Then in R:

```{r}
install.packages("Rcpp")

install.packages(path.expand("~/Downloads/Reservoir-main"),
                 repos = NULL, type = "source")
```

### Restart R

**Session: Restart R** (Ctrl/Cmd+Shift+F10) or go to the session tab at the top and select **Restart R**.


### Load Library

In R:
```{r}
library(Reservoir)
```







