© 2025 The Regents of the University of Michigan  
Carson Dudley, Reiden Magdaleno -- Michigan Public Health Integrated Center for Outbreak Analytics and Modeling

# Reservoir: A Large-Scale Simulated Dataset for Training and Evaluating Epidemiological Models

![License: PolyForm Noncommercial 1.0.0](https://img.shields.io/badge/license-PolyForm--Noncommercial%201.0.0-blue)
![R](https://img.shields.io/badge/r-4.4%2B-276DC3)

---

## Try It

[Run the Colab Tutorial](_______)  
_No installation or coding required._

---

## Paper

- **Reservoir Paper:** [Reservoir: A Large-Scale Simulated Dataset for Training and Evaluating Epidemiological Models]() <br>
  **
---

## Installation


Download `Reservoir-main` as a ZIP from GitHub (green **Code** button to **Download ZIP**) and unzip it. That folder **is** the package.

Reservoir contains C++ that has to be compiled on your machine, so you need a compiler:

| OS      | what to run                                       |
|---------|---------------------------------------------------|
| macOS   | `xcode-select --install` in Terminal              |
| Windows | install Rtools from CRAN, matching your R version |
| Linux   | `sudo apt install r-base-dev`                     |

Then in R:

```{r}
install.packages("Rcpp")

install.packages(path.expand("~/Downloads/Reservoir-main"),
                 repos = NULL, type = "source")
```





