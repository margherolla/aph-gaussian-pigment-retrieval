# Gaussian Decomposition and Pigment Retrieval from Phytoplankton Absorption Spectra

MATLAB workflow for Gaussian decomposition of phytoplankton absorption spectra (**a<sub>ph</sub>(λ)**), retrieval of phytoplankton pigment concentrations using calibrated Multiple Linear Regression (MLR) models, and hierarchical clustering of pigments based on their predicted pigment composition.

---

## Overview

This repository provides a complete workflow to:

* perform Gaussian decomposition of phytoplankton absorption spectra (**a<sub>ph</sub>(λ)**);
* estimate chlorophyll-*a* (TChla) from the 676 nm Gaussian peak;
* retrieve concentrations of 13 phytoplankton pigments using calibrated MLR models;
* apply a regime-dependent retrieval for Hex-fucoxanthin (Hex) and 19'-Butanoyloxyfucoxanthin (But);
* perform hierarchical clustering of pigments based on predicted pigment-to-TChla ratios.

The workflow is based on the methodology presented in:

> **Costanzo et al.** *(manuscript in preparation).*

---

## Prerequisite

This repository requires **phytoplankton absorption spectra** (**a<sub>ph</sub>(λ)**) as input.

If your starting point is raw particulate absorption measured with a **Sea-Bird AC-S**, **a<sub>ph</sub>(λ)** should first be obtained by estimating and removing non-algal particle absorption (**a<sub>NAP</sub>(λ)**) using the companion repository:

**Companion repository**

https://github.com/margherolla/anap-pigments-retrieval

That repository derives **a<sub>NAP</sub>(λ)** from **c<sub>p</sub>(660)** and reconstructs the full **a<sub>NAP</sub>(λ)** spectrum, producing the **a<sub>ph</sub>(λ)** spectra required by the present workflow.

---

## Workflow

```text
Raw AC-S measurements
        │
        ▼
Estimate aNAP(λ) from cp(660)
(companion repository)
        │
        ▼
Compute aph(λ)
        │
        ▼
Gaussian decomposition
        │
        ▼
Gaussian peak amplitudes
        │
        ▼
MLR pigment retrieval
        │
        ▼
Pigment concentrations
        │
        ▼
Pigment ratios
        │
        ▼
Hierarchical clustering of pigments
```

---

## Repository structure

```text
.
├── data/
│   └── example_aph_stations.xlsx
│
├── models/
│   └── Pigment_MLR_final_14vars_HexBut2reg.xlsx
│
├── functions/
│   ├── spectral_decomp_aph.m
│   └── ...
│
├── scripts/
│   └── run_gaussian_pigment_example.m
│
├── outputs/
│
├── figures/
│
├── README.md
└── LICENSE
```

---

## Input data

The example input file is provided as an Excel spreadsheet.

The script automatically reads all columns named

```
aph_(400)
aph_(401)
...
```

as phytoplankton absorption spectra, together with their associated uncertainties

```
aph_sd_(400)
aph_sd_(401)
...
```

The wavelength is automatically extracted from the column names.

---

## Running the workflow

Open MATLAB, navigate to the **scripts** folder and run

```matlab
run_gaussian_pigment_example
```

The script will:

1. read the input **a<sub>ph</sub>(λ)** spectra;
2. perform Gaussian decomposition;
4. retrieve TChla and other 13 pigment concentrations;
5. perform hierarchical clustering of pigments;
6. save figures and output tables.

---

## Outputs

The workflow generates:

* Gaussian peak amplitudes
* Predicted pigment concentrations
* Pigment ratios
* Hierarchical clustering dendrogram
* Output tables in Excel format

---

## Companion repository

This repository is part of a two-step workflow.

### Step 1

Estimate **a<sub>NAP</sub>(λ)** and derive **a<sub>ph</sub>(λ)** from AC-S measurements:

https://github.com/margherolla/anap-pigments-retrieval

### Step 2

Use the resulting **a<sub>ph</sub>(λ)** spectra with the present repository to estimate phytoplankton pigments and analyse pigment composition.

---

## Citation

If you use this repository, please cite:

**Costanzo, M., Brando, V. E., Boss, E., Chase, A. P., & Doxaran, D.**, 2026. *A novel approach to estimate non-algal particle absorption for improved retrieval of pigment concentrations in coastal waters.* 

and 

**Costanzo, M., Brando, V. E., Boss, E., Chase, A. P., Doxaran, D. & Organelli, E.**, under review. *Estimation of phytoplankton community composition from continuous spectrophotometric measurements in coastal waters.*