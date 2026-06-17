# Two-Phase Flow Characterization in a Rectangular Microchannel

![Python](https://img.shields.io/badge/Python-3.x-blue) ![Status](https://img.shields.io/badge/status-in%20progress-yellow)

Experimental and computational study of gas-liquid two-phase flow in an 8 mm rectangular microchannel. The work covers flow pattern identification, liquid film thickness measurement via confocal microscopy, and validation of void fraction and pressure drop correlations from the literature.

> Part of PhD research. Raw experimental data is not included in this repository.

---

## Overview

Air-water flows were generated across a 0.75 m test section under controlled superficial velocity conditions (164 experimental runs). Flow patterns were visually classified into seven regimes (bubbly, dispersed, plug, slug, wavy, churn, annular/stratified). Liquid film thickness was measured with a confocal chromatic sensor at three wall positions and processed to correct for optical effects. Void fraction was then computed from the film data and benchmarked against 10+ empirical correlations.

---

## Repository Layout

```
phd_paper/
├── notebooks/          # Jupyter notebooks — full analysis pipeline
│   ├── process_raw_experimental_data.ipynb
│   ├── process_raw_confocal_data.ipynb
│   ├── thickness.ipynb
│   ├── void_fraction_calculation.ipynb
│   ├── flow_pattern_maps.ipynb
│   ├── kde_plots.ipynb
│   ├── boxplots.ipynb
│   ├── profile_pic.ipynb
│   └── selected_conditions.ipynb
├── src/                # Python modules
│   ├── load_data.py                  # LabVIEW & confocal data loaders
│   ├── combined_properties.py        # Thermophysical properties (CoolProp)
│   ├── dimensionless_number.py       # Re, Fr, Bo, We, Martinelli parameter
│   ├── alpha_correlations.py         # 10+ void fraction correlations
│   ├── dpdz_correlations.py          # 13+ pressure drop correlations
│   ├── void_fraction_per_pattern.py  # Pattern-specific void fraction logic
│   ├── functions.py                  # Plotting & error metrics utilities
│   └── modules.py                    # Calibration routines
└── data/
    ├── bronze_layer/   # Raw LabVIEW & confocal CSV files
    ├── silver_layer/   # Cleaned run-level datasets (.pkl / .csv)
    └── gold_layer/     # Final figures & analysis outputs
```

---

## Data Pipeline

```
Raw measurements (LabVIEW logs + confocal sensor CSVs)
              │
              ▼
      [Bronze layer]  ── date-organized raw files
              │
              ▼ cleaning, calibration, optical correction
      [Silver layer]  ── per-run datasets, aggregated statistics
              │
              ▼ correlation validation, visualization
       [Gold layer]   ── flow maps, KDE plots, void fraction results
```

---

## Key Analyses

- **Flow pattern maps** — regime boundaries plotted on superficial velocity axes (jl vs jg), compared against Garimella et al. framework
- **Film thickness distributions** — KDE and boxplot analysis per flow regime; refractive-index correction applied to confocal readings
- **Void fraction** — experimental values derived from film geometry; benchmarked against Lockhart-Martinelli, Zivi, Premoli, Hughmark, Rouhani, Steiner, Kanizawa, Tibirica, and homogeneous models
- **Pressure drop** — 13 literature correlations implemented (Friedel, Chisholm, Sun-Mishima, Muller-Heck, and others)
- **Error metrics** — MAE, RMSE, MAPE, and η(n%) (fraction of predictions within ±n%)

---

## Tech Stack

| Purpose | Libraries |
|---|---|
| Data processing | `pandas`, `numpy` |
| Thermophysical properties | `CoolProp` |
| Statistics & curve fitting | `scipy`, `scikit-learn` |
| Visualization | `matplotlib`, `seaborn` |
| Notebooks | `Jupyter` |

---

## Status

Work in progress. Experimental campaign complete; analysis and manuscript preparation are ongoing.
