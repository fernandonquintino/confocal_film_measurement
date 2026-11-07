# Two-Phase Flow Characterization in a Square Channel with a Confocal Chromatic Sensor

![Python](https://img.shields.io/badge/Python-3.x-blue) ![Status](https://img.shields.io/badge/status-in%20progress-yellow)

Experimental study of horizontal air–water two-phase flow in a square channel with an 8 mm internal side (hydraulic diameter $D_h$ = 8 mm). A confocal chromatic sensor on the bottom wall measures the local liquid film and interface heights. The work covers the treatment of the confocal signal, an objective flow pattern criterion built from the confocal signal and the mean pressure gradient, flow pattern maps, and void fraction and pressure drop correlations from the literature.

> Part of PhD research. Raw and cleaned experimental data (bronze and silver layers) and the manuscripts are not included in this repository.

---

## Overview

- **Facility:** horizontal air–water flow; differential pressure measured over 0.75 m; high-speed visualization at a glass section. 164 runs; **157 analysed** ($j_l \geq$ 0.035 m/s).
- **Sensor:** Micro-Epsilon IFC2451 controller with an IFS2405-10 probe. Each sample gives the wall thickness `t1` ($\delta_1$) and the distance to the first interface `t2` ($\delta_2$). Most runs were recorded at 100 Hz, the others at 1.3 kHz to 59 kHz.
- **Visual flow patterns:** bubbly, plug, slug, dispersed, churn, wavy (including the visual stratified runs) and annular.
- **Objective criterion:** liquid-bridge statistics, the bottom film estimate and the KDE of $\delta_2$, combined with two mean pressure-gradient ratios ($R_h$, $\phi_l^2$). The current version is `objective_criteria_11`: agreement 87.9 % and balanced accuracy 88.1 % against the visual classes. The method is described in [FLOW_PATTERN_CRITERIA.md](FLOW_PATTERN_CRITERIA.md).

---

## Manuscripts (local only, never committed)

| Folder | Content |
|---|---|
| `artigo/` | Original complete manuscript `confocal_paper.tex` (Elsevier class, 21 pages, `objective_criteria_10`). Kept unchanged as the reference version; basis of a later paper on all flow patterns (*Flow Measurement and Instrumentation*). |
| `artigo_fluids/` | `confocal_full.tex`: full copy of the original with every dependency, length thresholds relative to $D_h$ (`objective_criteria_11`, figures from `paper_figures_1`). `confocal_fluids.tex`: working manuscript for *Fluids* (MDPI template in `Definitions/`) on the confocal method and slug flow. |

Both folders are in `.gitignore` and must never be pushed.

---

## Repository layout

```
phd_paper/
├── README.md                  # this file
├── FLOW_PATTERN_CRITERIA.md   # data treatment and objective criterion (scientific reference)
├── PROJECT_HISTORY.md         # revision rounds, results versions, settled decisions
├── docs/archive/              # original revision briefs, changelogs, validation notes (local only)
├── notebooks/                 # analysis pipeline (see below)
├── src/                       # Python modules
├── scripts/                   # reproducibility and verification scripts
└── data/
    ├── bronze_layer/          # raw LabVIEW and confocal CSV files (local only)
    ├── silver_layer/          # cleaned per-run datasets (local only)
    └── gold_layer/            # figures and results (tracked)
```

**Main modules (`src/`)**

| Module | Role |
|---|---|
| `load_data.py` | LabVIEW and confocal loaders (bronze → silver) |
| `confocal_treatment.py` | Mixer-wave removal (`processed_v2`) |
| `flow_pattern_criteria.py` | Criterion features and rules (`FEATURE_PARAMS`, `CRITERION_PARAMS`, `classify_frame`); `CHANNEL_SIDE` = 8.0 mm for the dimensionless lengths |
| `pressure_features.py`, `dpdz_correlations.py` | $R_h$, $\phi_l^2$ (Churchill/Jones friction factor) and 13+ pressure drop correlations |
| `flow_pattern_report.py`, `film_dynamics.py` | Figures and film-dynamics helpers of the criterion |
| `output_versions.py` | Numbered result folders that are never overwritten |
| `combined_properties.py`, `dimensionless_number.py` | Fluid properties (CoolProp) and dimensionless numbers |
| `alpha_correlations.py`, `void_fraction_per_pattern.py` | Void fraction correlations and per-pattern logic |
| `functions.py`, `modules.py`, `boxplots.py`, `plot_flow_pattern_maps.py` | Plotting, error metrics, calibration routines, maps |

---

## Data pipeline

```
[Bronze] raw LabVIEW logs + confocal CSVs (dated folders)
   │  process_raw_experimental_data.ipynb, process_raw_confocal_data.ipynb
   ▼
[Silver] experimental_results.csv; confocal_runs (regular), _processed (trimming,
         tube-thickness correction), _processed_v2 (mixer waves removed)
   │  kde_plots.ipynb (criterion)  →  criterion_verification.ipynb  →  paper_figures.ipynb
   ▼
[Gold]   objective_criteria_<N>/ (classes, features, verification, paper_figures_<N>), maps, plots
```

**Run order for the criterion results** (from `notebooks/`, headless):

```
jupyter-execute --inplace --timeout=-1 kde_plots.ipynb               # new objective_criteria_<N>
jupyter-execute --inplace --timeout=-1 criterion_verification.ipynb  # set CRITERIA_VERSION first
jupyter-execute --inplace --timeout=-1 paper_figures.ipynb           # set CRITERIA_VERSION first
```

- `kde_plots.ipynb` creates the next free `objective_criteria_<N>` folder with a `run_manifest.json` and appends a row to `data/gold_layer/objective_criteria_versions.md` (the version registry).
- `criterion_verification.ipynb` writes to `<version>/verification/`. **Set `CRITERIA_VERSION` before running**, otherwise it writes into the older version.
- `paper_figures.ipynb` creates the next free `<version>/paper_figures_<N>` folder.

Other notebooks: `thickness`, `boxplots`, `flow_pattern_maps`, `selected_conditions`, `profile_pic(_v2)`, `void_fraction_calculation`, `optional_flow_map` (a draft that is not used in the manuscript).

**Scripts (`scripts/`)**
- `convective_scaling_eval.py`: evaluation of convective time scaling (`objective_criteria_11/convective_scaling/`, results in `findings.md`).
- `minimum_film_thickness_table.py`, `flow_pattern_map*.py`: manuscript tables and maps.
- `verify_v3.py`, `verify_final_v6.py`, `refit_thresholds_rh1.py`: historical checks of `objective_criteria_10`. They use the former mm parameter keys, so run them at commit `22932c7`.

---

## Data and privacy rules

- `data/bronze_layer/` and `data/silver_layer/` hold unpublished data. They are ignored by Git and must never be pushed before publication.
- `artigo/`, `artigo_fluids/` and `docs/archive/` are ignored as well.
- `data/gold_layer/` (figures and result tables) and the code are tracked.
- Stage files by explicit path; do not use `git add -A`.

---

## Tech stack

| Purpose | Libraries |
|---|---|
| Data processing | `pandas`, `numpy` |
| Thermophysical properties | `CoolProp` |
| Statistics and curve fitting | `scipy`, `scikit-learn` |
| Visualization | `matplotlib`, `seaborn` |
| Notebooks | `Jupyter` |
| Manuscripts | LaTeX (`latexmk`; Elsevier and MDPI classes) |

---

## Status

The experimental campaign is complete. The original manuscript is finished (revision rounds in [PROJECT_HISTORY.md](PROJECT_HISTORY.md)). The *Fluids* paper is in preparation.
