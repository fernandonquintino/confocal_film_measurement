# Project history

Record of how the code, the results versions and the main decisions evolved. The per-version registry with the scores is `data/gold_layer/objective_criteria_versions.md`, and each version folder has its `run_manifest.json`. The original revision records (briefs, changelogs, validation notes) are in `docs/archive/` (local only; index in `docs/archive/README.md`). Unresolved items are in `OPEN_ITEMS.md` (local only).

Scores are given as agreement / balanced accuracy (BA) on `processed_v2`, 157 runs, seven classes.

---

## 1. Objective criterion development (versions 1–6, 2026-10-05 to 2026-10-19)

| Version | Change | Score |
|---|---|---|
| 1 | First criterion: 7 classes, $j_l \geq$ 0.035 m/s; dispersed by film jumps or ≥ 2 KDE peaks; objective stratified rule | 84.2 % / 75.2 % |
| 2 | Visual stratified evaluated as wavy (no stratified rule); film-dynamics features for churn vs wavy | 87.0 % / 88.2 % (6 classes) |
| 3 | Churn back as a class; agitated film only when $(dp/dz) \geq$ homogeneous prediction; dispersed vs churn by KDE peaks or $\phi_l^2 \leq 5$ | 86.6 % / 85.5 % |
| 4 | Liquid-rich branch (thin film seen, $R_h \geq 1$, $\phi_l^2 \leq 15$) | 86.6 % / 88.5 % |
| 5 | Same criterion; runs 10–12 relabelled wavy-annular → churn after reviewing the photos | 88.5 % / 89.4 % |
| 6 | Fluid properties at the water temperature (was the mean of water and air) | 88.5 % / 89.4 % |

The film-dynamics features of versions 2–3 were dropped from the rules. `PAPER_REVISION_PLAN.md` (now archived) planned the rewrite of the manuscript around the version 6 criterion.

## 2. Manuscript revision rounds (2026-10-20 to 2026-10-21)

**Round 1 (v6 → v7 → v8).**
- v7: reproducibility check of v6 (identical results).
- v8: laminar friction factor of a square duct in the pressure features, $f = 56.91/Re$ (Shah and London), instead of $64/Re$. No class changed: 88.5 % / 89.4 %.
- New `notebooks/criterion_verification.ipynb` (precision and recall, real sampling rates, record lengths, pressure-tolerance bounds, wall-reference sensitivity, 100 Hz downsampling, baseline comparison).
- `paper_figures.ipynb` became the generator of the manuscript figures. `scripts/minimum_film_thickness_table.py` recomputed the minimum-film table on the trimmed runs.

**Round 2 (v8 → v9 → v10).**
- v9, two corrections of the data reduction, thresholds unchanged:
  - gas mass flux with the density at the orifice meter, consistent with $j_g = Q_g/A$ at the meter (`two_phase_state`); the test-section density had underestimated $G_g$ by 15 % on average;
  - Churchill (1977) friction factor with the laminar-equivalent Reynolds number of Jones (1976), $Re^* = (64/56.91)\,Re$, continuous over all regimes (was 56.91/Re with Fang et al. above 2300).
  - Result: churn runs 27, 138, 140 became wavy, 86.6 % / 86.3 %.
- v10, threshold refit at the author's request with the $R_h$ limit fixed at 1 (`scripts/refit_thresholds_rh1.py`): dispersed $\phi_l^2$ 5 → 5.5 and liquid-rich $\phi_l^2$ 15 → 20. Runs 136 and 138 became correct: **87.9 % / 88.1 %**, 19 disagreements. This is the version of `artigo/confocal_paper.tex`.
- Verification extended: friction near the former switch, $P_{tp}$ ± 10 kPa, tolerance bounds, duration by pattern and rate, downsampling at every phase, undefined features.

**Round 3 (v10, text only).** Every number was reproduced from the code and the silver data (`scripts/verify_v3.py`); no code or result changed.

**Round 4, final pass (v10, text only).** Tables, thresholds and text numbers were checked against the outputs (`scripts/verify_final_v6.py`). One factual correction came from the data: the gap near 0.95 s in the wavy profile (run 44) is film below the detection limit, not a removed mixer wave.

## 3. Settled decisions

| Topic | Decision |
|---|---|
| Gas flow | $Q_g$ comes from the LabVIEW program (orifice meter; zero from run 0 of each day). All gas quantities refer to the gas line; $G_g$ uses the density at the meter |
| Fluid properties | At the test-section pressure and the water temperature |
| Test-section pressure | $P_{tp}$ = 93 kPa + 3.5 Δp. The coefficient comes from a linear regression of tests with an absolute transducer at the first tap, made after the runs |
| Pressure gradient | $(dp/dz)_{tp}$ = Δp / 0.75 m (tap distance) |
| Friction references | Churchill with Jones $Re^*$, $Po$ = 56.91, for both the homogeneous and the liquid-alone references |
| $R_h$ limit | Kept at 1 (homogeneous reference), although 0.75–0.95 would score up to 89.8 % BA |
| Reported thresholds | One value per threshold; the coarse-grid and leave-one-out results stay in the outputs only |
| Reference classes | Visual labels; transitional labels grouped (stratified and wavy-annular → wavy, bubbly-churn → dispersed, churn-annular → churn); runs 10–12 churn |
| Run selection | $j_l \geq$ 0.035 m/s (`JL_MIN`); no run discarded for missing data |
| Refractive index | $n_{BK7}$ = 1.519 (code value) |
| Wall reference | Fixed 0.985 mm, described but not discussed further |
| Calibration | Gravimetric points hard-coded in `old_notebooks/calibration_q1/q2.ipynb`; Coriolis data no longer exist; no full uncertainty propagation |
| ±0.15 kPa | Confirmed by the author as the transducer tolerance used for the bounds (its type is still open, see `OPEN_ITEMS.md`) |
| Optional transition map | Draft in `data/gold_layer/optional_flow_map/` (`notebooks/optional_flow_map.ipynb`, re-run with v10); not inserted in the manuscript |

## 4. Not reproducible with the project materials

- Orifice equation, plate dimensions, $C_o$ and ε: only in the LabVIEW program; only $Q_g$ is stored.
- $P_{tp}$ regression data: number of tests, conditions and residuals.
- Raw gravimetric calibration data (only the hard-coded points), Coriolis comparison, transducer calibration files.
- Camera frames are not synchronized with the confocal records, so bridges, mixer waves and individual disagreements cannot be checked frame by frame.
- The 100 Hz downsampling test reproduces the time resolution only, not the longer exposure of a real 100 Hz acquisition.

## 5. Fluids paper preparation (from 2026-11-06)

The original manuscript is split: a *Fluids* paper on the confocal method and slug flow (`artigo_fluids/`), and later a *Flow Measurement and Instrumentation* paper on all flow patterns.

- **Stage 1:** inventory of the dimensional criterion quantities. Lengths divided by $D_h$ change no class (checked down to the samples); convective time scaling evaluated.
- **Stage 2, v11 (commits `7bce457`, `82747be`):**
  - length thresholds of the criterion and the mixer-wave removal relative to $D_h$ = 8.0 mm (`CHANNEL_SIDE`, `*_rel` keys). Bit-identical features and classes to v10;
  - time thresholds kept dimensional; old mm keys rejected by `merge_params`;
  - `paper_figures.ipynb` writes to numbered `paper_figures_<N>` folders (`paper_figures_1` has the $D_h$ labels); `kde_plots.ipynb` writes its KDE figures inside the version folder;
  - `scripts/convective_scaling_eval.py` and `objective_criteria_11/convective_scaling/findings.md` (evaluation only, not in the manuscript);
  - `artigo_fluids/confocal_full.tex`: full copy of the original with $D_h$ notation.
- **Stage 3:** MDPI template cleaned into `artigo_fluids/confocal_fluids.tex`.
- **Stage 4:** documentation consolidated (this file, `OPEN_ITEMS.md`, updated `README.md` and `FLOW_PATTERN_CRITERIA.md`, archive). `notebooks/prompt.md` (an unused draft) was removed.
- **Carried forward from Stage 2:** bridge duration is not a slug length (slug units must merge bridges and aerated stretches); 98 of the 157 runs are at 100 Hz, which limits short-event resolution; the closest threshold margin is run 161, 0.0001 mm from the annular $\sigma_f$ limit.
