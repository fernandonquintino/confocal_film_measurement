# Confocal data treatment and objective flow pattern criteria

This document describes exactly how the confocal runs are treated (silver layer variants `processed` and `processed_v2`) and the objective criterion that assigns a flow pattern to each run. Values are those of results version **`objective_criteria_11`** (`data/gold_layer/objective_criteria_11/run_manifest.json`). Version 11 is version 10 with the length thresholds written relative to the hydraulic diameter $D_h$ (section 4.4); the classification is identical. Version 10 is the criterion of version 4 with the two $\phi_l^2$ limits refitted (dispersed 5.5, liquid-rich 20; $R_h$ limit kept at 1), evaluated against the revised visual labels (runs 10 to 12 relabelled churn, section 2), with the fluid properties at the water temperature, the gas mass flux with the density at the orifice meter and the Churchill friction factor with the laminar-equivalent Reynolds number of a square duct (section 3.3). Checks of this version are in `notebooks/criterion_verification.ipynb` (outputs in `objective_criteria_11/verification/`).

| Item | Location |
|---|---|
| Treatment code | `notebooks/process_raw_confocal_data.ipynb`, `src/load_data.py`, `src/confocal_treatment.py` |
| Criterion code | `src/flow_pattern_criteria.py`, `src/pressure_features.py` |
| Analysis notebook | `notebooks/kde_plots.ipynb`, section *Objective flow pattern criterion* |
| Verification notebook | `notebooks/criterion_verification.ipynb` (precision/recall, sampling rates, pressure tolerance, wall reference, 100 Hz downsampling, baseline comparison) |
| Results | `data/gold_layer/objective_criteria_<N>/` (one folder per version, see `data/gold_layer/objective_criteria_versions.md`) |
| Convective scaling (evaluation only) | `scripts/convective_scaling_eval.py` → `objective_criteria_11/convective_scaling/` (`findings.md`) |

Each confocal run holds about 10 000 samples of `time` and three thicknesses: `t1` (tube wall, nominally 0.98 mm), `t2` (liquid film, measured from the inner wall) and `t3`. Channel: square, horizontal, air-water, internal side 8 mm (nominal), equal to the hydraulic diameter $D_h$ = 8.0 mm (`CHANNEL_SIDE`). The inner side measured with the probe at the visualization section is 8.02 mm; it is used only in the sensor treatment (section 1).

---

## 1. Data variants

| Variant | Silver folder | Treatment |
|---|---|---|
| `regular` | `data/silver_layer/confocal_runs` | loading only (section 1.1) |
| `processed` | `data/silver_layer/confocal_runs_processed` | regular + trimming + tube-thickness correction (section 1.2) |
| `processed_v2` | `data/silver_layer/confocal_runs_processed_v2` | processed + mixer-wave removal (section 1.3) |

### 1.1 Regular (bronze → silver)

`load_confocal` (`src/load_data.py`) reads each bronze CSV, renames the columns to `time`, `t1`, `t2`, `t3`, and then:

1. Marks invalid readings as missing: `t1` > 8 mm, `t2` > 10 mm, `t3` > 8 mm, and the sensor's no-signal codes (-10).
2. Shifts `time` so that the file starts at 0 s.
3. Caps `t2` at 8.02 mm (channel height).

The runs are numbered sequentially (`confocal_run_1.csv` ... `confocal_run_164.csv`), in the same order as `data/silver_layer/experimental_results.csv`.

### 1.2 Processed

Applied to every regular run in this order (`process_raw_confocal_data.ipynb`, section *Process tube thickness*).

**1. Trimming to the first valid `t1`.** Many files begin with a long block in which the sensor was not yet acquiring: `t1` is missing in the first rows of 68 of the 164 runs, up to 96 % of the file. Rows before the first valid `t1` are removed, and `time` is shifted so that `time = 0` at the first valid `t1`. After the trimming, `t1` is never missing again.

**2. Tube-thickness correction.** When a liquid film thinner than the sensor's resolution lies on the wall, the sensor merges it with the wall reading, so `t1` grows above the tube thickness. With $\delta_1$ = `t1` and $\delta_2$ = `t2`:

$$d_1 = \max(\delta_1 - 0.985,\ 0), \qquad d_2 = d_1\, \frac{n_\mathrm{water}}{n_\mathrm{BK7}}, \qquad n_\mathrm{water} = 1.333,\ n_\mathrm{BK7} = 1.519$$

$$\delta_{1,\mathrm{adj}} = \delta_1 - d_1, \qquad \delta_{2,\mathrm{adj}} = \delta_2 + d_2 \quad (\text{if } \delta_2 \text{ is valid})$$

- If $\delta_2$ is missing and $d_2 > 0.373$ mm, then $\delta_{2,\mathrm{adj}} = d_2$. Below 0.373 mm (the 5th percentile of a stable run) the missing value is kept.
- $\delta_{2,\mathrm{adj}}$ is capped at 8.02 mm.

The run-level summary is saved to `data/silver_layer/confocal_results_processed.pkl/.csv`. It includes `n_samples` and `duration` because some trimmed runs are short.

### 1.3 Processed_v2: removal of mixer waves

In separated flows, liquid accumulated at the mixer is released as short waves that almost fill the channel. At the probe they appear as isolated tall readings (`t2` of 6 mm to 8 mm) lasting 1 ms to 3 ms inside a thin film. Real liquid slugs and plugs last 15 ms to 130 ms, and the liquid of bubbly flow is surrounded by other tall readings.

Starting from `processed` (`src/confocal_treatment.py`, `MIXER_WAVE_PARAMS`):

1. **Event detection.** Tall readings (`t2` ≥ 0.5 $D_h$ = 4.0 mm) closer than 20 ms are grouped into one event.
2. **Removal test.** An event is removed only if all three conditions hold:
   - **short**: its longest contiguous run of tall readings lasts less than 5 ms;
   - **inside a film**: the median valid `t2` within 0.2 s before and after the event is below 0.25 $D_h$ (2.0 mm);
   - **isolated**: fewer than 10 % of the readings within that window are tall.
3. **Extent.** The removed event is extended to the neighbouring samples with `t2` ≥ 0.25 $D_h$ (2.0 mm) (its rise and fall).
4. **Output.** Removed samples get `t2` = NaN and are flagged in a boolean `mixer_wave` column. The analysis then excludes them entirely instead of counting them as undetected.

At 100 Hz one sample lasts 10 ms, so the duration test cannot be met and those runs are unchanged: mixer waves and short slugs cannot be told apart at that rate.

The removal changes 21 runs, all recorded at 1.3 kHz or more: 9 visual wavy, 6 slug, 3 stratified, 2 plug and 1 churn (18 of them in the analysis set). The full list is in `data/silver_layer/mixer_wave_removal.csv`; the run-level summary is in `confocal_results_processed_v2.pkl/.csv`.

---

## 2. Analysis set and reference classes

- **Runs:** the 157 runs with $j_l \geq 0.035$ m/s. The 7 runs below it (38, 39, 40, 49, 104, 105, 164; mostly stratified) are excluded. The evaluated region is drawn in `data/gold_layer/flow_pattern_map_evaluated_region.pdf`.
- **Classes:** bubbly, plug, slug, churn, dispersed, wavy, annular.
- **Visual stratified → wavy:** the visual *stratified* runs are evaluated as **wavy**, because the probe sees waves in all of them.
- **Revised visual labels (2026-10-19):** after reviewing the photos, runs 10, 11 and 12 (0228_06 to 0228_08, $j_l$ = 1.24 to 1.32 m/s, $j_g$ = 8.2 to 10.0 m/s) were relabelled from wavy-annular to **churn** in `notebooks/process_raw_experimental_data.ipynb`. Run 139 (0325_07, $j_l$ = 0.80 m/s, $j_g$ = 12.0 m/s) stays wavy-annular (wavy).
- **No run is discarded for missing data.** A missing `t2` with a valid `t1` means the wall is seen but no interface is detected (film below the detection limit, or interface too steep). That is used as information.

---

## 3. Quantities used by the criterion

All quantities are computed per run, after excluding the `mixer_wave` samples (processed_v2).

### 3.1 Sample states (from `t2`)

| State | Definition | Physical meaning |
|---|---|---|
| liquid | `t2`/$D_h$ ≥ 0.9375 (7.5 mm) | liquid bridges the channel (liquid slug, bubbly liquid) |
| bubble | 0.25 ≤ `t2`/$D_h$ < 0.9375 | bubble, or nose/tail of an elongated bubble |
| film | `t2`/$D_h$ < 0.25 (2.0 mm) | wall film (separated flow or below an elongated bubble) |
| undetected | `t2` missing | film below detection or steep interface |

### 3.2 Features

| Feature | Definition |
|---|---|
| `f_liquid`, `f_undetected` | fraction of samples in the liquid / undetected state |
| `bubble_share` | f_bubble / (f_bubble + f_film): whether the interfaces seen are bubbles or a wall film |
| `bridge_freq_hz` | number of liquid bridges per second. A bridge is a contiguous run of liquid samples lasting at least 5 ms. |
| `bridge_mean_duration_s` | mean duration of those bridges |
| `n_peaks` | number of peaks of the KDE of the detected `t2`: Gaussian KDE (Scott bandwidth) on 0–1.0625 $D_h$ (0–8.5 mm), prominence ≥ 5 % of the maximum; needs ≥ 50 detected samples |
| bottom film estimate | per sample: `t2` if detected and below the liquid state; else, if `t1` > 0.985 mm, the film merged with the wall reading $(t_1 - 0.985)\,n_\mathrm{water}/n_\mathrm{BK7}$ (unclipped `t1` from the regular run); else no value. It includes films thinner than the detection limit. |
| `film_est_mean_mm`, `film_est_std_mm` | mean and standard deviation of the bottom film estimate (needs ≥ 5 values); used in the measurement figures |
| `film_est_mean_rel`, `film_est_std_rel` | the same divided by $D_h$ ($\bar\delta_f/D_h$, $\sigma_f/D_h$); used by the rules |
| `film_jump_frac` | among pairs of detected samples (below the liquid state) 10 ms apart, the fraction whose thickness changes by more than 0.0375 $D_h$ (0.3 mm) (needs ≥ 30 pairs). The fixed time lag keeps 100 Hz and kHz runs comparable. High values mean the probe sees either the film or a bubble passing at a random height; low values mean smooth waves. |
| `ratio_homog` | measured $(dp/dz)$ / homogeneous prediction (section 3.3) |
| `phi2_l` | measured $(dp/dz)$ / liquid-alone prediction: the two-phase multiplier $\phi_l^2$ (section 3.3) |

### 3.3 Mean pressure-drop features (`src/pressure_features.py`)

Only the **mean** pressure gradient of each run is available: $(dp/dz)_\mathrm{tp} = \Delta p / L$, with $L$ = 0.75 m. Properties come from CoolProp at the test-section pressure and at the **water temperature**, used as an approximation of the unmeasured test-section temperature: water carries over 83 % of the heat capacity rate in every run (median 99.8 %), so an adiabatic mixture stays within 0.6 °C of it. The former choice, the mean of the water and air temperatures, was 2.4 °C lower on average; with it, or with the gas at the air temperature, no run changes class and the features change by 1 % to 2 % (median). $D_h$ = 8 mm.

**Gas mass flux.** `Qg` (LabVIEW) is the volumetric flow at the orifice meter, so $j_g = Q_g/A$ refers to the meter and the gas mass flux uses the density at the meter: $G_g = P_g/(R_{air} T_g)\,j_g$ (`two_phase_state` in `src/pressure_features.py`; up to version 8 the density at the test section was used, which underestimated $G_g$ by 15 % on average).

**Friction factor.** The two reference gradients use the Darcy friction factor of Churchill (1977) for smooth ducts, continuous over the laminar, transitional and turbulent regimes, with the laminar-equivalent Reynolds number of Jones (1976), $Re^* = (64/Po)\,Re = 1.1246\,Re$, $Po$ = 56.91 for a square duct (Shah and London, 1978; `SQUARE_DUCT_PO`). In the laminar regime $f = 56.91/Re$. Up to version 7: 64/Re and Fang et al. (2011) above Re = 2300; version 8: 56.91/Re and Fang (discontinuous at 2300, f from 0.025 to 0.047).

$$f = 8\left[(8/Re^*)^{12} + (a+b)^{-3/2}\right]^{1/12}, \quad a = \left[2.457\ln\left(1/(7/Re^*)^{0.9}\right)\right]^{16}, \quad b = (37530/Re^*)^{16}$$

**Homogeneous (McAdams), no slip.** With $G = G_l + G_g$ and $x = G_g/G$:

$$\rho_h = \left[\frac{x}{\rho_g} + \frac{1-x}{\rho_l}\right]^{-1}, \qquad \mu_h = \left[\frac{x}{\mu_g} + \frac{1-x}{\mu_l}\right]^{-1}, \qquad Re = \frac{G D_h}{\mu_h}, \qquad \left(\frac{dp}{dz}\right)_h = \frac{f\,G^2}{2\,D_h\,\rho_h}$$

**Liquid alone:**

$$\left(\frac{dp}{dz}\right)_l = \frac{f(Re_l)\,G_l^2}{2\,D_h\,\rho_l}, \qquad Re_l = \frac{G_l\, D_h}{\mu_l}$$

**The two features:**
- $\mathrm{ratio\_homog} = (dp/dz)_\mathrm{tp} / (dp/dz)_h$ ($R_h$, homogeneous-normalized pressure-gradient ratio)
  - ≥ 1: observed for dispersed and churn flows, consistent with strongly mixed phases;
  - < 1: observed for wavy flow, consistent with gas slipping over a separated liquid layer.
- $\phi_l^2 = (dp/dz)_\mathrm{tp} / (dp/dz)_l$ (liquid-reference two-phase multiplier)
  - small: the friction is dominated by the liquid;
  - large: the gas adds much of the friction.

---

## 4. Objective criteria

The rules are applied in this order. The first rule that matches assigns the pattern.

### 4.1 Intermittent flows: the liquid bridges the channel

A run is intermittent when liquid bridges occur at least **0.1 times per second** (`bridge_freq_hz` ≥ 0.1).

| Pattern | Condition (within the intermittent family) | Physical reading |
|---|---|---|
| **Bubbly** | `bubble_share` ≥ 0.5 | the interfaces seen are mostly bubbles, not a wall film |
| **Plug** | `bridge_mean_duration_s` ≥ 0.02 s **and** `f_liquid` ≥ 0.10 | long liquid bridges occupying at least 10 % of the time |
| **Slug** | otherwise | short, aerated liquid slugs |

### 4.2 Separated flows: no liquid bridges

| Pattern | Condition | Physical reading |
|---|---|---|
| **Annular** | $\bar\delta_f/D_h$ < 0.0125 (0.1 mm) **and** $\sigma_f/D_h$ < 0.025 (0.2 mm) | thin and uniform bottom film, permanently below the detection limit (runs without any film information also fall here) |
| **Dispersed** | agitated or liquid-rich (definitions below) **and** (`n_peaks` ≥ 2 **or** $\phi_l^2$ ≤ 5.5) | liquid-dominated friction; film plus bubble layer |
| **Churn** | agitated or liquid-rich **and** not dispersed | mixed flow in which the gas adds much of the friction |
| **Wavy** | otherwise | film with smooth waves (includes the visual stratified runs) |

**Agitated film** (wavy-churn-dispersed transition). Both conditions hold:

1. **Film signature**, one of:
   - `n_peaks` ≥ 2 (film plus bubble layer), or
   - a liquid film present with abruptly varying thickness: $\bar\delta_f/D_h$ ≥ 0.05 (0.4 mm), `film_jump_frac` ≥ 0.4, and the second thickness detected in ≥ 15 % of the samples (1 − `f_undetected` ≥ 0.15).
2. **Pressure:** `ratio_homog` ≥ 1.0 (pressure gradient at least the homogeneous prediction).

Below the homogeneous prediction, a thick, jumpy film is a separated flow with the gas slipping over it (low $j_l$ runs), and the run is wavy.

**Liquid-rich film** (high-velocity runs at the churn-annular transition). The run is neither annular nor agitated, and all of the following hold:

- $\bar\delta_f/D_h$ < 0.0625 (0.5 mm): the probe sees only a thin, broken film, because the interface is lost to aeration;
- `ratio_homog` ≥ 1.0: the flow is mixed;
- $\phi_l^2$ ≤ 20: below every thin-film wavy run with `ratio_homog` ≥ 1 (≥ 24.8) and every annular run; the churn and dispersed thin-film runs reach 15.1.

### 4.3 Thresholds

| Parameter | Value | Rule |
|---|---|---|
| `f_bridge_hz` | 0.1 Hz | intermittent family |
| `bubble_share_min` | 0.5 | bubbly |
| `plug_min_duration` | 0.02 s | plug |
| `plug_min_f_liquid` | 0.10 | plug |
| `annular_max_mean_rel` | 0.0125 (0.1 mm) | annular |
| `annular_max_std_rel` | 0.025 (0.2 mm) | annular |
| `dispersed_min_peaks` | 2 | agitated film / dispersed |
| `dispersed_min_mean_rel` | 0.05 (0.4 mm) | agitated film |
| `dispersed_min_jump` | 0.4 | agitated film |
| `dispersed_min_detected` | 0.15 | agitated film |
| `agitated_min_homog_ratio` | 1.0 | agitated and liquid-rich films |
| `dispersed_max_phi2` | 5.5 | dispersed vs churn |
| `liquid_rich_max_phi2` | 20 | liquid-rich film |
| `liquid_rich_max_film_rel` | 0.0625 (0.5 mm) | liquid-rich film |

Feature parameters (`FEATURE_PARAMS`):

| Parameter | Key | Value |
|---|---|---|
| Liquid / film state heights | `h_liq_rel`, `h_film_rel` | 0.9375 / 0.25 $D_h$ (7.5 / 2.0 mm) |
| Minimum liquid-bridge duration | `min_dwell` | 5 ms |
| KDE peak prominence | `kde_prominence` | 5 % of the maximum |
| Jump lag / jump size | `jump_lag`, `jump_min_change_rel` | 10 ms / 0.0375 $D_h$ (0.3 mm) |

### 4.4 Dimensional and dimensionless quantities

- **Lengths of the criterion** are relative to the nominal hydraulic diameter $D_h$ = 8.0 mm (`CHANNEL_SIDE`; keys and columns ending in `_rel`). Because 8 is a power of two, value_rel = value_mm / 8 is exact in floating point: version 11 reproduces every feature and class of version 10 bit for bit. The mixer-wave removal uses the same convention (`MIXER_WAVE_PARAMS`: `h_tall_rel`, `h_film_rel`, `context_max_film_rel`).
- **Lengths of the sensor treatment stay in mm**, because they are properties of the wall and the instrument, not of the flow: wall reference 0.985 mm, fill limit 0.373 mm, cap 8.02 mm (measured inner side, `CHANNEL_HEIGHT`), loading limits (section 1.1).
- **Time thresholds stay dimensional** (0.1 Hz, 20 ms, 5 ms, 10 ms; mixer waves 20 ms, 5 ms, 0.2 s). The minimum bridge duration and the jump lag are tied to the sample interval; the family and plug limits to the record length.
- **Old mm keys are rejected:** `merge_params` raises a `KeyError` naming the new key, so an override with a former key (e.g. `h_liq`, `annular_max_mean`) is never silently ignored. The historical scripts `verify_v3.py`, `verify_final_v6.py` and `refit_thresholds_rh1.py` use the former keys and run at commit `22932c7`.
- **Convective scaling (evaluated, not adopted).** With $U = j_m$ (also $j_l$ and a drift-flux $U_t$), the Strouhal form $St = f_{bridge} D_h/U$ separates the families as well as $f_{bridge}$ (BA 88.1 %), while $L^* = \bar t_{bridge} U/D_h$ degrades the plug/slug split (BA 86.3 %, leave-one-out 84.6 %). $L^*$ is not a liquid-slug length: bridges are the unaerated stretches of a slug. Processing times scaled convectively would fall below one sample in 76 (minimum bridge) and 53 (jump lag) of the 98 runs at 100 Hz. Details: `objective_criteria_11/convective_scaling/findings.md`.

---

## 5. Agreement with the visual classification (versions 10 and 11)

157 runs, seven classes. Agreement is the fraction of runs whose objective pattern equals the visual (reference) one. Balanced accuracy is the mean recall per class, so the large classes don't dominate it.

| Variant | Agreement | Balanced accuracy |
|---|---|---|
| regular | 76.4 % | 76.8 % |
| processed | 87.9 % | 88.1 % |
| processed_v2 | 87.9 % | 88.1 % |

Per class (processed): bubbly 11/11, plug 30/33, slug 18/25, churn 12/14, dispersed 14/20, wavy 43/44, annular 10/10. Version 10 refits the two $\phi_l^2$ limits on the version 9 features (coarse grid with the `ratio_homog` limit fixed at 1; dispersed 5 → 5.5, liquid-rich 15 → 20): runs 136 (dispersed) and 138 (churn) become correct (version 9: 86.6 % / 86.3 %). Version 9 (gas mass flux at the meter density, Churchill/Jones friction factor) changes three classes with respect to version 8: churn runs 27, 138 and 140 become wavy (version 8: 88.5 % / 89.4 %). The objective patterns are the same as in version 4; only the reference labels of runs 10 to 12 changed (version 4, with them as wavy-annular: 86.6 % / 88.5 %). Version 6 (properties at the water temperature) gives the same objective patterns as version 5, and version 8 (square-duct laminar friction factor: $\phi_l^2$ +12.5 % in the 58 runs with laminar $Re_l$, $R_h$ up to +12.5 % in the 15 runs with laminar $Re_{tp}$) the same as version 6.

Precision per class (runs correctly identified / runs assigned to the class): bubbly 11/11, plug 30/31, slug 18/23, churn 12/16, dispersed 14/14, wavy 43/52, annular 10/10.

**Known limitations:**
- The probe sits on the **bottom wall** only. It cannot see the top wall, so a thin stratified layer and the bottom of a wavy film look alike.
- **Thin-film runs near the annular transition** are classified by the mean pressure drop rather than by the film.
- **Four dispersed runs** (8, 9, 133, 134) appear as churn: $\phi_l^2$ = 6.2 to 8.1, above the dispersed threshold of 5.5. Two churn runs (27, 140) appear as wavy: `ratio_homog` 0.98 and 0.94, just below 1.
- **Some disagreements are expected at transitions:**
  - the 0313 dispersed runs (54, 55) show liquid bridges and appear as slug;
  - slugs at $j_g$ > 3 m/s without liquid bridges at the probe appear as wavy.
- **The pressure gate at 1.0 is close to the highest agitated-film run below it:** slug run 91 (visual slug, classified wavy) has `ratio_homog` = 0.981; the low-$j_l$ wavy and slug runs with an agitated film, and churn runs 27 and 140, span 0.68 to 0.98. Within the ±0.15 kPa of the differential pressure transducer (±200 Pa/m) runs 91 and 116 can change between churn and wavy; with the flow-rate tolerances (±3 % liquid, ±5 % gas, worst combination) also runs 21, 27 and 140 (churn/wavy) and 95 and 136 (churn/dispersed) (`verification/pressure_bounds.csv`). A shift of ±10 kPa in the estimated test-section pressure changes `ratio_homog` by about ±7 % and moves four runs: 116 at −10 kPa; 27, 91 and 140 at +10 kPa (`verification/ptp_sensitivity.csv`). Using (Tl + Tg)/2 or Tg for the properties changes no class.
- **Short records:** 10 runs have fewer than 1 000 samples (flag `low_sample_flag`). Runs 97 and 158 are slug runs too short to catch slow slugs.
- **Sensitivity:** the bubble share, plug duration, liquid-rich film and liquid-rich $\phi_l^2$ (16 to 24) limits sit on plateaus; `agitated_min_homog_ratio` and `dispersed_max_phi2` are the most sensitive (a change of 0.05 or 0.5 moves one or two runs). `agitated_min_homog_ratio` of 0.75 to 0.95 would give up to 89.8 % balanced accuracy, but it is kept at 1 (homogeneous reference). The coarse grid (11 thresholds, 103 680 combinations) reaches 90.1 % balanced accuracy in-sample (24 tied combinations); leave-one-out re-tuning gives 86.6 % (ties broken by grid order) or 89.1 % (ties towards the chosen thresholds). These grid results are not reported in the paper (author decision).
- **Acquisition rate:** downsampling the 59 kHz runs to 100 Hz changes 1 to 13 classes depending on the sampling phase (7 at the first sample; 16 runs change at some phase, 2 at every phase), so the criterion depends on the rate (`verification/downsample_100hz_phases.csv`).
- **Wall reference:** with 0.995 mm instead of 0.985 mm no class changes; 0.975 mm (below the dry-wall readings) changes runs 1, 4 and 37 (`verification/wall_reference_sensitivity.csv`).
- **Friction model:** with the gas mass flux at the meter density, changing only the friction model (piecewise vs Churchill/Jones) changes run 116; 29 runs have $Re_{tp}$ and 31 $Re_l$ between 1800 and 4000 (`verification/friction_switch.csv`).
