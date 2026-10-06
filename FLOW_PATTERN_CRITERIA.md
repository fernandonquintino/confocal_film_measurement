# Confocal data treatment and objective flow pattern criteria

This document describes exactly how the confocal runs are treated (silver layer variants `processed` and `processed_v2`) and the objective criterion that assigns a flow pattern to each run. Values are those of results version **`objective_criteria_4`** (`data/gold_layer/objective_criteria_4/run_manifest.json`).

| Item | Location |
|---|---|
| Treatment code | `notebooks/process_raw_confocal_data.ipynb`, `src/load_data.py`, `src/confocal_treatment.py` |
| Criterion code | `src/flow_pattern_criteria.py`, `src/pressure_features.py` |
| Analysis notebook | `notebooks/kde_plots.ipynb`, section *Objective flow pattern criterion* |
| Results | `data/gold_layer/objective_criteria_<N>/` (one folder per version, see `data/gold_layer/objective_criteria_versions.md`) |

Each confocal run holds about 10 000 samples of `time` and three thicknesses: `t1` (tube wall, nominally 0.98 mm), `t2` (liquid film, measured from the inner wall) and `t3`. Channel: 8 mm square, horizontal, air-water.

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

1. **Event detection.** Tall readings (`t2` ≥ 4.0 mm) closer than 20 ms are grouped into one event.
2. **Removal test.** An event is removed only if all three conditions hold:
   - **short**: its longest contiguous run of tall readings lasts less than 5 ms;
   - **inside a film**: the median valid `t2` within 0.2 s before and after the event is below 2.0 mm;
   - **isolated**: fewer than 10 % of the readings within that window are tall.
3. **Extent.** The removed event is extended to the neighbouring samples with `t2` ≥ 2.0 mm (its rise and fall).
4. **Output.** Removed samples get `t2` = NaN and are flagged in a boolean `mixer_wave` column. The analysis then excludes them entirely instead of counting them as undetected.

At 100 Hz one sample lasts 10 ms, so the duration test cannot be met and those runs are unchanged: mixer waves and short slugs cannot be told apart at that rate.

The removal changes 21 runs, all recorded at 1.3 kHz or more: 9 visual wavy, 6 slug, 3 stratified, 2 plug and 1 churn (18 of them in the analysis set). The full list is in `data/silver_layer/mixer_wave_removal.csv`; the run-level summary is in `confocal_results_processed_v2.pkl/.csv`.

---

## 2. Analysis set and reference classes

- **Runs:** the 157 runs with $j_l \geq 0.035$ m/s. The 7 runs below it (38, 39, 40, 49, 104, 105, 164; mostly stratified) are excluded. The evaluated region is drawn in `data/gold_layer/flow_pattern_map_evaluated_region.pdf`.
- **Classes:** bubbly, plug, slug, churn, dispersed, wavy, annular.
- **Visual stratified → wavy:** the visual *stratified* runs are evaluated as **wavy**, because the probe sees waves in all of them.
- **No run is discarded for missing data.** A missing `t2` with a valid `t1` means the wall is seen but no interface is detected (film below the detection limit, or interface too steep). That is used as information.

---

## 3. Quantities used by the criterion

All quantities are computed per run, after excluding the `mixer_wave` samples (processed_v2).

### 3.1 Sample states (from `t2`)

| State | Definition | Physical meaning |
|---|---|---|
| liquid | `t2` ≥ 7.5 mm | liquid bridges the channel (liquid slug, bubbly liquid) |
| bubble | 2.0 ≤ `t2` < 7.5 mm | bubble, or nose/tail of an elongated bubble |
| film | `t2` < 2.0 mm | wall film (separated flow or below an elongated bubble) |
| undetected | `t2` missing | film below detection or steep interface |

### 3.2 Features

| Feature | Definition |
|---|---|
| `f_liquid`, `f_undetected` | fraction of samples in the liquid / undetected state |
| `bubble_share` | f_bubble / (f_bubble + f_film): whether the interfaces seen are bubbles or a wall film |
| `bridge_freq_hz` | number of liquid bridges per second. A bridge is a contiguous run of liquid samples lasting at least 5 ms. |
| `bridge_mean_duration_s` | mean duration of those bridges |
| `n_peaks` | number of peaks of the KDE of the detected `t2`: Gaussian KDE (Scott bandwidth) on 0–8.5 mm, prominence ≥ 5 % of the maximum; needs ≥ 50 detected samples |
| bottom film estimate | per sample: `t2` if detected and < 7.5 mm; else, if `t1` > 0.985 mm, the film merged with the wall reading $(t_1 - 0.985)\,n_\mathrm{water}/n_\mathrm{BK7}$ (unclipped `t1` from the regular run); else no value. It includes films thinner than the detection limit. |
| `film_est_mean_mm`, `film_est_std_mm` | mean and standard deviation of the bottom film estimate (needs ≥ 5 values) |
| `film_jump_frac` | among pairs of detected samples (`t2` < 7.5 mm) 10 ms apart, the fraction whose thickness changes by more than 0.3 mm (needs ≥ 30 pairs). The fixed time lag keeps 100 Hz and kHz runs comparable. High values mean the probe sees either the film or a bubble passing at a random height; low values mean smooth waves. |
| `ratio_homog` | measured $(dp/dz)$ / homogeneous prediction (section 3.3) |
| `phi2_l` | measured $(dp/dz)$ / liquid-alone prediction: the two-phase multiplier $\phi_l^2$ (section 3.3) |

### 3.3 Mean pressure-drop features (`src/pressure_features.py`)

Only the **mean** pressure gradient of each run is available: $(dp/dz)_\mathrm{tp} = \Delta p / L$, with $L$ = 0.75 m. Properties come from CoolProp at the test-section pressure and temperature; $d$ = 8 mm.

The two reference gradients use the same friction factor:

$$f = 64/Re \ (Re < 2300), \qquad f = 0.25\left[\log_{10}\!\left(\frac{150.39}{Re^{0.98865}} - \frac{152.66}{Re}\right)\right]^{-2} \ \text{otherwise}$$

**Homogeneous (McAdams), no slip.** With $G = G_l + G_g$ and $x = G_g/G$:

$$\rho_h = \left[\frac{x}{\rho_g} + \frac{1-x}{\rho_l}\right]^{-1}, \qquad \mu_h = \left[\frac{x}{\mu_g} + \frac{1-x}{\mu_l}\right]^{-1}, \qquad Re = \frac{G d}{\mu_h}, \qquad \left(\frac{dp}{dz}\right)_h = \frac{f\,G^2}{2\,d\,\rho_h}$$

**Liquid alone:**

$$\left(\frac{dp}{dz}\right)_l = \frac{f(Re_l)\,G_l^2}{2\,d\,\rho_l}, \qquad Re_l = \frac{G_l\, d}{\mu_l}$$

**The two features:**
- $\mathrm{ratio\_homog} = (dp/dz)_\mathrm{tp} / (dp/dz)_h$
  - ≥ 1: the phases are strongly mixed;
  - < 1: the gas slips over a separated liquid layer.
- $\phi_l^2 = (dp/dz)_\mathrm{tp} / (dp/dz)_l$
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
| **Annular** | `film_est_mean_mm` < 0.1 **and** `film_est_std_mm` < 0.2 mm | thin and uniform bottom film, permanently below the detection limit (runs without any film information also fall here) |
| **Dispersed** | agitated or liquid-rich (definitions below) **and** (`n_peaks` ≥ 2 **or** $\phi_l^2$ ≤ 5) | liquid-dominated friction; film plus bubble layer |
| **Churn** | agitated or liquid-rich **and** not dispersed | mixed flow in which the gas adds much of the friction |
| **Wavy** | otherwise | film with smooth waves (includes the visual stratified runs) |

**Agitated film** (wavy-churn-dispersed transition). Both conditions hold:

1. **Film signature**, one of:
   - `n_peaks` ≥ 2 (film plus bubble layer), or
   - a liquid film present with abruptly varying thickness: `film_est_mean_mm` ≥ 0.4 mm, `film_jump_frac` ≥ 0.4, and the film detected in ≥ 15 % of the samples (1 − `f_undetected` ≥ 0.15).
2. **Pressure:** `ratio_homog` ≥ 1.0 (pressure gradient at least the homogeneous prediction).

Below the homogeneous prediction, a thick, jumpy film is a separated flow with the gas slipping over it (low $j_l$ runs), and the run is wavy.

**Liquid-rich film** (high-velocity runs at the churn-annular transition). The run is neither annular nor agitated, and all of the following hold:

- `film_est_mean_mm` < 0.5 mm: the probe sees only a thin, broken film, because the interface is lost to aeration;
- `ratio_homog` ≥ 1.0: the flow is mixed;
- $\phi_l^2$ ≤ 15: liquid-dominated friction, below every thin-film wavy and annular run.

### 4.3 Thresholds

| Parameter | Value | Rule |
|---|---|---|
| `f_bridge_hz` | 0.1 Hz | intermittent family |
| `bubble_share_min` | 0.5 | bubbly |
| `plug_min_duration` | 0.02 s | plug |
| `plug_min_f_liquid` | 0.10 | plug |
| `annular_max_mean` | 0.1 mm | annular |
| `annular_max_std` | 0.2 mm | annular |
| `dispersed_min_peaks` | 2 | agitated film / dispersed |
| `dispersed_min_mean` | 0.4 mm | agitated film |
| `dispersed_min_jump` | 0.4 | agitated film |
| `dispersed_min_detected` | 0.15 | agitated film |
| `agitated_min_homog_ratio` | 1.0 | agitated and liquid-rich films |
| `dispersed_max_phi2` | 5 | dispersed vs churn |
| `liquid_rich_max_phi2` | 15 | liquid-rich film |
| `liquid_rich_max_film` | 0.5 mm | liquid-rich film |

Feature parameters:

| Parameter | Value |
|---|---|
| Liquid / film state heights | 7.5 mm / 2.0 mm |
| Minimum liquid-bridge duration | 5 ms |
| KDE peak prominence | 5 % of the maximum |
| Jump lag / jump size | 10 ms / 0.3 mm |

---

## 5. Agreement with the visual classification (version 4)

157 runs, seven classes. Agreement is the fraction of runs whose objective pattern equals the visual (reference) one. Balanced accuracy is the mean recall per class, so the large classes don't dominate it.

| Variant | Agreement | Balanced accuracy |
|---|---|---|
| regular | 75.2 % | 76.1 % |
| processed | 86.6 % | 88.5 % |
| processed_v2 | 86.6 % | 88.5 % |

Per class (processed): bubbly 11/11, plug 30/33, slug 18/25, churn 11/11, dispersed 13/20, wavy 43/47, annular 10/10.

**Known limitations:**
- The probe sits on the **bottom wall** only. It cannot see the top wall, so a thin stratified layer and the bottom of a wavy film look alike.
- **Thin-film runs near the annular transition** are classified by the mean pressure drop rather than by the film.
- **Five dispersed runs** (8, 9, 133, 134, 136) appear as churn: $\phi_l^2$ = 5.2 to 7.7, just above the dispersed threshold.
- **Some disagreements are expected at transitions:**
  - the 0313 dispersed runs (54, 55) show liquid bridges and appear as slug;
  - slugs at $j_g$ > 3 m/s without liquid bridges at the probe appear as wavy;
  - the 0228 wavy-annular runs 10 to 12 ($j_l$ > 1.2 m/s) appear as churn.
- **The pressure gate at 1.0 is close to the highest low-$j_l$ wavy run (0.99).** Propagating the uncertainty of the differential pressure transducer (±0.15 kPa) is recommended before quoting it.
- **Short records:** 10 runs have fewer than 1 000 samples (flag `low_sample_flag`). Runs 97 and 158 are slug runs too short to catch slow slugs.
- **Sensitivity:** the thresholds sit on plateaus; re-tuning all of them jointly does not improve the leave-one-out agreement.
