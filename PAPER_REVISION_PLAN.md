# Paper revision plan: adding the objective flow pattern criterion

Guide for rewriting `confocal_paper.pdf` around the objective flow pattern criterion. For each section it lists **what the section has to contain**: the facts, arguments, numbers and physical explanations, plus the figures, tables and equations it uses. The wording is left to the author.

| Source | Use |
|---|---|
| `confocal_paper.pdf` | current version (10 pages, 21 figures, 4 tables, 8 equations). "Old" numbers below refer to it. |
| `FLOW_PATTERN_CRITERIA.md` | treatment and criterion definitions (cited as **FPC §x**) |
| `data/gold_layer/objective_criteria_5/` | results of the criterion, version 5 (the version to report: same criterion as version 4, with runs 10–12 relabelled churn) |
| `data/gold_layer/objective_criteria_5/run_manifest.json` | every threshold and parameter |
| `data/gold_layer/objective_criteria_5/processed_v2/objective_flow_pattern.csv` | features and objective pattern of each of the 157 runs |
| `data/gold_layer/objective_criteria_5/variant_comparison_criterion.csv` | agreement per data variant |
| `notebooks/paper_figures.ipynb` → `data/gold_layer/paper_figures/` | every data figure of the revised paper (§4) |

Decisions taken:
- Section 2.5 describes the **processed_v2** treatment: trimming + tube-thickness correction + mixer-wave removal.
- The criterion gets its **own Section 3**. Results move to Section 4.
- Results come from `objective_criteria_5/processed_v2/`, figures from `data/gold_layer/paper_figures/`. processed and processed_v2 give identical results: the same confusion matrix and the same 18 disagreements.

Status labels used below:

| Label | Meaning |
|---|---|
| **keep** | unchanged |
| **edit** | same item, text or caption changes |
| **regenerate** | same kind of figure, re-plot with new data or classes |
| **ready** | new item, the file already exists |
| **to create** | new item, does not exist yet |
| **remove** | drop from the paper |

---

## 1. Number of valid runs (checked against the data)

**The analysed set is 157 runs. Removing stratified did not remove any run:** the 3 visual stratified runs that remain inside the evaluated region are evaluated as wavy. Checked against `data/silver_layer/experimental_results.csv` and `objective_flow_pattern.csv`.

| Step | Runs | Detail |
|---|---|---|
| Experimental campaign | **164** | all runs with confocal data |
| Removed: $j_l$ < 0.035 m/s | −7 | 38, 39, 49, 105 (stratified); 40 (stratified-wavy); 104 (wavy); 164 (annular) |
| **Evaluated set** | **157** | used by every result of the criterion |

In those 157 runs, the visual log labels map onto the 7 reference classes like this (the paper should state this mapping, Table 6):

| Reference class | Runs | Visual labels included |
|---|---|---|
| bubbly | 11 | bubbly (11) |
| plug | 33 | plug (32), bubble (1: run 124, $j_l$ = 0.24, $j_g$ = 0.26 m/s; **check this label**) |
| slug | 25 | slug (25) |
| dispersed | 20 | bubbly-churn (16), dispersed bubbles (4) |
| churn | 14 | churn (11: 8 logged as churn + runs 10, 11, 12, first logged as wavy-annular and relabelled churn after reviewing the photos), churn-annular (3: runs 1, 4, 138) |
| wavy | 44 | wavy (20), wavy-annular (12: incl. run 139, kept as wavy after reviewing the photos), annular-wavy (1: run 157), stratified-wavy (8), **stratified (3: runs 48, 108, 147)** |
| annular | 10 | annular (10) |
| **Total** | **157** | |

Relabelled after reviewing the photos (2026-10-19): runs 10, 11, 12 (0228_06 to 0228_08; $j_l$ = 1.24–1.32, $j_g$ = 8.2–10.0 m/s) changed from wavy-annular to **churn** in `notebooks/process_raw_experimental_data.ipynb`. Run 139 (0325_07, $j_l$ = 0.80, $j_g$ = 12.0 m/s) stays wavy. The criterion results are version **`objective_criteria_5`**: same rules and thresholds as version 4, new reference labels.

What happened to stratified:
- **Pure visual stratified:** 7 runs in the campaign.
  - 4 have $j_l$ < 0.035 m/s and are excluded.
  - 3 remain (48, 108, 147) and are evaluated as wavy. The criterion classifies all three as wavy.
- **Visual stratified-wavy:** 9 runs. One is excluded (40); the other 8 were already in the wavy class.

Alternative, if you prefer to drop the 3 remaining stratified runs instead of relabelling them: **154 runs**, agreement 88.3 %, balanced accuracy 89.4 %, wavy 40/41. The conclusions don't change. The rest of this plan uses 157 (the current gold results).

For the "53 valid runs" sentence of the old paper: after trimming, 69 runs have < 50 % missing `t2`. This number only matters to explain why the filter is dropped.

---

## 2. What changes, in one page

| Topic | Current paper | Revised paper |
|---|---|---|
| Runs analysed | 53 "valid" runs (< 50 % missing `t2`) | **157 runs** with $j_l \geq 0.035$ m/s |
| Missing `t2` | a reason to discard a run | **information**: the wall is seen but the interface is not detected |
| Classes | 8 visual (incl. stratified) | **7 reference classes**: bubbly, plug, slug, churn, dispersed, wavy, annular |
| Treatment | tube-thickness correction | + **trimming** to the first valid `t1` (68 runs) + **mixer-wave removal** (21 runs) |
| Objective identification | suggested from the PDF shapes | **rule-based criterion**: confocal signal + **mean pressure drop** |
| Result | – | **88.5 % agreement, 89.4 % balanced accuracy** (77.1 % / 78.2 % on untreated data) |

Key message: the variant `regular (trimmed only)` already reaches 88.5 % / 89.4 %. **Trimming is the treatment step that drives the gain.** Tube-thickness correction and mixer-wave removal are physically motivated but do not change the classification.

Phrases to find and fix everywhere:
- "53 valid runs", "missing rate lesser than 50 %", and the solid/open-marker explanation;
- "stratified" in the lists of observed patterns;
- "can be applied for objective flow pattern identification" (abstract) and "possible to use those as an objective criteria" (Sec. 3.3);
- footer "Preprint submitted to Nuclear Physics B" (template leftover).

---

## 3. Section-by-section content

### Abstract (edit)

Must contain, in this order:
1. **Context:** a liquid film is present in most two-phase flow patterns and controls heat transfer, pressure drop and void fraction. Measuring it is difficult because it is thin and changes quickly. (Keep the current text.)
2. **Technique:** confocal chromatic microscopy, non-invasive and high-speed, applied here to horizontal air-water flow in an 8 mm square channel. (Keep.)
3. **Dataset:** 164 runs; 157 evaluated over seven patterns (bubbly, plug, slug, churn, dispersed, wavy, annular).
4. **Treatment:** trimming, tube-thickness correction (recovers films below the detection limit from the wall reading), and removal of mixer waves. Missing film readings are used as information, not discarded.
5. **Criterion:** the run is classified by rules on physically interpretable quantities:
   - liquid-bridge statistics (intermittent flows);
   - the bottom-film mean, variability and abruptness (separated flows);
   - the mean pressure gradient, normalised by the homogeneous and liquid-alone predictions.
6. **Result:** 88.5 % agreement with the visual classification and 89.4 % balanced accuracy. Every bubbly, churn and annular run is identified correctly. The disagreements are at transitions.
7. **Limitation:** the probe sees one wall only, and films below about 0.3–0.5 mm are not detected directly.

Remove: "indicating that the same technique can be applied for objective flow pattern identification" (it is now done).

### 1. Introduction (edit)

Keep the paragraphs on film importance, horizontal asymmetry, non-circular channels, the review of measurement techniques, and the confocal studies [1]–[20].

Must add or change:
1. **New paragraph: why an objective criterion.**
   - Flow patterns are usually identified visually (high-speed camera). This is subjective, especially at transitions (wavy–annular, churn–dispersed, plug–slug), and depends on optical access.
   - Objective methods classify a run from a measured signal instead. Typical inputs:
     - the PDF of a local signal (void fraction, conductance, film thickness), whose shape differs by pattern;
     - statistics of the pressure signal.
   - Give 2–4 references (see §8).
2. **Gap:** no study used a confocal chromatic sensor to identify flow patterns objectively. Its high rate and µm resolution resolve liquid bridges and film dynamics. On the other hand, it has a detection limit and sees a single wall.
3. **Objective of the work:**
   - (i) assess the technique for all flow patterns of a horizontal square channel;
   - (ii) propose a treatment that keeps every run;
   - (iii) build an objective criterion from the confocal signal and the mean pressure drop, and evaluate it against the visual classification.
4. Change the list "stratified, wave, slug, churn and annular" to the patterns actually observed: bubbly, plug, slug, churn, dispersed, wavy, annular.

### 2. Experimental methodology

#### 2.1 Two-phase flow test facility (keep)
- Fig. 1 and Fig. 2; the current text.
- Add one sentence: the differential pressure transducer of the test section ($\Delta P_{tp}$, Table 1), between taps 0.75 m apart, gives the mean pressure gradient used by the criterion. Only its time average per run is available; there is no dynamic signal.

#### 2.2 Experimental data reduction (keep + one addition)
- Eqs. 1–4 and the current text.
- Add:
  - $(dp/dz)_{tp} = \Delta p_{tp}/L$ with $L$ = 0.75 m;
  - water and air properties come from CoolProp at the test-section pressure and temperature;
  - mass fluxes $G_l = \rho_l j_l$, $G_g = \rho_g j_g$ and quality $x = G_g/(G_l+G_g)$.
- These definitions are needed in Section 3.3.

#### 2.3 Uncertainty analysis (edit)
- Table 1 and the current calibration text.
- Add, if done: the uncertainty of `ratio_homog`, propagated from ±0.15 kPa in $\Delta P_{tp}$ and the flow-rate uncertainties. It matters because the criterion uses `ratio_homog` ≥ 1.0, and the highest low-$j_l$ wavy run is at 0.99. If not done, state it as a limitation in Section 4.6.

#### 2.4 Confocal measurement system (keep + move)
- Figs. 3–5, Table 2 and the current text: principle, refractive-index correction, 0.5 mm minimum from the manufacturer, dark reference, empty-tube check.
- Move here from Sec. 2.5:
  - the measurement rate was set per run: about 1.5 kHz for thick films and 0.1 kHz for thin films; some 0228 runs at about 3.6 kHz;
  - a higher rate means a shorter exposure, so less light returns from weak interfaces.
- In the 157 runs: 98 at 100 Hz, 50 at 1.3–2 kHz, 9 above 2 kHz.
- Output columns: `t1` (tube wall, nominally 0.98 mm), `t2` (liquid film, measured from the inner wall), `t3` (not used). Each file has about 10 000 samples.

#### 2.5 Confocal data treatment (rewrite)

The old section is built around discarding runs. The new one describes a pipeline that keeps every run. Five subsections.

**2.5.1 Loading.**
- Sensor no-signal codes (-10, 2147483640) become missing values.
- Physically impossible readings become missing: `t1` > 8 mm, `t2` > 10 mm, `t3` > 8 mm.
- Time starts at 0. `t2` is capped at 8.02 mm, the measured inner side of the channel (Table 2).

**2.5.2 Trimming to the first valid `t1` (new).**
- **Observation:** 68 of the 164 files start with a block in which the sensor was not yet acquiring. In that block both `t1` and `t2` are missing, for up to 96 % of the file.
- **Why it matters:** counted as missing film readings, these blocks inflated the missing fraction. They were the main reason the old < 50 % filter rejected most runs: 69 runs pass that filter after trimming, against 53 before.
- **Procedure:** rows before the first valid `t1` are removed, and time is reset to 0 there.
- **Consequence:**
  - After trimming, `t1` is never missing again. So every remaining missing `t2` means "wall seen, interface not detected", which is physically meaningful (2.5.5).
  - Some records become short: 10 runs have < 1 000 samples (flagged).
- **Effect on the criterion:** untreated 77.1 % → trimmed 88.5 % agreement (Table 11). Mention it here or in Section 4.5.

**2.5.3 Tube-thickness correction (keep the current argument, reorganised).**
- **Physical mechanism:** when a liquid film thinner than the sensor's resolution lies on the inner wall, the reflections from the glass–liquid and liquid–gas interfaces cannot be separated. The sensor reports one, thicker "wall". So `t1` rises above the tube thickness in patterns with thin films. New Fig. 6 (157 runs, medians): wavy 1.038 and annular 1.035 mm, dispersed and churn 0.998 mm, against 0.977–0.980 mm for bubbly, plug and slug; the highest run reaches 1.104 mm.
- **Fig. 6:** `t1` before the treatment, by pattern.
- **Table 3:** choice of the 0.985 mm threshold. The minimum is unaffected; at higher thresholds the film starts to be included in all statistics.
- **Correction:** the excess $d_1$ above 0.985 mm is removed from `t1`. It is added to `t2` after converting from the glass refractive index to the water one ($d_2 = d_1\,n_\mathrm{water}/n_\mathrm{BK7}$, Eqs. 5–8).
- **Minimum film:** when `t2` is missing, the excess is used as the film only if $d_2$ > 0.373 mm. That is the 5th percentile of a stable run (Table 4). Smaller excesses are kept as missing, to avoid adding noise; this fills 0.3 % of the missing values, against 25 % with no threshold.
- **Cap:** the corrected `t2` is capped at 8.02 mm.
- **Fig. 7:** `t1` after the treatment. Every run lies between 0.970 and 0.985 mm (medians 0.976–0.985 mm). Annular sits at the 0.985 mm threshold, because its `t1` is above it almost all the time.
- Point forward: Section 3 reuses the same excess, without the 0.373 mm threshold, to estimate the bottom film when it is below detection (Eq. 9).

**2.5.4 Mixer-wave removal (new).**
- **Physical mechanism:** in separated flows, liquid accumulates at the Y-shaped mixer and is released periodically as short waves that almost fill the channel. At the probe, 1.33 m downstream, they appear as isolated tall readings (`t2` of 6–8 mm) lasting 1–3 ms inside a thin film.
- **Why remove them:** they are not part of the developed flow pattern. Real liquid slugs and plugs last 15–130 ms, and the liquid of bubbly flow is surrounded by other tall readings. Left in, they look like liquid bridges or bubbles.
- **Detection:**
  - tall readings (`t2` ≥ 4 mm) closer than 20 ms form one event;
  - the event is removed only if it is short (longest contiguous tall stretch < 5 ms), inside a film (median `t2` within ±0.2 s < 2 mm) and isolated (< 10 % tall readings in that window);
  - the removed stretch includes the rise and fall (neighbouring samples ≥ 2 mm);
  - removed samples are flagged and excluded from every statistic, not counted as undetected.
- **Limitation:** at 100 Hz one sample lasts 10 ms, so the duration test can't be met. Those runs are unchanged: mixer waves and short slugs can't be told apart at that rate.
- **Result:** 21 runs changed, all at ≥ 1.3 kHz (9 wavy, 6 slug, 3 stratified, 2 plug, 1 churn; 18 in the evaluated set).
  - The objective classification doesn't change. A likely reason: the 5 ms minimum liquid-bridge duration (Section 3.2) already ignores such short events.
  - The value of the step is a cleaner signal, for the profiles and the PDFs.
- **Fig. 8** (ready), **Table 5** (parameters).

**2.5.5 Missing `t2` as information (new).**
- No run is discarded for missing data.
- A missing `t2` with a valid `t1` has a physical meaning: the wall is seen, but no interface is detected. Either the film is below the detection limit (about 0.25–0.5 mm, Table 4), or the interface is too steep or too aerated to reflect light back to the probe (the maximum interface angle is 14°).
- **The missing fraction can't be a criterion,** because it depends on the acquisition settings. Dispersed runs at similar $j_l$ and $j_g$ lose 74–82 % of `t2` at about 3.6 kHz (0228 runs), but about 0 % at 100 Hz (0318 runs).
- So the criterion uses three things instead:
  - whether the channel is bridged by liquid;
  - the detected thickness;
  - the film hidden in `t1` (Eq. 9).

#### 2.6 Experimental range (edit)

Must contain:
- **Fig. 9 = `paper_figures/fig09_evaluated_region.pdf`** (same as `data/gold_layer/flow_pattern_map_evaluated_region.pdf`; replaces old Fig. 8).
  - 164 runs; the range of $j_l$ and $j_g$.
  - The **evaluated region** holds the 157 runs with $j_l \geq 0.035$ m/s. The 7 excluded runs are mostly stratified (4 of 7).
  - The transition lines are those of Coleman and Garimella [22] (5.36 mm hydraulic diameter, aspect ratio 0.725, horizontal air-water).
  - The visual classification of each run appears later, next to the objective one (Fig. 24).
- **Visual classification:** high-speed camera, mostly 8000 fps. The log has compound labels (bubbly-churn, wavy-annular…). These were mapped to 7 reference classes: **Table 6** (the mapping in §1 of this plan).
- **Fig. 10 = old Fig. 9 without panel (d) stratified.** Seven panels, ordered by the two families of the criterion: (a) bubbly, (b) plug, (c) slug | (d) dispersed, (e) churn, (f) wavy, (g) annular.
- **Physical description of each pattern,** as seen in the images. It should match the criterion's reading:
  - **Bubbly:** small bubbles carried by a continuous liquid, mostly in the upper half because of buoyancy. No wall film below them.
  - **Plug:** long, smooth elongated bubbles in the upper part of the channel, with a thin film below them, separated by long liquid plugs without bubbles.
  - **Slug:** as plug, but the higher gas velocity aerates the liquid slugs and shortens them. The bubble tail is disturbed.
  - **Dispersed:** high liquid velocity. Gas broken into bubbles over a liquid film whose thickness varies abruptly. In this channel it is better described as the wavy–churn–dispersed transition (keep the current remark that velocities were not high enough for pure dispersed bubbles).
  - **Churn:** chaotic, oscillating flow. Neither a stable film nor stable liquid slugs. The gas carries much of the momentum.
  - **Wavy:** a continuous liquid layer at the bottom with waves on its surface. The gas flows above it and slips over it.
  - **Annular:** the gas flows in the core, and the liquid forms a thin film all around the wall, thicker at the bottom because of gravity.
- **Stratified:** the visual stratified runs are evaluated as wavy, because the bottom probe sees waves in all of them. A smooth stratified layer and the bottom of a wavy film are indistinguishable from one wall.

### 3. Objective flow pattern criterion (new section)

This section is the core of the revision. Everything is computed per run on processed_v2 data, with the mixer-wave samples excluded. For each quantity this plan gives the definition and then the **physical reading**, so the text can explain what the probe sees.

**Opening paragraph.**
- The criterion classifies each run from its confocal record and its mean pressure gradient. It uses a few physically interpretable quantities and ordered rules, not a black-box classifier.
- Superficial velocities are **not** inputs. The map (Fig. 24) is only a physical check.
- **What the probe measures:** `t2` is the distance from the inner bottom wall to the first interface met by the light. A tall reading means the light went far into the channel before meeting an interface. A short one means a film on the wall. A missing one means no interface returned light.

#### 3.1 Sample states

Each sample gets one of four states from `t2` (**Table 7**):

| State | Definition | Physical reading |
|---|---|---|
| liquid | `t2` ≥ 7.5 mm | No interface between the bottom and the top wall: the probe's vertical line is full of liquid. This is a **liquid bridge**: a liquid slug or plug, or the continuous liquid of bubbly flow. |
| bubble | 2.0 ≤ `t2` < 7.5 mm | The interface seen is the underside of a gas structure well above the wall: a dispersed bubble, or the nose or tail of an elongated bubble. |
| film | `t2` < 2.0 mm | The interface is close to the wall: the liquid film below an elongated bubble, or the liquid layer of a separated flow. |
| undetected | `t2` missing | Wall seen, interface not seen: film below the detection limit, or a steep or aerated interface (2.5.5). |

Why these thresholds:
- **7.5 mm:** close to the 8.02 mm channel height, leaving a margin for readings at the top wall.
- **2.0 mm:** separates films from bubbles. The film below elongated bubbles starts at about 2.2 mm and drains to about 0.36 mm (old Fig. 11a). The films of separated flows are below 1.5 mm.

#### 3.2 Features from the confocal signal

Present only the features the rules use (**Table 8**). The full list can go in an appendix. Each feature needs its definition and its physical meaning. Median values per class are in Table 10, from `objective_flow_pattern.csv`.

| Feature | Definition | Physical reading | Median per class |
|---|---|---|---|
| `bridge_freq_hz` | liquid bridges per second. A bridge is a contiguous run of liquid samples lasting ≥ 5 ms. | How often the channel is filled with liquid at the probe. This is the signature of intermittent flow. The 5 ms minimum rejects one-sample spikes at kHz rates; at 100 Hz any single sample counts. | bubbly 23 Hz, plug 6.0, slug 0.68; separated 0 |
| `bridge_mean_duration_s` | mean duration of the bridges | How long the liquid slug or plug takes to pass the probe. Long for unaerated plugs, short when bubbles break the slug. | plug 43 ms, slug 12 ms, bubbly 8 ms |
| `f_liquid` | fraction of samples in the liquid state | Share of time the channel is full of liquid (the liquid "hold-up" seen by the probe) | bubbly 0.38, plug 0.34, slug 0.008 |
| `bubble_share` | $f_\mathrm{bubble}/(f_\mathrm{bubble}+f_\mathrm{film})$ (Eq. 10) | Among the interfaces seen that are not bridges, whether they are bubbles at mid-height or a wall film | bubbly 0.83; plug 0.08, slug 0.10 |
| `n_peaks` | number of peaks of the KDE of the detected `t2` (Gaussian, Scott bandwidth, 0–8.5 mm, prominence ≥ 5 % of the maximum, ≥ 50 detected samples) | Number of distinct heights where interfaces are found. Two peaks = film + liquid bridges (intermittent), or film + a layer of bubbles (dispersed). One peak = a single film. | bubbly, plug, dispersed 2; slug, churn, wavy 1 |
| bottom film estimate (Eq. 9) | per sample: `t2` if detected and < 7.5 mm; otherwise $(t_1-0.985)\,n_\mathrm{water}/n_\mathrm{BK7}$ if $t_1$ > 0.985 mm; otherwise no value | Thickness of the liquid on the bottom wall, **including films below the detection limit**, which are recovered from the wall reading (same mechanism as 2.5.3). Liquid bridges are excluded. | – |
| `film_est_mean_mm` | mean of the bottom film estimate (≥ 5 values) | How much liquid lies on the bottom wall | bubbly 3.0 mm (bubble heights), plug 1.0, slug 1.1, dispersed 1.05, churn 0.45, wavy 0.27, annular 0.055 |
| `film_est_std_mm` | standard deviation of the same | How much the bottom film varies (waves, passing structures) | wavy 0.43, churn 0.42, dispersed 0.46, annular 0.036 |
| `film_jump_frac` | among pairs of detected samples (`t2` < 7.5 mm) 10 ms apart, the fraction whose thickness changes by > 0.3 mm (≥ 30 pairs) (Eq. 11) | **Smooth vs abrupt film.** A wave changes the thickness continuously, so readings 10 ms apart are close. In an agitated film, the probe sees either the wall film or the underside of a bubble passing at a random height, so the readings jump. The fixed 10 ms lag keeps 100 Hz and kHz runs comparable. | dispersed 0.48, churn 0.45, wavy 0.38, plug 0.19 |
| `f_undetected` | fraction of undetected samples | Used only as a reliability condition (the film must be seen in ≥ 15 % of the samples), never as a pattern signature | annular 0.997, wavy 0.87 |

**Eq. 9** (bottom film estimate) is a piecewise definition. **Eq. 10**: bubble share. **Eq. 11**: film jump fraction.

#### 3.3 Features from the mean pressure drop

- **Why pressure is needed:** a single bottom probe can't always tell a mixed flow from a separated one. A thick, abruptly varying film is seen both in dispersed/churn flow and in low-$j_l$ wavy flow with large waves. At high velocities the interface is lost to aeration and churn flow looks like a thin film. The mean pressure gradient adds an integral, cross-section-wide view.
- Only the mean gradient exists per run. It is compared with two reference flows that use the same friction factor:
  - **Eq. 12:** friction factor $f = 64/Re$ for $Re$ < 2300, otherwise $f = 0.25[\log_{10}(150.39/Re^{0.98865} - 152.66/Re)]^{-2}$ (smooth pipe; `friction_f` in `src/dpdz_correlations.py`).
  - **Eq. 13:** homogeneous mixture (McAdams): $\rho_h = [x/\rho_g + (1-x)/\rho_l]^{-1}$, $\mu_h = [x/\mu_g + (1-x)/\mu_l]^{-1}$, $Re = G d/\mu_h$, with $G = G_l + G_g$.
  - **Eq. 14:** $(dp/dz)_h = f\,G^2/(2 d \rho_h)$.
  - **Eq. 15:** liquid alone: $(dp/dz)_l = f(Re_l)\,G_l^2/(2 d \rho_l)$ with $Re_l = G_l d/\mu_l$.
  - **Eq. 16:** $\mathrm{ratio\_homog} = (dp/dz)_{tp}/(dp/dz)_h$ and $\phi_l^2 = (dp/dz)_{tp}/(dp/dz)_l$.
  - Square channel: $d$ = 8 mm is both the side and the hydraulic diameter.
- **Physical reading of `ratio_homog`.** The homogeneous model assumes both phases move at the same velocity, as a well-mixed fluid.
  - **≥ 1:** the measured gradient reaches the no-slip prediction, so the phases are strongly mixed and exchange momentum (churn, dispersed). Medians: dispersed 1.41, churn 1.32.
  - **< 1:** the gas slips over a separated liquid layer and friction is lower than if the phases were mixed (wavy). Median wavy 0.75; low-$j_l$ wavy and slug runs 0.70–0.99.
  - Annular runs have high values (median 2.9). The thin film and the fast gas core create more wall shear than the homogeneous model predicts. This is why annular is identified first, from the film (Section 3.4).
- **Physical reading of $\phi_l^2$** (Lockhart–Martinelli two-phase multiplier): how many times the gradient exceeds that of the liquid flowing alone.
  - **Small:** the liquid dominates the friction and the gas is carried as bubbles (dispersed, median 4.3).
  - **Large:** the gas adds much of the friction (churn 9.3, wavy 26, annular 254).

#### 3.4 Decision rules

The rules are applied in order, and the first rule that matches assigns the pattern. **Fig. 11** shows the decision tree (`paper_figures/fig11_decision_tree.tex`, TikZ, thresholds read from `run_manifest.json`). **Table 9** lists the rules and thresholds (FPC §4.1–4.3). For each rule, give the condition and the physical reading:

**Step 1. Family: is the channel bridged by liquid?**
- **Intermittent** if `bridge_freq_hz` ≥ 0.1 Hz.
- **Physical reading:** in intermittent flows, liquid slugs or plugs, or the continuous liquid of bubbly flow, fill the channel at the probe again and again. In separated flows the liquid stays at the bottom and the probe always finds an interface (film) or nothing (undetected).
- The threshold sits in a wide gap: separated runs have 0 Hz (a few spurs up to about 0.05 Hz), and the slowest slug runs have about 0.13 Hz (Fig. 21a).

**Step 2a. Intermittent family.**

| Pattern | Condition | Physical reading |
|---|---|---|
| **Bubbly** | `bubble_share` ≥ 0.5 | Between bridges, the interfaces seen are mostly **bubbles** at mid-height, not a wall film. The liquid is continuous and the gas is dispersed in it, so there are no elongated bubbles with a film below. |
| **Plug** | `bridge_mean_duration_s` ≥ 0.02 s **and** `f_liquid` ≥ 0.10 | **Long, unaerated liquid plugs** occupying at least 10 % of the time. At low $j_g$ the liquid between elongated bubbles stays continuous. |
| **Slug** | otherwise | **Short, aerated liquid slugs.** The higher $j_g$ entrains bubbles into the slug, which break the liquid bridge into short pieces. Liquid is seen < 1 % of the time (median). |

**Step 2b. Separated family: no liquid bridges.**

| Pattern | Condition | Physical reading |
|---|---|---|
| **Annular** | `film_est_mean_mm` < 0.1 mm **and** `film_est_std_mm` < 0.2 mm | **Thin and uniform bottom film,** permanently below the detection limit (seen only through `t1`). The fast gas core ($j_g$ 12.6–25 m/s) shears the film thin and keeps it free of large waves. Runs with no film information at all also fall here, because the film is always below detection. |
| **Dispersed** | agitated or liquid-rich film (below) **and** (`n_peaks` ≥ 2 **or** $\phi_l^2$ ≤ 5) | Mixed flow where the **liquid dominates the friction:** the gas is carried as bubbles by a fast liquid ($j_l$ 1.38–3.12 m/s). Two KDE peaks = the wall film plus a layer of bubbles above it. |
| **Churn** | agitated or liquid-rich film **and** not dispersed | Mixed flow where **the gas adds much of the friction** ($\phi_l^2$ > 5): chaotic, oscillating liquid without a stable film or stable slugs. |
| **Wavy** | otherwise | **Film with smooth waves** and gas slipping over it. Includes the visual stratified runs. |

**Agitated film** (the wavy–churn–dispersed transition). Both conditions must hold:
1. **Film signature**, one of:
   - `n_peaks` ≥ 2: interfaces at two distinct heights (wall film + bubble layer);
   - **a real liquid film with abruptly varying thickness:** `film_est_mean_mm` ≥ 0.4 mm, `film_jump_frac` ≥ 0.4, and the film detected in ≥ 15 % of samples.
     - The mean and detection conditions make sure the jumps belong to a real film, not to dropouts near the detection limit: thin wavy films can look jumpy with few detected points (Fig. 22b).
2. **Pressure gate:** `ratio_homog` ≥ 1.0. The flow is at least as dissipative as a well-mixed flow.
   - **Physical reading:** below the homogeneous prediction, a thick jumpy film is a separated flow with large waves and gas slipping over it. These are the low-$j_l$ wavy and slug runs (0.70–0.99), so the run is wavy.

**Liquid-rich film** (high-velocity runs at the churn–annular transition). The run is neither annular nor agitated, and all of the following hold:
- `film_est_mean_mm` < 0.5 mm: the probe sees only a thin, broken film, because at high $j_l$ and $j_g$ the interface is aerated and scatters the light;
- `ratio_homog` ≥ 1.0: the pressure says the flow is mixed, not separated;
- $\phi_l^2$ ≤ 15: friction still dominated by the liquid. Thin-film wavy runs with `ratio_homog` ≥ 1 have $\phi_l^2$ ≥ 17.9 (runs 7, 19, 20, 21, 22, 139), and annular runs are far higher. Runs 10–12 ($\phi_l^2$ = 9.3–11.9) fall in this branch: they were first logged as wavy-annular and are churn after reviewing the photos, in agreement with the criterion.
- **Physical reading:** the bottom probe is "blinded" by aeration, so the film looks thin, as in wavy flow. The pressure gradient reveals a mixed, liquid-rich flow, so the run is churn or dispersed. Examples: runs 1, 4, 10, 11, 12, 138 (Fig. 22d, Fig. 17).

**Why this order:**
- Intermittency is the clearest signature, so it is tested first.
- Annular is tested before the pressure gates, because annular runs also have `ratio_homog` > 1 (3.3).
- Wavy is the default for any separated flow without mixing evidence.

#### 3.5 Evaluation of the criterion
- **Reference:** the visual classification of the 157 runs mapped to 7 classes (Table 6).
- **Eq. 17:**
  - agreement $= N_\mathrm{correct}/N$;
  - balanced accuracy $= \frac{1}{K}\sum_k \mathrm{recall}_k$, the mean recall per class.
  - Why both: the classes are unbalanced (wavy 44 runs, annular 10). Balanced accuracy keeps the large classes from dominating.
- Honest statement: the thresholds were chosen on this same dataset. Robustness is shown by:
  - plateaus in the sensitivity analysis (Fig. 25);
  - re-tuning all thresholds on a grid does not improve the defaults: in-sample the best balanced accuracy equals the default (89.4 %), and with leave-one-out validation it drops to 85.6 % (`variant_comparison_sensitivity_grid.csv`).
- The comparison between data variants (regular, trimmed only, processed, processed_v2) uses the same rules and thresholds.

### 4. Results (old Section 3)

Intro paragraph (edit). Three items:
- (i) the film thickness profile per pattern;
- (ii) film thickness and PDF per pattern;
- (iii) the objective criterion: feature planes, agreement, disagreements and sensitivity.

#### 4.1 Liquid film thickness profile (edit / regenerate)

Must contain, for each pattern: what the profile shows, and **which feature of the criterion it illustrates**. Profiles are produced from processed_v2 by `notebooks/paper_figures.ipynb` (`fig12a_...` to `fig17c_...`), with the same runs and time windows as the current paper (0.5 s bubbly, 1 s otherwise).

| Fig. | Pattern | Points to make |
|---|---|---|
| 12 (old 10) | Bubbly | `t2` mostly at 8 mm (liquid state), interrupted by bubbles at 2–6 mm. Practically no wall film, hence high bubble share. At higher $j_l$ the bubbles are closer together and less regular (keep the current comparison). |
| 13 (old 11) | Plug | Long liquid bridges at 8 mm alternating with elongated bubbles. Below each bubble the film drains from about 2.2 mm (nose) to about 0.36 mm (tail). Bubbles take about 0.15 s to pass. Long bridges, so plug. At higher $j_l$ the bubbles are shorter, close to bubbly. |
| 14 (old 12) | Slug | Same structure, but the bridges are short and broken by entrained bubbles, and the film varies more. This is why the mean bridge duration falls below 20 ms. Note: run 67 (14a) is classified **plug** by the criterion (Table 12); run 56 (14b) is slug. |
| 15 (old 13) | Dispersed | No bridges. A film of about 0.5–1.5 mm with frequent abrupt jumps (high jump fraction). Describe it as the wavy–churn–dispersed transition, and keep the remark that no tube-height reading appears. Note: runs 8 and 9 are classified **churn** by the criterion ($\phi_l^2$ just above 5, Table 12); runs 74–77, 80, 81 are dispersed for both. |
| 16 (old 14) | Wavy | A continuous thin film with smooth waves, often below detection (gaps). The old text about the spur at 0.95 s changes: it was a mixer wave, now removed (Fig. 8). |
| 17 (new) | Churn and annular | Churn: few detected points, a broken film, and much of the film only visible through `t1`. Annular: almost everything below detection, film from `t1` about 0.05 mm, uniform. Files: churn run 93 ($j_l$ = 1.11, $j_g$ = 3.95 m/s), churn run 1 (liquid-rich branch, $j_l$ = 0.92, $j_g$ = 5.83 m/s), annular run 85 ($j_l$ = 0.14, $j_g$ = 20.95 m/s); 1 s window with the most detected samples. Both signals are drawn (detected `t2`, and the film from `t1` in orange) to show that "missing" data still carries the film. |

Optional: draw the 2.0 mm and 7.5 mm state thresholds (Table 7) as dashed lines on all profiles.

#### 4.2 Mean liquid film thickness by flow pattern (regenerate)
- **Fig. 18:** bottom film estimate per pattern, 7 classes, 157 runs. File: `paper_figures/fig18_bottom_film_by_pattern.pdf`.
  - Why this one: it excludes liquid bridges and includes films below detection, so it is a "film" thickness.
  - The old Fig. 15 mixed films with bridges: bubbly 6.65, plug 4.47, slug 2.47 mm.
- Points to make (medians):
  - **bubbly 3.0 mm:** not a film; these are the heights of the bubbles;
  - plug 1.0, slug 1.1 and dispersed 1.05 mm;
  - churn 0.45 mm;
  - wavy 0.27 mm;
  - annular 0.055 mm.
  - In the separated flows the film thins as $j_g$ grows (gas shear): dispersed and churn > wavy > annular.
- **Fig. 19** (old Fig. 16, zoom, `paper_figures/fig19_bottom_film_separated_zoom.pdf`): the separated patterns only (dispersed, churn, wavy, annular), same colours as Fig. 18.
  - Wavy has the widest spread because it covers the widest range of conditions.
  - Annular is below the 0.5 mm manufacturer minimum. This is why the `t1` correction is essential.
- Keep the remark that a sensor with a thinner range would be needed for annular films.

#### 4.3 Probability density function of the film thickness (regenerate)
- **Fig. 20:** one multi-panel figure replacing old Figs. 17–21, from processed_v2 (`paper_figures/fig20_kde_by_pattern.pdf`; single panels `fig20a_...` to `fig20f_...`). One panel per pattern with the runs of the current paper, **adding churn** (runs 26, 27, 94, 138, 1, 140). Run 40 (wavy) is now outside the evaluated region and is replaced by run 45. For annular, say that almost no `t2` is detected (99.7 % undetected, median), so there is no meaningful KDE.
- Points to make, each connected to a feature:
  - **Bubbly:** a dominant peak at 8 mm (liquid) and a broad one at 2–6 mm (bubbles). This motivates `bubble_share`.
  - **Plug and slug:** two peaks, film (0.5–1.5 mm) and channel height (8 mm), whose relative height depends on conditions. This motivates the bridge statistics; `n_peaks` = 2 for plug.
  - **Dispersed:** an irregular curve with several peaks between 1 and 2 mm (film plus bubbles). This motivates `n_peaks` ≥ 2 and the jump fraction.
  - **Wavy:** a single peak at the film thickness, lower and more regular. **Churn:** similar to wavy, single peak (median `n_peaks` = 1).
- Conclude: the PDF separates the intermittent patterns well, but **not churn, dispersed and wavy**. This motivates the film-dynamics and pressure features of Section 3.
- Remove "it is possible to use those as an objective criteria".

#### 4.4 Feature planes of the criterion (new, ready)
Files in `data/gold_layer/paper_figures/` (`fig21a_...` to `fig22d_...`, re-plotted from `objective_criteria_5`). Each figure has the thresholds as dashed lines, and the runs coloured by visual pattern.

**Fig. 21: intermittent family.**
- (a) `fig21a_criterion_family.pdf`, bridge frequency vs bottom-film CV:
  - all intermittent runs lie right of 0.1 Hz;
  - separated runs pile up at 0;
  - the exceptions are the disagreements: dispersed 54 and 55 and wavy 110 have bridges; slug 15, 16, 24, 91, 97, 158 and plug 99 have none or almost none.
- (b) `fig21b_criterion_intermittent.pdf`, bubble share vs mean bridge duration:
  - bubbly is isolated at bubble share > 0.6;
  - plug sits above 20 ms and slug below.

**Fig. 22: separated family.**
- (a) `fig22a_criterion_separated.pdf`, film mean vs std:
  - annular sits in the corner below 0.1 / 0.2 mm;
  - wavy has a larger std for the same mean than dispersed (waves vs a thick, mixed film).
- (b) `fig22b_criterion_dispersed.pdf`, film mean vs jump fraction:
  - dispersed and churn are in the "thick and jumpy" quadrant;
  - several wavy runs are there too. This is why a pressure gate is needed.
- (c) `fig22c_criterion_agitated.pdf`, `ratio_homog` vs $\phi_l^2$:
  - the gate at 1.0 removes the jumpy wavy and slug runs;
  - $\phi_l^2$ = 5 separates dispersed (below) from churn (above).
- (d) `fig22d_criterion_liquid_rich.pdf`, same axes for the thin-film runs (only the runs discussed in the text are labelled):
  - the high-velocity runs (visual churn 1, 4, 10, 11, 12, 138; visual dispersed 9, 134, 136) have ratio > 1 and $\phi_l^2$ < 15;
  - every thin-film wavy run has ratio < 1 or $\phi_l^2$ ≥ 17.9.

Optional **Table 10:** median feature values per class (the medians in Table 8 above), from `objective_flow_pattern.csv`, with stratified mapped to wavy.

#### 4.5 Agreement with the visual classification (new)
- **Fig. 23:** `fig23_confusion_matrix.pdf` (157 runs, same pattern order on both axes). Points to make:
  - bubbly 11/11, churn 14/14, annular 10/10;
  - wavy 43/44, plug 30/33, slug 18/25, dispersed 13/20;
  - the errors sit next to the diagonal, between physically neighbouring patterns: plug↔slug, slug→wavy, dispersed→churn or slug, wavy→slug (run 110).
- **Fig. 24:** visual vs objective map side by side.
  - (a) `fig24a_flow_pattern_map_visual.pdf`: visual classification of the 157 runs with the 7 reference classes (runs 10–12 as churn).
  - (b) `fig24b_flow_pattern_map_objective.pdf`, same markers.
  - Point to make: without using $j_l$ or $j_g$, the objective patterns fall in coherent regions consistent with Coleman and Garimella. Dispersed sits at high $j_l$, annular at the highest $j_g$, and plug at low $j_g$.
- **Table 11:** agreement, balanced accuracy and per-class recall per data variant (`variant_comparison_criterion.csv`):

| Variant | Agreement | Balanced accuracy |
|---|---|---|
| regular (untreated) | 77.1 % | 78.2 % |
| regular, trimmed only | 88.5 % | 89.4 % |
| processed | 88.5 % | 89.4 % |
| processed_v2 | 88.5 % | 89.4 % |

  Point to make: the gain comes from trimming. Without it, the leading blocks with no signal distort every fraction (bridges, film detection), and plug, slug, churn and dispersed are the most affected (plug recall 0.64 → 0.91, churn 0.71 → 1.0).

#### 4.6 Disagreements, sensitivity and limitations (new)
**Table 12:** the 18 disagreements grouped by cause (`objective_criteria_5/processed_v2/reclassified_runs.csv`). For each group, the physical explanation:

| Visual → objective | Runs | Physical explanation |
|---|---|---|
| dispersed → slug | 54, 55 | liquid bridges reach the probe at $j_g$ ≈ 1.2 m/s, at the bubbly–intermittent transition |
| dispersed → churn | 8, 9, 133, 134, 136 | $\phi_l^2$ = 5.2–7.7, just above the dispersed threshold of 5. At the churn–dispersed transition the gas still adds part of the friction. |
| plug → slug | 35, 36 | mean bridge duration 12–15 ms < 20 ms: plugs short enough to look like slugs at the plug–slug transition |
| plug → wavy | 99 | very rare bridges (0.02 Hz) in a 100 Hz record |
| slug → plug | 67 | long bridges (25 ms) at $j_l$ = 0.97 m/s |
| slug → wavy | 15, 16, 24, 91 | $j_g$ > 3 m/s. The slugs no longer reach the bottom probe as liquid bridges (aerated or bypassing) |
| slug → wavy | 97, 158 | short records (< 1 000 samples): slow slugs not caught |
| stratified-wavy → slug | 110 | occasional bridges (0.17 Hz) at 100 Hz: large waves touching the top wall, or unresolved mixer waves |

- **Fig. 25:** `fig25_sensitivity_thresholds.pdf`. Each threshold sits on a plateau of balanced accuracy (e.g. the liquid-rich $\phi_l^2$ limit gives the same result from 15 to 17, the liquid-rich film limit from 0.40 to 1.0 mm), and jointly re-tuning all thresholds does not improve the defaults (leave-one-out 85.6 % vs 89.4 %). Can go to an appendix if space is short.
- **Limitations** to state:
  - one probe on the bottom wall: a thin stratified layer and the bottom of a wavy film look alike, and the top film of annular flow is not seen;
  - mixer waves and short slugs can't be resolved at 100 Hz;
  - the pressure gate at 1.0 is close to the highest low-$j_l$ wavy run (0.99), so the uncertainty matters (2.3);
  - the thresholds were tuned and evaluated on the same 157 runs;
  - the visual reference itself is uncertain at transitions.

### 5. Conclusions (edit)

Must contain:
1. **Technique (keep, update):** confocal chromatic microscopy resolves the film and the liquid bridges in all intermittent patterns. It is most reliable for plug, slug and bubbly. Films below about 0.3–0.5 mm, or steep or aerated interfaces, aren't detected directly; a sensor with a thinner range is needed for annular films (keep the turndown-ratio remark).
2. **Treatment (new):**
   - trimming the leading no-signal blocks is the decisive step (77.1 % → 88.5 %);
   - the wall reading recovers films below detection;
   - mixer waves can be removed at ≥ 1.3 kHz;
   - missing readings carry information and no run needs to be discarded.
3. **Objective criterion (new):**
   - liquid-bridge statistics separate intermittent from separated flows, and bubbly, plug and slug;
   - the bottom-film mean, variability and abruptness identify annular and agitated films;
   - the mean pressure gradient (ratio to homogeneous, $\phi_l^2$) separates churn, dispersed and wavy where one bottom probe cannot.
   - Result: 88.5 % agreement, 89.4 % balanced accuracy, with disagreements at transitions.
4. **PDF (rewritten):** the PDF shapes support the intermittent rules, but they aren't enough for the separated patterns.
5. **Map (rewritten):** the objective map agrees with the visual one and with Coleman and Garimella; the disagreements concentrate at the transitions.
6. **Future work:** a second probe on the top wall, a dynamic pressure signal, and one sampling rate (≥ 1.3 kHz) for all runs.

### Appendix (optional)
- Full feature list (FPC §3.2) and the full threshold table.
- `objective_criteria_5/processed_v2/mutual_information.pdf`: feature relevance, supporting the choice of features. Replace the code names with symbols first.
- Fig. 25, if moved out of 4.6.

---

## 4. Figure map

All data figures of the revised paper are produced by `notebooks/paper_figures.ipynb` in **`data/gold_layer/paper_figures/`**. The notebook reads processed_v2 and the `objective_criteria_5` results, does not re-run the criterion, and checks that the classification is unchanged (88.5 % / 89.4 %).

| New | Old | Content | File in `data/gold_layer/paper_figures/` | Status |
|---|---|---|---|---|
| 1–5 | 1–5 | facility, photo, confocal schematic, principle, thickness indicatives | – | keep |
| 6 | 6 | tube thickness before treatment (157 runs, 7 classes) | `fig06_tube_thickness_before_treatment.pdf` | done |
| 7 | 7 | tube thickness after treatment (157 runs, 7 classes) | `fig07_tube_thickness_after_treatment.pdf` | done |
| 8 | – | mixer-wave removal, run 46 | `fig08_mixer_wave_removal_run46.pdf` | done |
| 9 | 8 | evaluated region | `fig09_evaluated_region.pdf` | done |
| 10 | 9 | visual definition | camera images | edit in LaTeX (remove (d) stratified) |
| 11 | – | decision-tree flowchart (TikZ, thresholds of version 5 = version 4) | `fig11_decision_tree.tex` / `.pdf` | done |
| 12–16 | 10–14 | bubbly, plug, slug, dispersed, wavy profiles (same runs: 65, 59; 33, 28; 67, 56; 8, 9; 45, 44) | `fig12a_...` to `fig16b_...` | done |
| 17 | – | churn (runs 93, 1) and annular (run 85) profiles, with the film from $t_1$ | `fig17a_...`, `fig17b_...`, `fig17c_...` | done |
| 18 | 15 | bottom film per pattern | `fig18_bottom_film_by_pattern.pdf` | done |
| 19 | 16 | zoom on separated patterns | `fig19_bottom_film_separated_zoom.pdf` | done |
| 20 | 17–21 | PDFs per pattern, incl. churn (multi-panel and single panels) | `fig20_kde_by_pattern.pdf`, `fig20a_...` to `fig20f_...` | done |
| 21 | – | intermittent family (a, b) | `fig21a_criterion_family.pdf`, `fig21b_criterion_intermittent.pdf` | done |
| 22 | – | separated family (a–d) | `fig22a_...` to `fig22d_...` | done |
| 23 | – | confusion matrix (same order on both axes) | `fig23_confusion_matrix.pdf` | done |
| 24 | – | visual vs objective map | `fig24a_flow_pattern_map_visual.pdf`, `fig24b_flow_pattern_map_objective.pdf` | done |
| 25 | – | threshold sensitivity | `fig25_sensitivity_thresholds.pdf` | done |
| – | 9(d) | stratified image | – | remove |

Note on the example runs: the dispersed profiles (runs 8, 9, Fig. 15) are classified as **churn** by the criterion, and the slug profile of Fig. 14a (run 67) as **plug** (Table 12). Either say so in the text or pick other examples (`PROFILE_RUNS` in the notebook). Runs classified dispersed include 74–77, 80, 81 (0318); slug runs classified slug include 13, 14, 25.

## 5. Table map

| New | Old | Content | Source | Status |
|---|---|---|---|---|
| 1–3 | 1–3 | sensors, empty tube, first-thickness threshold | – | keep |
| 4 | 4 | low-percentile film thickness | – | edit (caption: justification of 0.373 mm, not a filter) |
| 5 | – | mixer-wave parameters | `run_manifest.json` → `mixer_wave_params` | to create |
| 6 | – | visual labels → reference classes (§1 of this plan) | `experimental_results.csv` | to create |
| 7 | – | sample states | Section 3.1 above | to create |
| 8 | – | features, definitions, physical reading | Section 3.2 above | to create |
| 9 | – | rules and thresholds | Section 3.4 above; `run_manifest.json` | to create |
| 10 | – | median feature values per class (optional) | `objective_criteria_5/processed_v2/objective_flow_pattern.csv` | to create |
| 11 | – | agreement per data variant | `objective_criteria_5/variant_comparison_criterion.csv` | to create |
| 12 | – | disagreements by cause | `objective_criteria_5/processed_v2/reclassified_runs.csv` | to create |

## 6. Equation map

| New | Old | Content |
|---|---|---|
| 1–4 | 1–4 | orifice flow, air density, superficial velocities (keep) |
| 5–8 | 5–8 | tube-thickness correction (keep) |
| 9 | – | bottom film estimate |
| 10 | – | bubble share |
| 11 | – | film jump fraction |
| 12 | – | friction factor |
| 13 | – | homogeneous $x$, $\rho_h$, $\mu_h$, $Re$ |
| 14 | – | homogeneous gradient |
| 15 | – | liquid-alone gradient |
| 16 | – | `ratio_homog` and $\phi_l^2$ |
| 17 | – | agreement and balanced accuracy |
| (opt.) | – | uncertainty of `ratio_homog` |

## 7. Numbers cheat-sheet

| Quantity | Value |
|---|---|
| Runs: campaign / evaluated / alternative without stratified | 164 / **157** / 154 |
| Excluded ($j_l$ < 0.035 m/s) | 38, 39, 40, 49, 104, 105, 164 |
| Visual stratified in the evaluated set (evaluated as wavy) | 48, 108, 147 (all classified wavy) |
| Class sizes | bubbly 11, plug 33, slug 25, dispersed 20, churn 14, wavy 44, annular 10 |
| Sampling rates in the evaluated set | 98 at 100 Hz, 50 at 1.3–2 kHz, 9 > 2 kHz |
| Trimmed runs | 68 of 164 (up to 96 % of the file) |
| Runs < 50 % missing `t2`, before → after trimming | 53 → 69 |
| Mixer-wave runs changed | 21 (18 in the evaluated set) |
| Short records (< 1 000 samples) | 10 |
| Agreement / balanced accuracy | regular 77.1 / 78.2 %; trimmed, processed, processed_v2 88.5 / 89.4 %; without the 3 stratified runs 88.3 / 89.4 % |
| Per class | bubbly 11/11, plug 30/33, slug 18/25, churn 14/14, dispersed 13/20, wavy 43/44, annular 10/10 |
| Disagreements | 18 of 157 (Table 12) |
| Threshold grid (balanced accuracy) | default 89.4 %; best in-sample 89.4 %; leave-one-out 85.6 % |
| Version 4 (runs 10–12 still wavy-annular), for reference | 86.6 / 88.5 %, 21 disagreements |
| Tube thickness $t_1$ medians, before → after treatment [mm] | wavy 1.038 → 0.982; annular 1.035 → 0.985; dispersed 0.998 → 0.978; churn 0.998 → 0.976; bubbly, plug, slug 0.977–0.980 → 0.976–0.977 |
| Missing `t2`, dispersed at similar $j_l$, $j_g$ | 74–82 % at ≈ 3.6 kHz vs ≈ 0 % at 100 Hz |
| $j_l$ / $j_g$ range per class [m/s] | bubbly 0.96–2.00 / 0.26–0.80; plug 0.12–0.97 / 0.25–2.09; slug 0.12–1.50 / 0.80–6.50; dispersed 1.38–3.12 / 1.18–7.98; churn 0.60–1.40 / 3.38–10.46; wavy 0.05–0.80 / 2.71–15.76; annular 0.05–0.52 / 12.56–25.16 |

## 8. References to add

Suggestions; check each one before citing.

| Topic | Where | Suggested reference |
|---|---|---|
| Homogeneous viscosity | Eq. 13 | McAdams et al. (1942) |
| Smooth-pipe friction factor (Eq. 12) | Eq. 12 | Fang, Xu and Zhou (2011), Nucl. Eng. Des. The rough-pipe branch of `friction_f` is the same paper's rough-pipe correlation |
| Two-phase multiplier | Eq. 16 | Lockhart and Martinelli (1949) |
| KDE bandwidth | Section 3.2 | Scott (1992) |
| Objective identification from signal PDFs | Introduction | Jones and Zuber (1975) |
| Objective identification with local probes | Introduction | Barnea, Shoham and Taitel (1980) |
| Identification from pressure fluctuations | Introduction | Matsui (1984) |
| Balanced accuracy | Section 3.5 | Brodersen et al. (2010) |

## 9. To-do before submission

### Figures to produce
- [x] Figs. 6, 7: tube thickness, 157 trimmed runs, 7 classes.
- [ ] Fig. 10: remove the stratified panel, re-layout (LaTeX).
- [x] Fig. 11: decision-tree flowchart.
- [x] Figs. 12–16: profiles from processed_v2.
- [x] Fig. 17: churn + annular profiles.
- [x] Fig. 19: zoom on the separated patterns.
- [x] Fig. 20: multi-panel PDFs, including churn.
- [x] Fig. 24(a): visual map of the 157 runs, 7 classes.
- [ ] Decide whether to keep runs 8, 9 (dispersed → churn) and 67 (slug → plug) as profile examples.

### Fixes to existing gold figures (done in `paper_figures/`; the `objective_criteria_5` files are untouched)
- [x] Liquid-rich figure: only the discussed runs are labelled, and the axis is widened.
- [x] Confusion matrix: same order on both axes.
- [x] Objective map: the region labels are drawn above the markers.
- [x] Sensitivity figure: shorter axis labels.
- [ ] Once Tables 8–9 fix the symbols, align the axis labels of Figs. 21–22 and 25 with them.

### Text and data
- [ ] Check the visual label of run 124 ("bubble", counted as plug).
- [ ] Uncertainty of `ratio_homog` (Section 2.3).
- [ ] Footer journal name.
- [ ] Keep "stratified" only where it refers to the visual label.
