# Convective scaling of the time rules: findings (objective_criteria_11)

Produced by `scripts/convective_scaling_eval.py` (Fluids paper, Stage 1/2). Evaluation only: the reported criterion keeps its dimensional time thresholds, and nothing here is in the manuscript.

## Set-up
- Rules replaced one at a time, everything else as reported:
  - family split: f_bridge → St = f_bridge·D_h/U;
  - plug split: t̄_bridge → L* = t̄_bridge·U/D_h.
- U = j_m (jl + jg), j_l, or U_t = 1.2·j_m + 0.54·√(g·D_h). U_t is a Bendiksen-type assumption, not validated for square ducts.
- j is taken at the meter conditions; D_h = H = 8.0 mm.
- Scores: BA / AG of the whole criterion against the visual classes, both in-sample and with the threshold chosen by leave-one-out (LOO).

## Results (`summary.csv`)
| Rule | Variable | Best BA | LOO BA | Near-optimal band ratio | Runs changed |
|---|---|---|---|---|---|
| family | f_bridge (reported, 0.1 Hz) | 88.1 | 87.5 | 18.6 | 0 |
| family | St, U = j_m | 88.1 | 87.2 | 14.5 | 0 |
| family | St, U = j_l | 88.4 | 87.2 | 5.6 | 2 |
| family | St, U = U_t | 88.1 | 87.5 | 14.5 | 0 |
| plug | t̄_bridge (reported, 20 ms) | 88.1 | 88.1 | 2.1 | 0 |
| plug | L*, U = j_m | 86.3 | 84.6 | 3.9 | 4 |
| plug | L*, U = j_l | 85.8 | 85.8 | 18.5 | 6 |
| plug | L*, U = U_t | 87.0 | 85.6 | 1.4 | 2 |

- **Family rule:** St with j_m is an equivalent dimensionless form on these data, but no better. Its threshold (St ≈ 1e-4) marks "bridges recur within the record", not a slug-frequency transition.
- **Plug rule:** convective scaling degrades it. With both rules convective (j_m), runs 13, 55, 57 and 58 move from slug to plug (`per_run.csv`).

## Notes carried forward (Stages 6–8)
1. **Bridge duration is not slug length.** A liquid bridge is the unaerated stretch of a slug, and entrained bubbles split aerated slugs into short bridges.
   - Visual plug vs slug: t̄_bridge medians are 44.5 vs 12.0 ms, overlapping only between 12.5 and 25.3 ms.
   - L* with j_m: 2.24–8.40 vs 2.33–7.19, a near-complete overlap.
   - Slug lengths need slug-unit detection (merging bridges and bubbles), plus a stated velocity assumption (U_t, C0, U_d).
2. **Sampling resolution:**
   - 98/157 runs were recorded at 100 Hz, where mixer waves and short slugs cannot be resolved; kHz records last only 3–7 s.
   - Processing time scales in convective form would fall below one sample: min_dwell in 76/98 and jump_lag in 53/98 of the 100 Hz runs (`processing_time_scales.csv`).
3. **Classification sensitivity:**
   - Run 161 lies 0.0001 mm (1.25e-5 H) from the annular σ_f limit.
   - The R_h and dispersed φ_l² limits move 1–2 runs with small shifts.
   - The family and plug thresholds depend on record length.
