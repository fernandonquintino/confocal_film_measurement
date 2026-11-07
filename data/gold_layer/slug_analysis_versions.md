# slug_analysis versions

| Version | Date | Description | Results |
|---|---|---|---|
| slug_analysis_1 | 2026-11-07 | First slug-unit detection from the confocal states (tall readings incl. mixer-flagged grouped within 20 ms; min slug 5 ms; min bubble 20 ms); camera interface velocity (preliminary) | plug/slug median f_s 5.64/4.34 Hz; superseded: fixed 20 ms grouping merges short bubbles at high j_l (runs 58, 66) |
| slug_analysis_2 | 2026-11-07 | Bubble limit as a length: tall readings closer than max(2 D_h / j_m, 2 samples) form one slug; sensitivity grid L_B,min 1/2/4 D_h x min slug 0/5/10 ms; uncertainty flags | plug/slug median f_s 5.63/4.39 Hz; superseded: min slug 0 ms in the grid counts isolated readings |
| slug_analysis_3 | 2026-11-07 | As version 2 with the sensitivity grid restricted to plausible values (min slug 2.5/5/10 ms) | plug/slug median f_s 5.63/4.39 Hz, intermittency 0.62/0.91; 17 runs with uncertain frequency (p90/p10 > 2); camera velocity consistent in 25 of 58 runs |
