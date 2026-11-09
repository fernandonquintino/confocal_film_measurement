# slug_analysis versions

| Version | Date | Description | Results |
|---|---|---|---|
| slug_analysis_1 | 2026-11-07 | First slug-unit detection from the confocal states (tall readings incl. mixer-flagged grouped within 20 ms; min slug 5 ms; min bubble 20 ms); camera interface velocity (preliminary) | plug/slug median f_s 5.64/4.34 Hz; superseded: fixed 20 ms grouping merges short bubbles at high j_l (runs 58, 66) |
| slug_analysis_2 | 2026-11-07 | Bubble limit as a length: tall readings closer than max(2 D_h / j_m, 2 samples) form one slug; sensitivity grid L_B,min 1/2/4 D_h x min slug 0/5/10 ms; uncertainty flags | plug/slug median f_s 5.63/4.39 Hz; superseded: min slug 0 ms in the grid counts isolated readings |
| slug_analysis_3 | 2026-11-07 | As version 2 with the sensitivity grid restricted to plausible values (min slug 2.5/5/10 ms) | plug/slug median f_s 5.63/4.39 Hz, intermittency 0.62/0.91; 17 runs with uncertain frequency (p90/p10 > 2); camera velocity consistent in 25 of 58 runs |
| slug_analysis_4 | 2026-11-07 | As version 3; figures as PDF only, one panel per file (detection QC split per run in detection_qc/), runs identified by pattern and superficial velocities, confocal time axis labelled as confocal record; length-ratio panel flags runs uncertain in frequency or slug duration | Same numbers as version 3 (plug/slug median f_s 5.63/4.39 Hz, intermittency 0.62/0.91) |
| slug_analysis_5 | 2026-11-08 | As version 4 in the paper figure style (LaTeX text with txfonts, 252 pt column width, 2.2 in height, capitalized labels, intermittency labelled gamma); detection QC at column width with the legend below | Same numbers as versions 3 and 4; figures used in artigo_fluids/confocal_fluids.tex (second draft) |
