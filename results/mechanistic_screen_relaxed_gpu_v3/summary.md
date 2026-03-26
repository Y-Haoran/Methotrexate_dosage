# Mechanistic Screen v1

- candidate matrix: `data/preformulation/paracetamol_mechanistic_screen_demo.csv`
- checkpoint: `checkpoints/best_inference_ckpt.pt`
- current status: `relaxed local compatibility screen`
- confidence tier: `neutral candidates only`
- current use: `shortlist generation, not decision-making`
- next gate: `DFT spot-checks`
- do not over-interpret: `monoatomic ion cases use a zero-reference proxy for isolated ions`

## Neutral Ranking

| candidate_id   | api_name    | polymer_key   | cosolvent_system   |   pre_interaction_energy_eV |   post_interaction_energy_eV |   score_delta_relaxation |   post_api_polymer_min_distance_A |   post_api_api_min_distance_A |   post_solvent_polymer_min_distance_A | reactive_hotspot_risk   |   screening_score_proxy |
|:---------------|:------------|:--------------|:-------------------|----------------------------:|-----------------------------:|-------------------------:|----------------------------------:|------------------------------:|--------------------------------------:|:------------------------|------------------------:|
| screen_001     | paracetamol | pvp           | water              |                     86.7084 |                      2.92662 |                  83.2817 |                           1.39379 |                     nan       |                               8.25983 | low                     |                -1.42662 |
| screen_003     | paracetamol | hpc           | ethanol_water      |                    144.488  |                      9.3717  |                 133.917  |                           1.41984 |                     nan       |                               5.7167  | medium                  |                -7.0717  |
| screen_002     | paracetamol | peg           | water              |                    126.578  |                     13.4666  |                 112.111  |                           1.68123 |                     nan       |                               8.05768 | medium                  |               -12.1666  |
| screen_004     | paracetamol | pectin        | glycerol_water     |                    233.448  |                     24.5023  |                 207.246  |                           1.32266 |                     nan       |                               6.42648 | low                     |               -21.7023  |
| screen_006     | nan         | pvp           | water              |                    243.995  |                     34.6862  |                 208.559  |                         nan       |                     nan       |                               3.22959 | low                     |               -31.4362  |
| screen_005     | paracetamol | pvp           | water              |                    288.478  |                     47.3797  |                 239.298  |                           1.3289  |                       1.38711 |                               2.5132  | low                     |               -42.4797  |

## Ionic Ranking

| candidate_id   | api_name         | polymer_key   | cosolvent_system   |   pre_interaction_energy_eV |   post_interaction_energy_eV |   score_delta_relaxation |   post_api_polymer_min_distance_A |   post_ion_polymer_min_distance_A | reactive_hotspot_risk   |   screening_score_proxy |
|:---------------|:-----------------|:--------------|:-------------------|----------------------------:|-----------------------------:|-------------------------:|----------------------------------:|----------------------------------:|:------------------------|------------------------:|
| screen_007     | diclofenac_anion | pvp           | water              |                    -4237.69 |                     -4389.68 |                   150.39 |                           1.35863 |                            3.8753 | low                     |                 4393.88 |

## Score Components

| candidate_id   | ranking_group   |   score_interaction_term |   score_api_polymer_term |   score_api_api_term |   score_polymer_polymer_term |   score_hotspot_term |   screening_score_proxy |
|:---------------|:----------------|-------------------------:|-------------------------:|---------------------:|-----------------------------:|---------------------:|------------------------:|
| screen_007     | ionic           |               4389.68    |                      4.2 |                 -0   |                         0    |                   -0 |              4393.88    |
| screen_001     | neutral         |                 -2.92662 |                      1.5 |                 -0   |                         0    |                   -0 |                -1.42662 |
| screen_003     | neutral         |                 -9.3717  |                      2.3 |                 -0   |                         0    |                   -0 |                -7.0717  |
| screen_002     | neutral         |                -13.4666  |                      1.3 |                 -0   |                         0    |                   -0 |               -12.1666  |
| screen_004     | neutral         |                -24.5023  |                      2.8 |                 -0   |                         0    |                   -0 |               -21.7023  |
| screen_006     | neutral         |                -34.6862  |                      0   |                 -0   |                         3.25 |                   -0 |               -31.4362  |
| screen_005     | neutral         |                -47.3797  |                      7.7 |                 -2.8 |                         0    |                   -0 |               -42.4797  |

## Notes

- This is a local-cluster proxy screen, not a rigorous free-energy calculation.
- `screening_score_proxy` is now decomposed into explicit interaction, aggregation, cohesion, and hotspot terms.
- Neutral and ionic candidates are separated because the ionic branch is not yet cross-comparable to the neutral branch.
- HPMC and Carbopol-like proxies are not yet wired into this first runnable screen; current supported polymer keys come from the existing descriptor workflow.
