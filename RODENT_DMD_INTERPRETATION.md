# Rodent DMD Interpretation Notes

These notes summarize how the current rodent DMD analysis connects to the operator-response interpretation used in the human gait manuscript and to prior mouse neurodevelopmental phenotyping work.

## Framing

The human manuscript treats delay-embedded DMD models as compact kinematic response surrogates. The fitted operator is not interpreted as a literal neural controller or a direct perturbation-response measurement. Instead, its metrics describe comparative response geometry: normalized non-normality reflects non-orthogonal transient structure, transient growth reflects finite-horizon amplification, resolvent gain reflects model-implied sensitivity to periodic inputs, participation ratio reflects effective dimensionality, and resolvent input-output maps describe where modeled responses are routed.

The same interpretation is appropriate for the rodent data. Here the state space is centroid geometry: each paw is represented by radius from the paw centroid and angle about that centroid. This makes the analysis less about absolute position and more about the coordination geometry of the limbs around the body-centered paw configuration.

## Prior Phenotype Expectations

The main reference point is Klibaite et al., "Deep phenotyping reveals movement phenotypes in mouse neurodevelopmental models" (Molecular Autism, 2022; https://doi.org/10.1186/s13229-022-00492-8). That study used open-field pose tracking to compare Cntnap2 knockout mice and L7-Tsc1 mutant mice, both motivated as neurodevelopmental/autism-related models.

Their key gait results were:

- Both Cntnap2 knockout and L7-Tsc1 mutant mice showed gait defects, including forelimb lag.
- Wild-type mice transitioned from walking to trotting around 0.1 m/s.
- Cntnap2 knockout mice had a similar walk-to-trot transition speed to littermates, but showed larger opposite-limb phase differences during trotting.
- L7-Tsc1 mutants used walking/ambling more strongly and transitioned to trotting at a higher speed, around 0.2 m/s.
- L7-Tsc1 mutants showed the larger phase disruption, with opposite-limb phase differences around 1.0 rad at 0.2 m/s compared with about 0.5 rad in littermates.
- Cntnap2 and Tsc1 phenotypes were not identical: Cntnap2 was closer to a hyperactivity / fast-locomotion phenotype, whereas L7-Tsc1 showed slower, less coordinated locomotion and weaker behavioral adaptation.

Thus, before looking at DMD metrics, the expectation is not simply "mutants are worse." A more specific expectation is that both models may alter phase coordination, but TSC/L7-Tsc1 should show a stronger disruption of gait organization, speed regime, and inter-limb coordination.

## What The Current DMD Analysis Adds

The current results are broadly consistent with that expectation. The largest separation is not a uniform genotype effect inside each line; instead, the major operator-level difference is between the Cntnap and TSC families.

The TSC groups show higher model-implied frequency response. Peak log10 resolvent gain is higher in TSC groups than in Cntnap groups:

- Cntnap controls: 3.09 +/- 0.03 SEM
- Cntnap homozygotes: 3.09 +/- 0.03 SEM
- TSC controls: 3.30 +/- 0.04 SEM
- TSC heterozygotes: 3.37 +/- 0.04 SEM
- TSC homozygotes: 3.24 +/- 0.03 SEM
- new_preds: 3.11 +/- 0.01 SEM

The 0.75-2x stride-frequency band shows the same pattern: Cntnap groups and new_preds sit near 2.64-2.75, whereas TSC groups sit higher, around 2.86-2.90. Interpreted cautiously, TSC locomotion requires a fitted operator with stronger amplification near stride-related frequencies. That aligns with the prior report that TSC mutants show altered gait regime and stronger phase disruption.

The input-output structure also separates the lines. TSC groups have larger mean off-diagonal IO, fore-hind coupling, and left-right coupling than Cntnap groups:

- Mean off-diagonal IO is about 0.19 in Cntnap groups, 0.20-0.22 in TSC groups, and 0.20 in new_preds.
- Fore-hind coupling is about 0.205 in Cntnap groups, 0.214-0.231 in TSC groups, and 0.212 in new_preds.
- Left-right coupling is about 0.205-0.209 in Cntnap groups, 0.214-0.232 in TSC groups, and 0.218 in new_preds.

This suggests that the TSC operators are less self-contained and more distributed across limb blocks. In the language of the human manuscript, this is a reorganization of transient coordination pathways rather than a scalar stability loss.

The centroid-geometry decomposition gives an additional interpretation. Some models show visually clean angle-to-angle structure, while others show more angle/radius cross-coupling. In the group summaries, TSC has stronger angle-to-radius coupling than Cntnap:

- Cntnap angle-to-radius coupling: about 0.13-0.16
- TSC angle-to-radius coupling: about 0.26-0.29
- new_preds: about 0.17

Angle-to-angle coupling is also higher in TSC than Cntnap, but the stronger distinction is the increased coupling from angular coordination into radial excursions. This may be a DMD-level signature of less clean phase coordination: the fitted dynamics cannot describe limb cycling as mostly angular progression around the centroid, and instead recruit radial amplitude changes as part of the response pathway.

## Relation To Non-Normality And Dimensionality

The current normalized Henrici/non-normality summaries are not the strongest group discriminator. Cntnap controls and heterozygotes are around 0.55, while TSC and new_preds are closer to 0.53-0.54. This reinforces the manuscript's caution: non-normality alone is not a severity score. Its interpretation depends on what it does in the fitted operator.

In this dataset, the more interpretable pattern is that TSC has higher resolvent gain and more distributed IO despite similar or slightly lower normalized non-normality. That is similar in spirit to the human stroke result where the consequence of non-normal geometry depended on group context. The important biological claim is not "higher non-normality means worse gait"; it is that different groups organize transient response geometry differently.

Participation ratio is also informative but not definitive. Cntnap homozygotes have higher participation ratio than Cntnap controls, whereas TSC heterozygotes and controls are lower and TSC homozygotes partially rebound. This may reflect different speed distributions and bout composition as much as genotype. It should be interpreted with speed-matched or hierarchical modeling rather than as a simple rank ordering.

## Working Interpretation

The DMD analysis appears to recover a coordination-level distinction that is compatible with the deep-phenotyping paper:

- Cntnap models look closer to the `new_preds` aggregate and to their own controls in operator-response structure. Their known phenotype may be more about activity level, speed occupancy, and phase offsets than a large reorganization of the fitted centroid-geometry response operator.
- TSC models show stronger resolvent gain and more distributed IO coupling, especially angle-to-radius coupling. This is consistent with a locomotor pattern that is less cleanly organized as phase progression and more dependent on coupled angular/radial corrections.
- The TSC pattern fits the prior observation of more persistent slow/ambling locomotion, delayed walk-to-trot transition, and larger phase offsets.

The most useful narrative is therefore not a generic ASD-mouse impairment story. A better narrative is:

> DMD reveals that neurodevelopmental gait phenotypes differ in the response geometry required to reproduce limb coordination. Cntnap and TSC models both have reported gait-phase abnormalities, but TSC locomotion is represented by operators with stronger stride-frequency amplification and more distributed angle/radius input-output coupling. This suggests a shift away from clean angular limb cycling toward a response structure in which angular phase and radial limb excursions are more strongly mixed.

## Caveats

These are fitted-operator metrics, not direct measurements of neural coupling, muscle synergies, or externally evoked perturbation responses. The current models also combine bouts within subject-speed conditions and use centroid geometry with Savitzky-Golay smoothing. Absolute values may change with state space, delay length, rank, smoothing, and speed-binning choices.

The strongest next robustness checks would be:

- repeat Figure 3 across centroid_angles, centroid_geometry, and centroid_geometry_sincos;
- compare speed-matched subsets within each line;
- model group effects with subject-level hierarchy and speed covariates;
- test whether angle-to-radius coupling specifically tracks the phase-lag phenotype from the original paper.

## Sources

- Human gait interpretation: `../Gait_NonNormality_Paper/RSI_paper/main.tex`
- Klibaite U, Kislin M, Verpeut JL, Bergeler S, Sun X, Shaevitz JW, Wang SSH. Deep phenotyping reveals movement phenotypes in mouse neurodevelopmental models. Molecular Autism. 2022;13:12. https://doi.org/10.1186/s13229-022-00492-8
- Current local summaries: `results/rodent_group_summary.csv`, `results/rodent_io_summary_bootstrap.csv`, `results/rodent_resolvent_summary_bootstrap.csv`
