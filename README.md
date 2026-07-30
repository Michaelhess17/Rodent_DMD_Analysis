# Rodent DMD Analysis

This is a small, standalone analysis project for applying `DMDAnalysis` to the
rodent bout data stored in the original `Gait_NonNormality_Paper/data` folder.
It intentionally keeps the rodent analysis separate from the human paper repo.

## Key Design

The DMD training unit is a subject by speed-condition model. Each model pools
multiple bouts from the same subject and speed bin, but transitions are only
formed within each bout. No `X_t -> X_{t+1}` pair crosses a bout boundary.

Default data source:

```text
../Gait_NonNormality_Paper/data
```

Override it with:

```bash
RODENT_DATA_DIR=/path/to/data julia --project=. scripts/run_all.jl
```

## Run

Fast sampled run:

```bash
julia --project=. scripts/run_all.jl
```

`run_all.jl` uses compact defaults (`--n-delays 6 --n-angles 48`) so the whole
loader, pooled DMD fit, summaries, and direct figures can be checked quickly.
Pass arguments to override those defaults.

Full run over all bout CSVs:

```bash
julia --project=. scripts/00_fit_rodent_dmd_metrics.jl --full
julia --project=. scripts/01_make_rodent_figures.jl
```

Useful options:

```bash
julia --project=. scripts/00_fit_rodent_dmd_metrics.jl \
  --dataset all \
  --max-per-group 250 \
  --max-per-subject 30 \
  --state-space centroid_xy \
  --normalize center \
  --speed-condition-mode subject \
  --speed-conditions 3 \
  --min-bouts-per-model 5 \
  --min-strides-per-model 8 \
  --n-delays 26 \
  --rank 40 \
  --n-angles 240 \
  --checkpoint-interval 25 \
  --quiet-progressbar \
  --log-interval 1
```

Subject-speed conditions are excluded unless they pass the model-level
sampling thresholds. The default production thresholds require at least 5 bouts
and about 8 estimated stride cycles per fitted DMD operator. For exploratory
diagnostics, lower `--min-bouts-per-model` and `--min-strides-per-model`.
Use `--rank 40` or another positive rank cap for practical resolvent runs.
Without a rank cap, the default energy rule can keep the full delayed state
dimension, e.g. 8 state variables x 26 delays = 208, which makes the resolvent
SVD loop extremely slow.

Long runs periodically write:

```text
results/rodent_subject_speed_metrics.checkpoint.csv
results/rodent_dmd_results.checkpoint.jls
results/rodent_config.checkpoint.jls
```

For background or `tee` runs, use `--quiet-progressbar` so timestamped messages
stay readable:

```bash
julia --project=. scripts/00_fit_rodent_dmd_metrics.jl \
  --full \
  --state-space centroid_angle_sincos \
  --min-bouts-per-model 8 \
  --min-strides-per-model 8 \
  --n-delays 26 \
  --rank 40 \
  --n-angles 240 \
  --checkpoint-interval 25 \
  --quiet-progressbar \
  --force 2>&1 | tee rodent_full_run.log
```

State-space options:

- `raw_xy`: original paw coordinates, then optional bout normalization.
- `centroid_xy`: subtract the per-frame paw centroid before normalization.
- `body_aligned_xy`: subtract centroid, rotate to the fore-hind body axis, and scale by body length.
- `centroid_geometry`: paw radius and angle from the per-frame centroid.
- `centroid_angles`: unwrapped paw angles from the per-frame centroid.
- `centroid_angle_sincos`: paw angles from the centroid encoded as cosine/sine pairs; this is the preferred angle-only state for DMD because it avoids circular wrap jumps.
- `centroid_geometry_sincos`: alias for `centroid_angle_sincos`.

## Outputs

```text
results/rodent_subject_speed_metrics.csv
results/rodent_group_summary.csv
results/rodent_speed_bin_summary.csv
results/rodent_speed_and_henrici_trends.csv
results/rodent_io_structure_summary.csv
figures/Figure 1.svg
figures/Figure 2.svg
```

Figure 1 summarizes group-level trends across speed bins, non-normality,
transient growth, harmonics, and participation ratio. Figure 2 summarizes the
resolvent input-output structure as group heatmaps.
