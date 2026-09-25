# Simulation results and figure guide

R scripts for **Optimal mixture designs for alloy compositions**. All generated plots, numerical results, checkpoints, and session records go in [`figures/`](figures/).

## Run

Tested with R 4.5.2. Install missing packages once:

```r
install.packages(c("spacefillr", "SPlit", "cluster", "support", "rkriging",
                   "Ternary", "magick", "PlotTools", "MOFAT",
                   "ggplot2", "tidyr", "tidyselect"))
```

From the repository root:

```sh
Rscript --vanilla simulation_results/run_all.R
Rscript --vanilla simulation_results/verify_results.R
```

This runs Figures 2–4, D4 and D7 experiments, their prediction/runtime plots, and the D7 coordinate-maxima analysis. The full run uses 30 replications and takes tens of minutes on this machine. Each script also works from `simulation_results/` (e.g. `Rscript --vanilla Figure_2.R`) or via an absolute path from another directory. Rerunning regenerates its outputs.

For individual plots, run `D4.R` before `plot_D4.R`, `D7.R` before `plot_D7.R`, and both experiments before `plot_runtime.R`. `D7_max_level.R` runs independently. Partial experiment checkpoints are saved after each design size; plotting requires a completed run. Restart an interrupted experiment by rerunning its script.

## Paper → script → output

All filenames below are relative to `figures/`. Prediction plots contain **full**, **boundary**, and **interior** test-set panels, in that order.

| Paper figure | Script(s) | What is plotted | Figure output |
| --- | --- | --- | --- |
| 2(a–c) | [`Figure_2.R`](Figure_2.R) | Random, Sobol QMC, and K-medoids designs on the three-component simplex | `rand_X.png`, `sobol_X.png`, `k_med_X.png` |
| 3(a–b) | [`Figure_3.R`](Figure_3.R) | SP and SSP MSE surfaces at isotropic Gaussian lengthscale `theta = 0.45`, with a common color scale | `MSE_ternary_SP.png`, `MSE_ternary_SSP_kappa_1.png` |
| 4 | [`Figure_4.R`](Figure_4.R) | SP (red circles) → SSP (green diamonds), joined by arrows | `all_X.png` |
| 5 | [`D4.R`](D4.R), [`D7.R`](D7.R) → [`plot_runtime.R`](plot_runtime.R) | Median sampling-only elapsed time over 30 replications | `D4_D7_runtime.jpeg` |
| 6 / S1 | `D4.R` → [`plot_D4.R`](plot_D4.R) | D4 bulk-modulus NMAE (%) / RMSE (GPa) | `combined_Bulk_Modulus_D4_{NMAE,RMSE}.jpeg` |
| 7 / S2 | `D4.R` → `plot_D4.R` | D4 formation-energy NMAE (%) / RMSE (eV) | `combined_Formation_Energy_D4_{NMAE,RMSE}.jpeg` |
| 8 / S3 | `D7.R` → [`plot_D7.R`](plot_D7.R) | D7 bulk-modulus NMAE (%) / RMSE (GPa) | `combined_Bulk_Modulus_D7_{NMAE,RMSE}.jpeg` |
| 9 | [`D7_max_level.R`](D7_max_level.R) | Coordinate maxima at `n = 70`; each box contains 30 designs, with jittered points | `D7_max_level.jpeg` |


## Data and saved results

- **D4:** `pca_pspall4.csv`, 6,545 Al–Nb–Ti–Zr compositions; formation energy and bulk modulus; `n = 20, 30, 40, 50, 60`.
- **D7:** `RHEA7unique.csv`, 12,012 Mo–Nb–Ta–Ti–V–W–Zr compositions; bulk modulus; `n = 50, 60, 70, 80, 90`.
- Atom counts are divided by 128. Each experiment uses Random, QMC, K-medoids, SP, and SSP (`kappa = 1`), seeds 1–30, and Gaussian `rkriging` fits. Test sets exclude training indices; boundary points have at least one fraction equal to 0 or 1. NMAE is divided by the mean absolute response over the **whole candidate set** and multiplied by 100.
- `D4.Rdata` / `D7.Rdata` store all metrics, sampling runtimes, selected row indices, and run metadata; `D4_metrics.csv` / `D7_metrics.csv` provide 750 rows each. `D4_session.txt` / `D7_session.txt` record R/package versions.
- `Figure_2_results.rds`, `Figure_3_results.rds`, and `Figure_4_results.rds` save their designs and plot inputs. All three use 10,000 Sobol candidates; Figure 3 evaluates 40,000 Sobol points plus transformed corners. `D7_max_level.Rdata` and `.csv` save coordinate minima/maxima and selected indices (indices in Rdata); `D4_D7_runtime.csv` saves median runtimes.

Shared helpers: [`sampling_strategy.R`](sampling_strategy.R) implements the five designs and coordinate-wise SSP scaling/renormalization (Eqs. 14–17); [`MMSE_calculation.R`](MMSE_calculation.R) implements Gaussian MSE and its maximum; [`experiment_helpers.R`](experiment_helpers.R) runs GP evaluation and saves results; [`common.R`](common.R) handles paths, dependencies, and shared plotting settings. [`run_all.R`](run_all.R) executes the workflow; [`verify_results.R`](verify_results.R) checks all 1,500 prediction designs, 150 coordinate-maxima designs, saved metrics, and 14 image files.
