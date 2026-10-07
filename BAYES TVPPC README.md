# Did the Austrian Phillips Curve Flatten After the Euro?

A Bayesian time-varying-parameter (TVP) analysis of the relationship between inflation, unemployment and the output gap in Austria, 1997 to 2008.

Course project for *Bayesian Inference* (Master's level), Paris Lodron University of Salzburg, Summer Semester 2026. All code is in R.

## Summary

I estimate a New Keynesian-style Phillips curve for Austria with coefficients that drift over time, and ask whether inflation became less responsive to unemployment after euro adoption in 1999 Q1.

- **Unemployment slope:** there is moderate evidence of flattening. The posterior probability that the unemployment coefficient moved toward zero between 1999 Q1 and 2008 Q4 is about 0.72 to 0.75, but the posterior distributions overlap substantially.
- **Output gap coefficient:** no evidence of flattening. It moved away from zero over the same period, and the posterior probability of flattening is only 0.33 to 0.39.
- **Timing:** coefficient paths change gradually rather than breaking at 1999, so euro adoption alone is a weak explanation for the change.
- **Stochastic volatility:** in an extension adapted from published code (see *Credits* below), the apparent shift disappears and estimated inflation volatility falls over the sample.

These are associations. The design has no controls for other macroeconomic developments of the period.

## Data

| Variable | Source | Construction |
|---|---|---|
| Inflation | Eurostat HICP, all items, index 2015 = 100 (`prc_hicp_midx`) | Year-on-year rate from the monthly index, averaged to quarters |
| Unemployment | Eurostat monthly unemployment rate, seasonally adjusted, total (`une_rt_m`) | Quarterly average of monthly values |
| Output gap | OECD Economic Outlook No. 118, Austria | Quarterly series, read from `data/AT_output_gap.xlsx` |
| Lagged inflation | Constructed | One-quarter lag of inflation, used as a control |

The sample runs from 1997 Q2 to 2008 Q4 (47 quarterly observations after the lag). All series are used in levels. The data files are not included in this repository. The script downloads the Eurostat series directly. For the output gap, download the Austrian series from the OECD Economic Outlook and save it as `data/AT_output_gap.xlsx`, with the columns `Time period` (e.g. `1997-Q2`) and the gap series, as the script expects.

## Model

Observation equation, with time-varying coefficients:

```
pi_t = b0_t + b1_t * pi_{t-1} + b2_t * u_t + b3_t * gap_t + e_t,     e_t ~ N(0, sigma^2)
```

State equation, a random walk:

```
beta_t = beta_{t-1} + eta_t,     eta_t ~ N(0, Q)
```

Three specifications are estimated:

1. **Full model:** lagged inflation, unemployment and output gap
2. **Baseline:** lagged inflation and unemployment (used because of multicollinearity in the full model)
3. **Robustness:** lagged inflation and output gap

### Priors

Weakly informative and conjugate, so every conditional posterior has a closed form:

- Initial state: `beta_0 ~ N(0, 10 * I)`
- Innovation covariance: `Q ~ Inverse-Wishart(0.01 * I, K + 1)`
- Observation variance: `sigma^2 ~ Inverse-Gamma(0.01, 0.01)`

### Estimation

A Gibbs sampler cycles through three conditional draws:

1. The full coefficient path, drawn jointly in one block with the Carter and Kohn (1994) forward-filtering, backward-sampling algorithm (Kalman filter forward, simulation smoother backward)
2. `sigma^2` from its inverse-gamma posterior
3. `Q` from its inverse-Wishart posterior

I run 7,000 iterations and discard the first 2,000 as burn-in, which leaves 5,000 draws for inference.

### Identifying a change

I compare the posterior distribution of each coefficient at 1999 Q1 (euro adoption, a hypothesised break date) with its posterior at the final quarter, 2008 Q4. The statistic is `P(|beta_2008Q4| < |beta_1999Q1| | data)`, the posterior probability that the coefficient is closer to zero at the end of the sample. Kernel density plots of the two posteriors show the shift visually.

## Results

| Specification | Coefficient | P(flattened) |
|---|---|---|
| 1 (full) | Unemployment | 0.725 |
| 1 (full) | Output gap | 0.334 |
| 2 (baseline) | Unemployment | 0.751 |
| 3 | Output gap | 0.392 |

- In every specification the unemployment and output gap coefficients drift gradually, so flattening is a slow-moving process and not a structural break.
- The unemployment posteriors overlap substantially, so there is real uncertainty about the size of the shift.
- The output gap posteriors are clearly separated, but the coefficient moves away from zero, which contradicts the simple flattening hypothesis.
- The full specification shows the largest swings, probably because of multicollinearity.

## Extension: stochastic volatility

The baseline model assumes a constant error variance, so it may attribute changes in inflation volatility to coefficient drift. To check this, I re-estimated baseline specification 2 with stochastic volatility using an adapted version of the TVP-SVD sampler described below.

With stochastic volatility, the coefficient paths are nearly flat, the 1999 Q1 and 2008 Q4 posteriors almost coincide, and the earlier evidence of flattening disappears. The estimated volatility path declines over the sample, from about 2.6 early on to about 1.4 by 2008 Q4.

## Limitations

- **Association, not causation.** There is no identification strategy for euro adoption, and 1999 Q1 is a hypothesised break date.
- **Small sample.** There are 47 quarterly observations. Output gap data are only easily available quarterly.
- **Multicollinearity.** Unemployment and the output gap are correlated, which is why the baseline specification drops the output gap.
- **The extension changes more than volatility.** The TVP-SVD sampler differs from my baseline in more than stochastic volatility. By default it uses white-noise state evolution and a sparse-mixture shrinkage prior, which can flatten coefficient paths by itself. The disappearance of the shift therefore cannot be attributed to stochastic volatility alone. A cleaner test would add stochastic volatility to the baseline random-walk model.
- **Data revisions.** The output gap is an estimated, revised series.

## Repository structure

```
austrian-phillips-curve-tvp/
  README.md
  Application_TVPPC.R        # data download, baseline TVP Gibbs sampler, three specifications, posterior comparisons
  TVP_SV_extension.R         # stochastic volatility extension applied to the Austrian data
  external/                  # TVP-SVD files from the replication archive (see Credits)
  figures/                   # coefficient paths and kernel density plots
```

### Requirements

R with the packages `eurostat`, `dplyr`, `zoo`, `readxl`, `MASS`, `ggplot2`, `reshape2`, `coda`, `GIGrvg`, `Matrix`, `mvtnorm`, `shrinkTVP`, `stochvol` and `bayesm`.

### Running

1. Save the output gap data to `data/AT_output_gap.xlsx`.
2. Set the working directory to the repository root. The script currently uses a hard-coded `setwd()` path, so edit or remove it first.
3. Run `Application_TVPPC.R` for the baseline results and the posterior probabilities.
4. For the extension, download the TVP-SVD files into `external/` (see Credits) and run `[SV script name].R`.

## Credits and authorship

**Written by Ben Bell with assistance of Claude AI:** the data pipeline, the baseline TVP model, the Carter-Kohn Gibbs sampler, the three specifications, the posterior comparisons and all plots in `Application_TVPPC.R` and 'TVP_SV_extension.R '.

**Adapted from other work:** the stochastic volatility extension uses an adapted version of the TVP-SVD estimator from the replication archive of Hauzenberger, Huber, Koop and Onorante (2021), available at https://github.com/fhuber7/replication-archive/tree/main/TVPSVD_replication. The files I used (`main_replication.R` and `tvpsvd_estim.R`) state in their header that they are a repackaged, modified extract of the authors' original replication package. I did not write that estimator. I applied it to the Austrian data. Please check the archive for its licence and use the authors' original code if you need an official version.

**AI assistance:** [This project partially integrated Claude AI into the writing of Code]

## References

- Carter, C.K. and Kohn, R. (1994). On Gibbs sampling for state space models. *Biometrika*, 81(3), 541-553.
- Hauzenberger, N., Huber, F., Koop, G. and Onorante, L. (2021). Fast and flexible Bayesian inference in time-varying parameter regression models. *Journal of Business & Economic Statistics*, 39(4), 1096-1111.
- Huber, F. and Schreiner, J. (2025). *Are Phillips Curves in CESEE Still Alive and Well Behaved?* Oesterreichische Nationalbank.
- Koop, G. and Korobilis, D. (2009). Bayesian multivariate time series methods for empirical macroeconomics. *Foundations and Trends in Econometrics*, 3(4), 267-358.
- Musso, A., Stracca, L. and van Dijk, D. (2007). *Instability and Nonlinearity in the Euro Area Phillips Curve.* ECB Working Paper No. 811.
- Rumler, F. (2006). The New Keynesian Phillips Curve for Austria: an extension for the open economy. *Monetary Policy & the Economy*, Q4/06, 69-87.
- Götz, T.B. and Hauzenberger, K. (2021). Large mixed-frequency VARs with a parsimonious time-varying parameter structure. *The Econometrics Journal*, 24(3), 442-461.
- Data: Eurostat (HICP and unemployment); OECD (2025), *OECD Economic Outlook*, No. 118.
