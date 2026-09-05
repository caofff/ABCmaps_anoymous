# ABCmaps

`ABCmaps` is a lightweight R package containing three algorithms for approximate
Bayesian computation (ABC):

- `MAPS()` for multifidelity adaptive particle sampling.
- `PrefilterIS()` for prior-proposal prefiltered importance sampling (PF-IS).
- `ABC_ASMC()` for adaptive sequential Monte Carlo ABC.
- `estimate_pilot_fpr()` for estimating low-fidelity screening error on a
  shared pilot population.

The package also provides `maps_gaussian_proposal()`, the default proposal used
by `MAPS()` and `ABC_ASMC()`. The reproducible Toy model and comparison runner
are kept in `inst/examples/`, separate from the package algorithm code. The
package contains no other inference algorithms and has no external R package
dependencies.

## Download

Download the anonymous source archive from the [ABCmaps v0.10.0 GitHub release](https://github.com/caofff/ABCmaps/releases/download/v0.10.0/ABCmaps_0.10.0-anonymous.tar.gz).

## Package contents

```text
ABCmaps/
├── R/                 MAPS, PF-IS, ABC-ASMC, diagnostics, and utilities
├── man/               Function documentation
├── inst/examples/
│   └── toy-example.R   Complete 10-repetition Toy comparison
├── tests/smoke.R      Installation and algorithm smoke tests
├── DESCRIPTION        Package metadata
└── NAMESPACE          Exported functions
```

## Requirements

- A working R installation.
- Permission to write to an R library directory.
- No additional R packages are required.

Check that R is available from a terminal:

```sh
R --version
```

## Local installation

### Option 1: Install the local source archive from R

Replace the example path below with the absolute path to the downloaded or
built package archive:

```r
package_file <- normalizePath(
  "/path/to/MFSMC/ABCmaps_0.10.0.tar.gz",
  mustWork = TRUE
)

install.packages(
  package_file,
  repos = NULL,
  type = "source"
)
```

On Windows, either use forward slashes or escape backslashes:

```r
package_file <- "C:/path/to/MFSMC/ABCmaps_0.10.0.tar.gz"
```

### Option 2: Install the local source archive from a terminal

```sh
cd /path/to/MFSMC
R CMD INSTALL ABCmaps_0.10.0.tar.gz
```

### Option 3: Install directly from the unpacked package directory

```sh
R CMD INSTALL /path/to/MFSMC/ABCmaps
```

### Install into a personal library

If the default R library is not writable, create a personal library and pass it
to `install.packages()`:

```r
local_library <- path.expand("~/R/abcmaps-library")
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)

install.packages(
  "/path/to/MFSMC/ABCmaps_0.10.0.tar.gz",
  repos = NULL,
  type = "source",
  lib = local_library
)

.libPaths(c(local_library, .libPaths()))
library(ABCmaps)
```

### Verify the installation

```r
library(ABCmaps)

packageVersion("ABCmaps")
getNamespaceExports("ABCmaps")
```

The exported functions should be `ABC_ASMC`, `MAPS`, `PrefilterIS`,
`estimate_pilot_fpr`, and `maps_gaussian_proposal`.

## Model interface

All three algorithms use the same prior and discrepancy interfaces.

### Prior

`Prior` must be a list containing a sampling function. `MAPS()` also requires a
density function.

```r
Prior <- list(
  sample = function(n) {
    # Return a finite numeric n by d matrix.
  },
  density = function(theta) {
    # Return one finite, non-negative prior density for each row of theta.
  }
)
```

For a one-dimensional parameter, `sample(n)` may return either an `n` by `1`
matrix or a numeric vector of length `n`. For a multidimensional parameter, it
must return an `n` by `d` matrix. `density(theta)` must accept a particle matrix
and return a numeric vector of length `nrow(theta)`.

`PrefilterIS()` only calls `Prior$sample()`. A density function may still be
included so that the same `Prior` object can be used by both algorithms.

### Discrepancy function

`Dis_fun` must accept a particle matrix, a replication count, and a fidelity
label:

```r
Dis_fun <- function(theta, n_rep, type = "high") {
  # type is either "low" or "high".
  # Return a finite numeric nrow(theta) by n_rep matrix.
}
```

Each row corresponds to one parameter particle. Each column corresponds to one
independent simulator replicate. Smaller discrepancy values represent a closer
match to the observed data.

## `estimate_pilot_fpr()`

This function evaluates the low- and high-fidelity simulators on the same pilot
particles. It treats the high-fidelity ABC decision as the reference decision
and estimates the error introduced by low-fidelity screening.

The reported false-positive rate is the count-based quantity
`P(low pass | high reject)` in the pilot population. The function also reports
the cost-based upper bound `fpr_max` and a high-fidelity-posterior-weighted
false-negative rate,
`P(low reject | high pass)`, where the HF replicate acceptance fraction is used
as an ABC likelihood estimate.

### Usage

```r
estimate_pilot_fpr(
  Theta,
  Dis_fun,
  n_L,
  n_H,
  weights = NULL,
  kappa = 0.1,
  epsilon_target = NULL,
  epsilon = NULL,
  epsilon_L = NULL,
  a_L = 0.05,
  cost_mode = c("fixed", "measured"),
  c_L = 1,
  c_H = 10,
  verbose = FALSE
)
```

### Parameters

| Parameter | Description |
| --- | --- |
| `Theta` | Pilot particle vector for a one-dimensional parameter or an `N` by `d` numeric matrix. |
| `Dis_fun` | Low- and high-fidelity discrepancy function described above. |
| `n_L` | Positive integer number of low-fidelity replicates per pilot particle. |
| `n_H` | Positive integer number of high-fidelity replicates per pilot particle. |
| `weights` | Optional finite, non-negative pilot weights. Uniform weights are used when this is `NULL`. |
| `kappa` | Quantile in `(0, 1)` used to select `epsilon` when no fixed high-fidelity threshold is supplied. |
| `epsilon_target` | Optional finite lower bound for the selected high-fidelity threshold. |
| `epsilon` | Optional fixed high-fidelity threshold. If this is `NULL`, the weighted `kappa` quantile of the minimum HF discrepancies is used. |
| `epsilon_L` | Optional fixed low-fidelity threshold. If this is `NULL`, it is estimated from the HF-ABC posterior weights. |
| `a_L` | Target posterior-weighted false-negative rate in `[0, 1)` used for automatic `epsilon_L` selection. |
| `cost_mode` | Cost source used for `fpr_max`: either `"fixed"` or `"measured"`. The default is `"fixed"`. |
| `c_L` | Positive average cost of one LF simulation in fixed-cost mode. The default is `1`. |
| `c_H` | Positive average cost of one HF simulation in fixed-cost mode. The default is `10`. |
| `verbose` | If `TRUE`, reports the thresholds, FPR, and posterior-weighted false-negative rate. |

When `epsilon_L` is not supplied, the function calculates the weighted
`1 - a_L` quantile of the minimum LF discrepancies. The weights are proportional
to the pilot weight multiplied by the fraction of HF replicates accepted at
`epsilon`. This makes `a_L` directly relevant to the pilot threshold selection.

### Cost options

Use fixed costs when a benchmark, hardware specification, or known cost ratio is
available:

```r
fixed_cost_pilot <- estimate_pilot_fpr(
  Theta = pilot_theta,
  Dis_fun = Dis_fun,
  n_L = 2,
  n_H = 3,
  cost_mode = "fixed",
  c_L = 1,
  c_H = 10
)
```

Use measured costs to time the LF and HF simulator calls made by the pilot:

```r
measured_cost_pilot <- estimate_pilot_fpr(
  Theta = pilot_theta,
  Dis_fun = Dis_fun,
  n_L = 2,
  n_H = 3,
  cost_mode = "measured"
)
```

Measured mode does not run additional simulations. It records the wall-clock
time of the existing LF and HF pilot calls and calculates
`measured_c_L = low_elapsed_seconds / (N0*n_L)` and
`measured_c_H = high_elapsed_seconds / (N0*n_H)`. These measured values are
then used in the `fpr_max` formula. For very fast vectorized simulators, use a
larger pilot to obtain stable timings or use fixed-cost mode.

### Return value

The result contains four components:

| Component | Description |
| --- | --- |
| `summary` | One-row data frame containing thresholds, rates, confusion counts, and pilot settings. |
| `particles` | Particle-level table with discrepancies, match counts, weights, pass decisions, and classification flags. |
| `low_discrepancy` | Raw pilot LF discrepancy matrix. |
| `high_discrepancy` | Raw pilot HF discrepancy matrix. |

Important `summary` columns are:

| Column | Description |
| --- | --- |
| `epsilon`, `epsilon_L` | HF and LF thresholds used by the diagnostic. |
| `epsilon_selected`, `epsilon_L_selected` | Whether each threshold was selected automatically. |
| `p_H` | Pilot HF passing probability, calculated as the fraction of pilot particles with at least one accepted HF replicate. |
| `fpr` | Count-based fraction of HF-rejected particles that pass the LF screen. |
| `fpr_weighted` | Base-weighted FPR, included as a secondary diagnostic when non-uniform pilot weights are supplied. |
| `fpr_max` | Cost-based upper bound `1 - n_L*c_L / ((1-p_H)*n_H*c_H)`. |
| `fpr_below_max` | Whether the strict cost condition `fpr < fpr_max` holds. |
| `fnr_hf_posterior` | HF-ABC-posterior-weighted false-negative rate. |
| `fnr_unweighted` | Unweighted false-negative fraction among HF-passing particles. |
| `high_acceptance_rate`, `low_acceptance_rate` | Weighted pilot acceptance rates. |
| `high_posterior_ess` | Effective sample size of the normalized HF-ABC posterior weights. |
| `n_true_positive`, `n_false_positive` | Pilot classification counts for LF passes. |
| `n_true_negative`, `n_false_negative` | Pilot classification counts for LF rejections. |
| `cost_mode` | Cost mode used for the diagnostic. |
| `c_L`, `c_H` | Per-simulation costs actually used to calculate `fpr_max`. |
| `measured_c_L`, `measured_c_H` | Per-simulation wall-clock costs measured during the pilot, returned in both modes. |
| `low_elapsed_seconds`, `high_elapsed_seconds` | Total measured wall-clock times of the LF and HF pilot calls. |

If every pilot particle passes the HF screen, there are no HF-negative cases on
which to estimate an FPR or its cost bound, so `fpr`, `fpr_weighted`, and
`fpr_max` are returned as `NA`.

The cost comparison follows

```text
(1 - p_H) * (1 - fpr) * n_H * c_H > n_L * c_L
```

which is equivalent to `fpr < fpr_max`. A negative `fpr_max` is retained rather
than truncated; it indicates that the LF simulation cost is already too large
for the inequality to hold for any valid FPR.

### Pilot-to-algorithm workflow

```r
set.seed(100)
pilot_theta <- Prior$sample(1000)

pilot <- estimate_pilot_fpr(
  Theta = pilot_theta,
  Dis_fun = Dis_fun,
  n_L = 2,
  n_H = 3,
  kappa = 0.1,
  epsilon_target = 0.1,
  a_L = 0.05,
  cost_mode = "fixed",
  c_L = 1,
  c_H = 10,
  verbose = TRUE
)

pilot$summary[c(
  "epsilon", "epsilon_L", "p_H", "fpr", "fpr_max",
  "fpr_below_max", "fnr_hf_posterior"
)]

maps_fit <- MAPS(
  Prior = Prior,
  Dis_fun = Dis_fun,
  N = 1000,
  n_L = 2,
  n_H = 3,
  alpha = 0.7,
  alpha_L = 0.7,
  epsilon_T = 0.1,
  epsilon_L_T = pilot$summary$epsilon_L
)

pfis_fit <- PrefilterIS(
  Prior = Prior,
  Dis_fun = Dis_fun,
  N = 10000,
  n_L = 2,
  n_H = 3,
  epsilon = 0.1,
  epsilon_L = pilot$summary$epsilon_L
)
```

## `ABC_ASMC()`

`ABC_ASMC()` is the high-fidelity adaptive SMC baseline. It uses the same
adaptive quantile schedule, weighted covariance proposal, and ESS resampling
rule as MAPS, without an LF screening stage.

### Usage

```r
ABC_ASMC(
  Prior,
  Dis_fun,
  proposal_fun = maps_gaussian_proposal,
  N,
  n_H,
  alpha,
  epsilon_T,
  gamma = 0.5,
  type = "high",
  verbose = FALSE,
  max_iter = 100L
)
```

### Parameters

| Parameter | Description |
| --- | --- |
| `Prior` | Prior interface containing `sample(n)` and `density(theta)`. |
| `Dis_fun` | Simulator discrepancy function. |
| `proposal_fun` | Particle perturbation function. The default is `maps_gaussian_proposal`. |
| `N` | Positive integer number of particles. |
| `n_H` | Positive integer number of simulator replicates per evaluated particle. |
| `alpha` | Adaptive discrepancy-threshold quantile in `(0, 1)`. |
| `epsilon_T` | Finite terminal ABC threshold. |
| `gamma` | Resampling threshold as a fraction of `N`; it must lie in `[0, 1)`. |
| `type` | Fidelity label passed to `Dis_fun`; the default is `"high"`. |
| `verbose` | If `TRUE`, prints one progress message per iteration. |
| `max_iter` | Positive integer maximum number of iterations. |

### Return value

| Component | Description |
| --- | --- |
| `Theta` | Final `N` by `d` particle matrix. |
| `weights` | Final normalized particle weights. |
| `Eps` | Adaptive discrepancy threshold at each iteration. |
| `ESS` | Effective sample size after each iteration. |
| `n_high` | Cumulative number of particle rows evaluated by the simulator. Multiply by `n_H` for the replicate count. |
| `elapsed` | Cumulative elapsed time in seconds. |
| `smc_Data_ite` | Particle and weight history indexed by iteration. |

### Minimal call

```r
set.seed(1)

asmc_fit <- ABC_ASMC(
  Prior = Prior,
  Dis_fun = Dis_fun,
  N = 1000,
  n_H = 3,
  alpha = 0.7,
  epsilon_T = 0.1,
  gamma = 0.5,
  max_iter = 100
)

stopifnot(abs(sum(asmc_fit$weights) - 1) < 1e-12)
```

## `MAPS()`

### Usage

```r
MAPS(
  Prior,
  Dis_fun,
  proposal_fun = maps_gaussian_proposal,
  N,
  n_L,
  n_H,
  alpha,
  alpha_L,
  epsilon_T,
  epsilon_L_T = -Inf,
  gamma = 0.5,
  verbose = FALSE,
  max_iter = 100L
)
```

### Parameters

| Parameter | Description |
| --- | --- |
| `Prior` | Prior interface containing `sample(n)` and `density(theta)`. |
| `Dis_fun` | Low- and high-fidelity discrepancy function. |
| `proposal_fun` | Particle perturbation function. The default is `maps_gaussian_proposal`. |
| `N` | Positive integer number of particles maintained by MAPS. |
| `n_L` | Positive integer number of low-fidelity replicates per evaluated particle. |
| `n_H` | Positive integer number of high-fidelity replicates per evaluated particle. |
| `alpha` | High-fidelity threshold quantile in `(0, 1)`. |
| `alpha_L` | Low-fidelity threshold quantile in `(0, 1)`. |
| `epsilon_T` | Finite final high-fidelity ABC threshold. |
| `epsilon_L_T` | Lower bound for the adaptive low-fidelity threshold. Set it to `-Inf` to disable the bound. |
| `gamma` | Resampling threshold as a fraction of `N`; it must lie in `[0, 1)`. |
| `verbose` | If `TRUE`, prints one progress message per iteration. |
| `max_iter` | Positive integer maximum number of MAPS iterations. |

The custom proposal interface is:

```r
proposal_fun <- function(theta, Theta, weights) {
  # theta: active particles to perturb
  # Theta: complete particle population
  # weights: current normalized weights
  # Return a matrix with the same dimensions as theta.
}
```

### Return value

`MAPS()` returns a list with these components:

| Component | Description |
| --- | --- |
| `Theta` | Final `N` by `d` particle matrix. |
| `weights` | Final normalized particle weights. |
| `Eps` | High-fidelity threshold at each iteration. |
| `Eps_L` | Low-fidelity threshold at each iteration. |
| `ESS` | Effective sample size after each high-fidelity correction. Identical resampled particles are grouped when calculating this value. |
| `n_low` | Cumulative number of particle rows evaluated at low fidelity. Multiply by `n_L` for the corresponding replicate count. |
| `n_high` | Cumulative number of particle rows evaluated at high fidelity. Multiply by `n_H` for the corresponding replicate count. |
| `elapsed` | Cumulative elapsed time in seconds. |
| `MF_Data_ite` | Particle history containing parameter columns, weights, iteration, fidelity, and group labels. |
| `epsilon_L_T` | Low-fidelity threshold lower bound used in the run. |
| `epsilon_L_floor_fixed` | Whether a finite low-fidelity lower bound was used. |

### Minimal call

```r
set.seed(1)

maps_fit <- MAPS(
  Prior = Prior,
  Dis_fun = Dis_fun,
  N = 1000,
  n_L = 2,
  n_H = 3,
  alpha = 0.7,
  alpha_L = 0.7,
  epsilon_T = 0.1,
  epsilon_L_T = 0.2,
  gamma = 0.5,
  max_iter = 100
)

stopifnot(abs(sum(maps_fit$weights) - 1) < 1e-12)
weighted_parameter_mean <- colSums(maps_fit$Theta * maps_fit$weights)
```

## `PrefilterIS()`

### Usage

```r
PrefilterIS(
  Prior,
  Dis_fun,
  N,
  n_L,
  n_H,
  epsilon,
  epsilon_L,
  chunk_size = 5000L,
  verbose = FALSE
)
```

### Parameters

| Parameter | Description |
| --- | --- |
| `Prior` | Prior interface containing at least `sample(n)`. |
| `Dis_fun` | Low- and high-fidelity discrepancy function. |
| `N` | Positive integer number of candidates drawn from the prior. |
| `n_L` | Positive integer number of low-fidelity replicates per candidate. |
| `n_H` | Positive integer number of high-fidelity replicates per low-fidelity survivor. |
| `epsilon` | Finite high-fidelity ABC acceptance threshold. |
| `epsilon_L` | Finite low-fidelity screening threshold. |
| `chunk_size` | Positive integer number of prior candidates processed in one chunk. Smaller values reduce temporary memory use. |
| `verbose` | If `TRUE`, reports candidate and survivor counts for each chunk. |

PF-IS samples all candidates from the prior. A candidate is passed to the
high-fidelity simulator when at least one of its low-fidelity discrepancies is
at most `epsilon_L`. Its unnormalized importance weight is the number of
high-fidelity discrepancies at most `epsilon`.

### Return value

`PrefilterIS()` returns only particles with positive high-fidelity weight:

| Component | Description |
| --- | --- |
| `Theta` | Positive-weight particle matrix. |
| `weights` | Normalized importance weights. |
| `candidate_count` | Total number of prior candidates. |
| `lf_survivor_count` | Number of candidates passing the low-fidelity screen. |
| `lf_survival_rate` | Fraction of candidates passing the low-fidelity screen. |
| `positive_weight_count` | Number of returned particles. |
| `positive_weight_rate` | Fraction of all candidates receiving positive weight. |
| `low_simulations` | Total number of low-fidelity simulator replicates. |
| `high_simulations` | Total number of high-fidelity simulator replicates. |
| `n_L`, `n_H` | Replication counts supplied to the function. |
| `epsilon`, `epsilon_L` | Thresholds supplied to the function. |

### Minimal call

```r
set.seed(2)

pfis_fit <- PrefilterIS(
  Prior = Prior,
  Dis_fun = Dis_fun,
  N = 10000,
  n_L = 2,
  n_H = 3,
  epsilon = 0.1,
  epsilon_L = 0.2,
  chunk_size = 2000
)

stopifnot(abs(sum(pfis_fit$weights) - 1) < 1e-12)
weighted_parameter_mean <- colSums(pfis_fit$Theta * pfis_fit$weights)
```

## `maps_gaussian_proposal()`

### Usage

```r
maps_gaussian_proposal(theta, Theta, weights)
```

### Parameters and return value

| Parameter | Description |
| --- | --- |
| `theta` | Active particles to perturb, stored in an `m` by `d` matrix. |
| `Theta` | Complete particle population, stored in an `N` by `d` matrix. |
| `weights` | Finite, non-negative weights of length `N` with a positive sum. |

The function returns an `m` by `d` matrix. It uses a Gaussian random walk with
covariance equal to twice the weighted covariance of the complete population.
A small numerical floor is applied to degenerate covariance directions.

## Complete Toy example

The example uses `yobs = 1`, `N = 10000`, `n_H = n_L = 1`,
`epsilon_T = 0.05`, `alpha = alpha_L = 0.7`, and 10 independent repetitions.
The MAPS pilot uses `pilot_N = 5000`, `kappa = 0.05`, and `a_L = 0.01`.

Run the single example script:

```r
library(ABCmaps)

example_dir <- system.file("examples", package = "ABCmaps")
source(file.path(example_dir, "toy-example.R"))
```

The script calculates ISE against the exact HF-ABC density and writes the
following files beside `toy-example.R`:

- `toy-comparison.csv`
- `toy-summary.csv`
- `toy-comparison.pdf`
- `toy-comparison.png`

The complete English-commented example is available at
[`inst/examples/toy-example.R`](inst/examples/toy-example.R). It contains the
model settings, algorithm calls, ISE calculation, numerical summaries, and
plotting code in one file.

## Reproducibility

All algorithms use R's standard random-number generator. Call `set.seed()`
immediately before an algorithm call when reproducible output is required:

```r
set.seed(42)
fit <- PrefilterIS(
  Prior, Dis_fun,
  N = 5000, n_L = 2, n_H = 3,
  epsilon = 0.1, epsilon_L = 0.2
)
```

## Common errors

- **Unexpected prior sample size:** ensure `Prior$sample(n)` returns exactly
  `n` particle rows.
- **Unexpected discrepancy dimensions:** ensure `Dis_fun()` returns exactly
  `nrow(theta)` rows and `n_rep` columns.
- **Zero MAPS weights:** increase the relevant threshold, increase `N`, or
  inspect whether the low- and high-fidelity discrepancies overlap.
- **No positive PF-IS particles:** increase `N`, `epsilon_L`, or `epsilon`, then
  verify that the simulators can generate discrepancies below both thresholds.
- **MAPS reaches `max_iter`:** increase `max_iter` or use a less restrictive
  `epsilon_T`.
- **ABC-ASMC reaches `max_iter`:** increase `max_iter` or use a less restrictive
  `epsilon_T`.
- **Package not found after personal-library installation:** add the library
  path with `.libPaths()` before calling `library(ABCmaps)`.

## Build and check from source

Run these commands from the directory containing `ABCmaps/`:

```sh
R CMD build ABCmaps
R CMD check --no-manual ABCmaps_0.10.0.tar.gz
```

The smoke test can also be run through `R CMD check`; it verifies installation,
weight normalization, output dimensions, both pilot cost modes, and basic MAPS,
PF-IS, and ABC-ASMC execution.
