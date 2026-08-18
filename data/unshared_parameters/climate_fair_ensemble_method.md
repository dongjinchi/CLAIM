# FaIR constrained climate ensemble (B2)

This note documents `climate_fair_ensemble.csv`. The file is a parameter ensemble only; no SCC values were computed.

## Prior

The prior is built from the CMIP6-derived Leach et al. parameter tables supplied in `C:\Users\dongj\AppData\Local\Temp\claude\F--OneDrive-03researches-05SCC\40c8268c-0d5b-4b28-8a11-0a9dfae3ca35\scratchpad`. CO2 forcing, heat-response amplitudes, and heat-response timescales are drawn jointly with a Liu-West (2001) smoothed bootstrap over the 28 CMIP6 rows: `z' = z_bar + sqrt(1 - h^2) (z_i - z_bar) + h L eps`, with `L L' = Sigma_hat` and bandwidth `h = 0.05`. The shrinkage term makes the kernel variance-preserving, so the smoothed prior keeps the empirical mean and covariance exactly. Sampling the forcing and amplitude coefficients jointly is what preserves the CMIP6 `F2x`-`sum(q)` anticorrelation, which is part of the ECS distribution.

`q1`, `q2`, `q3` are carried in log space and `f1`, `f2`, `f3` in raw space (`f1` is negative for one member). Two earlier parameterisations were rejected on measurement, both because of the five near-two-box members whose `d3` is 1e4-1e9 years and whose `q3` is 2.4e-06 to 3.4e-02:

* Fitting one GLOBAL Gaussian to `log(q3)` collapses the slow box. That distribution's median is the geometric mean of `q3`, which those five members drag from an empirical median of 0.371 down to 0.049, shrinking `sum(q)` to 0.785 of its CMIP6 value and the prior ECS median to `3.071 * 0.785 = 2.41` -- reproduced exactly as 2.418 in the previous revision. TCR was left almost untouched (+0.013) because the slow box carries weight `v3 ~ 69.66/(2 d3) ~ 0.11` in the transient response, and that asymmetry is what identified the fault.
* Perturbing `q` in RAW space destroys the same five members from the other side: the cross-member sd of `q3` is ~0.4, so even a narrow kernel applies a +-0.02 absolute perturbation to values of order 1e-6 and pushes about half of them negative. The positivity filter in `draw_prior` then truncates only the low side, and the measured prior ECS median rose monotonically with bandwidth (rejection 11-14%).

A smoothed bootstrap in log space avoids both: the kernel is local and multiplicative, so degenerate members stay near their own tiny `q3` while the rest stay near 0.4, and positivity is automatic (measured rejection 0.0%). `d3` is winsorized at 10,000 years inside the log-timescale kernel only, so genuine near-two-box behaviour can still appear without unbounded numerical timescales dominating the covariance.

Carbon cycle (`r0`, `rU`, `rT`, `rA`) remains a separate Gaussian block over the 11 C4MIP members, sampled independently of the thermal block: only six CMIP6 members pair all blocks, too few to estimate cross-block covariance. This independence is an assumption, not a finding.

**Scope limitation.** This is a *smoothed CMIP6 empirical prior*, not an independent perturbed-parameter ensemble of the Leach kind. The bandwidth cannot be raised much above 0.05 without inflating the ECS 5-95% width past the empirical spread (measured width ratios 1.09, 1.34, 1.47, 1.48, 1.52, 1.48 at h = 0.05, 0.15, 0.25, 0.35, 0.45, 0.60), because the log-`q3` direction carries the degenerate members' enormous variance. The prior therefore cannot fill the space between CMIP6 models. What it adds over using the 28 members directly is observational reweighting, free recombination of carbon-cycle and thermal blocks rather than only the 6 that pair, and an ensemble of arbitrary size.

## Observational constraint

Likelihoods were computed by `fair_likelihood.jl` using the repository Julia FaIR chain. The Julia script builds one `get_model(:national; climate=:fair, scenario="ssp245", end_year=2024, bau_ref_year=2024)` model, checks that explicitly overriding default thermal parameters reproduces the default historical temperature path, and then loops through candidate members with in-place `update_param!` followed by `run(m)`.

The constraint compares full-forcing FaIR global temperature to Berkeley Earth `climate_T_global.csv`: mean 2010-2019 warming minus mean 1950-1979 warming, and the 1995-2024 linear trend in degC per decade. The level term is now an in-model-period increment on both sides, so the 1950 initial offset cancels rather than contaminating the likelihood. This is a pragmatic full-forcing constraint, not Leach et al.'s full GWI multi-product anthropogenic-warming pipeline.

The error model is `Sigma = diag(sigma_k^2) + tau^2 * 1 1'`, i.e. a shared model-discrepancy offset `delta ~ N(0, tau^2)` common to all four increments, marginalised analytically via Sherman-Morrison. Here tau = 0.080 degC.

The diagonal terms are 0.104, 0.104, 0.104, 0.109 degC, each the quadrature sum of an observational term (0.050) and an internal-variability term (0.044, 0.044, 0.044, 0.055). Only genuinely period-independent errors are kept there. The internal-variability term is derived, not asserted: residuals of the observed series about a 31-year centred moving average have sd 0.111 degC with lag-1 autocorrelation 0.396, an effective decorrelation time of 2.31 years, hence 0.031 degC on a 30-year mean and 0.044 degC on the difference of two such means.

Model structural error is what tau carries. It used to sit on the diagonal at 0.08 degC per increment, treated as independent between periods -- which is wrong, because a model that is structurally too slow in mid-century is systematically so rather than randomly so. tau is that same asserted 0.08 degC relocated, not enlarged. The correction matters because ECS and TCR variation moves all four increments together (member-wise correlation +0.99 among the first three), so the direction the metrics can see is exactly the common direction; treating four correlated-error metrics as independent shrank sigma on that direction by sqrt(4 / 1.52) = 1.62x. Pass `--shared-bias-tau 0` to recover the previous over-confident likelihood.

The scored prior is written to `F:\OneDrive\03researches\05SCC\03Claim\data\unshared_parameters\climate_fair_ensemble_scored_prior.csv` and contains all candidates, parameter values, per-increment diagnostics, and unnormalized log likelihoods. The posterior CSV is a resampled ensemble, so `weight` is the equal consumption weight for each row; `likelihood_weight` keeps the normalized pre-resampling likelihood weight of the source candidate.

## Diagnostics

- prior candidates scored: 2000
- posterior members written: 2000
- ESS: 484.9 (24.2% of scored prior)
- default-parameter residuals by increment (degC): -0.004, -0.110, -0.052, +0.129
- default-parameter residuals by increment (sigma): -0.04, -1.06, -0.50, +1.18
- observed increments (degC): 0.098, 0.226, 0.444, 0.476
- observed total warming 2011-2024 vs 1850-1900: 1.243 degC
- Earth energy inventory 1971-2018 (x1e22 J): observed 43.5, default-parameter residual +13.5 (+2.03 sigma)
- Earth energy inventory 1971-2018 (x1e22 J): prior median 40.3 -> posterior median 43.5
- Spearman correlation between likelihood weight and modelled heat uptake: 0.19
- prior TCR median / 5-95%: 1.887 / 1.508-2.701
- prior ECS median / 5-95%: 3.006 / 1.916-5.966
- posterior TCR median / 5-95%: 1.715 / 1.500-2.137
- posterior ECS median / 5-95%: 2.534 / 1.951-3.637
- corr(F2x, sum(q)): CMIP6 -0.689; prior -0.676

### ECS distribution

| sample | median | 5% | 95% |
|---|---:|---:|---:|
| CMIP6 empirical | 3.071 | 1.936 | 5.713 |
| prior | 3.006 | 1.916 | 5.966 |
| posterior | 2.534 | 1.951 | 3.637 |
| prior shift | -0.065 |  |  |
| likelihood shift | -0.472 |  |  |

### TCR distribution

| sample | median | 5% | 95% |
|---|---:|---:|---:|
| CMIP6 empirical | 1.890 | 1.537 | 2.613 |
| prior | 1.887 | 1.508 | 2.701 |
| posterior | 1.715 | 1.500 | 2.137 |
| prior shift | -0.003 |  |  |
| likelihood shift | -0.172 |  |  |

### Prior component medians

| component | CMIP6 median | prior median | ratio |
|---|---:|---:|---:|
| F2x | 3.5397 | 3.5042 | 0.990 |
| sum_q | 0.87706 | 0.87689 | 1.000 |
| q1 | 0.19869 | 0.19918 | 1.003 |
| q2 | 0.38592 | 0.37747 | 0.978 |
| q3 | 0.37112 | 0.35871 | 0.967 |
| d1 | 0.97469 | 0.97878 | 1.004 |
| d2 | 8.429 | 8.4452 | 1.002 |
| d3 | 318.03 | 322.43 | 1.014 |

## Comparison with AR6, and how to use this ensemble

AR6 assessed TCR at 1.8 degC (likely 1.4-2.2) and ECS at 3.0 degC (likely 2.5-4.0).

**TCR behaves as intended.** The posterior median is close to the AR6 best estimate and essentially the whole posterior falls inside the AR6 likely range. This is expected: the trend metric constrains the transient response directly.

**ECS is shifted low and stays wide.** Seventy-five years of record cannot pin the equilibrium response; ECS is constrained only indirectly, through its correlation with TCR inside the CMIP6 prior. It is therefore pulled down without being sharpened.

**Part of that downward shift is an artefact, not evidence about ECS.** The scored prior shows the likelihood penalising warming monotonically: the Spearman correlation between likelihood weight and TCR is -0.37, and the weight share by prior-ECS quartile runs 44% / 37% / 17% / 2% from coldest to hottest. Some of that is the constraint doing its job on the known CMIP6 hot-model problem. Whether any of it is instead a systematic bias no choice of sensitivity can fix is answered by the total-warming diagnostic: 50% of prior members overshoot the observed total warming (median overshoot -0.002 degC), and the coldest ECS quartile sits at 1.073 degC against an observed 1.243 degC. If the coldest quartile also overshoots, the shift is an artefact rather than evidence; if it straddles the observation, the constraint is discriminating properly.

Note the metric history here, because it mattered enormously. An earlier revision scored a 1995-2024 trend and diagnosed a 1.5 sigma common-mode warm bias, which suppressed the posterior ECS by an amount that was 98% artefact. That bias was an artefact of the window itself: on a common 1850-1900 baseline the model is slightly *cold* over 1850-2024 and 1900-2024, and only looks hot in short recent windows whose first half contains the 1998-2013 hiatus.

**The posterior is over-confident for TCR.** Its 5-95 percent width is materially narrower than the AR6 66 percent likely range. Two metrics treated as independent Gaussians, one scenario, one observational product, and unsampled forcing uncertainty all push in that direction.

**Therefore, when propagating to SCC:** report the interval from this ensemble as a *lower bound* on climate uncertainty, not as the climate uncertainty. The 28-member CMIP6 spread (unconstrained, hot-biased) is the natural upper reference. Quoting this posterior alone as the climate uncertainty range would understate it. Fixing this properly requires sampling aerosol and external forcing uncertainty, which the repository currently has no data for; inflating sigma to widen the posterior would be curve-fitting to AR6 and is deliberately not done.
