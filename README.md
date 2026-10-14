# Monte Carlo Option Pricing Engine

A parallel Monte Carlo engine in C++17 / OpenMP for European, Asian and American
options under geometric Brownian motion. Every estimator reports a standard
error and 95% confidence interval and is validated against a closed form, a
lattice, or published results. The output is bit-identical for any number of
threads.

Headline results (Intel i5-6200U, 2 cores / 4 threads, MSVC 19.44 `/O2`):

- **Variance reduction at equal cost:** control variates cut variance **6.9x**
  for an at-the-money European call and **~1,300x** for an arithmetic Asian call
  (geometric-average control, correlation 0.9996).
- **Accuracy:** 4M-path prices and Greeks land within 1.2 standard errors of
  Black-Scholes. Longstaff-Schwartz gives **4.473** on the standard American-put
  test case (paper: 4.472).
- **Throughput:** 65.6M paths/s on 4 threads (2.8x over 1 thread), with results
  bit-identical across thread counts.

## What it prices

| Product | Method | Validated against |
|---|---|---|
| European call/put | Plain MC, antithetic variates, control variate (discounted S_T) | Black-Scholes |
| Greeks: delta, vega, gamma | Pathwise (delta, vega); likelihood-ratio/pathwise hybrid (gamma) | Black-Scholes Greeks |
| Arithmetic-average Asian | Plain MC, geometric-average control variate | Kemna-Vorst closed form (control), plain MC |
| American / Bermudan put | Longstaff-Schwartz least squares, two-pass (out-of-sample) | CRR binomial tree, Longstaff & Schwartz (2001) Table 1 |
| Implied volatility | Safeguarded Newton-Raphson (Newton steps, bisection fallback) | Round trip |

## Results

All numbers come from `bench` (see [results/BENCHMARKS.md](results/BENCHMARKS.md); CSVs in
[results/](results/)). Base case: S = K = 100, r = 5%, sigma = 20%, T = 1.

### European call: estimators at equal cost (4M paths, 4 threads)

Black-Scholes = 10.450584.

| Method | Price | Std. error | Error / SE | Time (ms) | Variance reduction | Efficiency gain |
|---|---|---|---|---|---|---|
| Plain | 10.45139 | 0.00736 | +0.11 | 60.7 | 1.00x | 1.00x |
| Antithetic | 10.45644 | 0.00520 | +1.13 | 42.4 | 2.00x | 2.87x |
| Control variate | 10.45657 | 0.00281 | +2.13 | 68.8 | 6.88x | 6.07x |

*Variance reduction* compares SE² x paths, so an antithetic pair is charged for
both of its paths. *Efficiency gain* also charges for run time (variance x time),
which is the figure that matters in practice.

### Thread scaling (plain estimator, 16M paths)

| Threads | Time (ms) | Speedup | Paths / s | Price | Bit-identical to 1 thread |
|---|---|---|---|---|---|
| 1 | 687.2 | 1.00x | 23.3M | 10.455089 | yes |
| 2 | 392.5 | 1.75x | 40.8M | 10.455089 | yes |
| 4 | 243.9 | 2.82x | 65.6M | 10.455089 | yes |

The test machine has 2 physical cores, so the step from 2 to 4 threads only
gains from hyper-threading.

### Greeks from a single simulation (4M paths)

| Greek | Estimator | Monte Carlo | Std. error | Black-Scholes | Error / SE |
|---|---|---|---|---|---|
| Delta | pathwise | 0.636499 | 0.000288 | 0.636831 | -1.15 |
| Gamma | likelihood ratio / pathwise | 0.018771 | 0.000019 | 0.018762 | +0.48 |
| Vega | pathwise | 37.542286 | 0.037848 | 37.524035 | +0.48 |

### Arithmetic Asian call: geometric control variate (1M paths)

| Monitoring dates | Plain price | Plain SE | CV price | CV SE | Correlation | Variance reduction |
|---|---|---|---|---|---|---|
| 12 | 6.14960 | 0.00851 | 6.156109 | 0.000237 | 0.99961 | 1292x |
| 52 | 5.85172 | 0.00811 | 5.854182 | 0.000224 | 0.99962 | 1304x |
| 252 | 5.79274 | 0.00802 | 5.782067 | 0.000220 | 0.99962 | 1322x |

The control variate costs almost nothing extra per path, so a 1,300x variance
reduction means the same accuracy with roughly 1,300x fewer paths.

### American put: Longstaff-Schwartz vs CRR binomial

K = 40, r = 6%, T = 1, 50 exercise dates. The exercise rule is fitted on 100k
paths and the option priced on 1M independent paths.

| S | sigma | European (BS) | CRR American | LSM out-of-sample | SE | LSM in-sample | LSM - CRR |
|---|---|---|---|---|---|---|---|
| 36 | 0.2 | 3.8443 | 4.4867 | 4.4730 | 0.0029 | 4.4793 | -0.0137 |
| 40 | 0.2 | 2.0664 | 2.3195 | 2.3114 | 0.0027 | 2.3053 | -0.0081 |
| 44 | 0.2 | 1.0169 | 1.1130 | 1.1069 | 0.0021 | 1.1095 | -0.0061 |
| 36 | 0.4 | 6.7114 | 7.1092 | 7.1049 | 0.0061 | 7.0994 | -0.0043 |
| 40 | 0.4 | 5.0596 | 5.3181 | 5.3102 | 0.0057 | 5.3082 | -0.0080 |
| 44 | 0.4 | 3.7828 | 3.9531 | 3.9459 | 0.0052 | 3.9369 | -0.0072 |

The CRR tree allows exercise at any time, while LSM allows it on 50 dates and
its out-of-sample value is a lower bound. Small negative gaps are therefore
expected, and that is what the table shows.

## Method notes

**Exact simulation.** Under GBM, S_T = S_0 exp((r - sigma²/2)T + sigma sqrt(T) Z), so
European paths need one normal draw and have no discretisation error. The only
error left is statistical, which the standard error measures.

**Reproducible parallel random numbers.** Version 1 of this project shared a
single `std::mt19937` across OpenMP threads, which is a data race. Now paths are
split into fixed blocks of 16,384. Block *b* owns an `xoshiro256**` generator
seeded through SplitMix64 from (seed, *b*), and per-block statistics are merged
in block order. The result depends only on the seed, never on the thread count
or schedule, and the tests check this bit for bit. Normal variates use the
Marsaglia polar method instead of `std::normal_distribution`, whose algorithm
varies between standard libraries.

**Streaming statistics.** Means, variances and covariances are accumulated with
Welford's update and merged with Chan et al.'s pairwise formula. This avoids the
cancellation of naive sum-of-squares and makes per-block merging exact.

**Control variates.** With payoff Y and control X of known mean mu_X, the estimator is
Y_bar - beta (X_bar - mu_X), where beta = Cov(X,Y)/Var(X). Its variance is
Var(Y)(1 - rho²)/N. Estimating beta from the same paths adds an O(1/N) bias,
which is negligible at these path counts. The controls are:
- European: discounted S_T, with E[e^{-rT} S_T] = S_0.
- Asian: the discretely monitored geometric-average option. ln G is normal with
  mean ln S_0 + (r - sigma²/2) T (m+1)/(2m) and variance sigma² T (m+1)(2m+1)/(6m²),
  which gives a Black-Scholes-type closed form. With m = 1 it reduces to
  Black-Scholes exactly, which the tests check.

**Greeks.**
- Delta and vega use pathwise derivatives: d/dS_0 and d/dsigma of the discounted
  payoff along each path.
- Gamma cannot be computed pathwise, because the payoff's second derivative is a
  Dirac delta. Instead, the pathwise delta is differentiated once more with the
  likelihood-ratio method:
  Gamma = E[e^{-rT} 1{S_T > K} (S_T / S_0²)(Z / (sigma sqrt(T)) - 1)].
  This has far lower variance than the pure likelihood-ratio estimator.

**Longstaff-Schwartz.** Working backwards over the exercise dates, the code
regresses the realised discounted cash flows of in-the-money paths on
{1, x, x², x³} with x = S/K, and exercises where intrinsic value exceeds the
fitted continuation value. The normal equations are accumulated in parallel and
solved by Gaussian elimination with partial pivoting. The fitted rule is then
applied to *independent* paths. A fixed rule can only be sub-optimal, so the
reported price is a low-biased lower bound with no look-ahead bias. The
in-sample value is shown for comparison.

## Build and run

Requires a C++17 compiler with OpenMP (GCC, Clang with libomp, or MSVC).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure   # 54 checks
./build/bench results                                     # regenerates results/
./build/price --style american --type put --spot 36 --strike 40 --rate 0.06 --dates 50
```

Without CMake:

```bash
g++ -std=c++17 -O2 -fopenmp -Iinclude src/price.cpp -o price          # GCC / Clang
cl /std:c++17 /O2 /openmp:llvm /EHsc /Iinclude src\price.cpp          # MSVC
```

`price` options: `--style european|asian|american --type call|put --spot --strike --rate --vol --expiry --paths --dates --threads --seed`.

## Layout

```
include/mcpricer/
  rng.hpp            xoshiro256** + SplitMix64 seeding, Marsaglia polar normals
  stats.hpp          Welford / Chan streaming mean, variance, covariance; McResult
  parallel.hpp       deterministic block-parallel reduction (OpenMP)
  black_scholes.hpp  closed-form prices, Greeks, implied volatility
  european.hpp       plain / antithetic / control-variate pricing, MC Greeks
  asian.hpp          geometric closed form, arithmetic Asian with geometric control
  american.hpp       CRR binomial tree, two-pass Longstaff-Schwartz
src/price.cpp        command-line pricer
bench/bench.cpp      benchmarks -> results/*.csv, results/BENCHMARKS.md
tests/tests.cpp      54 checks: closed forms, finite differences, statistical bounds, determinism
```

## References

- P. Glasserman, *Monte Carlo Methods in Financial Engineering*, Springer, 2004.
- F. Longstaff, E. Schwartz, "Valuing American Options by Simulation: A Simple Least-Squares Approach", *Review of Financial Studies* 14(1), 2001.
- A. Kemna, A. Vorst, "A Pricing Method for Options Based on Average Asset Values", *Journal of Banking & Finance* 14, 1990.
- J. Cox, S. Ross, M. Rubinstein, "Option Pricing: A Simplified Approach", *Journal of Financial Economics* 7, 1979.
- D. Blackman, S. Vigna, "Scrambled Linear Pseudorandom Number Generators", *ACM TOMS* 47(4), 2021.
- T. Chan, G. Golub, R. LeVeque, "Algorithms for Computing the Sample Variance", *The American Statistician* 37(3), 1983.
