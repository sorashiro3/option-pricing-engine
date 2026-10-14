# Benchmarks

Compiler: MSVC 1944. OpenMP threads available: 4.

## European call: estimators at equal cost

S=K=100, r=5%, sigma=20%, T=1. Black-Scholes = 10.450584. 4,000,000 paths, 4 threads.

| Method | Price | Std. error | Error / SE | Time (ms) | Variance reduction | Efficiency gain |
|---|---|---|---|---|---|---|
| Plain | 10.45139 | 0.00736 | +0.11 | 60.7 | 1.00x | 1.00x |
| Antithetic | 10.45644 | 0.00520 | +1.13 | 42.4 | 2.00x | 2.87x |
| Control variate | 10.45657 | 0.00281 | +2.13 | 68.8 | 6.88x | 6.07x |

Variance reduction compares SE^2 x paths, so antithetic pairs are charged for both paths. Efficiency gain also charges for run time (variance x time).

## Thread scaling (plain estimator, 16,000,000 paths)

| Threads | Time (ms) | Speedup | Paths / s | Price | Bit-identical to 1 thread |
|---|---|---|---|---|---|
| 1 | 687.2 | 1.00x | 23.3M | 10.455089 | yes |
| 2 | 392.5 | 1.75x | 40.8M | 10.455089 | yes |
| 4 | 243.9 | 2.82x | 65.6M | 10.455089 | yes |

## Convergence: error shrinks as 1/sqrt(N)

| Paths | Plain abs. error | Plain SE | Control-variate abs. error | Control-variate SE |
|---|---|---|---|---|
| 1000 | 0.25481 | 0.45305 | 0.11253 | 0.17520 |
| 10000 | 0.11311 | 0.14755 | 0.01847 | 0.05654 |
| 100000 | 0.03174 | 0.04667 | 0.00631 | 0.01778 |
| 1000000 | 0.00173 | 0.01474 | 0.00907 | 0.00562 |
| 10000000 | 0.00607 | 0.00466 | 0.00269 | 0.00177 |

## Monte Carlo Greeks vs closed form (4,000,000 paths, one simulation)

| Greek | Estimator | Monte Carlo | Std. error | Black-Scholes | Error / SE |
|---|---|---|---|---|---|
| Price | discounted payoff | 10.451395 | 0.007361 | 10.450584 | +0.11 |
| Delta | pathwise | 0.636499 | 0.000288 | 0.636831 | -1.15 |
| Gamma | likelihood ratio / pathwise | 0.018771 | 0.000019 | 0.018762 | +0.48 |
| Vega | pathwise | 37.542286 | 0.037848 | 37.524035 | +0.48 |

## Arithmetic Asian call: geometric-average control variate

S=K=100, r=5%, sigma=20%, T=1, 1,000,000 paths.

| Monitoring dates | Plain price | Plain SE | CV price | CV SE | Correlation | Variance reduction | Time plain / CV (ms) |
|---|---|---|---|---|---|---|---|
| 12 | 6.14960 | 0.00851 | 6.156109 | 0.000237 | 0.99961 | 1292x | 168 / 175 |
| 52 | 5.85172 | 0.00811 | 5.854182 | 0.000224 | 0.99962 | 1304x | 909 / 842 |
| 252 | 5.79274 | 0.00802 | 5.782067 | 0.000220 | 0.99962 | 1322x | 3691 / 3898 |

## American put: Longstaff-Schwartz vs CRR binomial

K=40, r=6%, T=1, 50 exercise dates. LSM fits on 100,000 paths and prices on 1,000,000 independent paths. CRR uses 5,000 steps (continuous exercise).

| S | sigma | European (BS) | CRR American | LSM out-of-sample | SE | LSM in-sample | LSM - CRR | LSM time (ms) |
|---|---|---|---|---|---|---|---|---|
| 36 | 0.2 | 3.8443 | 4.4867 | 4.4730 | 0.0029 | 4.4793 | -0.0137 | 740 |
| 40 | 0.2 | 2.0664 | 2.3195 | 2.3114 | 0.0027 | 2.3053 | -0.0081 | 916 |
| 44 | 0.2 | 1.0169 | 1.1130 | 1.1069 | 0.0021 | 1.1095 | -0.0061 | 920 |
| 36 | 0.4 | 6.7114 | 7.1092 | 7.1049 | 0.0061 | 7.0994 | -0.0043 | 866 |
| 40 | 0.4 | 5.0596 | 5.3181 | 5.3102 | 0.0057 | 5.3082 | -0.0080 | 908 |
| 44 | 0.4 | 3.7828 | 3.9531 | 3.9459 | 0.0052 | 3.9369 | -0.0072 | 916 |

The 50-date Bermudan is worth slightly less than the continuously exercisable American option, and the out-of-sample LSM value is a lower bound, so small negative gaps are expected.
