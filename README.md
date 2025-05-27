# Option Pricing Engine

A high-performance simulator for pricing European call options using both **Black-Scholes** and **Monte Carlo** models—designed for real-time, low-latency environments.

---

## Overview

This project implements:

- 📈 **Black-Scholes Analytical Model**  
- 🎲 **Monte Carlo Simulation with Antithetic Variates**  
- ⚙️ **OpenMP-Based Parallelization** for fast, scalable simulation  
- 📊 **Performance Benchmarking** and convergence analysis

---

## Core Features

- Analytical and stochastic pricing of European call options
- Variance reduction using **Antithetic Variates**
- Speedup with **multithreaded Monte Carlo simulations**
- Execution time comparisons under real-time constraints
- Clean, modular C++ code with performance metrics

---

## Parameters (Default)

| Parameter | Description         | Value      |
|----------:|---------------------|------------|
| `S`       | Spot Price          | `100.0`    |
| `K`       | Strike Price        | `100.0`    |
| `r`       | Risk-Free Rate      | `0.05`     |
| `σ`       | Volatility          | `0.2`      |
| `T`       | Time to Maturity    | `1 year`   |
| `N`       | # of Simulations    | `1,000,000`|

---
