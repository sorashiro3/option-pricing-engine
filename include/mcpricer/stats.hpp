#pragma once
// Streaming estimators for Monte Carlo output.
//
// Both accumulators use Welford's update (no catastrophic cancellation from
// sum-of-squares) and Chan et al.'s pairwise merge, so per-block results can be
// combined exactly as if all samples had been added to one accumulator.

#include <cmath>
#include <cstdint>

namespace mcp {

// Mean and variance of one quantity.
struct RunningStats {
    std::int64_t n = 0;
    double mean = 0.0;
    double m2 = 0.0;  // sum of squared deviations from the mean

    void add(double x) {
        ++n;
        const double delta = x - mean;
        mean += delta / static_cast<double>(n);
        m2 += delta * (x - mean);
    }

    void merge(const RunningStats& other) {
        if (other.n == 0) return;
        if (n == 0) {
            *this = other;
            return;
        }
        const double total = static_cast<double>(n + other.n);
        const double delta = other.mean - mean;
        mean += delta * static_cast<double>(other.n) / total;
        m2 += other.m2 + delta * delta * static_cast<double>(n) * static_cast<double>(other.n) / total;
        n += other.n;
    }

    double variance() const { return n > 1 ? m2 / static_cast<double>(n - 1) : 0.0; }
    double std_error() const { return n > 0 ? std::sqrt(variance() / static_cast<double>(n)) : 0.0; }
};

// Means, variances and covariance of a pair (x, y). Used for control variates:
// y is the payoff, x the control whose expectation is known exactly.
struct RunningCovariance {
    std::int64_t n = 0;
    double mean_x = 0.0, mean_y = 0.0;
    double sxx = 0.0, syy = 0.0, sxy = 0.0;  // co-moments

    void add(double x, double y) {
        ++n;
        const double inv_n = 1.0 / static_cast<double>(n);
        const double dx = x - mean_x;
        const double dy = y - mean_y;
        mean_x += dx * inv_n;
        mean_y += dy * inv_n;
        sxx += dx * (x - mean_x);
        syy += dy * (y - mean_y);
        sxy += dx * (y - mean_y);
    }

    void merge(const RunningCovariance& other) {
        if (other.n == 0) return;
        if (n == 0) {
            *this = other;
            return;
        }
        const double total = static_cast<double>(n + other.n);
        const double weight = static_cast<double>(n) * static_cast<double>(other.n) / total;
        const double dx = other.mean_x - mean_x;
        const double dy = other.mean_y - mean_y;
        sxx += other.sxx + dx * dx * weight;
        syy += other.syy + dy * dy * weight;
        sxy += other.sxy + dx * dy * weight;
        mean_x += dx * static_cast<double>(other.n) / total;
        mean_y += dy * static_cast<double>(other.n) / total;
        n += other.n;
    }

    // Optimal control-variate coefficient beta* = Cov(x, y) / Var(x).
    double beta() const { return sxx > 0.0 ? sxy / sxx : 0.0; }

    double correlation() const { return (sxx > 0.0 && syy > 0.0) ? sxy / std::sqrt(sxx * syy) : 0.0; }
};

// Result of a Monte Carlo estimate.
struct McResult {
    double price = 0.0;
    double std_error = 0.0;
    std::int64_t paths = 0;  // simulated paths, i.e. the cost of the estimate

    double ci95_low() const { return price - 1.959963984540054 * std_error; }
    double ci95_high() const { return price + 1.959963984540054 * std_error; }

    // SE^2 x paths: the variance one simulated path contributes. Comparing this
    // across estimators measures variance reduction at equal simulation cost,
    // which is fair to antithetic sampling (two paths per sample).
    double variance_per_path() const { return std_error * std_error * static_cast<double>(paths); }
};

}  // namespace mcp
