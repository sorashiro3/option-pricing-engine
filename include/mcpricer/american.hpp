#pragma once
// American (Bermudan) options.
//
//  * binomial_price: Cox-Ross-Rubinstein lattice, the deterministic reference.
//  * price_american_lsm: Longstaff-Schwartz (2001) least-squares Monte Carlo.
//
// LSM runs in two passes so the estimate has no look-ahead bias:
//   1. Fit: simulate N_fit paths, step backwards through the exercise dates and
//      regress discounted future cash flows on polynomials in S (in-the-money
//      paths only) to estimate the continuation value at each date.
//   2. Price: simulate fresh, independent paths and exercise them with the
//      fitted rule. A fixed (possibly sub-optimal) rule can only under-exercise,
//      so this estimate is a low-biased, honest lower bound.
// The in-sample value from pass 1 is reported as well: it reuses the paths the
// rule was fitted on, so it tends to be biased upwards.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "black_scholes.hpp"
#include "european.hpp"
#include "parallel.hpp"
#include "rng.hpp"
#include "stats.hpp"

namespace mcp {

// CRR binomial tree. american = false prices the European option on the same
// tree, which the tests use to check the lattice against Black-Scholes.
inline double binomial_price(const Option& o, int steps, bool american = true) {
    const double dt = o.expiry / steps;
    const double up = std::exp(o.vol * std::sqrt(dt));
    const double down = 1.0 / up;
    const double p = (std::exp(o.rate * dt) - down) / (up - down);
    const double discount = std::exp(-o.rate * dt);

    std::vector<double> value(static_cast<std::size_t>(steps) + 1);
    for (int i = 0; i <= steps; ++i)
        value[static_cast<std::size_t>(i)] = payoff(o.type, o.spot * std::pow(up, 2 * i - steps), o.strike);

    for (int step = steps - 1; step >= 0; --step) {
        for (int i = 0; i <= step; ++i) {
            const auto idx = static_cast<std::size_t>(i);
            double v = discount * (p * value[idx + 1] + (1.0 - p) * value[idx]);
            if (american) v = std::max(v, payoff(o.type, o.spot * std::pow(up, 2 * i - step), o.strike));
            value[idx] = v;
        }
    }
    return value[0];
}

struct LsmConfig {
    int exercise_dates = 50;          // Bermudan exercise at t_k = k T / M, k = 1..M
    std::int64_t fit_paths = 100000;  // pass 1: regression
    std::int64_t price_paths = 1000000;  // pass 2: independent pricing paths
    std::uint64_t seed = 20260527;
    std::int64_t block_size = 1 << 12;
};

struct LsmResult {
    McResult price;         // out-of-sample (pass 2) estimate
    double in_sample = 0.0; // pass-1 estimate, for comparison
};

namespace detail {

constexpr int kBasis = 4;  // 1, x, x^2, x^3 with x = S / K
using Coefficients = std::array<double, kBasis>;

inline std::array<double, kBasis> basis(double moneyness) {
    return {1.0, moneyness, moneyness * moneyness, moneyness * moneyness * moneyness};
}

// Normal equations (Phi^T Phi) beta = Phi^T y, accumulated block by block.
struct NormalEquations {
    std::array<double, kBasis * kBasis> a{};
    std::array<double, kBasis> rhs{};
    std::int64_t n = 0;

    void add(const std::array<double, kBasis>& phi, double y) {
        for (int r = 0; r < kBasis; ++r) {
            rhs[r] += phi[r] * y;
            for (int c = 0; c < kBasis; ++c) a[r * kBasis + c] += phi[r] * phi[c];
        }
        ++n;
    }
    void merge(const NormalEquations& other) {
        for (int i = 0; i < kBasis * kBasis; ++i) a[i] += other.a[i];
        for (int i = 0; i < kBasis; ++i) rhs[i] += other.rhs[i];
        n += other.n;
    }

    // Gaussian elimination with partial pivoting. False if (nearly) singular.
    bool solve(Coefficients& beta) const {
        std::array<double, kBasis * kBasis> m = a;
        Coefficients b = rhs;
        for (int col = 0; col < kBasis; ++col) {
            int pivot = col;
            for (int r = col + 1; r < kBasis; ++r)
                if (std::fabs(m[r * kBasis + col]) > std::fabs(m[pivot * kBasis + col])) pivot = r;
            if (std::fabs(m[pivot * kBasis + col]) < 1e-12) return false;
            if (pivot != col) {
                for (int c = 0; c < kBasis; ++c) std::swap(m[col * kBasis + c], m[pivot * kBasis + c]);
                std::swap(b[col], b[pivot]);
            }
            for (int r = col + 1; r < kBasis; ++r) {
                const double f = m[r * kBasis + col] / m[col * kBasis + col];
                for (int c = col; c < kBasis; ++c) m[r * kBasis + c] -= f * m[col * kBasis + c];
                b[r] -= f * b[col];
            }
        }
        for (int r = kBasis - 1; r >= 0; --r) {
            double s = b[r];
            for (int c = r + 1; c < kBasis; ++c) s -= m[r * kBasis + c] * beta[c];
            beta[r] = s / m[r * kBasis + r];
        }
        return true;
    }
};

inline double continuation(const Coefficients& beta, double moneyness) {
    const auto phi = basis(moneyness);
    double v = 0.0;
    for (int i = 0; i < kBasis; ++i) v += beta[i] * phi[i];
    return v;
}

// Accumulator that holds nothing: lets parallel_blocks drive a pure "fill" loop.
struct NoAcc {
    void merge(const NoAcc&) {}
};

constexpr std::uint64_t kPricingStreamOffset = 1ULL << 40;  // keeps pass-2 streams disjoint from pass 1

}  // namespace detail

inline LsmResult price_american_lsm(const Option& o, const LsmConfig& cfg = LsmConfig()) {
    using namespace detail;
    const int m = cfg.exercise_dates;
    const double dt = o.expiry / m;
    const double step_drift = (o.rate - 0.5 * o.vol * o.vol) * dt;
    const double step_vol = o.vol * std::sqrt(dt);
    const double step_discount = std::exp(-o.rate * dt);
    const std::int64_t n = cfg.fit_paths;

    // ---- Pass 1: simulate and store the fitting paths (row = path, column = date).
    std::vector<double> spot(static_cast<std::size_t>(n) * static_cast<std::size_t>(m));
    parallel_blocks<NoAcc>(n, cfg.block_size, [&](std::int64_t b, std::int64_t begin, std::int64_t end, NoAcc&) {
        NormalGenerator z(cfg.seed, static_cast<std::uint64_t>(b));
        for (std::int64_t i = begin; i < end; ++i) {
            double s = o.spot;
            double* row = &spot[static_cast<std::size_t>(i) * static_cast<std::size_t>(m)];
            for (int k = 0; k < m; ++k) {
                s *= std::exp(step_drift + step_vol * z());
                row[k] = s;
            }
        }
    });
    auto at = [&](std::int64_t i, int k) {  // spot of path i at date k (1-based)
        return spot[static_cast<std::size_t>(i) * static_cast<std::size_t>(m) + static_cast<std::size_t>(k - 1)];
    };

    // cash[i]: value of path i's cash flow under the current rule, discounted to
    // the date being processed.
    std::vector<double> cash(static_cast<std::size_t>(n));
    for (std::int64_t i = 0; i < n; ++i) cash[static_cast<std::size_t>(i)] = payoff(o.type, at(i, m), o.strike);

    std::vector<Coefficients> beta(static_cast<std::size_t>(m) + 1);
    std::vector<char> has_rule(static_cast<std::size_t>(m) + 1, 0);

    for (int k = m - 1; k >= 1; --k) {
        for (double& c : cash) c *= step_discount;

        const NormalEquations eq = parallel_blocks<NormalEquations>(
            n, cfg.block_size, [&](std::int64_t, std::int64_t begin, std::int64_t end, NormalEquations& e) {
                for (std::int64_t i = begin; i < end; ++i) {
                    const double s = at(i, k);
                    if (payoff(o.type, s, o.strike) > 0.0) e.add(basis(s / o.strike), cash[static_cast<std::size_t>(i)]);
                }
            });
        if (eq.n < 2 * kBasis || !eq.solve(beta[static_cast<std::size_t>(k)])) continue;  // no exercise at this date
        has_rule[static_cast<std::size_t>(k)] = 1;

#pragma omp parallel for schedule(static)
        for (std::int64_t i = 0; i < n; ++i) {
            const double s = at(i, k);
            const double exercise = payoff(o.type, s, o.strike);
            if (exercise > 0.0 && exercise > continuation(beta[static_cast<std::size_t>(k)], s / o.strike))
                cash[static_cast<std::size_t>(i)] = exercise;
        }
    }
    double in_sample = 0.0;
    for (double c : cash) in_sample += c;
    in_sample = step_discount * in_sample / static_cast<double>(n);
    spot.clear();
    spot.shrink_to_fit();

    // ---- Pass 2: price with the fitted rule on independent paths.
    const RunningStats acc = parallel_blocks<RunningStats>(
        cfg.price_paths, cfg.block_size, [&](std::int64_t b, std::int64_t begin, std::int64_t end, RunningStats& s) {
            NormalGenerator z(cfg.seed, kPricingStreamOffset + static_cast<std::uint64_t>(b));
            for (std::int64_t i = begin; i < end; ++i) {
                double st = o.spot;
                double value = 0.0;
                for (int k = 1; k <= m; ++k) {
                    st *= std::exp(step_drift + step_vol * z());
                    const double exercise = payoff(o.type, st, o.strike);
                    if (exercise <= 0.0) continue;
                    const bool last = k == m;
                    const auto kk = static_cast<std::size_t>(k);
                    if (last || (has_rule[kk] && exercise > continuation(beta[kk], st / o.strike))) {
                        value = exercise * std::exp(-o.rate * dt * k);
                        break;
                    }
                }
                s.add(value);
            }
        });

    LsmResult result;
    result.price = {acc.mean, acc.std_error(), acc.n};
    result.in_sample = in_sample;
    // Exercising immediately is always allowed for an American option.
    const double intrinsic = payoff(o.type, o.spot, o.strike);
    if (intrinsic > result.price.price) result.price = {intrinsic, 0.0, acc.n};
    return result;
}

}  // namespace mcp
