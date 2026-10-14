#pragma once
// Arithmetic-average Asian options, discretely monitored at t_i = i T / m.
//
// The arithmetic average of lognormals has no closed form, but the geometric
// average G = (prod S_{t_i})^(1/m) is itself lognormal, so an option on G has an
// exact price (Kemna & Vorst, 1990). G tracks the arithmetic average A very
// closely path by path, which makes the geometric option an excellent control
// variate for the arithmetic one.

#include <cmath>
#include <cstdint>

#include "black_scholes.hpp"
#include "european.hpp"
#include "parallel.hpp"
#include "rng.hpp"
#include "stats.hpp"

namespace mcp {

// Exact price of a discretely monitored geometric-average Asian option.
// ln G ~ N(mu, v) with
//   mu = ln S0 + (r - sigma^2/2) T (m + 1) / (2m)
//   v  = sigma^2 T (m + 1)(2m + 1) / (6 m^2)
// With m = 1 this reduces to Black-Scholes, which the tests check.
inline double geometric_asian_price(const Option& o, int monitoring_dates) {
    const double m = static_cast<double>(monitoring_dates);
    const double mu = std::log(o.spot) + (o.rate - 0.5 * o.vol * o.vol) * o.expiry * (m + 1.0) / (2.0 * m);
    const double v = o.vol * o.vol * o.expiry * (m + 1.0) * (2.0 * m + 1.0) / (6.0 * m * m);
    const double sd = std::sqrt(v);
    const double d2 = (mu - std::log(o.strike)) / sd;
    const double d1 = d2 + sd;
    const double forward_g = std::exp(mu + 0.5 * v);  // E[G]
    const double discount = std::exp(-o.rate * o.expiry);
    return o.type == OptionType::Call ? discount * (forward_g * norm_cdf(d1) - o.strike * norm_cdf(d2))
                                      : discount * (o.strike * norm_cdf(-d2) - forward_g * norm_cdf(-d1));
}

enum class AsianMethod { Plain, GeometricControl };

inline McResult price_asian_arithmetic(const Option& o, int monitoring_dates, AsianMethod method,
                                       const McConfig& cfg = McConfig()) {
    const int m = monitoring_dates;
    const double dt = o.expiry / m;
    const double step_drift = (o.rate - 0.5 * o.vol * o.vol) * dt;
    const double step_vol = o.vol * std::sqrt(dt);
    const double discount = std::exp(-o.rate * o.expiry);
    const double log_spot = std::log(o.spot);

    // Simulates one path; returns discounted arithmetic and geometric payoffs.
    auto simulate = [&](NormalGenerator& z, double& arithmetic_payoff, double& geometric_payoff) {
        double log_s = log_spot;
        double sum_s = 0.0, sum_log_s = 0.0;
        for (int k = 0; k < m; ++k) {
            log_s += step_drift + step_vol * z();
            sum_s += std::exp(log_s);
            sum_log_s += log_s;
        }
        arithmetic_payoff = discount * payoff(o.type, sum_s / m, o.strike);
        geometric_payoff = discount * payoff(o.type, std::exp(sum_log_s / m), o.strike);
    };

    if (method == AsianMethod::Plain) {
        const RunningStats acc = parallel_blocks<RunningStats>(
            cfg.paths, cfg.block_size, [&](std::int64_t b, std::int64_t begin, std::int64_t end, RunningStats& s) {
                NormalGenerator z(cfg.seed, static_cast<std::uint64_t>(b));
                double a, g;
                for (std::int64_t i = begin; i < end; ++i) {
                    simulate(z, a, g);
                    s.add(a);
                }
            });
        return {acc.mean, acc.std_error(), acc.n};
    }

    const RunningCovariance acc = parallel_blocks<RunningCovariance>(
        cfg.paths, cfg.block_size, [&](std::int64_t b, std::int64_t begin, std::int64_t end, RunningCovariance& s) {
            NormalGenerator z(cfg.seed, static_cast<std::uint64_t>(b));
            double a, g;
            for (std::int64_t i = begin; i < end; ++i) {
                simulate(z, a, g);
                s.add(g, a);
            }
        });
    const double beta = acc.beta();
    const double price = acc.mean_y - beta * (acc.mean_x - geometric_asian_price(o, m));
    const double residual_var = (acc.syy - beta * acc.sxy) / static_cast<double>(acc.n - 1);
    return {price, std::sqrt(residual_var / static_cast<double>(acc.n)), acc.n};
}

}  // namespace mcp
