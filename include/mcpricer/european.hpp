#pragma once
// European options under geometric Brownian motion:
//   S_T = S_0 exp((r - sigma^2/2) T + sigma sqrt(T) Z),  Z ~ N(0, 1)
// One exact draw per path, so there is no time-discretisation error and the only
// error left is statistical, which the standard error measures.

#include <cmath>
#include <cstdint>

#include "black_scholes.hpp"
#include "parallel.hpp"
#include "rng.hpp"
#include "stats.hpp"

namespace mcp {

struct McConfig {
    std::int64_t paths = 1000000;
    std::uint64_t seed = 20260527;
    std::int64_t block_size = kDefaultBlockSize;
};

enum class EuroMethod {
    Plain,           // crude Monte Carlo
    Antithetic,      // pair Z with -Z
    ControlVariate,  // control: discounted S_T, whose risk-neutral mean is exactly S_0
};

inline const char* to_string(EuroMethod m) {
    switch (m) {
        case EuroMethod::Plain: return "Plain";
        case EuroMethod::Antithetic: return "Antithetic";
        case EuroMethod::ControlVariate: return "Control variate";
    }
    return "?";
}

namespace detail {
struct GbmTerminal {
    double drift, diffusion, discount;
    explicit GbmTerminal(const Option& o)
        : drift((o.rate - 0.5 * o.vol * o.vol) * o.expiry),
          diffusion(o.vol * std::sqrt(o.expiry)),
          discount(std::exp(-o.rate * o.expiry)) {}
};
}  // namespace detail

inline McResult price_european(const Option& o, EuroMethod method, const McConfig& cfg = McConfig()) {
    const detail::GbmTerminal g(o);

    if (method == EuroMethod::Plain) {
        const RunningStats acc = parallel_blocks<RunningStats>(
            cfg.paths, cfg.block_size, [&](std::int64_t b, std::int64_t begin, std::int64_t end, RunningStats& s) {
                NormalGenerator z(cfg.seed, static_cast<std::uint64_t>(b));
                for (std::int64_t i = begin; i < end; ++i) {
                    const double st = o.spot * std::exp(g.drift + g.diffusion * z());
                    s.add(g.discount * payoff(o.type, st, o.strike));
                }
            });
        return {acc.mean, acc.std_error(), acc.n};
    }

    if (method == EuroMethod::Antithetic) {
        // Each sample averages the payoffs of Z and -Z: two paths, one sample.
        const std::int64_t pairs = cfg.paths / 2;
        const RunningStats acc = parallel_blocks<RunningStats>(
            pairs, cfg.block_size / 2, [&](std::int64_t b, std::int64_t begin, std::int64_t end, RunningStats& s) {
                NormalGenerator z(cfg.seed, static_cast<std::uint64_t>(b));
                for (std::int64_t i = begin; i < end; ++i) {
                    const double shock = g.diffusion * z();
                    const double up = o.spot * std::exp(g.drift + shock);
                    const double down = o.spot * std::exp(g.drift - shock);
                    s.add(0.5 * g.discount * (payoff(o.type, up, o.strike) + payoff(o.type, down, o.strike)));
                }
            });
        return {acc.mean, acc.std_error(), 2 * acc.n};
    }

    // Control variate: Y = discounted payoff, X = discounted S_T with E[X] = S_0.
    // Estimator Y_bar - beta (X_bar - S_0) with beta = Cov(X, Y) / Var(X) fitted on
    // the same paths; its variance is Var(Y) (1 - rho^2).
    const RunningCovariance acc = parallel_blocks<RunningCovariance>(
        cfg.paths, cfg.block_size, [&](std::int64_t b, std::int64_t begin, std::int64_t end, RunningCovariance& s) {
            NormalGenerator z(cfg.seed, static_cast<std::uint64_t>(b));
            for (std::int64_t i = begin; i < end; ++i) {
                const double st = o.spot * std::exp(g.drift + g.diffusion * z());
                s.add(g.discount * st, g.discount * payoff(o.type, st, o.strike));
            }
        });
    const double beta = acc.beta();
    const double price = acc.mean_y - beta * (acc.mean_x - o.spot);
    const double residual_var = (acc.syy - beta * acc.sxy) / static_cast<double>(acc.n - 1);
    return {price, std::sqrt(residual_var / static_cast<double>(acc.n)), acc.n};
}

// Monte Carlo Greeks, all from one set of paths.
//   delta, vega: pathwise derivative (the payoff is Lipschitz in S_T).
//   gamma:       pathwise derivative is useless here (the payoff's second
//                derivative is a Dirac delta), so differentiate the pathwise
//                delta once more with the likelihood-ratio method (the "LR-PW"
//                estimator, Glasserman 2004, Sec. 7.3).
struct McGreeks {
    McResult price, delta, gamma, vega;
};

inline McGreeks greeks_european(const Option& o, const McConfig& cfg = McConfig()) {
    struct Acc {
        RunningStats price, delta, gamma, vega;
        void merge(const Acc& other) {
            price.merge(other.price);
            delta.merge(other.delta);
            gamma.merge(other.gamma);
            vega.merge(other.vega);
        }
    };

    const detail::GbmTerminal g(o);
    const double sqrt_t = std::sqrt(o.expiry);
    const double sign = o.type == OptionType::Call ? 1.0 : -1.0;

    const Acc acc = parallel_blocks<Acc>(
        cfg.paths, cfg.block_size, [&](std::int64_t b, std::int64_t begin, std::int64_t end, Acc& s) {
            NormalGenerator gen(cfg.seed, static_cast<std::uint64_t>(b));
            for (std::int64_t i = begin; i < end; ++i) {
                const double z = gen();
                const double st = o.spot * std::exp(g.drift + g.diffusion * z);
                const bool in_the_money = o.type == OptionType::Call ? st > o.strike : st < o.strike;
                // d payoff / d S_T: +1 / -1 in the money, 0 otherwise.
                const double dpayoff = in_the_money ? sign * g.discount : 0.0;

                s.price.add(g.discount * payoff(o.type, st, o.strike));
                s.delta.add(dpayoff * st / o.spot);
                s.vega.add(dpayoff * st * (sqrt_t * z - o.vol * o.expiry));
                s.gamma.add(dpayoff * st / (o.spot * o.spot) * (z / (o.vol * sqrt_t) - 1.0));
            }
        });

    auto result = [](const RunningStats& s) { return McResult{s.mean, s.std_error(), s.n}; };
    return {result(acc.price), result(acc.delta), result(acc.gamma), result(acc.vega)};
}

}  // namespace mcp
