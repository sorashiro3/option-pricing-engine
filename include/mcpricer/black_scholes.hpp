#pragma once
// Closed-form Black-Scholes prices and Greeks: the reference every Monte Carlo
// estimator in this project is validated against.

#include <algorithm>
#include <cmath>
#include <limits>

namespace mcp {

enum class OptionType { Call, Put };

struct Option {
    double spot;    // S0
    double strike;  // K
    double rate;    // r, continuously compounded risk-free rate
    double vol;     // sigma, annualised volatility
    double expiry;  // T, years
    OptionType type = OptionType::Call;
};

struct Greeks {
    double price = 0.0;
    double delta = 0.0;  // dV/dS
    double gamma = 0.0;  // d2V/dS2
    double vega = 0.0;   // dV/dsigma (per 1.00 of vol, not per 1%)
    double theta = 0.0;  // dV/dt (per year, calendar time passing)
    double rho = 0.0;    // dV/dr
};

inline double norm_pdf(double x) { return 0.3989422804014327 * std::exp(-0.5 * x * x); }

inline double norm_cdf(double x) { return 0.5 * std::erfc(-x * 0.7071067811865476); }

inline double payoff(OptionType type, double spot_at_expiry, double strike) {
    return type == OptionType::Call ? std::max(spot_at_expiry - strike, 0.0)
                                    : std::max(strike - spot_at_expiry, 0.0);
}

inline Greeks bs_greeks(const Option& o) {
    const double sqrt_t = std::sqrt(o.expiry);
    const double d1 = (std::log(o.spot / o.strike) + (o.rate + 0.5 * o.vol * o.vol) * o.expiry) / (o.vol * sqrt_t);
    const double d2 = d1 - o.vol * sqrt_t;
    const double discount = std::exp(-o.rate * o.expiry);
    const double pdf_d1 = norm_pdf(d1);

    Greeks g;
    g.gamma = pdf_d1 / (o.spot * o.vol * sqrt_t);
    g.vega = o.spot * pdf_d1 * sqrt_t;
    const double decay = -o.spot * pdf_d1 * o.vol / (2.0 * sqrt_t);
    if (o.type == OptionType::Call) {
        g.price = o.spot * norm_cdf(d1) - o.strike * discount * norm_cdf(d2);
        g.delta = norm_cdf(d1);
        g.theta = decay - o.rate * o.strike * discount * norm_cdf(d2);
        g.rho = o.strike * o.expiry * discount * norm_cdf(d2);
    } else {
        g.price = o.strike * discount * norm_cdf(-d2) - o.spot * norm_cdf(-d1);
        g.delta = norm_cdf(d1) - 1.0;
        g.theta = decay + o.rate * o.strike * discount * norm_cdf(-d2);
        g.rho = -o.strike * o.expiry * discount * norm_cdf(-d2);
    }
    return g;
}

inline double bs_price(const Option& o) { return bs_greeks(o).price; }

// Implied volatility by safeguarded Newton-Raphson: Newton steps on vega, with a
// bisection step whenever Newton would leave the bracket that still contains the
// root. The price is strictly increasing in vol, so the bracket always shrinks.
// Returns NaN when the target price lies outside the attainable range.
inline double implied_vol(Option o, double target_price, double tol = 1e-10, int max_iter = 200) {
    double lo = 1e-6, hi = 5.0;
    o.vol = lo;
    const double price_lo = bs_price(o);
    o.vol = hi;
    const double price_hi = bs_price(o);
    if (target_price < price_lo || target_price > price_hi) return std::numeric_limits<double>::quiet_NaN();

    double vol = 0.2;
    for (int i = 0; i < max_iter; ++i) {
        o.vol = vol;
        const Greeks g = bs_greeks(o);
        const double diff = g.price - target_price;
        if (std::fabs(diff) < tol) return vol;
        if (diff > 0.0)
            hi = vol;
        else
            lo = vol;
        double next = g.vega > 0.0 ? vol - diff / g.vega : 0.5 * (lo + hi);
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        vol = next;
    }
    return vol;
}

}  // namespace mcp
