// Command-line pricer.
//
//   price [--style european|asian|american] [--type call|put]
//         [--spot 100] [--strike 100] [--rate 0.05] [--vol 0.2] [--expiry 1]
//         [--paths 1000000] [--dates 12] [--threads N] [--seed 20260527]
//
// European: all three estimators, Monte Carlo Greeks and the closed form.
// Asian:    arithmetic average, plain vs geometric control variate.
// American: Longstaff-Schwartz vs a CRR binomial tree.

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "mcpricer/american.hpp"
#include "mcpricer/asian.hpp"
#include "mcpricer/black_scholes.hpp"
#include "mcpricer/european.hpp"

using namespace mcp;

namespace {

void usage() {
    std::printf(
        "usage: price [--style european|asian|american] [--type call|put]\n"
        "             [--spot S] [--strike K] [--rate r] [--vol sigma] [--expiry T]\n"
        "             [--paths N] [--dates M] [--threads N] [--seed N]\n");
}

void print_result(const char* label, const McResult& r, double ms) {
    std::printf("  %-18s %12.6f   SE %.6f   95%% CI [%.6f, %.6f]   %8.1f ms\n", label, r.price, r.std_error,
                r.ci95_low(), r.ci95_high(), ms);
}

template <class F>
double timed(F&& f) {
    const auto t0 = std::chrono::steady_clock::now();
    f();
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

int main(int argc, char** argv) {
    Option o{100.0, 100.0, 0.05, 0.20, 1.0, OptionType::Call};
    std::string style = "european";
    McConfig cfg;
    int dates = 12;

    for (int i = 1; i < argc; ++i) {
        const char* key = argv[i];
        if (std::strcmp(key, "--help") == 0 || std::strcmp(key, "-h") == 0) {
            usage();
            return 0;
        }
        if (i + 1 >= argc) {
            usage();
            return 2;
        }
        const char* value = argv[++i];
        if (std::strcmp(key, "--style") == 0) style = value;
        else if (std::strcmp(key, "--type") == 0) o.type = std::strcmp(value, "put") == 0 ? OptionType::Put : OptionType::Call;
        else if (std::strcmp(key, "--spot") == 0) o.spot = std::atof(value);
        else if (std::strcmp(key, "--strike") == 0) o.strike = std::atof(value);
        else if (std::strcmp(key, "--rate") == 0) o.rate = std::atof(value);
        else if (std::strcmp(key, "--vol") == 0) o.vol = std::atof(value);
        else if (std::strcmp(key, "--expiry") == 0) o.expiry = std::atof(value);
        else if (std::strcmp(key, "--paths") == 0) cfg.paths = std::atoll(value);
        else if (std::strcmp(key, "--dates") == 0) dates = std::atoi(value);
        else if (std::strcmp(key, "--threads") == 0) set_threads(std::atoi(value));
        else if (std::strcmp(key, "--seed") == 0) cfg.seed = std::strtoull(value, nullptr, 10);
        else {
            usage();
            return 2;
        }
    }
    if (o.spot <= 0 || o.strike <= 0 || o.vol <= 0 || o.expiry <= 0 || cfg.paths < 2 || dates < 1) {
        std::fprintf(stderr, "spot, strike, vol, expiry, paths and dates must be positive\n");
        return 2;
    }

    std::printf("%s %s  S=%g K=%g r=%g sigma=%g T=%g  paths=%lld threads=%d seed=%llu\n\n", style.c_str(),
                o.type == OptionType::Call ? "call" : "put", o.spot, o.strike, o.rate, o.vol, o.expiry,
                static_cast<long long>(cfg.paths), max_threads(), static_cast<unsigned long long>(cfg.seed));

    if (style == "european") {
        const Greeks bs = bs_greeks(o);
        std::printf("  %-18s %12.6f\n", "Black-Scholes", bs.price);
        for (EuroMethod m : {EuroMethod::Plain, EuroMethod::Antithetic, EuroMethod::ControlVariate}) {
            McResult r;
            const double ms = timed([&] { r = price_european(o, m, cfg); });
            print_result(to_string(m), r, ms);
        }
        const McGreeks g = greeks_european(o, cfg);
        std::printf("\n  %-8s %12s %10s %12s\n", "Greek", "Monte Carlo", "SE", "Closed form");
        std::printf("  %-8s %12.6f %10.6f %12.6f\n", "Delta", g.delta.price, g.delta.std_error, bs.delta);
        std::printf("  %-8s %12.6f %10.6f %12.6f\n", "Gamma", g.gamma.price, g.gamma.std_error, bs.gamma);
        std::printf("  %-8s %12.6f %10.6f %12.6f\n", "Vega", g.vega.price, g.vega.std_error, bs.vega);
        std::printf("  %-8s %12s %10s %12.6f\n", "Theta", "-", "-", bs.theta);
        std::printf("  %-8s %12s %10s %12.6f\n", "Rho", "-", "-", bs.rho);
    } else if (style == "asian") {
        std::printf("  %-18s %12.6f   (geometric average, closed form, %d dates)\n", "Geometric Asian",
                    geometric_asian_price(o, dates), dates);
        McResult r;
        double ms = timed([&] { r = price_asian_arithmetic(o, dates, AsianMethod::Plain, cfg); });
        print_result("Arithmetic, plain", r, ms);
        ms = timed([&] { r = price_asian_arithmetic(o, dates, AsianMethod::GeometricControl, cfg); });
        print_result("Arithmetic, CV", r, ms);
    } else if (style == "american") {
        LsmConfig lsm;
        lsm.exercise_dates = dates;
        lsm.price_paths = cfg.paths;
        lsm.seed = cfg.seed;
        std::printf("  %-18s %12.6f\n", "European (BS)", bs_price(o));
        std::printf("  %-18s %12.6f   (CRR, 5000 steps)\n", "American (CRR)", binomial_price(o, 5000, true));
        LsmResult r;
        const double ms = timed([&] { r = price_american_lsm(o, lsm); });
        print_result("LSM out-of-sample", r.price, ms);
        std::printf("  %-18s %12.6f\n", "LSM in-sample", r.in_sample);
    } else {
        usage();
        return 2;
    }
    return 0;
}
