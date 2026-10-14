// Self-contained test runner: no framework, so it builds anywhere the engine does.
// Statistical checks use a 4-sigma tolerance (false-failure rate ~6e-5 per check);
// every run uses fixed seeds, so a pass is reproducible rather than lucky.

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "mcpricer/american.hpp"
#include "mcpricer/asian.hpp"
#include "mcpricer/black_scholes.hpp"
#include "mcpricer/european.hpp"

using namespace mcp;

namespace {

int g_failures = 0;
int g_checks = 0;

void check(bool ok, const std::string& what) {
    ++g_checks;
    if (!ok) {
        ++g_failures;
        std::printf("    FAIL: %s\n", what.c_str());
    }
}

void check_close(double actual, double expected, double tol, const std::string& what) {
    char buf[256];
    std::snprintf(buf, sizeof buf, "%s: got %.8f, expected %.8f (tol %.2e)", what.c_str(), actual, expected, tol);
    check(std::fabs(actual - expected) <= tol, buf);
}

// |estimate - truth| within k standard errors.
void check_within_se(const McResult& r, double truth, double k, const std::string& what) {
    char buf[256];
    const double z = r.std_error > 0 ? (r.price - truth) / r.std_error : 0.0;
    std::snprintf(buf, sizeof buf, "%s: MC %.6f vs %.6f (SE %.2e, z = %+.2f)", what.c_str(), r.price, truth,
                  r.std_error, z);
    check(std::fabs(r.price - truth) <= k * r.std_error, buf);
}

const Option kAtmCall{100.0, 100.0, 0.05, 0.20, 1.0, OptionType::Call};
const Option kOtmPut{100.0, 90.0, 0.03, 0.30, 0.5, OptionType::Put};

void test_black_scholes() {
    // Put-call parity: C - P = S - K e^{-rT}.
    Option call = kAtmCall, put = kAtmCall;
    put.type = OptionType::Put;
    check_close(bs_price(call) - bs_price(put), call.spot - call.strike * std::exp(-call.rate * call.expiry), 1e-12,
                "put-call parity");
    // Hull's textbook value for S=K=100, r=5%, sigma=20%, T=1.
    check_close(bs_price(call), 10.450583572185565, 1e-9, "ATM call reference value");

    // Analytic Greeks against central finite differences.
    for (const Option& base : {kAtmCall, kOtmPut}) {
        const Greeks g = bs_greeks(base);
        const double h = 1e-4;
        Option up = base, dn = base;
        up.spot += h;
        dn.spot -= h;
        check_close(g.delta, (bs_price(up) - bs_price(dn)) / (2 * h), 1e-7, "delta vs finite difference");
        check_close(g.gamma, (bs_price(up) - 2 * bs_price(base) + bs_price(dn)) / (h * h), 1e-4,
                    "gamma vs finite difference");
        up = dn = base;
        up.vol += h;
        dn.vol -= h;
        check_close(g.vega, (bs_price(up) - bs_price(dn)) / (2 * h), 1e-6, "vega vs finite difference");
        up = dn = base;
        up.rate += h;
        dn.rate -= h;
        check_close(g.rho, (bs_price(up) - bs_price(dn)) / (2 * h), 1e-6, "rho vs finite difference");
        up = dn = base;
        up.expiry -= h;  // theta: value change as time passes, i.e. as T shrinks
        dn.expiry += h;
        check_close(g.theta, (bs_price(up) - bs_price(dn)) / (2 * h), 1e-5, "theta vs finite difference");
    }
}

void test_implied_vol() {
    for (double vol : {0.05, 0.2, 0.65, 1.5}) {
        for (const Option& base : {kAtmCall, kOtmPut}) {
            Option o = base;
            o.vol = vol;
            check_close(implied_vol(o, bs_price(o)), vol, 1e-8, "implied vol round trip");
        }
    }
    Option o = kAtmCall;
    check(std::isnan(implied_vol(o, 1000.0)), "implied vol of an impossible price is NaN");
}

void test_european_mc() {
    McConfig cfg;
    cfg.paths = 2000000;
    for (const Option& o : {kAtmCall, kOtmPut}) {
        const double truth = bs_price(o);
        for (EuroMethod m : {EuroMethod::Plain, EuroMethod::Antithetic, EuroMethod::ControlVariate})
            check_within_se(price_european(o, m, cfg), truth, 4.0, std::string("European ") + to_string(m));
    }
    // Both variance-reduction methods must actually reduce variance at equal cost.
    const double plain = price_european(kAtmCall, EuroMethod::Plain, cfg).variance_per_path();
    check(price_european(kAtmCall, EuroMethod::Antithetic, cfg).variance_per_path() < plain,
          "antithetic reduces variance per path");
    check(price_european(kAtmCall, EuroMethod::ControlVariate, cfg).variance_per_path() < plain,
          "control variate reduces variance per path");
}

void test_determinism_across_threads() {
    McConfig cfg;
    cfg.paths = 500000;
    const int saved = max_threads();
    set_threads(1);
    const McResult one = price_european(kAtmCall, EuroMethod::ControlVariate, cfg);
    set_threads(saved > 1 ? saved : 4);
    const McResult many = price_european(kAtmCall, EuroMethod::ControlVariate, cfg);
    set_threads(saved);
    check(one.price == many.price && one.std_error == many.std_error,
          "bit-identical price and SE for 1 thread and N threads");

    McConfig other = cfg;
    other.seed += 1;
    check(price_european(kAtmCall, EuroMethod::Plain, other).price !=
              price_european(kAtmCall, EuroMethod::Plain, cfg).price,
          "a different seed gives a different estimate");
}

void test_mc_greeks() {
    McConfig cfg;
    cfg.paths = 2000000;
    for (const Option& o : {kAtmCall, kOtmPut}) {
        const Greeks exact = bs_greeks(o);
        const McGreeks mc = greeks_european(o, cfg);
        check_within_se(mc.price, exact.price, 4.0, "MC price");
        check_within_se(mc.delta, exact.delta, 4.0, "pathwise delta");
        check_within_se(mc.vega, exact.vega, 4.0, "pathwise vega");
        check_within_se(mc.gamma, exact.gamma, 4.0, "LR-pathwise gamma");
    }
}

void test_asian() {
    // With one monitoring date the geometric Asian is a European option.
    for (const Option& o : {kAtmCall, kOtmPut})
        check_close(geometric_asian_price(o, 1), bs_price(o), 1e-12, "geometric Asian (m = 1) equals Black-Scholes");

    McConfig cfg;
    cfg.paths = 400000;
    for (const Option& o : {kAtmCall, kOtmPut}) {
        const McResult plain = price_asian_arithmetic(o, 12, AsianMethod::Plain, cfg);
        const McResult cv = price_asian_arithmetic(o, 12, AsianMethod::GeometricControl, cfg);
        const double combined_se = std::sqrt(plain.std_error * plain.std_error + cv.std_error * cv.std_error);
        char buf[160];
        std::snprintf(buf, sizeof buf, "Asian plain %.5f vs control variate %.5f", plain.price, cv.price);
        check(std::fabs(plain.price - cv.price) <= 4.0 * combined_se, buf);
        check(cv.std_error * 10.0 < plain.std_error, "geometric control variate cuts the SE at least 10x");
        // AM-GM: the arithmetic average dominates the geometric one path by path.
        if (o.type == OptionType::Call)
            check(cv.price > geometric_asian_price(o, 12), "arithmetic Asian call above geometric");
        else
            check(cv.price < geometric_asian_price(o, 12), "arithmetic Asian put below geometric");
    }
}

void test_binomial() {
    for (const Option& o : {kAtmCall, kOtmPut})
        check_close(binomial_price(o, 4000, false), bs_price(o), 5e-3, "CRR European converges to Black-Scholes");

    // Without dividends an American call is never exercised early.
    check_close(binomial_price(kAtmCall, 2000, true), binomial_price(kAtmCall, 2000, false), 1e-12,
                "American call (no dividends) = European call");
    Option put = kAtmCall;
    put.type = OptionType::Put;
    check(binomial_price(put, 2000, true) > binomial_price(put, 2000, false) + 0.1,
          "American put carries an early-exercise premium");
}

void test_lsm() {
    // Longstaff & Schwartz (2001), Table 1: S=36, K=40, r=6%, sigma=20%, T=1,
    // 50 exercise dates. Finite-difference American value 4.478, LSM 4.472.
    const Option o{36.0, 40.0, 0.06, 0.20, 1.0, OptionType::Put};
    LsmConfig cfg;
    cfg.fit_paths = 100000;
    cfg.price_paths = 400000;
    const LsmResult r = price_american_lsm(o, cfg);
    check_close(r.price.price, 4.472, 0.02, "LSM vs Longstaff-Schwartz Table 1");
    check(r.price.price > bs_price(o), "LSM American put above European put");
    check(r.price.price < binomial_price(o, 4000, true) + 3.0 * r.price.std_error,
          "out-of-sample LSM does not exceed the American value");
}

}  // namespace

int main() {
    const std::vector<std::pair<const char*, std::function<void()>>> tests = {
        {"black_scholes", test_black_scholes},
        {"implied_vol", test_implied_vol},
        {"european_mc", test_european_mc},
        {"determinism_across_threads", test_determinism_across_threads},
        {"mc_greeks", test_mc_greeks},
        {"asian", test_asian},
        {"binomial", test_binomial},
        {"lsm", test_lsm},
    };
    for (const auto& t : tests) {
        const int before = g_failures;
        t.second();
        std::printf("[%s] %s\n", g_failures == before ? " OK " : "FAIL", t.first);
    }
    std::printf("\n%d checks, %d failed\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
