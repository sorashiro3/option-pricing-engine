// Benchmarks: accuracy, variance reduction, Greeks, thread scaling, Asian and
// American pricing. Prints Markdown tables and writes CSVs to results/.
//
// Usage: bench [results_dir]

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "mcpricer/american.hpp"
#include "mcpricer/asian.hpp"
#include "mcpricer/black_scholes.hpp"
#include "mcpricer/european.hpp"

using namespace mcp;
using Clock = std::chrono::steady_clock;

namespace {

template <class F>
double time_ms(F&& f, int repeats = 3) {
    // Best of `repeats`: the least-disturbed run is the best estimate of the cost.
    double best = 1e300;
    for (int r = 0; r < repeats; ++r) {
        const auto t0 = Clock::now();
        f();
        const double ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
        if (ms < best) best = ms;
    }
    return best;
}

std::FILE* open_csv(const std::string& dir, const char* name, const char* header) {
    const std::string path = dir + "/" + name;
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (f) std::fprintf(f, "%s\n", header);
    return f;
}

void close_csv(std::FILE* f) {
    if (f) std::fclose(f);
}

std::string compiler() {
#if defined(_MSC_VER)
    return "MSVC " + std::to_string(_MSC_VER);
#elif defined(__clang__)
    return std::string("Clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("GCC ") + __VERSION__;
#else
    return "unknown compiler";
#endif
}

const Option kAtmCall{100.0, 100.0, 0.05, 0.20, 1.0, OptionType::Call};

void european_methods(const std::string& dir) {
    std::printf("\n## European call: estimators at equal cost\n\n");
    std::printf("S=K=100, r=5%%, sigma=20%%, T=1. Black-Scholes = %.6f. 4,000,000 paths, %d threads.\n\n",
                bs_price(kAtmCall), max_threads());
    std::printf("| Method | Price | Std. error | Error / SE | Time (ms) | Variance reduction | Efficiency gain |\n");
    std::printf("|---|---|---|---|---|---|---|\n");
    std::FILE* csv = open_csv(dir, "european_methods.csv",
                              "method,price,std_error,error_over_se,time_ms,variance_reduction,efficiency_gain");

    McConfig cfg;
    cfg.paths = 4000000;
    const double truth = bs_price(kAtmCall);
    double plain_var = 0.0, plain_ms = 0.0;
    for (EuroMethod m : {EuroMethod::Plain, EuroMethod::Antithetic, EuroMethod::ControlVariate}) {
        McResult r;
        const double ms = time_ms([&] { r = price_european(kAtmCall, m, cfg); });
        if (m == EuroMethod::Plain) {
            plain_var = r.variance_per_path();
            plain_ms = ms;
        }
        const double vrf = plain_var / r.variance_per_path();
        const double efficiency = vrf * plain_ms / ms;  // variance x time, relative to plain
        const double z = (r.price - truth) / r.std_error;
        std::printf("| %s | %.5f | %.5f | %+.2f | %.1f | %.2fx | %.2fx |\n", to_string(m), r.price, r.std_error, z, ms,
                    vrf, efficiency);
        if (csv) std::fprintf(csv, "%s,%.8f,%.8f,%.4f,%.3f,%.4f,%.4f\n", to_string(m), r.price, r.std_error, z, ms, vrf,
                              efficiency);
    }
    close_csv(csv);
    std::printf("\nVariance reduction compares SE^2 x paths, so antithetic pairs are charged for both paths. "
                "Efficiency gain also charges for run time (variance x time).\n");
}

void thread_scaling(const std::string& dir) {
    std::printf("\n## Thread scaling (plain estimator, 16,000,000 paths)\n\n");
    std::printf("| Threads | Time (ms) | Speedup | Paths / s | Price | Bit-identical to 1 thread |\n");
    std::printf("|---|---|---|---|---|---|\n");
    std::FILE* csv = open_csv(dir, "thread_scaling.csv", "threads,time_ms,speedup,paths_per_second,price");

    McConfig cfg;
    cfg.paths = 16000000;
    const int hardware = max_threads();
    double base_ms = 0.0, base_price = 0.0;
    for (int threads = 1; threads <= hardware; threads *= 2) {
        set_threads(threads);
        McResult r;
        const double ms = time_ms([&] { r = price_european(kAtmCall, EuroMethod::Plain, cfg); });
        if (threads == 1) {
            base_ms = ms;
            base_price = r.price;
        }
        const double rate = static_cast<double>(cfg.paths) / (ms / 1000.0);
        std::printf("| %d | %.1f | %.2fx | %.1fM | %.6f | %s |\n", threads, ms, base_ms / ms, rate / 1e6, r.price,
                    r.price == base_price ? "yes" : "NO");
        if (csv) std::fprintf(csv, "%d,%.3f,%.4f,%.0f,%.10f\n", threads, ms, base_ms / ms, rate, r.price);
        if (threads * 2 > hardware && threads != hardware) threads = hardware / 2;  // always finish on `hardware`
    }
    set_threads(hardware);
    close_csv(csv);
}

void convergence(const std::string& dir) {
    std::printf("\n## Convergence: error shrinks as 1/sqrt(N)\n\n");
    std::printf("| Paths | Plain abs. error | Plain SE | Control-variate abs. error | Control-variate SE |\n");
    std::printf("|---|---|---|---|---|\n");
    std::FILE* csv = open_csv(dir, "convergence.csv", "paths,plain_price,plain_se,cv_price,cv_se,black_scholes");

    const double truth = bs_price(kAtmCall);
    for (std::int64_t n = 1000; n <= 10000000; n *= 10) {
        McConfig cfg;
        cfg.paths = n;
        const McResult plain = price_european(kAtmCall, EuroMethod::Plain, cfg);
        const McResult cv = price_european(kAtmCall, EuroMethod::ControlVariate, cfg);
        std::printf("| %lld | %.5f | %.5f | %.5f | %.5f |\n", static_cast<long long>(n), std::fabs(plain.price - truth),
                    plain.std_error, std::fabs(cv.price - truth), cv.std_error);
        if (csv) std::fprintf(csv, "%lld,%.8f,%.8f,%.8f,%.8f,%.8f\n", static_cast<long long>(n), plain.price,
                              plain.std_error, cv.price, cv.std_error, truth);
    }
    close_csv(csv);
}

void greeks(const std::string& dir) {
    std::printf("\n## Monte Carlo Greeks vs closed form (4,000,000 paths, one simulation)\n\n");
    std::printf("| Greek | Estimator | Monte Carlo | Std. error | Black-Scholes | Error / SE |\n");
    std::printf("|---|---|---|---|---|---|\n");
    std::FILE* csv = open_csv(dir, "greeks.csv", "greek,estimator,mc,std_error,black_scholes");

    McConfig cfg;
    cfg.paths = 4000000;
    const McGreeks mc = greeks_european(kAtmCall, cfg);
    const Greeks bs = bs_greeks(kAtmCall);
    struct Row {
        const char* name;
        const char* estimator;
        McResult mc;
        double exact;
    };
    for (const Row& row : {Row{"Price", "discounted payoff", mc.price, bs.price},
                           Row{"Delta", "pathwise", mc.delta, bs.delta},
                           Row{"Gamma", "likelihood ratio / pathwise", mc.gamma, bs.gamma},
                           Row{"Vega", "pathwise", mc.vega, bs.vega}}) {
        std::printf("| %s | %s | %.6f | %.6f | %.6f | %+.2f |\n", row.name, row.estimator, row.mc.price,
                    row.mc.std_error, row.exact, (row.mc.price - row.exact) / row.mc.std_error);
        if (csv) std::fprintf(csv, "%s,%s,%.8f,%.8f,%.8f\n", row.name, row.estimator, row.mc.price, row.mc.std_error,
                              row.exact);
    }
    close_csv(csv);
}

void asian(const std::string& dir) {
    std::printf("\n## Arithmetic Asian call: geometric-average control variate\n\n");
    std::printf("S=K=100, r=5%%, sigma=20%%, T=1, 1,000,000 paths.\n\n");
    std::printf("| Monitoring dates | Plain price | Plain SE | CV price | CV SE | Correlation | Variance reduction | Time plain / CV (ms) |\n");
    std::printf("|---|---|---|---|---|---|---|---|\n");
    std::FILE* csv = open_csv(dir, "asian.csv",
                              "monitoring_dates,plain_price,plain_se,cv_price,cv_se,variance_reduction,plain_ms,cv_ms");

    McConfig cfg;
    cfg.paths = 1000000;
    for (int m : {12, 52, 252}) {
        McResult plain, cv;
        const double plain_ms = time_ms([&] { plain = price_asian_arithmetic(kAtmCall, m, AsianMethod::Plain, cfg); }, 1);
        const double cv_ms =
            time_ms([&] { cv = price_asian_arithmetic(kAtmCall, m, AsianMethod::GeometricControl, cfg); }, 1);
        const double vrf = plain.variance_per_path() / cv.variance_per_path();
        const double rho = std::sqrt(1.0 - 1.0 / vrf);  // VRF = 1 / (1 - rho^2)
        std::printf("| %d | %.5f | %.5f | %.6f | %.6f | %.5f | %.0fx | %.0f / %.0f |\n", m, plain.price,
                    plain.std_error, cv.price, cv.std_error, rho, vrf, plain_ms, cv_ms);
        if (csv) std::fprintf(csv, "%d,%.8f,%.8f,%.8f,%.8f,%.3f,%.3f,%.3f\n", m, plain.price, plain.std_error, cv.price,
                              cv.std_error, vrf, plain_ms, cv_ms);
    }
    close_csv(csv);
}

void american(const std::string& dir) {
    std::printf("\n## American put: Longstaff-Schwartz vs CRR binomial\n\n");
    std::printf("K=40, r=6%%, T=1, 50 exercise dates. LSM fits on 100,000 paths and prices on 1,000,000 "
                "independent paths. CRR uses 5,000 steps (continuous exercise).\n\n");
    std::printf("| S | sigma | European (BS) | CRR American | LSM out-of-sample | SE | LSM in-sample | LSM - CRR | LSM time (ms) |\n");
    std::printf("|---|---|---|---|---|---|---|---|---|\n");
    std::FILE* csv = open_csv(dir, "american.csv",
                              "spot,vol,european,crr_american,lsm,lsm_se,lsm_in_sample,lsm_ms");

    LsmConfig cfg;
    for (double vol : {0.2, 0.4}) {
        for (double spot : {36.0, 40.0, 44.0}) {
            const Option o{spot, 40.0, 0.06, vol, 1.0, OptionType::Put};
            const double crr = binomial_price(o, 5000, true);
            LsmResult r;
            const double ms = time_ms([&] { r = price_american_lsm(o, cfg); }, 1);
            std::printf("| %.0f | %.1f | %.4f | %.4f | %.4f | %.4f | %.4f | %+.4f | %.0f |\n", spot, vol, bs_price(o),
                        crr, r.price.price, r.price.std_error, r.in_sample, r.price.price - crr, ms);
            if (csv) std::fprintf(csv, "%.1f,%.2f,%.8f,%.8f,%.8f,%.8f,%.8f,%.3f\n", spot, vol, bs_price(o), crr,
                                  r.price.price, r.price.std_error, r.in_sample, ms);
        }
    }
    close_csv(csv);
    std::printf("\nThe 50-date Bermudan is worth slightly less than the continuously exercisable American option, "
                "and the out-of-sample LSM value is a lower bound, so small negative gaps are expected.\n");
}

}  // namespace

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : "results";
    std::printf("# Benchmarks\n\nCompiler: %s. OpenMP threads available: %d.\n", compiler().c_str(), max_threads());
    european_methods(dir);
    thread_scaling(dir);
    convergence(dir);
    greeks(dir);
    asian(dir);
    american(dir);
    return 0;
}
