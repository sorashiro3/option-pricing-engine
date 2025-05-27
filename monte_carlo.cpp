#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <omp.h>

using namespace std;

// Black-Scholes analytical solution for European Call Option
double black_scholes_call(double S, double K, double r, double sigma, double T) {
    double d1 = (log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    auto norm_cdf = [](double x) {
        return 0.5 * erfc(-x / sqrt(2));
    };
    return S * norm_cdf(d1) - K * exp(-r * T) * norm_cdf(d2);
}

// Monte Carlo simulation for European Call Option
double monte_carlo_call(double S, double K, double r, double sigma, double T, int num_paths) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> dist(0.0, 1.0);

    double sum_payoff = 0.0;
    #pragma omp parallel for reduction(+:sum_payoff)
    for (int i = 0; i < num_paths; ++i) {
        double Z = dist(gen);
        double ST = S * exp((r - 0.5 * sigma * sigma) * T + sigma * sqrt(T) * Z);
        sum_payoff += max(ST - K, 0.0);
    }
    return exp(-r * T) * (sum_payoff / num_paths);
}

// Antithetic Variates for variance reduction
double monte_carlo_call_antithetic(double S, double K, double r, double sigma, double T, int num_paths) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> dist(0.0, 1.0);

    double sum_payoff = 0.0;
    #pragma omp parallel for reduction(+:sum_payoff)
    for (int i = 0; i < num_paths / 2; ++i) {
        double Z = dist(gen);
        double ST1 = S * exp((r - 0.5 * sigma * sigma) * T + sigma * sqrt(T) * Z);
        double ST2 = S * exp((r - 0.5 * sigma * sigma) * T - sigma * sqrt(T) * Z);
        sum_payoff += 0.5 * (max(ST1 - K, 0.0) + max(ST2 - K, 0.0));
    }
    return exp(-r * T) * (sum_payoff / (num_paths / 2));
}

int main() {
    double S = 100.0;    // Spot price
    double K = 100.0;    // Strike price
    double r = 0.05;     // Risk-free rate
    double sigma = 0.2;  // Volatility
    double T = 1.0;      // Time to maturity
    int num_paths = 1000000;

    auto start = chrono::high_resolution_clock::now();
    double bs_price = black_scholes_call(S, K, r, sigma, T);
    auto end = chrono::high_resolution_clock::now();
    cout << "Black-Scholes Price: " << bs_price << ", Time: " << chrono::duration<double>(end - start).count() << "s\n";

    start = chrono::high_resolution_clock::now();
    double mc_price = monte_carlo_call(S, K, r, sigma, T, num_paths);
    end = chrono::high_resolution_clock::now();
    cout << "Monte Carlo Price:   " << mc_price << ", Time: " << chrono::duration<double>(end - start).count() << "s\n";

    start = chrono::high_resolution_clock::now();
    double mc_av_price = monte_carlo_call_antithetic(S, K, r, sigma, T, num_paths);
    end = chrono::high_resolution_clock::now();
    cout << "MC Antithetic Price: " << mc_av_price << ", Time: " << chrono::duration<double>(end - start).count() << "s\n";

    return 0;
}
