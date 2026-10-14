#pragma once
// Random number generation for parallel Monte Carlo.
//
// The original engine shared one std::mt19937 between OpenMP threads, which is a
// data race: results changed run to run and the draws were not independent.
// Here every block of paths owns its own generator, seeded from (seed, stream),
// so a block produces the same numbers no matter which thread runs it.

#include <cmath>
#include <cstdint>

namespace mcp {

// SplitMix64 (Steele, Lea & Flood). Used only to expand a 64-bit seed into the
// 256-bit state of xoshiro256**, as its authors recommend.
inline std::uint64_t splitmix64(std::uint64_t& x) {
    std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// xoshiro256** (Blackman & Vigna, 2018): fast, 2^256 - 1 period, passes BigCrush.
class Xoshiro256 {
public:
    Xoshiro256(std::uint64_t seed, std::uint64_t stream) {
        std::uint64_t x = seed ^ (stream * 0xD1B54A32D192ED03ULL);
        for (auto& word : s_) word = splitmix64(x);
    }

    std::uint64_t next() {
        const std::uint64_t result = rotl(s_[1] * 5, 7) * 9;
        const std::uint64_t t = s_[1] << 17;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);
        return result;
    }

    // Uniform on the open interval (0, 1), built from the top 53 bits.
    double uniform() { return (static_cast<double>(next() >> 11) + 0.5) * (1.0 / 9007199254740992.0); }

private:
    static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    std::uint64_t s_[4];
};

// Standard normal variates by the Marsaglia polar method. Written out rather than
// using std::normal_distribution, whose algorithm differs between standard
// libraries, so a given seed gives the same prices with GCC, Clang and MSVC.
class NormalGenerator {
public:
    NormalGenerator(std::uint64_t seed, std::uint64_t stream) : rng_(seed, stream) {}

    double operator()() {
        if (has_spare_) {
            has_spare_ = false;
            return spare_;
        }
        double u, v, s;
        do {
            u = 2.0 * rng_.uniform() - 1.0;
            v = 2.0 * rng_.uniform() - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);
        const double scale = std::sqrt(-2.0 * std::log(s) / s);
        spare_ = v * scale;
        has_spare_ = true;
        return u * scale;
    }

private:
    Xoshiro256 rng_;
    double spare_ = 0.0;
    bool has_spare_ = false;
};

}  // namespace mcp
