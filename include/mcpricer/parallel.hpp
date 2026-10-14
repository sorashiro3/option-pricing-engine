#pragma once
// Deterministic parallel reduction over blocks of Monte Carlo paths.
//
// Paths are split into fixed-size blocks. Block b always draws from RNG stream b
// and fills its own accumulator; the accumulators are then merged in block order
// on one thread. The floating-point result is therefore bit-identical for any
// number of threads and any OpenMP schedule, which makes runs reproducible and
// lets the tests compare 1-thread and N-thread output exactly.

#include <algorithm>
#include <cstdint>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace mcp {

constexpr std::int64_t kDefaultBlockSize = 1 << 14;  // 16,384 paths per block

// fn(block_index, begin, end, accumulator) simulates items [begin, end).
// Acc must be default-constructible and provide merge(const Acc&).
template <class Acc, class BlockFn>
Acc parallel_blocks(std::int64_t n_items, std::int64_t block_size, BlockFn fn) {
    if (block_size < 1) block_size = 1;
    const std::int64_t n_blocks = (n_items + block_size - 1) / block_size;
    std::vector<Acc> partial(static_cast<std::size_t>(n_blocks));

#pragma omp parallel for schedule(dynamic, 1)
    for (std::int64_t b = 0; b < n_blocks; ++b) {
        const std::int64_t begin = b * block_size;
        const std::int64_t end = std::min(begin + block_size, n_items);
        fn(b, begin, end, partial[static_cast<std::size_t>(b)]);
    }

    Acc total;
    for (const Acc& p : partial) total.merge(p);
    return total;
}

inline void set_threads(int n) {
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)n;
#endif
}

inline int max_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

}  // namespace mcp
