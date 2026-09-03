// Compares every popcount kernel this build contains, across row lengths, so
// the dispatch thresholds in hamming.cpp can be checked on the CPU at hand.
// Includes hamming.cpp so it measures the shipped kernels, not a copy.
//
// Build: make bench      Usage: ./bench [threads] [working set MiB]

#define main hamming_main
#include "hamming.cpp"
#undef main

#include <random>

static size_t g_T, g_iters;
static uint64_t *g_db, *g_x;
static uint32_t* g_out;

template <uint32_t Dist(const uint64_t*, const uint64_t*, size_t)>
static double gibs(size_t n, size_t w) {
    double best = 0;
    for (int rep = 0; rep < 3; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        scan<Dist>(g_db, g_x, n, w, g_out, g_T, g_iters);
        double s = std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - t0).count();
        double g = double(n) * w * 8 * g_iters / s / 1073741824.0;
        if (g > best) best = g;
    }
    return best;
}

int main(int argc, char** argv) {
    g_T = argc > 1 ? strtoull(argv[1], nullptr, 10) : 0;
    if (!g_T) g_T = std::thread::hardware_concurrency();
    size_t mib = argc > 2 ? strtoull(argv[2], nullptr, 10) : 8;

    printf("threads=%zu, working set ~%zu MiB (sized to stay cache-resident so the\n"
           "kernels, not DRAM, are what is being compared)\n\n", g_T, mib);
    printf("%-6s %-4s %12s", "bits", "w", "scalar");
#if defined(__AVX2__)
    printf(" %12s", "avx2");
#endif
#if defined(__AVX512BW__)
    printf(" %12s", "avx512 shuf");
#endif
#if defined(__AVX512VPOPCNTDQ__) && defined(__AVX512F__)
    printf(" %12s", "vpopcntq");
#endif
    printf("   GiB/s\n");

    std::mt19937_64 rng(1);
    for (size_t w : {1u, 2u, 4u, 6u, 8u, 12u, 16u, 24u, 32u, 64u}) {
        size_t n = mib * 1024 * 1024 / (w * 8);
        g_iters = 200;
        g_db = (uint64_t*)xmalloc(n * w * 8);
        g_x = (uint64_t*)xmalloc(w * 8);
        g_out = (uint32_t*)xmalloc(n * 4);
        for (size_t i = 0; i < n * w; ++i) g_db[i] = rng();
        for (size_t i = 0; i < w; ++i) g_x[i] = rng();

        printf("%-6zu %-4zu %12.2f", w * 64, w, gibs<dist_scalar>(n, w));
#if defined(__AVX2__)
        printf(" %12.2f", gibs<dist_avx2>(n, w));
#endif
#if defined(__AVX512BW__)
        printf(" %12.2f", gibs<dist_avx512>(n, w));
#endif
#if defined(__AVX512VPOPCNTDQ__) && defined(__AVX512F__)
        printf(" %12.2f", gibs<dist_vpopcnt>(n, w));
#endif
        printf("\n");
        free(g_db); free(g_x); free(g_out);
    }
    printf("\nSet W_VPOPCNT / W_AVX512 / W_AVX2 in hamming.cpp to the smallest w\n"
           "at which each kernel wins its column.\n");
    return 0;
}
