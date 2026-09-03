// Checks every compiled-in popcount kernel against the scalar reference.
// Includes hamming.cpp directly so it tests the shipped code, not a copy.
//
// Build: g++ -O3 -march=native -std=c++17 -pthread verify.cpp -o verify

#define main hamming_main
#include "hamming.cpp"
#undef main

#include <random>

static int failures = 0;

template <uint32_t Dist(const uint64_t*, const uint64_t*, size_t)>
static void check(const char* name, const uint64_t* a, const uint64_t* q, size_t w,
                  const char* what) {
    uint32_t want = dist_scalar(a, q, w), got = Dist(a, q, w);
    if (want != got) {
        printf("FAIL %-18s w=%-3zu %-12s scalar=%u kernel=%u\n", name, w, what, want, got);
        ++failures;
    }
}

static void check_all(const uint64_t* a, const uint64_t* q, size_t w, const char* what) {
#if defined(__AVX512VPOPCNTDQ__) && defined(__AVX512F__)
    check<dist_vpopcnt>("AVX-512 VPOPCNTQ", a, q, w, what);
#endif
#if defined(__AVX512BW__)
    check<dist_avx512>("AVX-512 vpshufb", a, q, w, what);
#endif
#if defined(__AVX2__)
    check<dist_avx2>("AVX2 vpshufb", a, q, w, what);
#endif
    (void)a; (void)q; (void)w; (void)what;
}

int main() {
    const size_t maxw = 40;
    std::vector<uint64_t> a(maxw), q(maxw);
    std::mt19937_64 rng(7);

    // Every width up to 40 words covers each kernel's main loop plus every
    // possible tail remainder (1..7 words for the 512-bit paths).
    for (size_t w = 1; w <= maxw; ++w) {
        for (int rep = 0; rep < 200; ++rep) {
            for (size_t i = 0; i < w; ++i) { a[i] = rng(); q[i] = rng(); }
            check_all(a.data(), q.data(), w, "random");
        }
        for (size_t i = 0; i < w; ++i) { a[i] = 0; q[i] = 0; }
        check_all(a.data(), q.data(), w, "zeros");            // distance 0
        for (size_t i = 0; i < w; ++i) { a[i] = ~0ull; q[i] = 0; }
        check_all(a.data(), q.data(), w, "all-ones");         // distance 64*w
        for (size_t i = 0; i < w; ++i) { a[i] = 0; q[i] = 0; }
        for (size_t bit = 0; bit < 64 * w; ++bit) {           // one bit set
            a[bit / 64] = uint64_t(1) << (bit % 64);
            check_all(a.data(), q.data(), w, "one-bit");
            a[bit / 64] = 0;
        }
    }

    // The byte counters in the nibble-lookup kernels must not overflow: with
    // every bit set, each byte holds 8 and vpsadbw sums 8 of them per lane.
    for (size_t i = 0; i < maxw; ++i) { a[i] = ~0ull; q[i] = 0; }
    for (size_t w = 1; w <= maxw; ++w) check_all(a.data(), q.data(), w, "saturated");

    printf(failures ? "%d FAILURES\n" : "all kernels match scalar (%d failures)\n", failures);
    return failures != 0;
}
