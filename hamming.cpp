// hamming.cpp
//
// Hamming distance between a query x (b bits wide) and every row of a DB
// with n rows. Rows are bit-packed into 64-bit words; the DB is one flat
// malloc'd buffer of n * ceil(b/64) words.
//
// Two compute paths, chosen at compile time by ISA availability:
//   - scalar : XOR + POPCNT (__builtin_popcountll)
//   - AVX2   : XOR + vpshufb nibble-lookup popcount, accumulated with vpsadbw
//   - AVX-512: XOR + native VPOPCNTQ when compiled with -mavx512vpopcntdq
//
// Build: g++ -O3 -march=native -std=c++17 hamming.cpp -o hamming
// Usage: ./hamming <b:data width in bits> <n:db rows> [iters]

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <random>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

static inline size_t words_per_row(size_t b) { return (b + 63) / 64; }

// ---------------------------------------------------------------------------
// DB allocation (plain malloc; vector loads below use unaligned loads, so no
// special alignment is required)
// ---------------------------------------------------------------------------
static uint64_t* db_alloc(size_t n, size_t b) {
    size_t w = words_per_row(b);
    uint64_t* db = static_cast<uint64_t*>(malloc(n * w * sizeof(uint64_t)));
    if (!db) {
        fprintf(stderr, "malloc failed (%zu bytes)\n", n * w * sizeof(uint64_t));
        exit(1);
    }
    return db;
}

// ---------------------------------------------------------------------------
// Scalar path: one POPCNT per 64-bit word
// ---------------------------------------------------------------------------
static inline uint32_t hamming_row_scalar(const uint64_t* row, const uint64_t* x, size_t w) {
    uint32_t d = 0;
    for (size_t i = 0; i < w; ++i)
        d += static_cast<uint32_t>(__builtin_popcountll(row[i] ^ x[i]));
    return d;
}

static void hamming_all_scalar(const uint64_t* db, const uint64_t* x,
                               size_t n, size_t b, uint32_t* out) {
    size_t w = words_per_row(b);
    for (size_t r = 0; r < n; ++r)
        out[r] = hamming_row_scalar(db + r * w, x, w);
}

#if defined(__AVX512VPOPCNTDQ__) && defined(__AVX512F__)
// ---------------------------------------------------------------------------
// AVX-512 path: native 64-bit lane popcount (VPOPCNTQ), 8 words per step
// ---------------------------------------------------------------------------
static inline uint32_t hamming_row_vec(const uint64_t* row, const uint64_t* x, size_t w) {
    __m512i acc = _mm512_setzero_si512();
    size_t i = 0;
    for (; i + 8 <= w; i += 8) {
        __m512i a = _mm512_loadu_si512(row + i);
        __m512i q = _mm512_loadu_si512(x + i);
        acc = _mm512_add_epi64(acc, _mm512_popcnt_epi64(_mm512_xor_si512(a, q)));
    }
    uint32_t d = static_cast<uint32_t>(_mm512_reduce_add_epi64(acc));
    for (; i < w; ++i)
        d += static_cast<uint32_t>(__builtin_popcountll(row[i] ^ x[i]));
    return d;
}
#define HAVE_VEC_PATH 1
#define VEC_NAME "AVX-512 VPOPCNTQ"

#elif defined(__AVX2__)
// ---------------------------------------------------------------------------
// AVX2 path: vpshufb nibble-lookup popcount, 4 words (256 bits) per step.
// Per-byte counts are summed into 64-bit lanes with vpsadbw against zero.
// ---------------------------------------------------------------------------
static inline uint32_t hamming_row_vec(const uint64_t* row, const uint64_t* x, size_t w) {
    const __m256i lut = _mm256_setr_epi8(
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
    const __m256i low_mask = _mm256_set1_epi8(0x0f);
    __m256i acc = _mm256_setzero_si256();

    size_t i = 0;
    for (; i + 4 <= w; i += 4) {
        __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(row + i));
        __m256i q = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x + i));
        __m256i v = _mm256_xor_si256(a, q);

        __m256i lo = _mm256_and_si256(v, low_mask);
        __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), low_mask);
        __m256i cnt = _mm256_add_epi8(_mm256_shuffle_epi8(lut, lo),
                                      _mm256_shuffle_epi8(lut, hi));
        acc = _mm256_add_epi64(acc, _mm256_sad_epu8(cnt, _mm256_setzero_si256()));
    }

    uint64_t lanes[4];
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(lanes), acc);
    uint32_t d = static_cast<uint32_t>(lanes[0] + lanes[1] + lanes[2] + lanes[3]);

    for (; i < w; ++i)
        d += static_cast<uint32_t>(__builtin_popcountll(row[i] ^ x[i]));
    return d;
}
#define HAVE_VEC_PATH 1
#define VEC_NAME "AVX2 vpshufb"
#endif

#ifdef HAVE_VEC_PATH
static void hamming_all_vec(const uint64_t* db, const uint64_t* x,
                            size_t n, size_t b, uint32_t* out) {
    size_t w = words_per_row(b);
    for (size_t r = 0; r < n; ++r)
        out[r] = hamming_row_vec(db + r * w, x, w);
}
#endif

// ---------------------------------------------------------------------------
// Test / benchmark driver
// ---------------------------------------------------------------------------
static void fill_random(uint64_t* buf, size_t nwords, size_t b, size_t w, std::mt19937_64& rng) {
    // Mask so bits beyond width b in the last word of each row stay zero.
    uint64_t tail_mask = (b % 64) ? ((uint64_t(1) << (b % 64)) - 1) : ~uint64_t(0);
    for (size_t i = 0; i < nwords; ++i) {
        buf[i] = rng();
        if (i % w == w - 1) buf[i] &= tail_mask;
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <b:data width in bits> <n:db rows> [iters]\n", argv[0]);
        return 1;
    }
    size_t b = strtoull(argv[1], nullptr, 10);
    size_t n = strtoull(argv[2], nullptr, 10);
    size_t iters = (argc > 3) ? strtoull(argv[3], nullptr, 10) : 10;
    if (b == 0 || n == 0) {
        fprintf(stderr, "b and n must be > 0\n");
        return 1;
    }
    size_t w = words_per_row(b);

    uint64_t* db = db_alloc(n, b);
    uint64_t* x  = static_cast<uint64_t*>(malloc(w * sizeof(uint64_t)));
    uint32_t* out_scalar = static_cast<uint32_t*>(malloc(n * sizeof(uint32_t)));
    uint32_t* out_vec    = static_cast<uint32_t*>(malloc(n * sizeof(uint32_t)));
    if (!x || !out_scalar || !out_vec) {
        fprintf(stderr, "malloc failed\n");
        return 1;
    }

    std::mt19937_64 rng(42);
    fill_random(db, n * w, b, w, rng);
    fill_random(x, w, b, w, rng);

    printf("b=%zu bits, n=%zu rows, %zu words/row, db=%.2f MiB\n",
           b, n, w, double(n * w * sizeof(uint64_t)) / (1024.0 * 1024.0));

    using clk = std::chrono::steady_clock;

    auto t0 = clk::now();
    for (size_t it = 0; it < iters; ++it)
        hamming_all_scalar(db, x, n, b, out_scalar);
    double scalar_ms =
        std::chrono::duration<double, std::milli>(clk::now() - t0).count() / iters;
    printf("scalar POPCNT   : %.3f ms/pass\n", scalar_ms);

#ifdef HAVE_VEC_PATH
    auto t1 = clk::now();
    for (size_t it = 0; it < iters; ++it)
        hamming_all_vec(db, x, n, b, out_vec);
    double vec_ms =
        std::chrono::duration<double, std::milli>(clk::now() - t1).count() / iters;
    printf("%-16s: %.3f ms/pass (%.2fx)\n", VEC_NAME, vec_ms, scalar_ms / vec_ms);

    for (size_t r = 0; r < n; ++r) {
        if (out_scalar[r] != out_vec[r]) {
            fprintf(stderr, "MISMATCH row %zu: scalar=%u vec=%u\n",
                    r, out_scalar[r], out_vec[r]);
            return 1;
        }
    }
    printf("verify          : scalar == vector for all %zu rows\n", n);
#else
    printf("(no AVX2/AVX-512 at compile time; scalar path only)\n");
#endif

    // Show a few results.
    size_t show = n < 5 ? n : 5;
    for (size_t r = 0; r < show; ++r)
        printf("dist(x, row[%zu]) = %u\n", r, out_scalar[r]);

    free(db);
    free(x);
    free(out_scalar);
    free(out_vec);
    return 0;
}
