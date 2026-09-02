// hamming.cpp
//
// Hamming distance between a query x (b bits wide) and every row of a DB
// with n rows. Rows are bit-packed into 64-bit words; the DB is one flat
// malloc'd buffer of n * ceil(b/64) words.
//
// Rows are sharded across a persistent worker pool: each thread owns a
// contiguous row range and writes only its own slice of the output, so the
// compute path needs no locking or atomics.
//
// Three compute paths, chosen at compile time by ISA availability:
//   - scalar : XOR + POPCNT (__builtin_popcountll)
//   - AVX2   : XOR + vpshufb nibble-lookup popcount, accumulated with vpsadbw
//   - AVX-512: XOR + native VPOPCNTQ when compiled with -mavx512vpopcntdq
//
// Build: g++ -O3 -march=native -std=c++17 -pthread hamming.cpp -o hamming
// Usage: ./hamming <b:data width in bits> <n:db rows> [iters] [threads]
//        b and n accept size suffixes: 1K/1M/1G (SI, 1000-based)
//                                      1Ki/1Mi/1Gi (IEC, 1024-based)

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

static inline size_t words_per_row(size_t b) { return (b + 63) / 64; }

// ---------------------------------------------------------------------------
// Size parsing: "4096", "1K", "16M", "2Gi", "512KiB"
//   K/M/G  -> 1000, 1000^2, 1000^3   (SI)
//   Ki/Mi/Gi -> 1024, 1024^2, 1024^3 (IEC)
// A trailing 'B' is accepted and ignored.
// ---------------------------------------------------------------------------
static bool parse_size(const char* s, size_t* out) {
    if (!s || !*s) return false;
    char* end = nullptr;
    errno = 0;
    unsigned long long v = strtoull(s, &end, 10);
    if (end == s || errno == ERANGE) return false;

    unsigned long long mult = 1;
    if (*end) {
        char unit = *end++;
        unsigned long long base = 1000;
        if (*end == 'i' || *end == 'I') { base = 1024; ++end; }
        switch (unit) {
            case 'k': case 'K': mult = base;                   break;
            case 'm': case 'M': mult = base * base;             break;
            case 'g': case 'G': mult = base * base * base;      break;
            default: return false;
        }
        if (*end == 'b' || *end == 'B') ++end;
        if (*end) return false;
    }
    if (v != 0 && mult > (~0ULL) / v) return false;  // overflow
    *out = static_cast<size_t>(v * mult);
    return true;
}

static void human_bytes(double bytes, char* buf, size_t buflen) {
    const char* unit[] = {"B", "KiB", "MiB", "GiB", "TiB"};
    int u = 0;
    while (bytes >= 1024.0 && u < 4) { bytes /= 1024.0; ++u; }
    snprintf(buf, buflen, "%.2f %s", bytes, unit[u]);
}

// ---------------------------------------------------------------------------
// Persistent worker pool. run() executes fn(t) for t in [0, size()) and
// returns once every shard is done. The calling thread runs shard 0, so a
// pool of size 1 costs nothing.
// ---------------------------------------------------------------------------
class ThreadPool {
public:
    explicit ThreadPool(size_t nthreads) : nthreads_(nthreads ? nthreads : 1) {
        for (size_t i = 1; i < nthreads_; ++i)
            workers_.emplace_back([this, i] { worker(i); });
    }

    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lk(m_);
            stop_ = true;
            ++gen_;
        }
        cv_.notify_all();
        for (auto& t : workers_) t.join();
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    size_t size() const { return nthreads_; }

    void run(std::function<void(size_t)> fn) {
        {
            std::lock_guard<std::mutex> lk(m_);
            fn_ = std::move(fn);
            done_ = 0;
            ++gen_;
        }
        cv_.notify_all();
        fn_(0);  // caller thread takes shard 0
        std::unique_lock<std::mutex> lk(m_);
        done_cv_.wait(lk, [this] { return done_ == nthreads_ - 1; });
    }

private:
    void worker(size_t id) {
        size_t local_gen = 0;
        for (;;) {
            std::function<void(size_t)> fn;
            {
                std::unique_lock<std::mutex> lk(m_);
                cv_.wait(lk, [this, local_gen] { return gen_ != local_gen; });
                local_gen = gen_;
                if (stop_) return;
                fn = fn_;
            }
            fn(id);
            {
                std::lock_guard<std::mutex> lk(m_);
                ++done_;
            }
            done_cv_.notify_one();
        }
    }

    size_t nthreads_;
    std::vector<std::thread> workers_;
    std::mutex m_;
    std::condition_variable cv_, done_cv_;
    std::function<void(size_t)> fn_;
    size_t gen_ = 0;
    size_t done_ = 0;
    bool stop_ = false;
};

// Contiguous row range owned by shard t of T (remainder spread over the
// first `n % T` shards).
static inline void shard_range(size_t n, size_t t, size_t T, size_t* r0, size_t* r1) {
    size_t base = n / T, rem = n % T;
    *r0 = t * base + (t < rem ? t : rem);
    *r1 = *r0 + base + (t < rem ? 1 : 0);
}

// ---------------------------------------------------------------------------
// DB allocation (plain malloc; vector loads below are unaligned, so no
// special alignment is required)
// ---------------------------------------------------------------------------
static void* xmalloc(size_t bytes, const char* what) {
    void* p = malloc(bytes);
    if (!p) {
        char hb[32];
        human_bytes(static_cast<double>(bytes), hb, sizeof hb);
        fprintf(stderr, "malloc failed for %s (%s)\n", what, hb);
        exit(1);
    }
    return p;
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

static void hamming_rows_scalar(const uint64_t* db, const uint64_t* x, size_t w,
                                size_t r0, size_t r1, uint32_t* out) {
    for (size_t r = r0; r < r1; ++r)
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
static void hamming_rows_vec(const uint64_t* db, const uint64_t* x, size_t w,
                             size_t r0, size_t r1, uint32_t* out) {
    for (size_t r = r0; r < r1; ++r)
        out[r] = hamming_row_vec(db + r * w, x, w);
}
#endif

// ---------------------------------------------------------------------------
// Multi-threaded entry point: distance from x to every row of db.
// ---------------------------------------------------------------------------
static void hamming_all(const uint64_t* db, const uint64_t* x, size_t n, size_t b,
                        uint32_t* out, ThreadPool& pool, bool use_vec) {
    const size_t w = words_per_row(b);
    const size_t T = pool.size();
    pool.run([&](size_t t) {
        size_t r0, r1;
        shard_range(n, t, T, &r0, &r1);
#ifdef HAVE_VEC_PATH
        if (use_vec) { hamming_rows_vec(db, x, w, r0, r1, out); return; }
#else
        (void)use_vec;
#endif
        hamming_rows_scalar(db, x, w, r0, r1, out);
    });
}

// ---------------------------------------------------------------------------
// Test / benchmark driver
// ---------------------------------------------------------------------------

// Fill rows [r0, r1) with random bits, masking off bits past width b in each
// row's last word so distances stay within [0, b].
static void fill_rows(uint64_t* db, size_t w, size_t b, size_t r0, size_t r1, uint64_t seed) {
    std::mt19937_64 rng(seed);
    uint64_t tail_mask = (b % 64) ? ((uint64_t(1) << (b % 64)) - 1) : ~uint64_t(0);
    for (size_t r = r0; r < r1; ++r) {
        uint64_t* row = db + r * w;
        for (size_t i = 0; i < w; ++i) row[i] = rng();
        row[w - 1] &= tail_mask;
    }
}

// Populated in parallel with the same sharding the compute pass uses, so on a
// NUMA box each thread first-touches the pages it will later read.
static void fill_db(uint64_t* db, size_t n, size_t b, ThreadPool& pool, uint64_t seed) {
    const size_t w = words_per_row(b);
    const size_t T = pool.size();
    pool.run([&](size_t t) {
        size_t r0, r1;
        shard_range(n, t, T, &r0, &r1);
        fill_rows(db, w, b, r0, r1, seed + t);
    });
}

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr,
                "usage: %s <b:data width in bits> <n:db rows> [iters] [threads]\n"
                "  b, n accept size suffixes: 1K/1M/1G (SI), 1Ki/1Mi/1Gi (IEC)\n"
                "  threads: 0 or omitted = hardware concurrency\n"
                "example: %s 1024 1M 10 8\n",
                argv[0], argv[0]);
        return 1;
    }

    size_t b = 0, n = 0;
    if (!parse_size(argv[1], &b)) { fprintf(stderr, "bad data width: %s\n", argv[1]); return 1; }
    if (!parse_size(argv[2], &n)) { fprintf(stderr, "bad row count: %s\n", argv[2]); return 1; }
    if (b == 0 || n == 0) { fprintf(stderr, "b and n must be > 0\n"); return 1; }

    size_t iters = (argc > 3) ? strtoull(argv[3], nullptr, 10) : 10;
    if (iters == 0) iters = 1;

    size_t nthreads = (argc > 4) ? strtoull(argv[4], nullptr, 10) : 0;
    if (nthreads == 0) {
        nthreads = std::thread::hardware_concurrency();
        if (nthreads == 0) nthreads = 1;
    }
    if (nthreads > n) nthreads = n;  // no point in empty shards

    const size_t w = words_per_row(b);
    if (n > (~size_t(0)) / (w * sizeof(uint64_t))) {
        fprintf(stderr, "db size (b=%zu, n=%zu) overflows size_t\n", b, n);
        return 1;
    }
    const size_t db_bytes = n * w * sizeof(uint64_t);

    ThreadPool pool(nthreads);

    uint64_t* db = static_cast<uint64_t*>(xmalloc(db_bytes, "db"));
    uint64_t* x  = static_cast<uint64_t*>(xmalloc(w * sizeof(uint64_t), "query"));
    uint32_t* out_scalar = static_cast<uint32_t*>(xmalloc(n * sizeof(uint32_t), "out_scalar"));
    uint32_t* out_vec    = static_cast<uint32_t*>(xmalloc(n * sizeof(uint32_t), "out_vec"));

    char hb[32];
    human_bytes(static_cast<double>(db_bytes), hb, sizeof hb);
    printf("b=%zu bits, n=%zu rows, %zu words/row, db=%s, threads=%zu\n",
           b, n, w, hb, nthreads);

    fill_db(db, n, b, pool, 42);
    fill_rows(x, w, b, 0, 1, 12345);

    using clk = std::chrono::steady_clock;
    auto gbps = [&](double ms) {
        return double(db_bytes) / (ms * 1e-3) / (1024.0 * 1024.0 * 1024.0);
    };

    auto t0 = clk::now();
    for (size_t it = 0; it < iters; ++it)
        hamming_all(db, x, n, b, out_scalar, pool, /*use_vec=*/false);
    double scalar_ms =
        std::chrono::duration<double, std::milli>(clk::now() - t0).count() / iters;
    printf("scalar POPCNT   : %8.3f ms/pass  %6.2f GiB/s\n", scalar_ms, gbps(scalar_ms));

#ifdef HAVE_VEC_PATH
    auto t1 = clk::now();
    for (size_t it = 0; it < iters; ++it)
        hamming_all(db, x, n, b, out_vec, pool, /*use_vec=*/true);
    double vec_ms =
        std::chrono::duration<double, std::milli>(clk::now() - t1).count() / iters;
    printf("%-16s: %8.3f ms/pass  %6.2f GiB/s  (%.2fx)\n",
           VEC_NAME, vec_ms, gbps(vec_ms), scalar_ms / vec_ms);

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

    size_t show = n < 5 ? n : 5;
    for (size_t r = 0; r < show; ++r)
        printf("dist(x, row[%zu]) = %u\n", r, out_scalar[r]);

    free(db);
    free(x);
    free(out_scalar);
    free(out_vec);
    return 0;
}
