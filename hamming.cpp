// Hamming distance from a query x to every row of a bit-packed DB.
// Rows are ceil(b/64) uint64_t words in one malloc'd buffer of n rows.
//
// Build: g++ -O3 -march=native -std=c++17 -pthread hamming.cpp -o hamming
// Usage: ./hamming <bits> <rows> [iters] [threads]
//        bits/rows take size suffixes: 1K/1M/1G (SI), 1Ki/1Mi/1Gi (IEC)
// Package and DRAM energy come from Intel RAPL via powercap sysfs; set
// RAPL_PATH to override /sys/class/powercap.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <cstring>
#include <chrono>
#include <random>
#include <string>
#include <thread>
#include <vector>
#include <dirent.h>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

static size_t words(size_t b) { return (b + 63) / 64; }

#if defined(__AVX512VPOPCNTDQ__) && defined(__AVX512F__)
static const char* kVecName = "AVX-512 VPOPCNTQ";
static uint32_t dist_vec(const uint64_t* a, const uint64_t* q, size_t w) {
    __m512i acc = _mm512_setzero_si512();
    size_t i = 0;
    for (; i + 8 <= w; i += 8) {
        __m512i v = _mm512_xor_si512(_mm512_loadu_si512(a + i), _mm512_loadu_si512(q + i));
        acc = _mm512_add_epi64(acc, _mm512_popcnt_epi64(v));
    }
    uint32_t d = _mm512_reduce_add_epi64(acc);
    for (; i < w; ++i) d += __builtin_popcountll(a[i] ^ q[i]);
    return d;
}

#elif defined(__AVX2__)
static const char* kVecName = "AVX2 vpshufb";
static uint32_t dist_vec(const uint64_t* a, const uint64_t* q, size_t w) {
    // Byte popcount by nibble lookup, then vpsadbw sums bytes into 64-bit lanes.
    const __m256i lut = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
                                         0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
    const __m256i mask = _mm256_set1_epi8(0x0f);
    __m256i acc = _mm256_setzero_si256();
    size_t i = 0;
    for (; i + 4 <= w; i += 4) {
        __m256i v = _mm256_xor_si256(_mm256_loadu_si256((const __m256i*)(a + i)),
                                     _mm256_loadu_si256((const __m256i*)(q + i)));
        __m256i lo = _mm256_shuffle_epi8(lut, _mm256_and_si256(v, mask));
        __m256i hi = _mm256_shuffle_epi8(lut, _mm256_and_si256(_mm256_srli_epi16(v, 4), mask));
        acc = _mm256_add_epi64(acc, _mm256_sad_epu8(_mm256_add_epi8(lo, hi),
                                                    _mm256_setzero_si256()));
    }
    uint64_t l[4];
    _mm256_storeu_si256((__m256i*)l, acc);
    uint32_t d = l[0] + l[1] + l[2] + l[3];
    for (; i < w; ++i) d += __builtin_popcountll(a[i] ^ q[i]);
    return d;
}

#else
static const char* kVecName = "scalar POPCNT";
static uint32_t dist_vec(const uint64_t* a, const uint64_t* q, size_t w) {
    uint32_t d = 0;
    for (size_t i = 0; i < w; ++i) d += __builtin_popcountll(a[i] ^ q[i]);
    return d;
}
#endif

static bool read_u64(const std::string& path, uint64_t* v) {
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return false;
    bool ok = fscanf(f, "%" SCNu64, v) == 1;
    fclose(f);
    return ok;
}

static std::string read_line(const std::string& path) {
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return "";
    char buf[64] = {0};
    if (!fgets(buf, sizeof buf, f)) buf[0] = 0;
    fclose(f);
    std::string s(buf);
    while (!s.empty() && (s.back() == '\n' || s.back() == ' ')) s.pop_back();
    return s;
}

// Intel RAPL package and DRAM energy counters. Both levels of the powercap
// hierarchy (intel-rapl:N and intel-rapl:N:M) are listed in one flat
// directory, so a single scan finds packages and their dram subdomains.
struct Rapl {
    struct Domain { std::string file; uint64_t range, prev; bool dram; };
    std::vector<Domain> domains;
    size_t npkg = 0, ndram = 0;

    Rapl() {
        const char* root = getenv("RAPL_PATH");
        if (!root) root = "/sys/class/powercap";
        DIR* dir = opendir(root);
        if (!dir) return;
        while (dirent* e = readdir(dir)) {
            if (strncmp(e->d_name, "intel-rapl:", 11) != 0) continue;
            std::string p = std::string(root) + "/" + e->d_name;
            std::string name = read_line(p + "/name");
            bool dram = (name == "dram");
            if (!dram && name.compare(0, 8, "package-") != 0) continue;
            uint64_t range, now;
            if (!read_u64(p + "/max_energy_range_uj", &range)) continue;
            if (!read_u64(p + "/energy_uj", &now)) continue;
            domains.push_back({p + "/energy_uj", range, now, dram});
            (dram ? ndram : npkg)++;
        }
        closedir(dir);
    }

    void start() {
        for (Domain& d : domains) read_u64(d.file, &d.prev);
    }

    // Counters are free-running and wrap at max_energy_range_uj.
    void stop(double* pkg_j, double* dram_j) {
        *pkg_j = *dram_j = 0;
        for (Domain& d : domains) {
            uint64_t now;
            if (!read_u64(d.file, &now)) continue;
            uint64_t delta = now >= d.prev ? now - d.prev : now + d.range - d.prev;
            (d.dram ? *dram_j : *pkg_j) += delta * 1e-6;
        }
    }
};

// Splits [0, n) into one contiguous range per thread; each range is disjoint,
// so callers need no locking.
template <class F>
static void parallel_for(size_t n, size_t nthreads, F fn) {
    std::vector<std::thread> ts;
    for (size_t t = 0; t < nthreads; ++t)
        ts.emplace_back(fn, t * n / nthreads, (t + 1) * n / nthreads);
    for (auto& t : ts) t.join();
}

// iters repeats the scan inside the parallel region, so a benchmark measures
// the kernel rather than thread creation.
template <uint32_t Dist(const uint64_t*, const uint64_t*, size_t)>
static void scan(const uint64_t* db, const uint64_t* x, size_t n, size_t w,
                 uint32_t* out, size_t nthreads, size_t iters = 1) {
    parallel_for(n, nthreads, [=](size_t r0, size_t r1) {
        for (size_t i = 0; i < iters; ++i)
            for (size_t r = r0; r < r1; ++r) out[r] = Dist(db + r * w, x, w);
    });
}

static void fill(uint64_t* rows, size_t r0, size_t r1, size_t w, size_t b, uint64_t seed) {
    std::mt19937_64 rng(seed);
    // Bits past width b in the last word stay 0, so distances stay in [0, b].
    uint64_t tail = (b % 64) ? ((uint64_t(1) << (b % 64)) - 1) : ~uint64_t(0);
    for (size_t r = r0; r < r1; ++r) {
        uint64_t* row = rows + r * w;
        for (size_t i = 0; i < w; ++i) row[i] = rng();
        row[w - 1] &= tail;
    }
}

static size_t parse_size(const char* s) {
    char* end;
    size_t v = strtoull(s, &end, 10);
    size_t base = (end[0] && (end[1] == 'i' || end[1] == 'I')) ? 1024 : 1000;
    switch (*end) {
        case 0: return v;
        case 'k': case 'K': return v * base;
        case 'm': case 'M': return v * base * base;
        case 'g': case 'G': return v * base * base * base;
        default: return 0;
    }
}

static void* xmalloc(size_t bytes) {
    void* p = malloc(bytes);
    if (!p) { fprintf(stderr, "malloc failed (%.2f MiB)\n", bytes / 1048576.0); exit(1); }
    return p;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <bits> <rows> [iters] [threads]\n"
                        "example: %s 1024 1M 10 8\n", argv[0], argv[0]);
        return 1;
    }
    size_t b = parse_size(argv[1]), n = parse_size(argv[2]);
    size_t iters = argc > 3 ? strtoull(argv[3], nullptr, 10) : 10;
    size_t nthreads = argc > 4 ? strtoull(argv[4], nullptr, 10) : 0;
    if (!nthreads) nthreads = std::thread::hardware_concurrency();
    if (!b || !n) { fprintf(stderr, "bits and rows must be > 0\n"); return 1; }
    if (!iters) iters = 1;
    if (nthreads > n) nthreads = n;

    size_t w = words(b);
    if (n > SIZE_MAX / (w * sizeof(uint64_t))) { fprintf(stderr, "db too large\n"); return 1; }
    size_t bytes = n * w * sizeof(uint64_t);

    uint64_t* db = (uint64_t*)xmalloc(bytes);
    uint64_t* x = (uint64_t*)xmalloc(w * sizeof(uint64_t));
    uint32_t* out = (uint32_t*)xmalloc(n * sizeof(uint32_t));

    printf("b=%zu bits, n=%zu rows, %zu words/row, db=%.2f MiB, threads=%zu\n",
           b, n, w, bytes / 1048576.0, nthreads);

    // Filled with the same sharding the scan uses, so each thread first-touches
    // the pages it will later read.
    parallel_for(n, nthreads, [=](size_t r0, size_t r1) { fill(db, r0, r1, w, b, 42 + r0); });
    fill(x, 0, 1, w, b, 12345);

    Rapl rapl;
    if (rapl.domains.empty())
        printf("rapl            : unavailable (no readable powercap counters; "
               "needs an Intel host and usually root)\n");

    double pkg_j, dram_j;
    rapl.start();
    auto t0 = std::chrono::steady_clock::now();
    scan<dist_vec>(db, x, n, w, out, nthreads, iters);
    double ms = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - t0).count() / iters;
    rapl.stop(&pkg_j, &dram_j);
    pkg_j /= iters;
    dram_j /= iters;

    printf("%-16s: %8.3f ms/iter  %6.2f GiB/s\n", kVecName, ms,
           bytes / (ms * 1e-3) / 1073741824.0);
    if (!rapl.domains.empty()) {
        printf("%16s  pkg %8.4f J/iter  %6.2f W", "", pkg_j, pkg_j / (ms * 1e-3));
        if (rapl.ndram)
            printf("   dram %8.4f J/iter  %6.2f W", dram_j, dram_j / (ms * 1e-3));
        else
            printf("   dram n/a");
        printf("\n");
    }
    for (size_t r = 0; r < n && r < 5; ++r) printf("dist(x, row[%zu]) = %u\n", r, out[r]);

    free(db); free(x); free(out);
    return 0;
}
