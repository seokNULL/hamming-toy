# hamming-toy

Data width가 `b` 비트, DB row 수가 `n`일 때, 입력 `x`와 DB의 모든 row 간
Hamming distance를 계산하는 C++ 토이 프로젝트.

## 설계

- 각 row는 `ceil(b/64)`개의 `uint64_t` 워드로 bit-packing되고, DB 전체는
  `malloc`으로 할당한 하나의 평탄한 버퍼(`n * ceil(b/64)` 워드)에 저장된다.
- 거리 계산은 `XOR` 후 popcount. 컴파일 타임에 가능한 가장 빠른 경로가 선택된다:
  - **scalar**: `__builtin_popcountll` (CPU `POPCNT` 명령)
  - **AVX2**: `vpshufb` 니블 룩업 popcount + `vpsadbw` 누적 (256비트/스텝)
  - **AVX-512**: `VPOPCNTQ` 네이티브 popcount (512비트/스텝,
    `-mavx512vpopcntdq` 지원 시)
- `b`가 64의 배수가 아니면 각 row 마지막 워드의 상위 비트는 0으로 유지되어
  거리가 항상 `[0, b]` 범위에 들어간다.
- 벡터 로드는 unaligned load(`loadu`)를 사용하므로 `malloc` 정렬로 충분하다.

## 빌드 및 실행

```sh
make                    # g++ -O3 -march=native
./hamming <b> <n> [iters]

./hamming 1024 100000   # 1024-bit rows, 100K rows
```

실행하면 scalar/vector 두 경로를 모두 돌려 결과를 교차 검증하고
pass당 소요 시간을 출력한다.
