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
- **멀티스레딩**: row 단위로 샤딩한다. 각 스레드는 연속된 row 구간을 맡아
  자기 구간의 출력만 쓰므로 계산 경로에 락도 atomic도 없다. 스레드는 pass마다
  새로 만들지 않고 상주 워커 풀(`ThreadPool`)을 재사용하므로, 작은 DB에서
  스레드 생성 비용이 측정을 왜곡하지 않는다.
- DB 초기화도 계산과 **동일한 샤딩**으로 병렬 수행한다. NUMA 장비에서 각
  스레드가 나중에 읽을 페이지를 자기가 first-touch 하게 되어 지역성이 좋다.
- `b`가 64의 배수가 아니면 각 row 마지막 워드의 상위 비트는 0으로 유지되어
  거리가 항상 `[0, b]` 범위에 들어간다.
- 벡터 로드는 unaligned load(`loadu`)를 사용하므로 `malloc` 정렬로 충분하다.

## 빌드 및 실행

```sh
make                       # g++ -O3 -march=native -pthread
./hamming <b> <n> [iters] [threads]
```

`b`와 `n`은 크기 접미사를 받는다:

| 표기 | 배수 | 예 |
|---|---|---|
| `K` / `M` / `G` | 1000 / 1000² / 1000³ (SI) | `1M` = 1,000,000 |
| `Ki` / `Mi` / `Gi` | 1024 / 1024² / 1024³ (IEC) | `1Mi` = 1,048,576 |

끝의 `B`는 무시된다(`512KiB` = 524288). `threads`를 생략하거나 `0`으로 주면
`hardware_concurrency()`를 쓴다.

```sh
./hamming 1024 1K          # 1024-bit rows, 1,000 rows
./hamming 1024 1M 10       # 1M rows, 10회 반복
./hamming 1024 8Mi 5 4     # ~1 GiB DB, 4 스레드
```

실행하면 scalar/vector 두 경로를 모두 돌려 전 row 결과를 교차 검증하고,
pass당 소요 시간과 메모리 처리량(GiB/s)을 출력한다.

## 측정 예 (4코어, AVX-512 VPOPCNTQ, b=1024, n=1M)

| 스레드 | scalar | vector |
|---|---|---|
| 1 | 16.95 ms | 14.67 ms |
| 2 | 7.34 ms | 4.55 ms |
| 4 | 2.67 ms | 1.60 ms |

큰 DB에서는 메모리 대역폭에 묶이므로 벡터 경로의 이득이 줄어든다.
