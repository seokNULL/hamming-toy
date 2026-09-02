# hamming-toy

Data width가 `b` 비트, DB row 수가 `n`일 때, 입력 `x`와 DB의 모든 row 간
Hamming distance를 계산하는 C++ 토이 프로젝트.

## 설계

- 각 row는 `ceil(b/64)`개의 `uint64_t` 워드로 bit-packing되고, DB 전체는
  `malloc`으로 할당한 하나의 평탄한 버퍼에 저장된다.
- 거리 계산은 `XOR` 후 popcount. 컴파일 타임에 가장 빠른 경로 하나가 선택되어
  그것만 실행된다:
  - **AVX-512**: `VPOPCNTQ` (512비트/스텝, `-mavx512vpopcntdq` 지원 시)
  - **AVX2**: `vpshufb` 니블 룩업 + `vpsadbw` 누적 (256비트/스텝)
  - **scalar**: `__builtin_popcountll` (CPU `POPCNT`), SIMD가 없을 때의 폴백
- 멀티스레딩은 `parallel_for` 하나로 처리한다. `[0, n)`을 스레드마다 연속된
  구간으로 나누고 각자 자기 구간의 출력만 쓰므로 락도 atomic도 없다. DB
  초기화도 같은 분할로 돌려 각 스레드가 나중에 읽을 페이지를 first-touch 한다.
- `b`가 64의 배수가 아니면 각 row 마지막 워드의 남는 비트는 0으로 유지되어
  거리가 항상 `[0, b]`에 들어간다.
- 벡터 로드는 unaligned(`loadu`)라서 `malloc` 정렬로 충분하다.

## 빌드 및 실행

```sh
make                       # g++ -O3 -march=native -pthread
./hamming <bits> <rows> [iters] [threads]
```

`bits`와 `rows`는 크기 접미사를 받는다: `K`/`M`/`G`는 1000 기반(SI),
`Ki`/`Mi`/`Gi`는 1024 기반(IEC). `threads`를 생략하면
`hardware_concurrency()`를 쓴다.

```sh
./hamming 1024 1K          # 1024-bit rows, 1,000 rows
./hamming 1024 1M 10       # 1M rows, 10회 반복
./hamming 1024 8Mi 5 4     # ~1 GiB DB, 4 스레드
```

선택된 경로만 `iters`번 돌려 iter당 시간과 메모리 처리량을 출력한다. 반복은
병렬 구간 **안**에서 돌기 때문에 측정값에 스레드 생성 비용이 섞이지 않는다.

## 스윕 스크립트

`sweep.sh`는 bit width, DB row 수, 스레드 수를 리스트로 받아 모든 조합을
실행하고 실행마다 시간과 에너지를 한 줄씩 출력한다.

```sh
./sweep.sh
B_LIST="256 1024" N_LIST="1M 8Mi" T_LIST="1 2 4 8" ITERS=20 ./sweep.sh
CSV=out.csv ./sweep.sh          # 표와 함께 CSV로도 저장
```

| 변수 | 기본값 | 뜻 |
|---|---|---|
| `B_LIST` | `256 512 1024` | bit width 목록 |
| `N_LIST` | `100K 1M` | row 수 목록 (바이너리와 같은 접미사) |
| `T_LIST` | `1 2 4` | 스레드 수 목록 |
| `ITERS` | `10` | 조합당 반복 횟수 |
| `CSV` | (없음) | 지정하면 결과를 CSV로도 기록 |
| `BIN` | `./hamming` | 실행할 바이너리 |

```
bits   rows     T   kernel              ms/iter    GiB/s     pkg_J    pkg_W    dram_J   dram_W
1024   200K     1   AVX2 vpshufb          2.782     8.57    0.2774    99.71    0.0445    16.00
1024   200K     2   AVX2 vpshufb          0.965    24.71    0.0964    99.96    0.0154    15.99
1024   200K     4   AVX2 vpshufb          0.449    53.13    0.0450   100.30    0.0072    16.05
```

T 컬럼은 요청한 값이 아니라 바이너리가 실제로 쓴 스레드 수다(`n`보다 크면
줄어든다). RAPL을 읽을 수 없으면 에너지 컬럼은 `-`, DRAM 도메인이 없으면
`n/a`가 찍힌다. 어떤 조합이 실패해도 스윕은 계속되고 그 행은 `FAILED`로
표시되며 이유는 stderr로 나간다.

## 에너지 측정 (Intel RAPL)

powercap sysfs(`/sys/class/powercap/intel-rapl:*`)에서 package와 DRAM 에너지를
읽어 iter당 J와 평균 W를 함께 출력한다. 소켓이 여러 개면 전부 합산하고,
`core`/`psys` 도메인은 제외한다. `RAPL_PATH`로 sysfs 경로를 바꿀 수 있다.

```
AVX2 vpshufb    :    0.401 ms/iter   59.42 GiB/s
                  pkg   0.0400 J/iter   99.78 W   dram   0.0064 J/iter   15.96 W
```

주의할 점:

- 카운터를 읽으려면 보통 root 권한이 필요하고, Intel 호스트여야 한다. 읽을 수
  없으면 `rapl: unavailable`을 출력하고 시간 측정만 진행한다. VM/컨테이너에는
  대개 노출되지 않는다.
- 값은 **패키지 전체** 소비량이라 같은 소켓에서 도는 다른 프로세스도 포함된다.
- 카운터는 `max_energy_range_uj`에서 순환한다. 한 번의 순환은 복구하지만,
  측정 구간의 실제 증가량이 이 범위를 넘으면 구분할 방법이 없어 과소 보고된다.
  실제 범위는 수십 kJ(100 W에서 수십 분)라 통상적인 실행에서는 문제되지 않는다.
- DRAM 도메인이 없는 CPU(주로 데스크톱)에서는 `dram n/a`로 표시된다.
