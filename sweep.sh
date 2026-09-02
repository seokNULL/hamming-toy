#!/usr/bin/env bash
# Sweeps ./hamming over bit widths, DB row counts, and thread counts,
# printing runtime and RAPL energy for every run.
#
#   ./sweep.sh
#   B_LIST="256 1024" N_LIST="1M 8Mi" T_LIST="1 2 4 8" ./sweep.sh
#   CSV=out.csv ./sweep.sh
#
# Row counts take the binary's size suffixes (1K/1M/1G, 1Ki/1Mi/1Gi).
# RAPL usually needs root; without it the energy columns show "-".

set -u

BIN=${BIN:-./hamming}
B_LIST=${B_LIST:-"256 512 1024"}
N_LIST=${N_LIST:-"100K 1M"}
T_LIST=${T_LIST:-"1 2 4"}
ITERS=${ITERS:-10}
CSV=${CSV:-}

[ -x "$BIN" ] || make || exit 1

hdr="%-6s %-8s %-3s %-17s %9s %8s %9s %8s %9s %8s\n"
row="%-6s %-8s %-3s %-17s %9s %8s %9s %8s %9s %8s\n"
# shellcheck disable=SC2059
printf "$hdr" bits rows T kernel ms/iter GiB/s pkg_J pkg_W dram_J dram_W
if [ -n "$CSV" ]; then
  echo "bits,rows,threads,kernel,ms_per_iter,gib_per_s,pkg_j,pkg_w,dram_j,dram_w" > "$CSV"
fi

for b in $B_LIST; do
  for n in $N_LIST; do
    for t in $T_LIST; do
      if ! out=$("$BIN" "$b" "$n" "$ITERS" "$t" 2>&1); then
        # shellcheck disable=SC2059
        printf "$row" "$b" "$n" "$t" FAILED - - - - - -
        echo "  -> $(echo "$out" | tail -1)" >&2
        continue
      fi
      IFS=$'\t' read -r kern ms gib pj pw dj dw <<<"$(echo "$out" | awk '
        /ms\/iter/ {
          split($0, a, ":"); kern = a[1]; sub(/ +$/, "", kern)
          for (i = 1; i <= NF; i++) {
            if ($i == "ms/iter") ms  = $(i-1)
            if ($i == "GiB/s")   gib = $(i-1)
          }
        }
        / pkg / {
          for (i = 1; i <= NF; i++) {
            if ($i == "pkg")  { pj = $(i+1); pw = $(i+3) }
            if ($i == "dram") { if ($(i+1) == "n/a") { dj = "n/a"; dw = "n/a" }
                                else                 { dj = $(i+1); dw = $(i+3) } }
          }
        }
        END {
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", kern, ms, gib,
                 (pj == "" ? "-" : pj), (pw == "" ? "-" : pw),
                 (dj == "" ? "-" : dj), (dw == "" ? "-" : dw)
        }')"
      # Report the thread count the binary actually used (it clamps T to n).
      at=$(echo "$out" | sed -n 's/.*threads=\([0-9]*\).*/\1/p')
      # shellcheck disable=SC2059
      printf "$row" "$b" "$n" "${at:-$t}" "$kern" "$ms" "$gib" "$pj" "$pw" "$dj" "$dw"
      if [ -n "$CSV" ]; then
        echo "$b,$n,${at:-$t},$kern,$ms,$gib,$pj,$pw,$dj,$dw" >> "$CSV"
      fi
    done
  done
done
