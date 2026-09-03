#!/usr/bin/env bash
# Sweeps ./hamming over bit widths, DB row counts, and thread counts,
# printing runtime and RAPL energy for every run.
#
#   ./sweep.sh
#   B_LIST="256 1024" N_LIST="1M 8Mi" T_LIST="1 2 4 8" ./sweep.sh
#   CSV=out.csv ./sweep.sh
#
# Row counts take the binary's size suffixes (1K/1M/1G, 1Ki/1Mi/1Gi).
# The binary times every kernel it was built with; the sweep reports the fastest
# of them. Each combination runs REPS times and the fastest run is reported;
# single runs on a shared machine vary by tens of percent.
# RAPL usually needs root; without it the energy columns show "-".

set -u

BIN=${BIN:-./hamming}
B_LIST=${B_LIST:-"256 512 1024"}
N_LIST=${N_LIST:-"100K 1M"}
T_LIST=${T_LIST:-"1 2 4"}
ITERS=${ITERS:-10}
REPS=${REPS:-3}
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
      best= bestout= bestgib=-1 failed=
      for _ in $(seq "$REPS"); do
        if ! out=$("$BIN" "$b" "$n" "$ITERS" "$t" 2>&1); then failed=$out; break; fi
        parsed=$(echo "$out" | awk '
        function flush() {
          if (!pending) return
          if (cg + 0 > bg + 0) { bg = cg; bk = ck; bm = cm
                                 bpj = cpj; bpw = cpw; bdj = cdj; bdw = cdw }
          pending = 0
        }
        # One block per kernel: a "ms/iter" line, optionally an energy line.
        /ms\/iter/ {
          flush()
          split($0, a, ":"); ck = a[1]; sub(/ +$/, "", ck)
          for (i = 1; i <= NF; i++) {
            if ($i == "ms/iter") cm = $(i-1)
            if ($i == "GiB/s")   cg = $(i-1)
          }
          cpj = "-"; cpw = "-"; cdj = "-"; cdw = "-"; pending = 1
        }
        / pkg / {
          for (i = 1; i <= NF; i++) {
            if ($i == "pkg")  { cpj = $(i+1); cpw = $(i+3) }
            if ($i == "dram") { if ($(i+1) == "n/a") { cdj = "n/a"; cdw = "n/a" }
                                else                 { cdj = $(i+1); cdw = $(i+3) } }
          }
        }
        END {
          flush()
          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", bk, bm, bg, bpj, bpw, bdj, bdw
        }')
        g=$(echo "$parsed" | cut -f3)
        if awk -v a="$g" -v c="$bestgib" 'BEGIN{exit !(a>c)}'; then
          bestgib=$g; best=$parsed; bestout=$out
        fi
      done
      if [ -n "$failed" ]; then
        # shellcheck disable=SC2059
        printf "$row" "$b" "$n" "$t" FAILED - - - - - -
        echo "  -> $(echo "$failed" | tail -1)" >&2
        continue
      fi
      IFS=$'\t' read -r kern ms gib pj pw dj dw <<<"$best"
      # Report the thread count the binary actually used (it clamps T to n).
      at=$(echo "$bestout" | sed -n 's/.*threads=\([0-9]*\).*/\1/p')
      # shellcheck disable=SC2059
      printf "$row" "$b" "$n" "${at:-$t}" "$kern" "$ms" "$gib" "$pj" "$pw" "$dj" "$dw"
      if [ -n "$CSV" ]; then
        echo "$b,$n,${at:-$t},$kern,$ms,$gib,$pj,$pw,$dj,$dw" >> "$CSV"
      fi
    done
  done
done
