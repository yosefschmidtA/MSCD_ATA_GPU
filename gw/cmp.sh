#!/bin/bash
# uso: gw/cmp.sh a.out b.out  -> compara chical (coluna 4) e r-factor
chi() { awk 'NF==5 && $1+0==$1 && $2 ~ /[eE]/ {print $4}' "$1"; }
paste <(chi $1) <(chi $2) | awk '{d=$2-$1; if(d<0)d=-d; s+=d*d; if(d>m)m=d; if(d>1e-4)b++} END{printf "n=%d max|dchi|=%.3e rms=%.3e acima1e-4=%d\n", NR, m, sqrt(s/NR), b+0}'
echo "  A: $(grep -h 'r-factora' $1)"; echo "  B: $(grep -h 'r-factora' $2)"
