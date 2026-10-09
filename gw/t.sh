#!/bin/bash
# uso: gw/t.sh <tag> [binario] [casos]   compara saida e log com as referencias
cd /home/madarame/Projetos/MSCDATA || exit 1
tag=$1; bin=${2:-./randmscd_gpu}; casos=${3:-"cov iron c963"}
filt() { grep -vE "calculated by|took|starting|ending|gw/|Job report|result saved|seconds|minutes"; }
ok=1
for f in $casos; do
  [ $f = c963 ] && [ ! -f gw/c963.in ] && continue
  gw/r.sh $bin $f $tag 1 gpu
  cmp -s <(filt < gw/${f}_refgpu.out) <(filt < gw/${f}_$tag.out) || { echo "  $f SAIDA DIFERE"; ok=0; }
  ref=gw/${f}_refgpu.log; [ $f = c963 ] && ref=gw/c963_ref.log
  diff -q <(filt < $ref | grep -v "^ *[0-9.]*% of") <(filt < gw/${f}_$tag.log | grep -v "^ *[0-9.]*% of") >/dev/null || { echo "  $f LOG DIFERE"; ok=0; }
done
[ $ok = 1 ] && echo "TUDO BIT A BIT ($tag)"
