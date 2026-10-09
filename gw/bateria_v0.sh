#!/bin/bash
# Bateria do executavel ORIGINAL (commit d4c8408, so' o limite 300->1250) no
# 1x2iron.in, -np dado na linha de comando. Grava tempo em gw/bateria_v0.csv.
cd /home/madarame/Projetos/MSCDATA || exit 1
bin=gw/v0/1250/randmscd_parallel
for np in "$@"; do
  sed "s|gw/iron.out|gw/iron_v0np$np.out|" gw/iron.in > gw/iron_v0np$np.in
  carga=$(cut -d' ' -f1 /proc/loadavg)
  t0=$(date +%s.%N)
  mpirun --use-hwthread-cpus -np $np $bin gw/iron_v0np$np.in > gw/iron_v0np$np.log 2> gw/iron_v0np$np.err
  rc=$?
  t1=$(date +%s.%N)
  w=$(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.1f", b-a}')
  rf=$(grep -h "r-factora =" gw/iron_v0np$np.out | awk '{print $3, $6}')
  echo "$np,$w,$rc,$carga,$rf" >> gw/bateria_v0.csv
  rm -f gw/iron_v0np$np.in
done
