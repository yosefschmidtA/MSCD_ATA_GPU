#!/bin/bash
# uso: gw/run.sh <binario> <cov|iron> <tag> <np> [gpu]
# roda da raiz do projeto, grava gw/<input>_<tag>.out e .log
cd /home/madarame/Projetos/MSCDATA || exit 1
bin=$1; in=$2; tag=$3; np=$4; g=$5
sed "s|gw/$in.out|gw/${in}_$tag.out|" gw/$in.in > gw/${in}_$tag.in
t0=$(date +%s.%N)
if [ "$g" = gpu ]; then export MSCD_GPU=1; else unset MSCD_GPU; fi
mpirun --use-hwthread-cpus --bind-to none -np $np $bin gw/${in}_$tag.in > gw/${in}_$tag.log 2> gw/${in}_$tag.err
rc=$?
t1=$(date +%s.%N)
rm -f gw/${in}_$tag.in
echo "$in $tag np=$np rc=$rc wall=$(echo "$t1 - $t0" | bc) $(grep -h 'r-factora' gw/${in}_$tag.out)"
