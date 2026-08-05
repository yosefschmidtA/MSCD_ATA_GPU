#!/bin/bash
# Campanha pareada V5 (CPU) x GPU, na mesma janela, alternada.
#
# POR QUE EXISTE: vigiar medição à mão foi o que mais gastou orçamento na
# sessão de 05/08/2026. Rode isto, vá fazer outra coisa, leia o CSV no fim.
#
# Uso:  ./baseline/campanha-gpu.sh [reps] [np_cpu]
#       ./baseline/campanha-gpu.sh 3 12      # padrão
#
# O braço de CPU roda no np que você pedir (12 é o melhor da CPU hoje); o de
# GPU roda sempre np=1, porque com uma placa np>1 não faz sentido para o laço.
# ALTERNADO de propósito: a deriva térmica desta máquina chega a 11% entre
# repetições, e pareamento é a única coisa que separa efeito de deriva.
set -u
REPS=${1:-3}
NPCPU=${2:-12}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
CSV=$ROOT/baseline/gpu-campanha-$(date +%Y%m%d-%H%M).csv
cd "$ROOT" || exit 1

for b in randmscd_parallel randmscd_gpu; do
  [ -x "$b" ] || { echo "falta $b -- veja o CLAUDE.md, seção Compilar e rodar"; exit 1; }
done
n=$(ps -eo args | grep -cE "^mpirun|^randmscd")
[ "$n" -gt 0 ] && { echo "ja ha $n mpirun/randmscd rodando -- aborte ou espere"; exit 1; }
echo "load antes: $(cat /proc/loadavg)"

echo "rep,versao,np,wall_segundos" > "$CSV"
run() { /usr/bin/time -f "%e" -o "$ROOT/.t.$$" env ${4:-X=1} \
  mpirun --use-hwthread-cpus --bind-to none -np "$2" "$1" Cov0.txt \
  >/dev/null 2>&1; cat "$ROOT/.t.$$"; }

echo "aquecimento (descartado)..."
run ./randmscd_parallel "$NPCPU" >/dev/null

for r in $(seq 1 "$REPS"); do
  a=$(run ./randmscd_parallel "$NPCPU"); sleep 15
  b=$(run ./randmscd_gpu 1 "" "MSCD_GPU=1"); sleep 15
  echo "$r,V5,$NPCPU,$a" >> "$CSV"
  echo "$r,GPU,1,$b"     >> "$CSV"
  echo "rep$r  V5(np=$NPCPU)=${a}s  GPU(np=1)=${b}s"
done
rm -f "$ROOT/.t.$$"

echo "--- minimos (a leitura oficial deste projeto) ---"
awk -F, 'NR>1 {if (!(($2) in m) || $4+0<m[$2]) m[$2]=$4+0}
  END { printf "  V5  = %.2f s\n  GPU = %.2f s\n", m["V5"], m["GPU"];
        printf "  razao GPU/V5 = %.2fx\n", m["V5"]/m["GPU"] }' "$CSV"
echo "bruto em $CSV"
echo "load depois: $(cat /proc/loadavg)"
