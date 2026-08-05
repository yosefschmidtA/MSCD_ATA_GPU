#!/bin/bash
# Escalabilidade do V5, para a figura V0 x V5 do OTIMIZACAO.md.
# Grava em baseline/escala.csv com versao=V5, mesmo esquema das campanhas
# anteriores.
#
# Binding: np=1 com --bind-to none (e' onde o OpenMP do V5 age, e sem a flag
# as 12 threads ficam empilhadas num nucleo -- ver "V5" no OTIMIZACAO.md);
# np>=2 com binding padrao, que e' a invocacao de producao. Medido em
# 05/08/2026, np=12: padrao 39,00 s contra bindnone 42,08 s -- desamarrar os
# ranks nao paga quando nao ha nucleo livre para o OpenMP pegar.
#
# Uso: ./baseline/campanha-v5-escala.sh [repeticoes]
# Requer binario COM -fopenmp.
set -u
REP=${1:-2}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
CSV=$ROOT/baseline/escala.csv
OUT=$ROOT/saida1Co-alterado-alexandre.txt
# Etiqueta com hora, nao so' o dia: em 05/08/2026 duas campanhas do mesmo dia
# (01:00 e 12:50) diferiram 5% no mesmo binario. Janela e' o intervalo em que a
# maquina esta no mesmo estado, nao a data.
JANELA=$(date +%Y-%m-%d-%H%M)

cd "$ROOT" || exit 1

n=$(ps -eo args | grep -cE "^mpirun|^randmscd") || true
if [ "$n" -ne 0 ]; then
  echo "ABORTADO: ha $n processo(s) mpirun/randmscd rodando."
  exit 1
fi
echo "load inicial: $(cat /proc/loadavg)"

corpo() { grep -v -E "This calculation took|starting on|and ending on|calculated by" "$1"; }

echo "--- aquecimento (descartado) ---"
mpirun --use-hwthread-cpus -np 4 randmscd_parallel Cov0.txt >/dev/null 2>&1

for r in $(seq 1 "$REP"); do
  for NP in 1 2 4 6 8 12; do
    if [ "$NP" -eq 1 ]; then flags="--bind-to none"; else flags=""; fi
    sleep 40
    rm -f "$OUT"
    out=$( { /usr/bin/time -f "%e" mpirun --use-hwthread-cpus $flags \
             -np "$NP" randmscd_parallel Cov0.txt >/dev/null; } 2>&1 | tail -1 )
    rc=$?
    if [ $rc -ne 0 ] || [ ! -s "$OUT" ]; then
      echo "FALHOU np=$NP rep=$r"; continue
    fi
    corpo "$ROOT/baseline/saida.txt" > /tmp/e5.a.$$
    corpo "$OUT" > /tmp/e5.b.$$
    if cmp -s /tmp/e5.a.$$ /tmp/e5.b.$$; then curva=identica; else curva=DIFERE; fi
    rm -f /tmp/e5.a.$$ /tmp/e5.b.$$
    printf "np=%-3s rep%-2s  %8s s  curva %s\n" "$NP" "$r" "$out" "$curva"
    echo "$JANELA,V5,$NP,$r,$out,$curva" >> "$CSV"
  done
done

echo "gravado em $CSV"
echo "load final: $(cat /proc/loadavg)"
