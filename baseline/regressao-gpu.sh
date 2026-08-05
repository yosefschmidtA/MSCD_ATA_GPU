#!/bin/bash
# Regressao do caminho de GPU. Uso: ./baseline/regressao-gpu.sh [np]
#
# POR QUE NAO E' O regressao.sh: a validacao byte a byte nao sobrevive ao
# port -- a reducao soma em outra ordem e em float isso mexe no ultimo bit.
# O criterio abaixo foi fixado ANTES de existir kernel (PLANO_CUDA.md, Fase 0),
# medido invertendo a ordem da soma sobre al em evenelem, que e' fisicamente
# identico e muda so' o arredondamento:
#
#     piso de ruido medido:  max|dchi| = 1,0e-5   rms = 3,9e-7
#     criterio:              max|dchi| <= 1e-4    rms|dchi| <= 1e-5
#
# 1e-4 e' 10x o piso e 0,011% da amplitude da curva (chi vai de -0,49043 a
# +0,42295). Acima disso nao e' reordenacao de soma: e' bug.
#
# O regressao.sh byte a byte continua valendo para mudancas de CPU e NAO deve
# ser relaxado -- perder o teste forte no CPU seria pagar duas vezes.
set -u
NP=${1:-1}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BASE=$ROOT/baseline
OUT=$ROOT/saida1Co-alterado-alexandre.txt
TMP=${TMPDIR:-/tmp}/regrgpu.$$

# chical e' a coluna 4 das linhas de dados (phi intensity background chical chiexp)
chi() { awk 'NF==5 && $1+0==$1 && $2 ~ /[eE]/ {print $4}' "$1"; }

cd "$ROOT" || exit 1
[ -f "$BASE/saida.txt" ] || { echo "sem linha de base em $BASE"; exit 1; }
[ -x "$ROOT/randmscd_gpu" ] || { echo "randmscd_gpu nao construido"; exit 1; }

MSCD_GPU=1 /usr/bin/time -f "WALL %e s  MAXRSS %M kB" -o "$TMP.time" \
  mpirun --use-hwthread-cpus --bind-to none -np "$NP" \
  randmscd_gpu Cov0.txt > "$TMP.stdout" 2>&1
rc=$?
echo "--- tempo (np=$NP, GPU) ---"; cat "$TMP.time"
[ $rc -ne 0 ] && { echo "FALHOU: exit=$rc"; tail -20 "$TMP.stdout"; exit 1; }

chi "$BASE/saida.txt" > "$TMP.a"; chi "$OUT" > "$TMP.b"
na=$(wc -l < "$TMP.a"); nb=$(wc -l < "$TMP.b")
if [ "$na" != "$nb" ] || [ "$na" -lt 100 ]; then
  echo "FALHOU: $na pontos na base, $nb agora"; exit 1
fi

paste "$TMP.a" "$TMP.b" | awk -v n="$na" '
{ d=$2-$1; if (d<0) d=-d; s+=d*d; if (d>m) m=d; if (d>1e-4) big++ }
END { rms=sqrt(s/n);
  printf "--- chi (%d pontos) ---\n", n;
  printf "  max|dchi| = %.3e   (criterio 1e-4)\n", m;
  printf "  rms|dchi| = %.3e   (criterio 1e-5)\n", rms;
  printf "  pontos acima de 1e-4: %d\n", big+0;
  if (m<=1e-4 && rms<=1e-5) { print "APROVADO"; exit 0 }
  print "REPROVADO"; exit 1 }'
st=$?
echo "--- factors (sanidade, NAO e' criterio) ---"
echo "  agora: $(grep -h 'factors =' "$TMP.stdout" | tail -1)"
echo "  base:  $(grep -h 'factors =' "$BASE/stdout.txt" | tail -1)"
[ $st -eq 0 ] && rm -f "$TMP".*
exit $st
